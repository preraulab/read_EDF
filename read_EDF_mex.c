/**
 * READ_EDF_MEX  High-performance EDF / EDF+ file reader for MATLAB
 *
 * Supports plain .edf, gzip-compressed .edf.gz, and zstd-compressed
 * .edf.zst inputs. Compressed files are decoded on the fly via zlib /
 * libzstd — no temp file is created. Reads stream record-by-record, so
 * peak memory is one record's worth of raw bytes plus the per-signal
 * output arrays.
 *
 * MATLAB usage:
 *   [header, sigheader, data, annotations] =
 *       read_EDF_mex(filename, channels, epochs, verbose, repair, debug)
 *
 * Compilation:
 *   mex -O -largeArrayDims read_EDF_mex.c -lz -lzstd
 *
 *   -largeArrayDims is REQUIRED.
 *   -lz links against system zlib (present on macOS via the Xcode SDK and
 *   on most Linux distributions). -lzstd links against libzstd (Homebrew
 *   on macOS, libzstd-dev on Linux). compile_edf_mex.m wires up the
 *   include/library paths automatically.
 *
 * Notes on compressed inputs:
 *   • RepairHeader is silently skipped for .gz / .zst inputs — we cannot
 *     rewrite a single byte of a compressed archive in place. The header
 *     is still corrected in the returned struct.
 *   • num_data_records reported in the EDF main header is trusted for
 *     header-only calls (nargout == 1). When data is requested, the actual
 *     record count is determined by streaming and the header is updated.
 */

#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>
#include <zlib.h>
#include <zstd.h>

#if defined(MX_COMPAT_32)
#error "This MEX requires -largeArrayDims"
#endif

#define EDF_HEADER_SIZE 256
#define MAX_LABEL_LEN   17   /* 16 chars + NUL */
#define MAX_CHANNELS    512

/* ============================================================
 *                         DEBUG MACROS
 * ============================================================ */
static int g_debug = 0;
#define DBG(...) do { if (g_debug) mexPrintf(__VA_ARGS__); } while (0)
#define ASSERT(cond) \
    do { if (g_debug && !(cond)) mexErrMsgIdAndTxt("read_EDF_mex:Assert", #cond); } while (0)

/* ============================================================
 *                         IO ABSTRACTION
 * ============================================================
 *
 * Wraps a plain FILE*, a zlib gzFile, or a zstd streaming decoder so
 * the rest of the code reads identically. Compression is detected from
 * the file extension on open (.gz -> zlib, .zst -> zstd, else plain).
 */

typedef enum { COMP_NONE = 0, COMP_GZ = 1, COMP_ZSTD = 2 } CompType;

typedef struct {
    CompType comp;
    FILE   *fp;            /* used for plain reads, and as raw source for zstd */
    gzFile  gz;

    /* zstd streaming state (only used when comp == COMP_ZSTD) */
    ZSTD_DStream *zds;
    unsigned char *z_in_buf;     /* compressed bytes read from fp */
    size_t        z_in_cap;
    ZSTD_inBuffer z_in;
    unsigned char *z_out_buf;    /* decompressed bytes ready for caller */
    size_t        z_out_cap;
    ZSTD_outBuffer z_out;
    size_t        z_out_consumed; /* bytes already returned from z_out_buf */
    int           z_eof;
} EDFReader;

static CompType detect_comp(const char *fname)
{
    size_t n = strlen(fname);
    if (n >= 3 &&
        fname[n-3] == '.' &&
        tolower((unsigned char)fname[n-2]) == 'g' &&
        tolower((unsigned char)fname[n-1]) == 'z')
        return COMP_GZ;
    if (n >= 4 &&
        fname[n-4] == '.' &&
        tolower((unsigned char)fname[n-3]) == 'z' &&
        tolower((unsigned char)fname[n-2]) == 's' &&
        tolower((unsigned char)fname[n-1]) == 't')
        return COMP_ZSTD;
    return COMP_NONE;
}

static int edfr_open(EDFReader *r, const char *fname)
{
    memset(r, 0, sizeof(*r));
    r->comp = detect_comp(fname);

    if (r->comp == COMP_GZ) {
        r->gz = gzopen(fname, "rb");
        return (r->gz != NULL);
    }
    if (r->comp == COMP_ZSTD) {
        r->fp = fopen(fname, "rb");
        if (!r->fp) return 0;
        r->zds = ZSTD_createDStream();
        if (!r->zds) { fclose(r->fp); r->fp = NULL; return 0; }
        ZSTD_initDStream(r->zds);
        r->z_in_cap  = ZSTD_DStreamInSize();
        r->z_out_cap = ZSTD_DStreamOutSize();
        r->z_in_buf  = (unsigned char*)malloc(r->z_in_cap);
        r->z_out_buf = (unsigned char*)malloc(r->z_out_cap);
        if (!r->z_in_buf || !r->z_out_buf) return 0;
        r->z_in.src  = r->z_in_buf; r->z_in.size = 0; r->z_in.pos = 0;
        r->z_out.dst = r->z_out_buf; r->z_out.size = r->z_out_cap; r->z_out.pos = 0;
        r->z_out_consumed = 0;
        r->z_eof = 0;
        return 1;
    }
    r->fp = fopen(fname, "rb");
    return (r->fp != NULL);
}

static void edfr_close(EDFReader *r)
{
    if (r->comp == COMP_GZ) {
        if (r->gz) { gzclose(r->gz); r->gz = NULL; }
    } else if (r->comp == COMP_ZSTD) {
        if (r->zds) { ZSTD_freeDStream(r->zds); r->zds = NULL; }
        if (r->z_in_buf)  { free(r->z_in_buf);  r->z_in_buf  = NULL; }
        if (r->z_out_buf) { free(r->z_out_buf); r->z_out_buf = NULL; }
        if (r->fp) { fclose(r->fp); r->fp = NULL; }
    } else {
        if (r->fp) { fclose(r->fp); r->fp = NULL; }
    }
}

/* Pull more decompressed bytes into r->z_out_buf. Returns 1 if any are
 * available afterward, 0 on clean EOF, -1 on error. */
static int zstd_refill(EDFReader *r)
{
    /* If output buffer fully consumed, reset. */
    if (r->z_out_consumed >= r->z_out.pos) {
        r->z_out.pos = 0;
        r->z_out_consumed = 0;
    }
    while (r->z_out_consumed >= r->z_out.pos) {
        /* Need to decompress. Refill input if exhausted. */
        if (r->z_in.pos >= r->z_in.size && !r->z_eof) {
            size_t got = fread(r->z_in_buf, 1, r->z_in_cap, r->fp);
            r->z_in.src  = r->z_in_buf;
            r->z_in.size = got;
            r->z_in.pos  = 0;
            if (got == 0) r->z_eof = 1;
        }
        if (r->z_in.pos >= r->z_in.size && r->z_eof) {
            return 0;  /* clean EOF */
        }
        /* Decompress some of input into out buffer. */
        r->z_out.pos = 0;
        r->z_out_consumed = 0;
        size_t ret = ZSTD_decompressStream(r->zds, &r->z_out, &r->z_in);
        if (ZSTD_isError(ret)) {
            mexPrintf("zstd decompress error: %s\n", ZSTD_getErrorName(ret));
            return -1;
        }
        if (r->z_out.pos > 0) return 1;
        /* No output produced this iteration; loop refills input or hits EOF. */
    }
    return 1;
}

/* Returns number of bytes actually read (0 on EOF or error). */
static size_t edfr_read(EDFReader *r, void *buf, size_t n)
{
    if (r->comp == COMP_GZ) {
        int got = gzread(r->gz, buf, (unsigned)n);
        if (got <= 0) return 0;
        return (size_t)got;
    }
    if (r->comp == COMP_ZSTD) {
        size_t total = 0;
        while (total < n) {
            int rc = zstd_refill(r);
            if (rc <= 0) break;
            size_t avail = r->z_out.pos - r->z_out_consumed;
            size_t want  = n - total;
            size_t take  = avail < want ? avail : want;
            memcpy((unsigned char*)buf + total,
                   r->z_out_buf + r->z_out_consumed, take);
            r->z_out_consumed += take;
            total += take;
        }
        return total;
    }
    return fread(buf, 1, n, r->fp);
}

/* Absolute seek from start of (uncompressed) stream. Plain only. */
static int edfr_seek_set(EDFReader *r, mwSize off)
{
    if (r->comp == COMP_GZ) {
        return (gzseek(r->gz, (z_off_t)off, SEEK_SET) >= 0);
    }
    if (r->comp == COMP_ZSTD) {
        return 0;  /* no seek support for zstd streams */
    }
    return (fseek(r->fp, (long)off, SEEK_SET) == 0);
}

/* ============================================================
 *                         STRUCTS
 * ============================================================ */

typedef struct {
    char edf_ver[9];
    char patient_id[81];
    char local_rec_id[81];
    char recording_startdate[9];
    char recording_starttime[9];
    int  num_header_bytes;
    int  num_data_records;
    double data_record_duration;
    int  num_signals;
} EDF_Header;

typedef struct {
    char signal_labels[MAX_LABEL_LEN];
    char transducer_type[81];
    char physical_dimension[9];
    double physical_min, physical_max;
    double digital_min, digital_max;
    char prefiltering[81];
    int  samples_in_record;
    double sampling_frequency;
} Signal_Header;

typedef struct {
    double onset;
    mxArray *texts;
} Annotation;

/* ============================================================
 *                  CHANNEL PLAN STRUCTS
 * ============================================================ */

typedef enum { CH_PLAIN, CH_REREF } ChannelKind;

typedef struct {
    ChannelKind kind;
    char label[64];
    int  raw_idx_a;
    int  raw_idx_b;
} ChannelPlan;

/* ============================================================
 *                         UTILITIES
 * ============================================================ */

static void trim_string(char *s)
{
    char *p = s;
    while (*p == ' ') p++;
    memmove(s, p, strlen(p) + 1);
    size_t n = strlen(s);
    while (n && s[n-1] == ' ') s[--n] = '\0';
}

static int ci_strcmp(const char *a, const char *b)
{
    while (*a && *b) {
        if (tolower((unsigned char)*a) != tolower((unsigned char)*b)) return 1;
        a++; b++;
    }
    return (*a != *b);
}

static int find_signal(const Signal_Header *sig, int nsig, const char *name)
{
    for (int i = 0; i < nsig; i++) {
        if (ci_strcmp(sig[i].signal_labels, name) == 0)
            return i;
    }
    return -1;
}

/* ============================================================
 *                         HEADER IO
 * ============================================================ */

static int read_edf_header(EDFReader *rd, EDF_Header *h)
{
    unsigned char buf[EDF_HEADER_SIZE];
    char tmp[16];

    if (edfr_read(rd, buf, EDF_HEADER_SIZE) != EDF_HEADER_SIZE) return 0;

    memcpy(h->edf_ver,             buf,     8); h->edf_ver[8]             = '\0'; trim_string(h->edf_ver);
    memcpy(h->patient_id,          buf+8,  80); h->patient_id[80]         = '\0'; trim_string(h->patient_id);
    memcpy(h->local_rec_id,        buf+88, 80); h->local_rec_id[80]       = '\0'; trim_string(h->local_rec_id);
    memcpy(h->recording_startdate, buf+168, 8); h->recording_startdate[8] = '\0';
    memcpy(h->recording_starttime, buf+176, 8); h->recording_starttime[8] = '\0';

    memcpy(tmp, buf+184, 8); tmp[8]='\0'; h->num_header_bytes      = atoi(tmp);
    memcpy(tmp, buf+236, 8); tmp[8]='\0'; h->num_data_records      = atoi(tmp);
    memcpy(tmp, buf+244, 8); tmp[8]='\0'; h->data_record_duration  = atof(tmp);
    memcpy(tmp, buf+252, 4); tmp[4]='\0'; h->num_signals           = atoi(tmp);

    return 1;
}

static int read_signal_headers(EDFReader *rd, Signal_Header *sh, int nsig)
{
    char tmp[16];
    int i;

    for (i=0;i<nsig;i++) { edfr_read(rd, sh[i].signal_labels,   16); sh[i].signal_labels[16]  ='\0'; trim_string(sh[i].signal_labels);   }
    for (i=0;i<nsig;i++) { edfr_read(rd, sh[i].transducer_type,  80); sh[i].transducer_type[80] ='\0'; trim_string(sh[i].transducer_type); }
    for (i=0;i<nsig;i++) { edfr_read(rd, sh[i].physical_dimension,8); sh[i].physical_dimension[8]='\0'; }
    for (i=0;i<nsig;i++) { edfr_read(rd, tmp,8); tmp[8]='\0'; sh[i].physical_min=atof(tmp); }
    for (i=0;i<nsig;i++) { edfr_read(rd, tmp,8); tmp[8]='\0'; sh[i].physical_max=atof(tmp); }
    for (i=0;i<nsig;i++) { edfr_read(rd, tmp,8); tmp[8]='\0'; sh[i].digital_min =atof(tmp); }
    for (i=0;i<nsig;i++) { edfr_read(rd, tmp,8); tmp[8]='\0'; sh[i].digital_max =atof(tmp); }
    for (i=0;i<nsig;i++) { edfr_read(rd, sh[i].prefiltering, 80); sh[i].prefiltering[80]='\0'; }
    for (i=0;i<nsig;i++) { edfr_read(rd, tmp,8); tmp[8]='\0'; sh[i].samples_in_record=atoi(tmp); }
    for (i=0;i<nsig;i++) { edfr_read(rd, tmp,32); } /* reserved */

    return 1;
}

/* ============================================================
 *                 EDF+ num_data_records FIX (plain only)
 * ============================================================
 *
 * Plain files: derive record count from file size and optionally rewrite
 * byte 236 of the header. Leaves stream positioned at hdr.num_header_bytes.
 *
 * Gz files: skipped entirely. The streaming data loop discovers the true
 * record count and updates hdr.num_data_records after the read. RepairHeader
 * cannot be honored on a .gz — we warn if verbose.
 */
static void fix_num_records_plain(EDFReader *rd, EDF_Header *hdr,
                                  Signal_Header *sig,
                                  const char *fname,
                                  int verbose,
                                  int repair)
{
    mwSize total_samp_per_rec = 0;
    for (int i = 0; i < hdr->num_signals; i++)
        total_samp_per_rec += sig[i].samples_in_record;

    DBG("\n===== DERIVED RECORD INFO =====\n");
    DBG("Total samples per record: %llu\n", (unsigned long long)total_samp_per_rec);
    for (int i = 0; i < hdr->num_signals; i++)
        DBG("Signal %d samples_in_record: %d\n", i, sig[i].samples_in_record);
    DBG("===============================\n\n");

    ASSERT(total_samp_per_rec > 0);

    mwSize bytes_per_rec = total_samp_per_rec * 2;
    if (bytes_per_rec == 0) return;

    fseek(rd->fp, 0, SEEK_END);
    mwSize fsize = (mwSize)ftell(rd->fp);
    mwSize data_bytes = fsize - hdr->num_header_bytes;

    mwSize actual_records = data_bytes / bytes_per_rec;

    if (hdr->num_data_records <= 0 ||
        hdr->num_data_records != (int)actual_records)
    {
        if (verbose)
            mexPrintf("Updating num_data_records from %d to %llu\n",
                      hdr->num_data_records, (unsigned long long)actual_records);

        if (repair) {
            FILE *fw = fopen(fname, "r+b");
            if (fw) {
                char rec_str[9];
                snprintf(rec_str, 9, "%-8llu", (unsigned long long)actual_records);
                fseek(fw, 236, SEEK_SET);
                fwrite(rec_str, 1, 8, fw);
                fclose(fw);
                if (verbose) mexPrintf("Header repaired on disk.\n");
            }
        }

        hdr->num_data_records = (int)actual_records;
    }

    fseek(rd->fp, hdr->num_header_bytes, SEEK_SET);
}

/* ============================================================
 *                   EDF+ ANNOTATION PARSER
 * ============================================================ */

static void parse_tal_block_fast(
    const unsigned char *buf, mwSize len,
    Annotation **out, mwSize *count)
{
    mwSize i = 0;

    while (i < len) {
        while (i < len && buf[i] == 0) i++;
        if (i >= len) break;

        if (buf[i] != '+' && buf[i] != '-') { i++; continue; }

        char onset_str[64];
        mwSize j = 0;
        while (i < len && buf[i] != 20 && j < sizeof(onset_str)-1)
            onset_str[j++] = buf[i++];
        onset_str[j] = '\0';
        if (i >= len) break;
        i++; /* skip 0x14 */

        char **texts = NULL;
        mwSize ntext = 0;

        while (i < len && buf[i] != 0) {
            char txt[256];
            mwSize k = 0;
            while (i < len && buf[i] != 20 && buf[i] != 0 && k < sizeof(txt)-1)
                txt[k++] = buf[i++];
            txt[k] = '\0';

            if (k > 0) {
                texts = mxRealloc(texts, (ntext+1)*sizeof(char*));
                texts[ntext] = mxMalloc(k+1);
                memcpy(texts[ntext], txt, k+1);
                ntext++;
            }
            if (i < len && buf[i] == 20) i++;
        }

        if (ntext > 0) {
            *out = mxRealloc(*out, (*count + 1) * sizeof(Annotation));
            (*out)[*count].onset = atof(onset_str);

            mxArray *cell = mxCreateCellMatrix(1, ntext);
            for (mwSize t = 0; t < ntext; t++) {
                mxSetCell(cell, t, mxCreateString(texts[t]));
                mxFree(texts[t]);
            }
            mxFree(texts);

            (*out)[*count].texts = cell;
            (*count)++;
        } else {
            mxFree(texts);
        }

        i++;
    }
}

/* ============================================================
 *                   CHANNEL PLAN BUILDER
 * ============================================================ */

static int build_channel_plan(const mxArray *ch_cell,
                              const Signal_Header *sig, int nsig,
                              ChannelPlan *plan,
                              int *need_raw)
{
    mwSize nchan = mxGetNumberOfElements(ch_cell);
    int nplan = 0;

    for (mwSize c = 0; c < nchan; c++) {
        mxArray *elem = mxGetCell(ch_cell, c);
        if (!elem) continue;

        char req[64];
        mxGetString(elem, req, sizeof(req));
        {
            char *p = req;
            while (*p == ' ') p++;
            memmove(req, p, strlen(p)+1);
            size_t n = strlen(req);
            while (n && req[n-1]==' ') req[--n]='\0';
        }
        if (req[0] == '\0') continue;

        int idx = find_signal(sig, nsig, req);
        if (idx >= 0) {
            plan[nplan].kind       = CH_PLAIN;
            plan[nplan].raw_idx_a  = idx;
            plan[nplan].raw_idx_b  = -1;
            strncpy(plan[nplan].label, req, sizeof(plan[nplan].label)-1);
            plan[nplan].label[sizeof(plan[nplan].label)-1] = '\0';
            need_raw[idx] = 1;
            nplan++;
            continue;
        }

        int found_reref = 0;
        size_t req_len = strlen(req);
        for (size_t d = 1; d < req_len; d++) {
            if (req[d] != '-') continue;

            char chA[64], chB[64];
            size_t la = d;
            size_t lb = req_len - d - 1;
            if (la == 0 || lb == 0) continue;

            strncpy(chA, req,     la); chA[la] = '\0';
            strncpy(chB, req+d+1, lb); chB[lb] = '\0';

            {
                char *p = chA; while(*p==' ')p++; memmove(chA,p,strlen(p)+1);
                size_t n=strlen(chA); while(n&&chA[n-1]==' ')chA[--n]='\0';
                p = chB; while(*p==' ')p++; memmove(chB,p,strlen(p)+1);
                n=strlen(chB); while(n&&chB[n-1]==' ')chB[--n]='\0';
            }

            int ia = find_signal(sig, nsig, chA);
            int ib = find_signal(sig, nsig, chB);

            if (ia >= 0 && ib >= 0) {
                plan[nplan].kind       = CH_REREF;
                plan[nplan].raw_idx_a  = ia;
                plan[nplan].raw_idx_b  = ib;
                strncpy(plan[nplan].label, req, sizeof(plan[nplan].label)-1);
                plan[nplan].label[sizeof(plan[nplan].label)-1] = '\0';
                need_raw[ia] = 1;
                need_raw[ib] = 1;
                nplan++;
                found_reref = 1;
                break;
            }
        }

        if (!found_reref) {
            plan[nplan].kind      = CH_PLAIN;
            plan[nplan].raw_idx_a = -1;
            plan[nplan].raw_idx_b = -1;
            strncpy(plan[nplan].label, req, sizeof(plan[nplan].label)-1);
            plan[nplan].label[sizeof(plan[nplan].label)-1] = '\0';
            nplan++;
        }
    }

    return nplan;
}

/* ============================================================
 *                         MEX ENTRY
 * ============================================================ */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs < 1)
        mexErrMsgIdAndTxt("read_EDF_mex:Input", "Filename required");

    int verbose = (nrhs > 3 && mxGetScalar(prhs[3]) != 0);
    g_debug     = (nrhs > 5 && mxGetScalar(prhs[5]) != 0);

    char fname[4096];
    mxGetString(prhs[0], fname, sizeof(fname));

    int use_channel_plan = 0;
    if (nrhs > 1 && mxIsCell(prhs[1]) && mxGetNumberOfElements(prhs[1]) > 0)
        use_channel_plan = 1;

    EDFReader rd;
    if (!edfr_open(&rd, fname))
        mexErrMsgIdAndTxt("read_EDF_mex:File", "Cannot open file");

    DBG("\n===== INPUT =====\n");
    DBG("File: %s\n", fname);
    DBG("Compressed: %s\n",
        rd.comp == COMP_GZ   ? "yes (.gz)" :
        rd.comp == COMP_ZSTD ? "yes (.zst)" : "no");
    DBG("=================\n\n");

    /* ============================================================
     *                       READ MAIN HEADER
     * ============================================================ */

    EDF_Header hdr;
    if (!read_edf_header(&rd, &hdr)) {
        edfr_close(&rd);
        mexErrMsgIdAndTxt("read_EDF_mex:Header", "Invalid EDF header");
    }

    DBG("\n===== EDF MAIN HEADER =====\n");
    DBG("Header bytes: %d\n",           hdr.num_header_bytes);
    DBG("Num data records (header): %d\n", hdr.num_data_records);
    DBG("Record duration: %.6f\n",      hdr.data_record_duration);
    DBG("Num signals: %d\n",            hdr.num_signals);
    DBG("===========================\n\n");

    ASSERT(hdr.num_signals > 0);
    ASSERT(hdr.num_header_bytes >= 256);
    ASSERT(hdr.data_record_duration > 0);

    /* ============================================================
     *                     READ SIGNAL HEADERS
     * ============================================================ */

    Signal_Header *sig = mxCalloc(hdr.num_signals, sizeof(Signal_Header));
    read_signal_headers(&rd, sig, hdr.num_signals);

    mwSize total_samp_per_rec = 0;
    int annot_idx = -1;

    DBG("\n===== SIGNAL HEADERS =====\n");
    for (int i = 0; i < hdr.num_signals; i++) {
        sig[i].sampling_frequency =
            sig[i].samples_in_record / hdr.data_record_duration;
        total_samp_per_rec += sig[i].samples_in_record;

        DBG("Signal %d: %s\n",           i, sig[i].signal_labels);
        DBG("  samples/record: %d\n",    sig[i].samples_in_record);
        DBG("  digital min/max: %f / %f\n", sig[i].digital_min, sig[i].digital_max);

        if (strcmp(sig[i].signal_labels, "EDF Annotations") == 0)
            annot_idx = i;
    }
    DBG("Total samples per record: %llu\n", (unsigned long long)total_samp_per_rec);
    DBG("===========================\n\n");

    ASSERT(total_samp_per_rec > 0);

    /* ============================================================
     *           FIX EDF+ num_data_records (plain files only)
     * ============================================================ */

    int repair = (nrhs > 4 && mxGetScalar(prhs[4]) != 0);

    if (rd.comp != COMP_NONE) {
        if (repair && verbose)
            mexPrintf("RepairHeader requested but file is compressed (%s) — "
                      "skipping on-disk repair.\n",
                      rd.comp == COMP_GZ ? ".gz" : ".zst");
    } else {
        fix_num_records_plain(&rd, &hdr, sig, fname, verbose, repair);
    }

    mwSize bytes_per_rec = total_samp_per_rec * 2;

    /* ============================================================
     *                   BUILD CHANNEL PLAN
     * ============================================================ */

    int *need_raw = mxCalloc(hdr.num_signals, sizeof(int));

    ChannelPlan *plan = NULL;
    int nplan = 0;

    if (use_channel_plan) {
        plan = mxMalloc(mxGetNumberOfElements(prhs[1]) * sizeof(ChannelPlan));
        nplan = build_channel_plan(prhs[1], sig, hdr.num_signals, plan, need_raw);

        if (nplan == 0) {
            edfr_close(&rd);
            mexErrMsgIdAndTxt("read_EDF_mex:NoChannels",
                              "No valid channels found in the channel list.");
        }

        /* Warn about unmatched plain channels at plan-build time so the
         * caller hears about a missing label even on header-only calls
         * (nargout < 3). The data loop below silently emits an empty cell
         * for these; without this pre-warn, header-only callers got no
         * indication at all. */
        for (int p = 0; p < nplan; p++) {
            if (plan[p].kind == CH_PLAIN && plan[p].raw_idx_a < 0) {
                mexWarnMsgIdAndTxt("read_EDF_mex:UnknownChannel",
                    "Channel '%s' not found in file.", plan[p].label);
            }
        }

        DBG("\n===== CHANNEL PLAN =====\n");
        for (int p = 0; p < nplan; p++) {
            if (plan[p].kind == CH_PLAIN)
                DBG("  PLAIN  [%d]  %s\n", plan[p].raw_idx_a, plan[p].label);
            else
                DBG("  REREF  [%d]-[%d]  %s\n",
                    plan[p].raw_idx_a, plan[p].raw_idx_b, plan[p].label);
        }
        DBG("========================\n\n");
    } else {
        for (int i = 0; i < hdr.num_signals; i++)
            need_raw[i] = 1;
    }

    /* ============================================================
     *               STREAMING DATA + ANNOTATION READ
     * ============================================================
     *
     * Read records one at a time. For each record:
     *   • decode samples for any signal in need_raw[]
     *   • parse annotations (if requested and present)
     *
     * Output arrays are pre-sized from hdr.num_data_records (or 64 if the
     * header's count is missing/invalid) and grown by doubling. After the
     * loop, hdr.num_data_records is set to the actual count.
     */

    int want_data   = (nlhs >= 3);
    int want_annots = (nlhs >= 4) && (annot_idx >= 0);

    double **raw_signals = NULL;
    mwSize  *sig_lengths = NULL;          /* total samples per signal (= spr * records_read) */
    mwSize  *sig_offset  = NULL;
    Annotation *alist = NULL;
    mwSize acount = 0;
    mwSize records_read = 0;

    if (want_data || want_annots) {

        /* Per-signal byte/sample offsets within a record. */
        sig_offset = mxMalloc(hdr.num_signals * sizeof(mwSize));
        {
            mwSize acc = 0;
            for (int s = 0; s < hdr.num_signals; s++) {
                sig_offset[s] = acc;
                acc += sig[s].samples_in_record;
            }
        }

        /* Initial output capacity in records. */
        mwSize records_capacity =
            (hdr.num_data_records > 0) ? (mwSize)hdr.num_data_records : 64;

        if (want_data) {
            raw_signals = mxCalloc(hdr.num_signals, sizeof(double *));
            sig_lengths = mxCalloc(hdr.num_signals, sizeof(mwSize));

            for (int s = 0; s < hdr.num_signals; s++) {
                if (!need_raw[s]) continue;
                ASSERT(sig[s].digital_max != sig[s].digital_min);
                mwSize len = (mwSize)sig[s].samples_in_record * records_capacity;
                raw_signals[s] = mxMalloc(len * sizeof(double));
            }
        }

        /* One record's worth of raw bytes. */
        unsigned char *rec_buf = mxMalloc(bytes_per_rec);

        DBG("\n===== STREAMING DATA =====\n");
        DBG("Initial records_capacity: %llu\n",
            (unsigned long long)records_capacity);

        for (;;) {
            size_t got = edfr_read(&rd, rec_buf, bytes_per_rec);
            if (got == 0) break;
            if (got < bytes_per_rec) {
                /* Trailing partial record — discard, matches plain-file
                 * behavior which truncates by integer division of file size. */
                break;
            }

            /* Grow output arrays if we exceeded the planned capacity. */
            if (want_data && records_read >= records_capacity) {
                mwSize new_cap = records_capacity * 2;
                for (int s = 0; s < hdr.num_signals; s++) {
                    if (!need_raw[s]) continue;
                    mwSize new_len = (mwSize)sig[s].samples_in_record * new_cap;
                    raw_signals[s] = mxRealloc(raw_signals[s],
                                               new_len * sizeof(double));
                }
                records_capacity = new_cap;
            }

            /* Decode each requested signal in this record. */
            if (want_data) {
                for (int s = 0; s < hdr.num_signals; s++) {
                    if (!need_raw[s]) continue;
                    mwSize spr = (mwSize)sig[s].samples_in_record;
                    mwSize byte_off = sig_offset[s] * 2;

                    double scale = (sig[s].physical_max - sig[s].physical_min) /
                                   (sig[s].digital_max  - sig[s].digital_min);
                    double offs  =  sig[s].physical_min - sig[s].digital_min * scale;

                    double *dst = raw_signals[s] + records_read * spr;
                    const unsigned char *src = rec_buf + byte_off;

                    for (mwSize k = 0; k < spr; k++) {
                        int16_t vraw = (int16_t)(src[2*k] | (src[2*k+1] << 8));
                        dst[k] = vraw * scale + offs;
                    }

                    if (g_debug && s == 0 && records_read == 0) {
                        for (mwSize k = 0; k < 5 && k < spr; k++) {
                            int16_t vraw = (int16_t)(src[2*k] | (src[2*k+1] << 8));
                            DBG("Raw[%llu]=%d\n", (unsigned long long)k, vraw);
                        }
                    }
                }
            }

            /* Parse this record's annotation slice in place. */
            if (want_annots) {
                mwSize annot_byte_off = sig_offset[annot_idx] * 2;
                mwSize annot_bytes = (mwSize)sig[annot_idx].samples_in_record * 2;
                parse_tal_block_fast(rec_buf + annot_byte_off, annot_bytes,
                                     &alist, &acount);
            }

            records_read++;
        }

        DBG("Records read: %llu\n", (unsigned long long)records_read);
        DBG("===========================\n\n");

        mxFree(rec_buf);

        /* Reconcile the header with the actual record count. */
        if ((mwSize)hdr.num_data_records != records_read) {
            if (verbose)
                mexPrintf("num_data_records: header=%d, actual=%llu\n",
                          hdr.num_data_records,
                          (unsigned long long)records_read);
            hdr.num_data_records = (int)records_read;
        }

        /* Trim per-signal output arrays to actual length. */
        if (want_data) {
            for (int s = 0; s < hdr.num_signals; s++) {
                if (!need_raw[s]) continue;
                mwSize spr = (mwSize)sig[s].samples_in_record;
                mwSize len = spr * records_read;
                if (records_read != records_capacity) {
                    raw_signals[s] = mxRealloc(raw_signals[s],
                                               (len ? len : 1) * sizeof(double));
                }
                sig_lengths[s] = len;
            }
        }
    }

    /* ============================================================
     *                       HEADER OUTPUT
     * ============================================================ */

    if (nlhs >= 1) {
        const char *f[] = {"edf_ver","patient_id","local_rec_id",
                           "recording_startdate","recording_starttime",
                           "num_header_bytes","num_data_records",
                           "data_record_duration","num_signals"};

        plhs[0] = mxCreateStructMatrix(1, 1, 9, f);

        mxSetField(plhs[0], 0, "edf_ver",               mxCreateString(hdr.edf_ver));
        mxSetField(plhs[0], 0, "patient_id",             mxCreateString(hdr.patient_id));
        mxSetField(plhs[0], 0, "local_rec_id",           mxCreateString(hdr.local_rec_id));
        mxSetField(plhs[0], 0, "recording_startdate",    mxCreateString(hdr.recording_startdate));
        mxSetField(plhs[0], 0, "recording_starttime",    mxCreateString(hdr.recording_starttime));
        mxSetField(plhs[0], 0, "num_header_bytes",       mxCreateDoubleScalar(hdr.num_header_bytes));
        mxSetField(plhs[0], 0, "num_data_records",       mxCreateDoubleScalar((double)hdr.num_data_records));
        mxSetField(plhs[0], 0, "data_record_duration",   mxCreateDoubleScalar(hdr.data_record_duration));
        mxSetField(plhs[0], 0, "num_signals",            mxCreateDoubleScalar(hdr.num_signals));
    }

    /* ============================================================
     *         SIGNAL HEADER + DATA OUTPUT (channel-plan aware)
     * ============================================================ */

    int out_count = use_channel_plan ? nplan : hdr.num_signals;

    if (nlhs >= 2) {
        const char *f[] = {"signal_labels","transducer_type",
                           "physical_dimension","physical_min",
                           "physical_max","digital_min",
                           "digital_max","prefiltering",
                           "samples_in_record","sampling_frequency"};

        plhs[1] = mxCreateStructMatrix(1, out_count, 10, f);

        for (int o = 0; o < out_count; o++) {
            int s;
            const char *lbl;

            if (use_channel_plan) {
                s   = (plan[o].raw_idx_a >= 0) ? plan[o].raw_idx_a : 0;
                lbl = plan[o].label;
            } else {
                s   = o;
                lbl = sig[o].signal_labels;
            }

            mxSetField(plhs[1], o, "signal_labels",
                       mxCreateString(lbl));
            mxSetField(plhs[1], o, "transducer_type",
                       mxCreateString(sig[s].transducer_type));
            mxSetField(plhs[1], o, "physical_dimension",
                       mxCreateString(sig[s].physical_dimension));
            mxSetField(plhs[1], o, "physical_min",
                       mxCreateDoubleScalar(sig[s].physical_min));
            mxSetField(plhs[1], o, "physical_max",
                       mxCreateDoubleScalar(sig[s].physical_max));
            mxSetField(plhs[1], o, "digital_min",
                       mxCreateDoubleScalar(sig[s].digital_min));
            mxSetField(plhs[1], o, "digital_max",
                       mxCreateDoubleScalar(sig[s].digital_max));
            mxSetField(plhs[1], o, "prefiltering",
                       mxCreateString(sig[s].prefiltering));
            mxSetField(plhs[1], o, "samples_in_record",
                       mxCreateDoubleScalar(sig[s].samples_in_record));
            mxSetField(plhs[1], o, "sampling_frequency",
                       mxCreateDoubleScalar(sig[s].sampling_frequency));
        }
    }

    if (nlhs >= 3 && raw_signals != NULL) {
        ASSERT(hdr.num_data_records >= 0);

        plhs[2] = mxCreateCellMatrix(1, out_count);

        for (int o = 0; o < out_count; o++) {
            if (use_channel_plan) {
                int ia = plan[o].raw_idx_a;
                int ib = plan[o].raw_idx_b;

                if (ia < 0) {
                    /* Already warned at plan-build time (above). Emit an
                     * empty cell so the output array stays index-aligned
                     * with the caller's Channels list. */
                    mxSetCell(plhs[2], o, mxCreateDoubleMatrix(1, 0, mxREAL));
                    continue;
                }

                mwSize len = sig_lengths[ia];

                mxArray *v = mxCreateDoubleMatrix(1, len, mxREAL);
                double  *d = mxGetPr(v);

                if (plan[o].kind == CH_PLAIN) {
                    memcpy(d, raw_signals[ia], len * sizeof(double));
                } else {
                    mwSize len_b = sig_lengths[ib];
                    if (len != len_b) {
                        mxDestroyArray(v);
                        mexErrMsgIdAndTxt("read_EDF_mex:RerefLength",
                            "Rereferenced channels '%s' have different lengths "
                            "(%llu vs %llu).",
                            plan[o].label,
                            (unsigned long long)len,
                            (unsigned long long)len_b);
                    }
                    double *da = raw_signals[ia];
                    double *db = raw_signals[ib];
                    for (mwSize k = 0; k < len; k++)
                        d[k] = da[k] - db[k];
                }

                mxSetCell(plhs[2], o, v);

            } else {
                mwSize len = sig_lengths[o];
                mxArray *v = mxCreateDoubleMatrix(1, len, mxREAL);
                if (len > 0)
                    memcpy(mxGetPr(v), raw_signals[o], len * sizeof(double));
                mxSetCell(plhs[2], o, v);
            }
        }
    }

    if (raw_signals) {
        for (int s = 0; s < hdr.num_signals; s++)
            if (raw_signals[s]) mxFree(raw_signals[s]);
        mxFree(raw_signals);
    }
    if (sig_lengths) mxFree(sig_lengths);
    if (sig_offset)  mxFree(sig_offset);

    /* ============================================================
     *                     ANNOTATION OUTPUT
     * ============================================================ */

    if (nlhs >= 4) {
        const char *f[] = {"onset", "text"};
        mxArray *A;

        if (acount > 0) {
            A = mxCreateStructMatrix(acount, 1, 2, f);
            for (mwSize i = 0; i < acount; i++) {
                mxSetField(A, i, "onset", mxCreateDoubleScalar(alist[i].onset));
                mxSetField(A, i, "text",  alist[i].texts);
            }
        } else {
            A = mxCreateStructMatrix(0, 0, 2, f);
        }

        if (alist) mxFree(alist);
        plhs[3] = A;
    }

    /* ---- Cleanup ---- */
    if (plan)     mxFree(plan);
    if (need_raw) mxFree(need_raw);
    mxFree(sig);
    edfr_close(&rd);
}
