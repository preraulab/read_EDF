/**
 * READ_EDF_MEX  High-performance EDF / EDF+ file reader for MATLAB
 *
 * -------------------------------------------------------------------------
 * DESCRIPTION
 * -------------------------------------------------------------------------
 *   READ_EDF_MEX is a compiled MEX implementation of a fast and robust
 *   European Data Format (EDF / EDF+) reader. It performs:
 *
 *       • Full main header parsing (256-byte EDF header)
 *       • Per-signal header parsing
 *       • Automatic correction of invalid num_data_records
 *       • Digital → physical unit conversion
 *       • Multi-channel extraction
 *       • EDF+ annotation (TAL) parsing
 *       • Optional on-disk header repair
 *       • Extensive debug instrumentation (when enabled)
 *
 *   The reader is designed for:
 *
 *       • Large PSG / EEG datasets
 *       • Clinical EDF+ files with malformed record counts
 *       • Efficient memory usage
 *       • Deterministic, release-build stability
 *
 * -------------------------------------------------------------------------
 * MATLAB USAGE
 * -------------------------------------------------------------------------
 *
 *   [header, sigheader, data, annotations] =
 *       read_EDF_mex(filename, channels, epochs, verbose, repair, debug)
 *
 *   Required:
 *       filename  : string
 *
 *   Optional:
 *       channels  : (currently unused placeholder)
 *       epochs    : (currently unused placeholder)
 *       verbose   : logical (prints header repair info)
 *       repair    : logical (rewrite corrected record count to disk)
 *       debug     : logical (enables internal assertions + tracing)
 *
 * -------------------------------------------------------------------------
 * OUTPUTS
 * -------------------------------------------------------------------------
 *
 *   header      : 1x1 struct
 *       edf_ver
 *       patient_id
 *       local_rec_id
 *       recording_startdate
 *       recording_starttime
 *       num_header_bytes
 *       num_data_records
 *       data_record_duration
 *       num_signals
 *
 *   sigheader   : 1xN struct array (per signal)
 *       signal_labels
 *       transducer_type
 *       physical_dimension
 *       physical_min
 *       physical_max
 *       digital_min
 *       digital_max
 *       prefiltering
 *       samples_in_record
 *       sampling_frequency
 *
 *   data        : 1xN cell array
 *       Each cell contains a 1x(total_samples) double vector
 *       converted to physical units.
 *
 *   annotations : struct array (EDF+ only)
 *       onset   : double (seconds)
 *       text    : 1xM cell array of annotation strings
 *
 * -------------------------------------------------------------------------
 * DIGITAL → PHYSICAL CONVERSION
 * -------------------------------------------------------------------------
 *
 *   For each signal:
 *
 *       scale = (phys_max - phys_min) / (dig_max - dig_min)
 *       offset = phys_min - dig_min * scale
 *
 *       physical_value = digital_value * scale + offset
 *
 * -------------------------------------------------------------------------
 * EDF+ RECORD COUNT CORRECTION
 * -------------------------------------------------------------------------
 *
 *   Many EDF+ files contain:
 *       num_data_records = -1
 *       or incorrect values.
 *
 *   This implementation derives the true record count from:
 *
 *       file_size
 *       header.num_header_bytes
 *       samples_in_record (all channels)
 *
 *   If mismatch is detected:
 *       • header is corrected in memory
 *       • optionally repaired on disk (repair = true)
 *
 * -------------------------------------------------------------------------
 * DEBUG MODE
 * -------------------------------------------------------------------------
 *
 *   When debug = true:
 *       • Internal assertions are enabled
 *       • Detailed header + signal diagnostics printed
 *       • Raw sample previews printed (first channel only)
 *
 *   Intended for development, not production.
 *
 * -------------------------------------------------------------------------
 * COMPILATION
 * -------------------------------------------------------------------------
 *
 *   mex -O -largeArrayDims read_EDF_mex.c
 *
 *   -largeArrayDims is REQUIRED.
 *
 * -------------------------------------------------------------------------
 * NOTES
 * -------------------------------------------------------------------------
 *
 *   • Assumes little-endian EDF (standard compliant).
 *   • Digital samples are 16-bit signed integers.
 *   • Designed for large datasets (uses mwSize indexing).
 *   • Safe for multi-GB EDF files.
 *
 * -------------------------------------------------------------------------
 */

#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#if defined(MX_COMPAT_32)
#error "This MEX requires -largeArrayDims"
#endif

#define EDF_HEADER_SIZE 256

/* ============================================================
*                         DEBUG MACROS
* ============================================================ */
static int g_debug = 0;
#define DBG(...) do { if (g_debug) mexPrintf(__VA_ARGS__); } while (0)
#define ASSERT(cond) \
do { if (g_debug && !(cond)) mexErrMsgIdAndTxt("read_EDF_mex:Assert", #cond); } while (0)

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
    char signal_labels[17];
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
*                         UTILITIES
* ============================================================ */

static void trim_string(char *s) {
    char *p = s;
    while (*p == ' ') p++;
    memmove(s, p, strlen(p) + 1);
    size_t n = strlen(s);
    while (n && s[n-1] == ' ') s[--n] = 0;
}

/* ============================================================
*                         HEADER IO
* ============================================================ */

static int read_edf_header(FILE *fid, EDF_Header *h) {
    unsigned char buf[EDF_HEADER_SIZE];
    char tmp[16];

    if (fread(buf, 1, EDF_HEADER_SIZE, fid) != EDF_HEADER_SIZE) return 0;

    memcpy(h->edf_ver, buf, 8); h->edf_ver[8] = 0; trim_string(h->edf_ver);
    memcpy(h->patient_id, buf+8, 80); h->patient_id[80] = 0; trim_string(h->patient_id);
    memcpy(h->local_rec_id, buf+88, 80); h->local_rec_id[80] = 0; trim_string(h->local_rec_id);
    memcpy(h->recording_startdate, buf+168, 8); h->recording_startdate[8] = 0;
    memcpy(h->recording_starttime, buf+176, 8); h->recording_starttime[8] = 0;

    memcpy(tmp, buf+184, 8); tmp[8]=0; h->num_header_bytes = atoi(tmp);
    memcpy(tmp, buf+236, 8); tmp[8]=0; h->num_data_records = atoi(tmp);
    memcpy(tmp, buf+244, 8); tmp[8]=0; h->data_record_duration = atof(tmp);
    memcpy(tmp, buf+252, 4); tmp[4]=0; h->num_signals = atoi(tmp);

    return 1;
}

static int read_signal_headers(FILE *fid, Signal_Header *sh, int nsig) {
    char tmp[16];
    int i;

    for (i=0;i<nsig;i++) fread(sh[i].signal_labels,16,1,fid), sh[i].signal_labels[16]=0, trim_string(sh[i].signal_labels);
    for (i=0;i<nsig;i++) fread(sh[i].transducer_type,80,1,fid), sh[i].transducer_type[80]=0, trim_string(sh[i].transducer_type);
    for (i=0;i<nsig;i++) fread(sh[i].physical_dimension,8,1,fid), sh[i].physical_dimension[8]=0;
    for (i=0;i<nsig;i++) fread(tmp,8,1,fid), tmp[8]=0, sh[i].physical_min=atof(tmp);
    for (i=0;i<nsig;i++) fread(tmp,8,1,fid), tmp[8]=0, sh[i].physical_max=atof(tmp);
    for (i=0;i<nsig;i++) fread(tmp,8,1,fid), tmp[8]=0, sh[i].digital_min=atof(tmp);
    for (i=0;i<nsig;i++) fread(tmp,8,1,fid), tmp[8]=0, sh[i].digital_max=atof(tmp);
    for (i=0;i<nsig;i++) fread(sh[i].prefiltering,80,1,fid), sh[i].prefiltering[80]=0;
    for (i=0;i<nsig;i++) fread(tmp,8,1,fid), tmp[8]=0, sh[i].samples_in_record=atoi(tmp);
    for (i=0;i<nsig;i++) fread(tmp,32,1,fid);

    return 1;
}

/* ============================================================
*                 EDF+ num_data_records FIX
* ============================================================ */

static void fix_num_records(FILE *fid, EDF_Header *hdr,
                            Signal_Header *sig,
                            const char *fname,
                            int verbose,
                            int repair)
{
    mwSize total_samp_per_rec = 0;
    for (int i=0;i<hdr->num_signals;i++)
        total_samp_per_rec += sig[i].samples_in_record;

    DBG("\n===== DERIVED RECORD INFO =====\n");
    DBG("Total samples per record: %llu\n",
        (unsigned long long)total_samp_per_rec);

    for (int i=0;i<hdr->num_signals;i++) {
        DBG("Signal %d samples_in_record: %d\n",
            i, sig[i].samples_in_record);
    }
    DBG("===============================\n\n");

    ASSERT(total_samp_per_rec > 0);

    mwSize bytes_per_rec = total_samp_per_rec * 2;

    fseek(fid, 0, SEEK_END);
    mwSize fsize = ftell(fid);
    mwSize data_bytes = fsize - hdr->num_header_bytes;

    if (bytes_per_rec == 0) return;

    mwSize actual_records = data_bytes / bytes_per_rec;

    if (hdr->num_data_records <= 0 ||
        hdr->num_data_records != (int)actual_records)
    {
        if (verbose)
            mexPrintf("Updating num_data_records from %d to %llu\n",
                      hdr->num_data_records,
                      (unsigned long long)actual_records);

        if (repair) {
            FILE *fw = fopen(fname, "r+b");
            if (fw) {
                char rec_str[9];
                snprintf(rec_str, 9, "%-8llu",
                         (unsigned long long)actual_records);
                fseek(fw, 236, SEEK_SET);
                fwrite(rec_str, 1, 8, fw);
                fclose(fw);

                if (verbose)
                    mexPrintf("Header repaired on disk.\n");
            }
        }

        hdr->num_data_records = (int)actual_records;
    }

    fseek(fid, hdr->num_header_bytes, SEEK_SET);
}

/* ============================================================
*                   EDF+ ANNOTATION PARSER
* ============================================================ */

static void parse_tal_block_fast(
    const unsigned char *buf, mwSize len,
    Annotation **out, mwSize *count
    ) {
    mwSize i = 0;

    while (i < len) {
        while (i < len && buf[i] == 0) i++;
        if (i >= len) break;

        if (buf[i] != '+' && buf[i] != '-') { i++; continue; }

        char onset_str[64];
        mwSize j = 0;
        while (i < len && buf[i] != 20 && j < sizeof(onset_str)-1)
            onset_str[j++] = buf[i++];
        onset_str[j] = 0;
        if (i >= len) break;
        i++; /* skip 0x14 */

        char **texts = NULL;
        mwSize ntext = 0;

        while (i < len && buf[i] != 0) {
            char txt[256];
            mwSize k = 0;
            while (i < len && buf[i] != 20 && buf[i] != 0 && k < sizeof(txt)-1)
                txt[k++] = buf[i++];
            txt[k] = 0;

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
*                         MEX ENTRY
* ============================================================ */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs < 1)
        mexErrMsgIdAndTxt("read_EDF_mex:Input","Filename required");

    int verbose = (nrhs>3 && mxGetScalar(prhs[3])!=0);
    g_debug     = (nrhs>5 && mxGetScalar(prhs[5])!=0);

    char fname[4096];
    mxGetString(prhs[0], fname, sizeof(fname));

    FILE *fid = fopen(fname,"rb");
    if (!fid) mexErrMsgIdAndTxt("read_EDF_mex:File","Cannot open file");

    /* ============================================================
    *                         READ MAIN HEADER
    * ============================================================ */

    EDF_Header hdr;
    if (!read_edf_header(fid,&hdr))
        mexErrMsgIdAndTxt("read_EDF_mex:Header","Invalid EDF header");

    DBG("\n===== EDF MAIN HEADER =====\n");
    DBG("Header bytes: %d\n", hdr.num_header_bytes);
    DBG("Num data records (header): %d\n", hdr.num_data_records);
    DBG("Record duration: %.6f\n", hdr.data_record_duration);
    DBG("Num signals: %d\n", hdr.num_signals);
    DBG("===========================\n\n");

    ASSERT(hdr.num_signals > 0);
    ASSERT(hdr.num_header_bytes >= 256);
    ASSERT(hdr.data_record_duration > 0);

    /* ============================================================
    *                       READ SIGNAL HEADERS
    * ============================================================ */

    Signal_Header *sig = mxCalloc(hdr.num_signals,sizeof(Signal_Header));
    read_signal_headers(fid,sig,hdr.num_signals);

    mwSize total_samp_per_rec = 0;
    int annot_idx = -1;

    DBG("\n===== SIGNAL HEADERS =====\n");

    for (int i=0;i<hdr.num_signals;i++) {
        sig[i].sampling_frequency =
            sig[i].samples_in_record / hdr.data_record_duration;

        total_samp_per_rec += sig[i].samples_in_record;

        DBG("Signal %d: %s\n", i, sig[i].signal_labels);
        DBG("  samples/record: %d\n", sig[i].samples_in_record);
        DBG("  digital min/max: %f / %f\n",
            sig[i].digital_min, sig[i].digital_max);

        if (strcmp(sig[i].signal_labels,"EDF Annotations")==0)
            annot_idx = i;
    }

    DBG("Total samples per record: %llu\n",
        (unsigned long long)total_samp_per_rec);
    DBG("===========================\n\n");

    ASSERT(total_samp_per_rec > 0);

    /* ============================================================
    *                 FIX EDF+ num_data_records
    * ============================================================ */

    int repair = (nrhs>4 && mxGetScalar(prhs[4])!=0);
    fix_num_records(fid, &hdr, sig, fname, verbose, repair);

    mwSize bytes_per_rec = total_samp_per_rec * 2;

    /* ============================================================
    *                 READ RAW DATA (if requested)
    * ============================================================ */

    unsigned char *raw = NULL;
    mwSize data_bytes = 0;

    if (nlhs >= 3)
    {
        fseek(fid,0,SEEK_END);
        mwSize fsize = ftell(fid);
        data_bytes = fsize - hdr.num_header_bytes;

        DBG("\n===== FILE SIZE CHECK =====\n");
        DBG("File size: %llu\n", (unsigned long long)fsize);
        DBG("Header bytes: %d\n", hdr.num_header_bytes);
        DBG("Data bytes: %llu\n", (unsigned long long)data_bytes);
        DBG("Bytes per record: %llu\n",
            (unsigned long long)bytes_per_rec);
        DBG("Expected data bytes: %llu\n",
            (unsigned long long)
            ((mwSize)hdr.num_data_records * bytes_per_rec));
        DBG("===========================\n\n");

        ASSERT(bytes_per_rec > 0);
        ASSERT(data_bytes > 0);

        raw = mxMalloc(data_bytes);

        fseek(fid,hdr.num_header_bytes,SEEK_SET);
        size_t nread = fread(raw,1,data_bytes,fid);

        DBG("fread requested: %llu\n",
            (unsigned long long)data_bytes);
        DBG("fread returned:  %llu\n",
            (unsigned long long)nread);

        ASSERT(nread == data_bytes);
    }

    /* ============================================================
    *                         HEADER OUTPUT
    * ============================================================ */

    if (nlhs>=1) {
        const char *f[]={"edf_ver","patient_id","local_rec_id",
                         "recording_startdate","recording_starttime",
                         "num_header_bytes","num_data_records",
                         "data_record_duration","num_signals"};

        plhs[0]=mxCreateStructMatrix(1,1,9,f);

        mxSetField(plhs[0],0,"edf_ver",
                   mxCreateString(hdr.edf_ver));
        mxSetField(plhs[0],0,"patient_id",
                   mxCreateString(hdr.patient_id));
        mxSetField(plhs[0],0,"local_rec_id",
                   mxCreateString(hdr.local_rec_id));
        mxSetField(plhs[0],0,"recording_startdate",
                   mxCreateString(hdr.recording_startdate));
        mxSetField(plhs[0],0,"recording_starttime",
                   mxCreateString(hdr.recording_starttime));
        mxSetField(plhs[0],0,"num_header_bytes",
                   mxCreateDoubleScalar(hdr.num_header_bytes));
        mxSetField(plhs[0],0,"num_data_records",
                   mxCreateDoubleScalar((double)hdr.num_data_records));
        mxSetField(plhs[0],0,"data_record_duration",
                   mxCreateDoubleScalar(hdr.data_record_duration));
        mxSetField(plhs[0],0,"num_signals",
                   mxCreateDoubleScalar(hdr.num_signals));
    }

    /* ============================================================
    *                     SIGNAL HEADER OUTPUT
    * ============================================================ */

    if (nlhs>=2) {
        const char *f[]={"signal_labels","transducer_type",
                         "physical_dimension","physical_min",
                         "physical_max","digital_min",
                         "digital_max","prefiltering",
                         "samples_in_record","sampling_frequency"};

        plhs[1]=mxCreateStructMatrix(1,hdr.num_signals,10,f);

        for (int i=0;i<hdr.num_signals;i++) {
            mxSetField(plhs[1],i,"signal_labels",
                       mxCreateString(sig[i].signal_labels));
            mxSetField(plhs[1],i,"transducer_type",
                       mxCreateString(sig[i].transducer_type));
            mxSetField(plhs[1],i,"physical_dimension",
                       mxCreateString(sig[i].physical_dimension));
            mxSetField(plhs[1],i,"physical_min",
                       mxCreateDoubleScalar(sig[i].physical_min));
            mxSetField(plhs[1],i,"physical_max",
                       mxCreateDoubleScalar(sig[i].physical_max));
            mxSetField(plhs[1],i,"digital_min",
                       mxCreateDoubleScalar(sig[i].digital_min));
            mxSetField(plhs[1],i,"digital_max",
                       mxCreateDoubleScalar(sig[i].digital_max));
            mxSetField(plhs[1],i,"samples_in_record",
                       mxCreateDoubleScalar(sig[i].samples_in_record));
            mxSetField(plhs[1],i,"sampling_frequency",
                       mxCreateDoubleScalar(sig[i].sampling_frequency));
        }
    }

    /* ============================================================
    *                       SIGNAL DATA OUTPUT
    * ============================================================ */

    if (nlhs>=3 && raw!=NULL)
    {
        DBG("\n===== OUTPUT LOOP =====\n");
        DBG("num_data_records: %d\n", hdr.num_data_records);
        DBG("=======================\n\n");

        ASSERT(hdr.num_data_records > 0);

        plhs[2]=mxCreateCellMatrix(1,hdr.num_signals);

        mwSize offset=0;

        for (int s=0;s<hdr.num_signals;s++)
        {
            mwSize spr=sig[s].samples_in_record;
            mwSize len=spr*hdr.num_data_records;

            DBG("Signal %d output length: %llu\n",
                s,(unsigned long long)len);

            ASSERT(len > 0);
            ASSERT(sig[s].digital_max != sig[s].digital_min);

            mxArray *v=mxCreateDoubleMatrix(1,len,mxREAL);
            double *out=mxGetPr(v);

            double scale=(sig[s].physical_max-sig[s].physical_min) /
                (sig[s].digital_max-sig[s].digital_min);

            double offs=sig[s].physical_min -
                sig[s].digital_min*scale;

            mwSize out_i=0;

            for (mwSize r=0;r<hdr.num_data_records;r++)
            {
                mwSize base=r*total_samp_per_rec + offset;

                for (mwSize k=0;k<spr;k++)
                {
                    mwSize idx=2*(base+k);
                    int16_t vraw=
                        (int16_t)(raw[idx] |
                        (raw[idx+1]<<8));

                    if (g_debug && s==0 && r==0 && k<5)
                        DBG("Raw[%llu]=%d\n",
                            (unsigned long long)k,
                            vraw);

                    out[out_i++]=vraw*scale+offs;
                }
            }

            mxSetCell(plhs[2],s,v);
            offset+=spr;
        }

        mxFree(raw);
    }

    /* ============================================================
    *                     ANNOTATION OUTPUT
    * ============================================================ */

    if (nlhs>=4 && annot_idx>=0)
    {
        const char *f[]={"onset","text"};
        mxArray *A = mxCreateStructMatrix(0,0,2,f);

        Annotation *alist = NULL;
        mwSize acount = 0;

        for (mwSize r=0;r<hdr.num_data_records;r++)
        {
            mwSize byte_off = 0;
            for (int i=0;i<annot_idx;i++)
                byte_off += sig[i].samples_in_record*2;

            mwSize annot_bytes =
                sig[annot_idx].samples_in_record*2;

            unsigned char *blk = mxMalloc(annot_bytes);

            fseek(fid,
                  hdr.num_header_bytes +
                  r*bytes_per_rec +
                  byte_off,
                  SEEK_SET);

            fread(blk, 1, annot_bytes, fid);

            parse_tal_block_fast(blk, annot_bytes,
                                 &alist, &acount);

            mxFree(blk);
        }

        if (acount > 0) {
            mxDestroyArray(A);
            A = mxCreateStructMatrix(acount,1,2,f);
            for (mwSize i=0;i<acount;i++) {
                mxSetField(A,i,"onset",
                           mxCreateDoubleScalar(alist[i].onset));
                mxSetField(A,i,"text",
                           alist[i].texts);
            }
        }

        mxFree(alist);
        plhs[3] = A;
    }

    mxFree(sig);
    fclose(fid);
}