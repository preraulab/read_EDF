/**
* READ_EDF_MEX - Fast, robust EDF/EDF+ reader (release build)
*
* MATLAB usage:
*   [header, sigheader, data, annotations] = read_EDF_mex(filename, channels, epochs, verbose, repair, debug)
*
* Compile:
*   mex -largeArrayDims read_EDF_mex.c
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
#define DBG(...) do { if (debug) mexPrintf(__VA_ARGS__); } while (0)
#define ASSERT(cond) \
do { if (debug && !(cond)) mexErrMsgIdAndTxt("read_EDF_mex:Assert", #cond); } while (0)

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
    for (i=0;i<nsig;i++) fread(tmp,32,1,fid); /* reserved */

    return 1;
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

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    if (nrhs < 1)
        mexErrMsgIdAndTxt("read_EDF_mex:Input","Filename required");

    int verbose = (nrhs>3 && mxGetScalar(prhs[3])!=0);
    int debug   = (nrhs>5 && mxGetScalar(prhs[5])!=0);

    char fname[4096];
    mxGetString(prhs[0], fname, sizeof(fname));

    FILE *fid = fopen(fname,"rb");
    if (!fid) mexErrMsgIdAndTxt("read_EDF_mex:File","Cannot open file");

    EDF_Header hdr;
    if (!read_edf_header(fid,&hdr))
        mexErrMsgIdAndTxt("read_EDF_mex:Header","Invalid EDF header");

    Signal_Header *sig = mxCalloc(hdr.num_signals,sizeof(Signal_Header));
    read_signal_headers(fid,sig,hdr.num_signals);

    mwSize total_samp_per_rec = 0;
    int annot_idx = -1;
    for (int i=0;i<hdr.num_signals;i++) {
        sig[i].sampling_frequency = sig[i].samples_in_record / hdr.data_record_duration;
        total_samp_per_rec += sig[i].samples_in_record;
        if (strcmp(sig[i].signal_labels,"EDF Annotations")==0)
            annot_idx = i;
    }

    /* ===================== Signal Data (conditional) ===================== */
    unsigned char *raw = NULL;
    mwSize nrec = 0, bytes_per_rec = total_samp_per_rec*2;
    if (nlhs >= 3) {  // only read raw data if user requested it
        fseek(fid,0,SEEK_END);
        mwSize fsize = ftell(fid);
        mwSize data_bytes = fsize - hdr.num_header_bytes;
        nrec = data_bytes / bytes_per_rec;

        raw = mxMalloc(data_bytes);
        fseek(fid,hdr.num_header_bytes,SEEK_SET);
        fread(raw,1,data_bytes,fid);
    }

    /* ===================== Header Output ===================== */
    if (nlhs>=1) {
        const char *f[]={"edf_ver","patient_id","local_rec_id","recording_startdate",
                         "recording_starttime","num_header_bytes","num_data_records",
                         "data_record_duration","num_signals"};
        plhs[0]=mxCreateStructMatrix(1,1,9,f);
        mxSetField(plhs[0],0,"edf_ver",mxCreateString(hdr.edf_ver));
        mxSetField(plhs[0],0,"patient_id",mxCreateString(hdr.patient_id));
        mxSetField(plhs[0],0,"local_rec_id",mxCreateString(hdr.local_rec_id));
        mxSetField(plhs[0],0,"recording_startdate",mxCreateString(hdr.recording_startdate));
        mxSetField(plhs[0],0,"recording_starttime",mxCreateString(hdr.recording_starttime));
        mxSetField(plhs[0],0,"num_header_bytes",mxCreateDoubleScalar(hdr.num_header_bytes));
        mxSetField(plhs[0],0,"num_data_records",mxCreateDoubleScalar((double)hdr.num_data_records));
        mxSetField(plhs[0],0,"data_record_duration",mxCreateDoubleScalar(hdr.data_record_duration));
        mxSetField(plhs[0],0,"num_signals",mxCreateDoubleScalar(hdr.num_signals));
    }

    /* ===================== Signal Header Output ===================== */
    if (nlhs>=2) {
        const char *f[]={"signal_labels","transducer_type","physical_dimension",
                         "physical_min","physical_max","digital_min","digital_max",
                         "prefiltering","samples_in_record","sampling_frequency"};
        plhs[1]=mxCreateStructMatrix(1,hdr.num_signals,10,f);
        for (int i=0;i<hdr.num_signals;i++) {
            mxSetField(plhs[1],i,"signal_labels",mxCreateString(sig[i].signal_labels));
            mxSetField(plhs[1],i,"transducer_type",mxCreateString(sig[i].transducer_type));
            mxSetField(plhs[1],i,"physical_dimension",mxCreateString(sig[i].physical_dimension));
            mxSetField(plhs[1],i,"physical_min",mxCreateDoubleScalar(sig[i].physical_min));
            mxSetField(plhs[1],i,"physical_max",mxCreateDoubleScalar(sig[i].physical_max));
            mxSetField(plhs[1],i,"digital_min",mxCreateDoubleScalar(sig[i].digital_min));
            mxSetField(plhs[1],i,"digital_max",mxCreateDoubleScalar(sig[i].digital_max));
            mxSetField(plhs[1],i,"samples_in_record",mxCreateDoubleScalar(sig[i].samples_in_record));
            mxSetField(plhs[1],i,"sampling_frequency",mxCreateDoubleScalar(sig[i].sampling_frequency));
        }
    }

    /* ===================== Signal Data Output ===================== */
    if (nlhs>=3 && raw!=NULL) {
        plhs[2]=mxCreateCellMatrix(1,hdr.num_signals);
        mwSize offset=0;
        for (int s=0;s<hdr.num_signals;s++) {
            mwSize spr=sig[s].samples_in_record;
            mwSize len=spr*hdr.num_data_records;
            mxArray *v=mxCreateDoubleMatrix(1,len,mxREAL);
            double *out=mxGetPr(v);

            double scale=(sig[s].physical_max-sig[s].physical_min) /
                (sig[s].digital_max-sig[s].digital_min);
            double offs=sig[s].physical_min - sig[s].digital_min*scale;

            mwSize out_i=0;
            for (mwSize r=0;r<hdr.num_data_records;r++) {
                mwSize base=r*total_samp_per_rec + offset;
                for (mwSize k=0;k<spr;k++) {
                    mwSize idx=2*(base+k);
                    int16_t vraw=(int16_t)(raw[idx] | (raw[idx+1]<<8));
                    out[out_i++]=vraw*scale+offs;
                }
            }
            mxSetCell(plhs[2],s,v);
            offset+=spr;
        }
        mxFree(raw);
    }

    /* ===================== Annotation Output ===================== */
    if (nlhs>=4 && annot_idx>=0) {
        const char *f[]={"onset","text"};
        mxArray *A = mxCreateStructMatrix(0,0,2,f);

        Annotation *alist = NULL;
        mwSize acount = 0;
        for (mwSize r=0;r<hdr.num_data_records;r++) {
            mwSize byte_off = 0;
            for (int i=0;i<annot_idx;i++)
                byte_off += sig[i].samples_in_record*2;

            mwSize annot_bytes = sig[annot_idx].samples_in_record*2;
            unsigned char *blk = mxMalloc(annot_bytes);
            fseek(fid, hdr.num_header_bytes + r*bytes_per_rec + byte_off, SEEK_SET);
            fread(blk, 1, annot_bytes, fid);
            parse_tal_block_fast(blk, annot_bytes, &alist, &acount);
            mxFree(blk);
        }

        if (acount > 0) {
            mxDestroyArray(A);
            A = mxCreateStructMatrix(acount,1,2,f);
            for (mwSize i=0;i<acount;i++) {
                mxSetField(A,i,"onset",mxCreateDoubleScalar(alist[i].onset));
                mxSetField(A,i,"text",alist[i].texts);
            }
        }
        mxFree(alist);
        plhs[3] = A;
    }

    mxFree(sig);
    fclose(fid);
}
