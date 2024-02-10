//
// Created by daran on 1/12/2017 to be used in ECE420 Sp17 for the first time.
// Modified by dwang49 on 1/1/2018 to adapt to Android 7.0 and Shield Tablet updates.
//

#include "ece420_main.h"
#include "ece420_lib.h"
#include "kiss_fft/kiss_fft.h"
#include <complex>
#include <vector>

// JNI Function
extern "C" {
JNIEXPORT float JNICALL
Java_com_ece420_lab4_MainActivity_getFreqUpdate(JNIEnv *env, jclass);
}

// Student Variables
#define F_S 48000
#define FRAME_SIZE 1024
#define VOICED_THRESHOLD 50000000  // Find your own threshold
float lastFreqDetected = -1;

void ece420ProcessFrame(sample_buf *dataBuf) {
    // Keep in mind, we only have 20ms to process each buffer!
    struct timeval start;
    struct timeval end;
    gettimeofday(&start, NULL);

    // Data is encoded in signed PCM-16, little-endian, mono
    float bufferIn[FRAME_SIZE];
    for (int i = 0; i < FRAME_SIZE; i++) {
        int16_t val = ((uint16_t) dataBuf->buf_[2 * i]) | (((uint16_t) dataBuf->buf_[2 * i + 1]) << 8);
        bufferIn[i] = (float) val;
    }

    // ********************** PITCH DETECTION ************************ //
    // In this section, you will be computing the autocorrelation of bufferIn
    // and picking the delay corresponding to the best match. Naively computing the
    // autocorrelation in the time domain is an O(N^2) operation and will not fit
    // in your timing window.
    //
    // First, you will have to detect whether or not a signal is voiced.
    // We will implement a simple voiced/unvoiced detector by thresholding
    // the power of the signal.
    //
    // Next, you will have to compute autocorrelation in its O(N logN) form.
    // Autocorrelation using the frequency domain is given as:
    //
    //  autoc = ifft(fft(x) * conj(fft(x)))
    //
    // where the fft multiplication is element-wise.
    //
    // You will then have to find the index corresponding to the maximum
    // of the autocorrelation. Consider that the signal is a maximum at idx = 0,
    // where there is zero delay and the signal matches perfectly.
    //
    // Finally, write the variable "lastFreqDetected" on completion. If voiced,
    // write your determined frequency. If unvoiced, write -1.
    // ********************* START YOUR CODE HERE *********************** //
    unsigned sum = 0;
    std::complex<float> mult;

    kiss_fft_cfg cfg = kiss_fft_alloc( FRAME_SIZE , false ,0,0 );
    kiss_fft_cfg cfg_i = kiss_fft_alloc( FRAME_SIZE , true ,0,0 );
    kiss_fft_cpx cx_in[FRAME_SIZE];
    kiss_fft_cpx cx_out[FRAME_SIZE];
    kiss_fft_cpx cx_out_conj[FRAME_SIZE];

    kiss_fft_cpx cx_in_i[FRAME_SIZE];
    kiss_fft_cpx cx_out_i[FRAME_SIZE];

    bool search = false;
    int maxVal = -1000;
    int maxIdx = -1;
    std::vector<float> peaks;
    float thresh = 0.5;
    float a[FRAME_SIZE];

    for(unsigned i = 0; i < FRAME_SIZE; i++){
        sum += pow(bufferIn[i],2);
        cx_in[i].r = bufferIn[i];
        cx_in[i].i = 0;

    }


    if(sum > VOICED_THRESHOLD){
        kiss_fft(cfg , cx_in, cx_out);
        kiss_fft(cfg , cx_in, cx_out_conj);
        for(unsigned i = 0; i < FRAME_SIZE; i++) {
            std::conj(cx_out_conj[i].i);
            mult = (cx_out[i].r + cx_out[i].i) * (cx_out_conj[i].r + cx_out_conj[i].i);
//            cx_in[i].r = cx_out[i].r * cx_out_conj[i].r;
//            cx_in[i].i = cx_out[i].i * cx_out_conj[i].i;
            cx_in_i[i].r = std::real(mult);
            cx_in_i[i].i = std::imag(mult);
        }

        kiss_fft(cfg_i , cx_in_i, cx_out_i);

        for(unsigned i = 0; i < sizeof(a); i++) {
            a[i] = (cx_out_i[i].r) / sum;
        }
        //sqrt(pow(cx_out[i].i, 2) + pow(cx_out[i].r, 2));
        //a[i] = (a[i] * b[i]) / sum;
        for(unsigned i = 0; i < sizeof(a); i++){
            if(a[i] >= thresh){
                search = true;
            }
            if(a[i] <= thresh && search == true){
                search = false;
                peaks.push_back(maxIdx);
                maxVal = -1000;
            }
            if(search == true){
                if(a[i]> maxVal){
                    maxIdx = i;
                    maxVal = a[i];
                }
            }
        }

        int L = findMaxArrayIdx(a,50,FRAME_SIZE/2);
//            for(unsigned i = 0; i < peaks.size(); i++){
//                if(peaks[i] < 170){
//                    continue;
//                }
//                else{
        lastFreqDetected = F_S/L;
        free(cfg);
        free(cfg_i);
//                    break;
//                }
//            }
    }
    else{
        lastFreqDetected = -1;
        free(cfg);
        free(cfg_i);
    }


    // ********************* END YOUR CODE HERE ************************* //
    gettimeofday(&end, NULL);
    LOGD("Time delay: %ld us",  ((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec)));
}

JNIEXPORT float JNICALL
Java_com_ece420_lab4_MainActivity_getFreqUpdate(JNIEnv *env, jclass) {
return lastFreqDetected;
}