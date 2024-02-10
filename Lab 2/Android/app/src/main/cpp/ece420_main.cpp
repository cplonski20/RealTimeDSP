//
// Created by daran on 1/12/2017 to be used in ECE420 Sp17 for the first time.
// Modified by dwang49 on 1/1/2018 to adapt to Android 7.0 and Shield Tablet updates.
//

#include "ece420_main.h"

// Student Variables
#define FRAME_SIZE 128

// FIR Filter Function Defined here located at the bottom
int16_t firFilter(int16_t sample);

void ece420ProcessFrame(sample_buf *dataBuf) {
    // Keep in mind, we only have a small amount of time to process each buffer!
    struct timeval start;
    gettimeofday(&start, NULL);

    // Using {} initializes all values in the array to zero
    int16_t bufferIn[FRAME_SIZE] = {};
    int16_t bufferOut[FRAME_SIZE] = {};

    // Your buffer conversion (unpacking) here
    // Fetch data sample from dataBuf->buf_[], unpack and put into bufferIn[]
    // ******************** START YOUR CODE HERE ******************** //
    for (int i=0; i<FRAME_SIZE; i++){
        bufferIn[i] = ((int16_t) ((dataBuf->buf_)[i*2+1])<<8)|(dataBuf->buf_[i*2]);
    }

    // ********************* END YOUR CODE HERE ********************* //
    // Loop code provided as a suggestion. This loop simulates sample-by-sample processing.
    for (int sampleIdx = 0; sampleIdx < FRAME_SIZE; sampleIdx++) {
        // Grab one sample from bufferIn[]
        int16_t sample = bufferIn[sampleIdx];
        // Call your filFilter funcion
        int16_t output = firFilter(sample);
        // Grab result and put into bufferOut[]
        bufferOut[sampleIdx] = output;
    }

    // Your buffer conversion (packing) here
    // Fetch data sample from bufferOut[], pack them and put back into dataBuf->buf_[]
    // ******************** START YOUR CODE HERE ******************** //

    for(int i=0; i<FRAME_SIZE; i++){
        (dataBuf->buf_)[2*i+1] = (uint8_t) (bufferOut[i] >>8);
        (dataBuf->buf_)[2*i] = (uint8_t) (bufferOut[i] & 0x00FF);
    }

    // ********************* END YOUR CODE HERE ********************* //

	// Log the processing time to Android Monitor or Logcat window at the bottom
    struct timeval end;
    gettimeofday(&end, NULL);
    LOGD("Loop timer: %ld us",  ((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec)));

}

// TODO: Change N_TAPS to match your filter design
#define N_TAPS 121
// TODO: Change myfilter to contain the coefficients of your designed filter.
double myfilter[N_TAPS] = {0.004053998819520238, 0.003413106430807618, 0.002469691024492706, 0.0012563675023683602, -0.00017122962902093442, -0.0017374156407246664, -0.003351252886173712, -0.004912883889730735, -0.006320964212783423, -0.007480643762514491, -0.00831146937132273, -0.008754553337058215, -0.008778376845193754, -0.008382674389659266, -0.007599971768654506, -0.006494518077016472, -0.005158549755240681, -0.0037060376635522363, -0.0022642799233249, -0.0009638968651225759, 7.205649750731406e-05, 0.0007390347898339504, 0.0009611273204637104, 0.0006995040272624597, -4.159278714064318e-05, -0.001211936366034262, -0.0027159025717696413, -0.004416726943234054, -0.006144328428580043, -0.007706316917130975, -0.008901537898370834, -0.009535280368041224, -0.00943510229175957, -0.008466126816130676, -0.006544642470635163, -0.0036489065529520426, 0.00017379884629305757, 0.004804582108314828, 0.010055263536722711, 0.015675705472397352, 0.021365900597841908, 0.026792219529816638, 0.03160690036717322, 0.035469598413421115, 0.038069624529441605, 0.039147401290951934, 0.038513666953476904, 0.036065060774267535, 0.031794924762343295, 0.025798444529720986, 0.018271607286175374, 0.009503854644748944, -0.0001352753052209069, -0.010214817010703531, -0.020265292699431132, -0.029804509731133272, -0.038364664220009476, -0.04551906664185049, -0.050906808537359745, -0.054253805477118465, 0.9446111269270208, -0.054253805477118465, -0.050906808537359745, -0.04551906664185049, -0.038364664220009476, -0.029804509731133272, -0.020265292699431132, -0.010214817010703531, -0.0001352753052209069, 0.009503854644748944, 0.018271607286175374, 0.025798444529720986, 0.031794924762343295, 0.036065060774267535, 0.038513666953476904, 0.039147401290951934, 0.038069624529441605, 0.035469598413421115, 0.03160690036717322, 0.026792219529816638, 0.021365900597841908, 0.015675705472397352, 0.010055263536722711, 0.004804582108314828, 0.00017379884629305757, -0.0036489065529520426, -0.006544642470635163, -0.008466126816130676, -0.00943510229175957, -0.009535280368041224, -0.008901537898370834, -0.007706316917130975, -0.006144328428580043, -0.004416726943234054, -0.0027159025717696413, -0.001211936366034262, -4.159278714064318e-05, 0.0006995040272624597, 0.0009611273204637104, 0.0007390347898339504, 7.205649750731406e-05, -0.0009638968651225759, -0.0022642799233249, -0.0037060376635522363, -0.005158549755240681, -0.006494518077016472, -0.007599971768654506, -0.008382674389659266, -0.008778376845193754, -0.008754553337058215, -0.00831146937132273, -0.007480643762514491, -0.006320964212783423, -0.004912883889730735, -0.003351252886173712, -0.0017374156407246664, -0.00017122962902093442, 0.0012563675023683602, 0.002469691024492706, 0.003413106430807618, 0.004053998819520238};


// Circular Buffer
int16_t circBuf[N_TAPS] = {};
int16_t circBufIdx = 0;

// FirFilter Function
int16_t firFilter(int16_t sample) {
    // This function simulates sample-by-sample processing. Here you will
    // implement an FIR filter such as:
    //
    // y[n] = a x[n] + b x[n-1] + c x[n-2] + ...
    //
    // You will maintain a circular buffer to store your prior samples
    // x[n-1], x[n-2], ..., x[n-k]. Suggested initializations circBuf
    // and circBufIdx are given.
    //
    // Input 'sample' is the current sample x[n].
    // ******************** START YOUR CODE HERE ******************** //
    int16_t output = 0;
    circBuf[circBufIdx] = sample;
    for(int n=0; n<N_TAPS; n++){
        int i;
        if(circBufIdx-n>=0){
            i = circBufIdx-n;
        }
        else {
            i =  circBufIdx-n+N_TAPS;
        }
        output += myfilter[n] * circBuf[i];
    }
    circBufIdx++;
    if(circBufIdx==N_TAPS){circBufIdx=0;}
    // ********************* END YOUR CODE HERE ********************* //
    return output;
}
