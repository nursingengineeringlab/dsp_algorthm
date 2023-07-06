#include <stddef.h>
#include <stdint.h>
#include <complex.h>
#include "fft.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define EMG_SIGNAL_SIZE	64
#define MA_FILTER_SIZE	8
#define ENVELOPE_SIZE 16


float emg_value_raw;
float emg_mean_absolute_value;
float emg_integrated;
float emg_ssi;
float emg_variance;
float emg_rms;
int16_t emg_myopulse_percent;


int16_t emg_array_raw[MA_FILTER_SIZE];
uint8_t emg_array_raw_index;
float emg_array_out[EMG_SIGNAL_SIZE];
uint8_t emg_array_out_index;
uint8_t emg_array_out_index_init;
int16_t emg_envelope[ENVELOPE_SIZE];
uint8_t emg_envelope_index;

uint8_t ndx;

float complex vector[EMG_SIGNAL_SIZE];
float samplerate;
float meanfreq;
float medianfreq;
float peakfrequ;
float totalpower;

static int ctz(size_t N)
{
	int ctz1 = 0;

	while( N ) {
		ctz1++;
		N >>= 1;
	}

	return ctz1-1;
}

static void nop_split(const float complex *x, float complex *X, size_t N)
{
	for(size_t n = 0; n < N/2; n++) {
		X[2*n+0] = x[0/2+n];
		X[2*n+1] = x[N/2+n];
	}
}

static void fft_split(const float complex *x, float complex *X, size_t N, float complex phi)
{
	for(size_t n = 0; n < N/2; n++) {
		X[2*n+0] = (x[0/2+n] + x[N/2+n]);
		X[2*n+1] = (x[0/2+n] - x[N/2+n]) * cexp(-2*(float)M_PI*I*phi);
	}
}

static size_t revbits(size_t v, int J)
{
	size_t r = 0;

	for(int j = 0; j < J; j++) {
		r |= ( (v>>j)&1 ) << (J-1-j);
	}

	return r;
}

static int nop_reverse(int b, float complex *buffers[2], size_t N)
{
	int J = ctz(N);

	for(int j = J-2; j >= 0; j--, b++) {
		size_t delta = N>>j;

		for(size_t n = 0; n < N; n += delta) {
			nop_split(buffers[b&1]+n, buffers[~b&1]+n, delta);
		}
	}

	return b;
}

static int fft_reverse(int b, float complex *buffers[2], size_t N)
{
	int J = ctz(N);

	for(int j = J-1; j >= 0; j--, b++) {
		size_t delta = N>>j;

		for(size_t n = 0; n < N; n += delta) {
			float complex phi = (float)revbits( n/delta, j) / (float)(2<<j);
			fft_split(buffers[b&1]+n, buffers[~b&1]+n, delta, phi);
		}
	}

	return b;
}

int fft(float complex *vector, size_t N)
{
	if( !N ) return 0;

	if( N & (N-1) ) return 1;

	float complex *buffers[2] = { vector, malloc(N*sizeof(float complex)) };

	if( !buffers[1] ) return -1;

	int b = 0;

	b = nop_reverse(b, buffers, N);
	b = fft_reverse(b, buffers, N);
	b = nop_reverse(b, buffers, N);

	memmove(vector, buffers[b&1], N*sizeof(float complex));

	free( buffers[1] );

	return 0;
}

float mag(float complex c)
{
    return cabs(c);
}
float powe(float complex* x, size_t n)
{
	float sum = 0.0;
	printf("Power of Frequencies:\n");
	for (int i = 0;i < n; i++){
		float power= x[i]*conj(x[i]);
		printf("%f\n",power);
		sum+=power;
	}
	return sum;
}
	
float mean(float complex* x, size_t n, float samplerate)
{
    float sum = 0.0;
    for (int i = 0; i < n; i++) {
        float magnitudevalue = mag(x[i]);
        float frequency = i * samplerate / n;
        sum += magnitudevalue * frequency;
    }
    return sum / n;
}
float median(float complex* x, size_t n, float samplerate)
{
    float frequencies[n];
    for (size_t i = 0; i < n; i++) {
        float magnitudevalue = mag(x[i]);
        float frequency = i * samplerate / n;
        frequencies[i] = magnitudevalue * frequency;
    }

    size_t middle = n / 2;
    if (n % 2 == 0) {
        float median = (frequencies[middle - 1] + frequencies[middle]) / 2.0;
        return median;
    }
    else {
        return frequencies[middle];
    }
}
float peakfreq(float complex* x, float n, float samplerate)
{
	float maxmag = 0.0;
	float peakfr = 0.0;
	for (size_t i = 0; i < n; i++){
		float magn = mag(x[i]);
		float frequency = i * samplerate / n;
		
		if (magn > maxmag){
			maxmag = magn;
			peakfr = frequency;
		}
	}
	return peakfr;
}

int main(){
	float MA_SUM = 0;
	
	for (;;)
	  {
	  
	  	// get raw signal input
	  	emg_value_raw = rand() % 97 * -1 ;
	  		  	
	  	// fill buffer with raw signals and increase index
	  	emg_array_raw[emg_array_raw_index++] = emg_value_raw;
	  	
	  	// reset array raw index
	  	if (emg_array_raw_index >= MA_FILTER_SIZE - 1)
	  		emg_array_raw_index = 0;

		// sum of filter size in array raw
		for(uint8_t i=0; i < MA_FILTER_SIZE; i++)
			MA_SUM += emg_array_raw[i];
		
		// fill out buffer with Moving average values
		
		emg_array_out[emg_array_out_index] = MA_SUM / MA_FILTER_SIZE;
		
		//Assigning complex vector for FFT
	  	
	  	vector[emg_array_out_index++] = MA_SUM / MA_FILTER_SIZE; 
		
		// a pointer to emg_array_out buffer
		float * c = emg_array_out;
		
		// print contents of emg_array_out ie MAV's
		while( *c != '\0'){
			printf("\n%f", *c);
			emg_integrated += (*c < 0) ? (*c * -1) : *c;
			emg_ssi += *c * *c;
			c++;
		}
		
		emg_mean_absolute_value = emg_integrated / EMG_SIGNAL_SIZE;
		emg_variance = emg_ssi / (EMG_SIGNAL_SIZE - 1);
		emg_rms = sqrt(emg_ssi / EMG_SIGNAL_SIZE);
		
		printf("\n");
		printf("\nemg_integrated: %f", emg_integrated);
		printf("\nemg_ssi: %f", emg_ssi);
		printf("\nemg_variance: %f", emg_variance);
		printf("\nemg_rms: %f", emg_rms);
		printf("\n\n");
		
		

		//printf("MA_SUM:%f, ", MA_SUM);
		//printf("MAV %f\n", MA_SUM / MA_FILTER_SIZE);
		
		MA_SUM = 0;
		emg_integrated = 0.0;
		emg_mean_absolute_value = 0.0;
		emg_ssi = 0.0;
		emg_rms = 0.0;
		emg_variance = 0.0;
		
		if (emg_array_out_index == EMG_SIGNAL_SIZE)
		{
			emg_array_out_index_init = emg_array_out_index;
			fft(vector, emg_array_out_index_init);
	
			printf("in frequency domain:\n");
	
			for(size_t x = 0; x < emg_array_out_index_init; x++) {
				printf("%f%+fi\n", creal(vector[x]), cimag(vector[x]));
			}
		
			meanfreq = mean(vector, emg_array_out_index_init, samplerate);
	    		medianfreq = median(vector, emg_array_out_index_init, samplerate);
	    		peakfrequ = peakfreq(vector, emg_array_out_index_init, samplerate);
			totalpower = powe(vector,emg_array_out_index_init);
			printf("EMG VALUE RAW: %lf\n", emg_value_raw);
	    		printf("Mean Frequency: %lf\n", meanfreq);
	    		printf("Median Frequency: %lf\n", medianfreq);
	    		printf("Peak Frequency: %lf\n", peakfrequ);	
			printf("Total Power: %lf\n", totalpower);
			printf("EMG Array Out Index %u\n",  emg_array_out_index);
		
		}
		
		// reset array out index
	  	if (emg_array_out_index >= EMG_SIGNAL_SIZE )
	  		emg_array_out_index = 0;
	}
	
	return 0;
	
}

