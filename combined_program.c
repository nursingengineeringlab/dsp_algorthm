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
#define FFT_SIZE 50


double emg_value_raw;
double emg_mean_absolute_value;
double emg_integrated;
double emg_ssi;
double emg_variance;
double emg_rms;
int16_t emg_myopulse_percent;


int16_t emg_array_raw[MA_FILTER_SIZE];
uint8_t emg_array_raw_index;
double emg_array_out[EMG_SIGNAL_SIZE];
uint8_t emg_array_out_index;
int16_t emg_envelope[ENVELOPE_SIZE];
uint8_t emg_envelope_index;

uint8_t ndx;

double complex vector[FFT_SIZE];
double samplerate;
double meanfreq;
double medianfreq;
double peakfrequ;
double totalpower;

static int ctz(size_t N)
{
	int ctz1 = 0;

	while( N ) {
		ctz1++;
		N >>= 1;
	}

	return ctz1-1;
}

static void nop_split(const double complex *x, double complex *X, size_t N)
{
	for(size_t n = 0; n < N/2; n++) {
		X[2*n+0] = x[0/2+n];
		X[2*n+1] = x[N/2+n];
	}
}

static void fft_split(const double complex *x, double complex *X, size_t N, double complex phi)
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

static int nop_reverse(int b, double complex *buffers[2], size_t N)
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

static int fft_reverse(int b, double complex *buffers[2], size_t N)
{
	int J = ctz(N);

	for(int j = J-1; j >= 0; j--, b++) {
		size_t delta = N>>j;

		for(size_t n = 0; n < N; n += delta) {
			double complex phi = (double)revbits( n/delta, j) / (double)(2<<j);
			fft_split(buffers[b&1]+n, buffers[~b&1]+n, delta, phi);
		}
	}

	return b;
}

int fft(double complex *vector, size_t N)
{
	if( !N ) return 0;

	if( N & (N-1) ) return 1;

	double complex *buffers[2] = { vector, malloc(N*sizeof(double complex)) };

	if( !buffers[1] ) return -1;

	int b = 0;

	b = nop_reverse(b, buffers, N);
	b = fft_reverse(b, buffers, N);
	b = nop_reverse(b, buffers, N);

	memmove(vector, buffers[b&1], N*sizeof(double complex));

	free( buffers[1] );

	return 0;
}

double mag(double complex c)
{
    return cabs(c);
}
double powe(double complex* x, size_t n)
{
	double sum = 0.0;
	printf("Power of Frequencies:\n");
	for (int i = 0;i < n; i++){
		double power= x[i]*conj(x[i]);
		printf("%f\n",power);
		sum+=power;
	}
	return sum;
}
	
double mean(double complex* x, size_t n, double samplerate)
{
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        double magnitudevalue = mag(x[i]);
        double frequency = i * samplerate / n;
        sum += magnitudevalue * frequency;
    }
    return sum / n;
}
double median(double complex* x, size_t n, double samplerate)
{
    double frequencies[n];
    for (size_t i = 0; i < n; i++) {
        double magnitudevalue = mag(x[i]);
        double frequency = i * samplerate / n;
        frequencies[i] = magnitudevalue * frequency;
    }

    size_t middle = n / 2;
    if (n % 2 == 0) {
        double median = (frequencies[middle - 1] + frequencies[middle]) / 2.0;
        return median;
    } else {
        return frequencies[middle];
    }
}
double peakfreq(double complex* x, double n, double samplerate)
{
	double maxmag = 0.0;
	double peakfr = 0.0;
	for (size_t i = 0; i < n; i++){
		double magn = mag(x[i]);
		double frequency = i * samplerate / n;
		
		if (magn > maxmag){
			maxmag = magn;
			peakfr = frequency;
		}
	}
	return peakfr;
}

int main(){
	double MA_SUM = 0;
	samplerate=64;



	for (size_t n = 0; n < FFT_SIZE; n++)
	  {
	  
	  
	  	// get raw signal input
	  	emg_value_raw = rand() % 97 * -1 ;
	  	
	  	//Assigning complex vector for FFT
	  	
	  	vector[emg_array_raw_index++] = emg_value_raw; //works for now but will need separated variable for index
	  	
	  	// fill buffer with raw signals and increase index
	  	emg_array_raw[emg_array_raw_index++] = emg_value_raw;
	  	
	  	// reset array raw index
	  	if (emg_array_raw_index >= MA_FILTER_SIZE - 1)
	  		emg_array_raw_index = 0;

		// sum of filter size in array raw
		for(uint8_t i=0; i < MA_FILTER_SIZE; i++)
			MA_SUM += emg_array_raw[i];
		
		// fill out buffer with Moving average values
		
		emg_array_out[emg_array_out_index++] = MA_SUM / MA_FILTER_SIZE;

		// reset array out index
	  	if (emg_array_out_index >= EMG_SIGNAL_SIZE - 1)
	  		emg_array_out_index = 0;

		// a pointer to emg_array_out buffer
		double * c = emg_array_out;
		
		// print contents of emg_array_out ie MAV's
		while( *c != '\0'){
			printf("out: %f ", *c);
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
		
		}
		
		size_t n=50;

		fft(vector, n);

		printf("in frequency domain:\n");

		for(size_t x = 0; x < n; x++) {
			printf("%lf%+lfi\n", creal(vector[x]), cimag(vector[x]));
		}
		meanfreq = mean(vector, n, samplerate);
    		medianfreq = median(vector, n, samplerate);
    		peakfrequ = peakfreq(vector, n, samplerate);
		totalpower = powe(vector,n);

    		printf("Mean Frequency: %lf\n", meanfreq);
    		printf("Median Frequency: %lf\n", medianfreq);
    		printf("Peak Frequency: %lf\n", peakfrequ);	
		printf("Total Power: %lf\n", totalpower);
	
	
		return 0;



	  

	

}
