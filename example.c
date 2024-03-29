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

double emg_value_raw;
double emg_mean_absolute_value;
double emg_integrated;
double emg_ssi;
double emg_variance;
double emg_rms;
int16_t emg_myopulse_percent;

uint8_t ndx;


int16_t emg_array_raw[MA_FILTER_SIZE];
uint8_t emg_array_raw_index;
double emg_array_out[EMG_SIGNAL_SIZE];
uint8_t emg_array_out_index;
int16_t emg_envelope[ENVELOPE_SIZE];
uint8_t emg_envelope_index;


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
int main()
{
	size_t N = 1<<3;

	double complex vector[N];

	for(size_t n = 0; n < N; n++) {
		vector[n] = n;
	}

	
	double samplerate;
    printf("Enter the sample rate: ");
    scanf("%lf", &samplerate);


	fft(vector, N);

	printf("in frequency domain:\n");

	for(size_t n = 0; n < N; n++) {
		printf("%lf%+lfi\n", creal(vector[n]), cimag(vector[n]));
	}
	meanfreq = mean(vector, N, samplerate);
    medianfreq = median(vector, N, samplerate);
    peakfrequ = peakfreq(vector, N, samplerate);
	totalpower = powe(vector,N);

    	printf("Mean Frequency: %f\n", meanfreq);
    	printf("Median Frequency: %f\n", medianfreq);
    	printf("Peak Frequency: %f\n", peakfrequ);	
	printf("Total Power: %f\n", totalpower);
	
	
	return 0;
}
