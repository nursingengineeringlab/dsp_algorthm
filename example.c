#include <stddef.h>
#include <complex.h>
#include "fft.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
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
    } else {
        return frequencies[middle];
    }
}
float peakfreq(float complex* x, size_t n, float samplerate)
{
	float maxmag = 0.0;
	float peakfr = 0.0;
	for (size_t i = 0; i < n; i++){
		float magn = mag(x[i]);
		printf("Magnitude: %f\n", magn);
		float frequency = i * samplerate / n;
		printf("Frequency: %f\n", frequency);
		
		if (magn > maxmag){
			maxmag = magn;
			peakfr = frequency;
			printf("Peak Frequency: %f\n", peakfr);
		}
	}
	return peakfr;
}
int main()
{
	size_t N = 1<<3;

	float complex vector[N];

	for(size_t n = 0; n < N; n++) {
		vector[n] = n;
	}

	
	float samplerate;
    printf("Enter the sample rate: ");
    scanf("%f", &samplerate);


	fft(vector, N);

	printf("in frequency domain:\n");

	for(size_t n = 0; n < N; n++) {
		printf("%f%+fi\n", creal(vector[n]), cimag(vector[n]));
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
