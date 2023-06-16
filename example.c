#include <stddef.h>
#include <complex.h>
#include "fft.h"
#include "ift.h"
#include <stdio.h>
float mag(float complex c)
{
    return cabs(c);
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

int main()
{
	size_t N = 1<<3;

	float complex vector[N];

	for(size_t n = 0; n < N; n++) {
		vector[n] = n;
	}

	printf("in time domain:\n");

	for(size_t n = 0; n < N; n++) {
		printf("%f%+fi\n", creal(vector[n]), cimag(vector[n]));
	}
	float samplerate;
    printf("Enter the sample rate: ");
    scanf("%lf", &samplerate);


	fft(vector, N);

	printf("in frequency domain:\n");

	for(size_t n = 0; n < N; n++) {
		printf("%f%+fi\n", creal(vector[n]), cimag(vector[n]));
	}
	float meanfreq = mean(vector, N, samplerate);
    float medianfreq = median(vector, N, samplerate);

    printf("\nMean Frequency: \n", meanfreq);
    printf("Median Frequency: \n", medianfreq);

	ift(vector, N);

	printf("in time domain:\n");

	for(size_t n = 0; n < N; n++) {
		printf("%f%+fi\n", creal(vector[n]), cimag(vector[n]));
	}

	return 0;
}
