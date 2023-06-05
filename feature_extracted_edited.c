#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
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


int16_t emg_array_raw[MA_FILTER_SIZE];
uint8_t emg_array_raw_index;
double emg_array_out[EMG_SIGNAL_SIZE];
uint8_t emg_array_out_index;
int16_t emg_envelope[ENVELOPE_SIZE];
uint8_t emg_envelope_index;

uint8_t ndx;




int main(){
	double MA_SUM = 0;



	for (ndx=0; ndx < 50; ndx++)
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

	

}

