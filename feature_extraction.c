#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define MA_filter_size 50
#define emg_signal_size 64

int32_t emg_value_raw;                           
int32_t emg_mean_absolute_value;                                                                                                
int32_t emg_integrated;                                                              
int32_t emg_variance;
int64_t emg_ssi;
int64_t emg_ssi_squared;
int64_t emg_root_mean_value;
int64_t emg_root_mean_square;

uint8_t emg_array_raw_index;
uint8_t emg_array_out_index;

int32_t emg_array_raw[MA_filter_size];
int32_t emg_array_out[emg_signal_size];
int8_t test_array[200];


int main(void){
    emg_array_raw_index = 0;
    emg_array_out_index = 0;
    double MA_sum = 0;

    uint8_t i;
    uint8_t j;
    uint8_t k;
    uint8_t s;
    uint8_t t=0;
    uint8_t z;

    for(z=0; z < 200; z++){
        test_array[z] = -z;
    }

    for(s=0; s<200; s++){
        emg_value_raw = test_array[t];
        t++;
        emg_array_raw[emg_array_raw_index] = emg_value_raw;

        printf("EMG: %d\n", emg_value_raw);

        if(emg_array_raw_index < MA_filter_size - 1){
            emg_array_raw_index++;
        }
        else{
            emg_array_raw_index = 0;
        }

        for(i = 0; i < MA_filter_size; i++){
            MA_sum += emg_array_raw[i];
        }
        emg_array_out[emg_array_out_index] = MA_sum / MA_filter_size;

        printf("EMG Out: %d\n", emg_array_out[emg_array_out_index]);
        if(emg_array_out_index < emg_signal_size - 1) {
            emg_array_out_index++;
        }
        else{
            emg_array_out_index = 0;
        }
        for(j = 0; j < emg_signal_size; j++){
            if(emg_array_out[j] < 0){
                emg_array_out[j] = emg_array_out[j] * (-1);
            }
            emg_integrated += emg_array_out[j];
        }

        emg_mean_absolute_value = emg_integrated / emg_signal_size;

        for(k = 0; k < emg_signal_size; k++){
            emg_ssi_squared = emg_array_out[k] * emg_array_out[k];
            emg_ssi += emg_ssi_squared;     
        }

        emg_variance = emg_ssi / (emg_signal_size-1);
        emg_root_mean_value = emg_ssi / emg_signal_size;
        emg_root_mean_square = sqrt(emg_root_mean_value);

        printf("iEMG: %d\n", emg_integrated);
        printf("MAV: %d\n", emg_mean_absolute_value);
        printf("SSI: %d\n", emg_ssi);
        printf("EMG Variance: %d\n", emg_variance);
        printf("RMS: %d\n", emg_root_mean_square);


        emg_mean_absolute_value = 0;                                                                                                  
        emg_integrated = 0;                                                              
        emg_variance = 0;
        emg_ssi = 0;
        emg_ssi_squared = 0;
        emg_root_mean_value = 0;
        emg_root_mean_square = 0;

    }
}