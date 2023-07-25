/**
 * @file       main.c
 * @version    1.0.0
 * @date       2022-09-25
 * @author     Abu Bony Amin, Ebenezer Asabre, Akshat Sahay 
 * @brief      EMG Signal Collection and Analysis 
 * @note       None
 * @example    None
 */






/* Includes ----------------------------------------------------------- */
#include <stdint.h>
#include <string.h>
#include <stdarg.h>
#include <stddef.h>
#include <math.h>

// for frequency
#include <stddef.h>
#include <stdint.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>



#include "nordic_common.h"
#include "app_error.h"
#include "boards.h"
#include "nrf.h"
#include "ble_hci.h"
#include "ble_advdata.h"
#include "ble_advertising.h"
#include "ble_conn_params.h"
#include "nrf_sdh.h"
#include "nrf_sdh_soc.h"
#include "nrf_sdh_ble.h"
#include "nrf_ble_gatt.h"
#include "nrf_ble_qwr.h"
#include "app_timer.h"
#include "app_util_platform.h"
#include "bsp_btn_ble.h"
#include "nrf_pwr_mgmt.h"
#include "nrf_delay.h"
#include "nrf_drv_gpiote.h"

#include "sys_bm.h"
#include "ble_bas.h"
#include "ble_dis.h"
#include "ble_acs.h"
#include "ble_mgs.h"
#include "ble_gys.h"
#include "bsp_hw.h"
#include "bsp_imu.h"
#include "bsp_afe.h"
#include "bsp_nand_flash.h"
#include "nrf52832_peripherals.h"
#include "uart.h"
#include "bno085.h"
#include "fft.h"

#if defined(UART_PRESENT)
#include "nrf_uart.h"
#endif
#if defined(UARTE_PRESENT)
#include "nrf_uarte.h"
#endif

/* Private defines ---------------------------------------------------- */
#define APP_BLE_CONN_CFG_TAG            1                                          /**< A tag identifying the SoftDevice BLE configuration. */

#define SENSORS_MEAS_INTERVAL           APP_TIMER_TICKS(2000)                      /**< Sensors measurement interval (ticks). */
#define BATT_LEVEL_MEAS_INTERVAL        APP_TIMER_TICKS(20000)                     /**< Battery level measurement interval (ticks). */

#define DEVICE_NAME                     "imu-lcd"                                  /**< Name of device. Will be included in the advertising data. */

#define MANUFACTURER_NAME               "miBEAT"                                   /**< Manufacturer. Will be passed to Device Information Service. */

#define ACS_SERVICE_UUID_TYPE           BLE_UUID_TYPE_VENDOR_BEGIN                  /**< UUID type for the Nordic UART Service (vendor specific). */

#define APP_BLE_OBSERVER_PRIO           3                                           /**< Application's BLE observer priority. You shouldn't need to modify this value. */

#define APP_ADV_INTERVAL                64                                          /**< The advertising interval (in units of 0.625 ms. This value corresponds to 40 ms). */

#define APP_ADV_DURATION                18000                                       /**< The advertising duration (180 seconds) in units of 10 milliseconds. */

#define MIN_CONN_INTERVAL               MSEC_TO_UNITS(20, UNIT_1_25_MS)             /**< Minimum acceptable connection interval (20 ms), Connection interval uses 1.25 ms units. */
#define MAX_CONN_INTERVAL               MSEC_TO_UNITS(75, UNIT_1_25_MS)             /**< Maximum acceptable connection interval (75 ms), Connection interval uses 1.25 ms units. */
#define SLAVE_LATENCY                   0                                           /**< Slave latency. */
#define CONN_SUP_TIMEOUT                MSEC_TO_UNITS(4000, UNIT_10_MS)             /**< Connection supervisory timeout (4 seconds), Supervision Timeout uses 10 ms units. */
#define FIRST_CONN_PARAMS_UPDATE_DELAY  APP_TIMER_TICKS(5000)                       /**< Time from initiating event (connect or start of notification) to first time sd_ble_gap_conn_param_update is called (5 seconds). */
#define NEXT_CONN_PARAMS_UPDATE_DELAY   APP_TIMER_TICKS(30000)                      /**< Time between each call to sd_ble_gap_conn_param_update after the first call (30 seconds). */
#define MAX_CONN_PARAMS_UPDATE_COUNT    3                                           /**< Number of attempts before giving up the connection parameter negotiation. */

// custom
#define LED_STATE_INTERVAL              APP_TIMER_TICKS(2000)                       /** led interval **/



#define DEAD_BEEF                       0xDEADBEEF                                  /**< Value used as error code on stack dump, can be used to identify stack location on stack unwind. */
#define DEVICE_ID                       0x55AA

#define SCALE_FACTOR_ACC                1000 
#define SCALE_FACTOR_GYR                10000
#define SCALE_FACTOR_MAG                1000 

#define TX_BYTE(x) do {} while (app_uart_put((x)) != NRF_SUCCESS)                   /**< Macro for UART TX. */ 


/* Private defines for time */
#define   EMG_SIGNAL_SIZE     64
#define   MA_FILTER_SIZE      8
#define   ENVELOPE_SIZE       16





/* Private macros ----------------------------------------------------- */                                                            /**< BLE HRNS service instance. */
BLE_ACS_DEF(m_acs);                                                                 /**< BLE ACS service instance. */
BLE_MGS_DEF(m_mgs);                                                                 /**< BLE MGS service instance. */
BLE_GYS_DEF(m_gys);                                                                 /**< BLE GYS service instance. */
BLE_BAS_DEF(m_bas);                                                                 /**< Structure used to identify the battery service. */
NRF_BLE_GATT_DEF(m_gatt);                                                           /**< GATT module instance. */
NRF_BLE_QWR_DEF(m_qwr);                                                             /**< Context for the Queued Write module.*/
BLE_ADVERTISING_DEF(m_advertising);                                                 /**< Advertising module instance. */
APP_TIMER_DEF(m_sensors_timer_id);                                                  /**< Sensor measurement timer. */
APP_TIMER_DEF(m_battery_timer_id);                                                  /**< Battery timer. */
APP_TIMER_DEF(m_led_timer_id);                                                      /**< LED timer **/


/* Private variables -------------------------------------------------- */
static uint16_t   m_conn_handle          = BLE_CONN_HANDLE_INVALID;                 /**< Handle of the current connection. */
static ble_uuid_t m_adv_uuids[]          =                                          /**< Universally unique service identifier. */
{
  {BLE_UUID_BATTERY_SERVICE,            BLE_UUID_TYPE_BLE},
  {BLE_UUID_DEVICE_INFORMATION_SERVICE, BLE_UUID_TYPE_BLE}
};

uint32_t app_time;                                                                  /**< Elapsed time in the app. */                                                                                                                                           
int16_t emg_value;                                                                  /**< Current raw ECG/EMG sample from AFE. */  
int16_t imu[9];                                                                     /**< IMU data. */
uint8_t imu_size = 0;                                                               /**< Count of valid IMU data. */

uint8_t msec = 0;                                                                   /**< Millisecond counter synchronyzed to AFE interrupts. */

sh2_SensorValue_t sensor_value;

/* Private function prototypes ---------------------------------------- */
static void nrf_qwr_error_handler(uint32_t nrf_error);
static void sleep_mode_enter(void);
static void log_init(void);
static void power_management_init(void);
static void idle_state_handle(void);

static void battery_level_meas_timeout_handler(void * p_context);
static void sensors_meas_timeout_handler(void * p_context);

static void battery_level_update(void);
static void u8p_send(uint8_t field, uint16_t value);
void setReports(void);  


// custom prototypes
static void timers_init(void);
static void application_timers_start(void);


// frequency domain prototypes
float peakfreq(float complex* x, float n, float samplerate);
float median(float complex* x, size_t n, float samplerate);
float mean(float complex* x, size_t n, float samplerate);
float powe(float complex* x, size_t n);
float mag(float complex c);
static int fft_reverse(int b, float complex *buffers[2], size_t N);
static int nop_reverse(int b, float complex *buffers[2], size_t N);
static size_t revbits(size_t v, int J);
static void fft_split(const float complex *x, float complex *X, size_t N, float complex phi);
static void nop_split(const float complex *x, float complex *X, size_t N);
static int ctz(size_t N);



// get LED
#define LED1  17
#define PIN_24  24

/** values for time domain features **/
float emg_value_raw;
float emg_mean_absolute_value;
float emg_integrated;
float emg_ssi;
float emg_variance;
float emg_rms;
int16_t emg_myopulse_percent;


float emg_array_raw[MA_FILTER_SIZE];
uint8_t emg_array_raw_index;
float emg_array_out[EMG_SIGNAL_SIZE];
uint8_t emg_array_out_index;
int16_t emg_envelope[ENVELOPE_SIZE];
uint8_t emg_envelope_index;
float MA_SUM;



float complex vector[EMG_SIGNAL_SIZE];
float samplerate;
float meanfreq;
float medianfreq;
float peakfrequ;
float totalpower;
uint8_t emg_array_out_index_init;






  void imu_show(int16_t * arr){

    printf("\n");
    while(*arr != '\0'){
      printf("%d, ", *arr);
      arr++;
    }  
  }



/* Function definitions ----------------------------------------------- */

/**
 * @brief Application main function.
 */
int main(void)
{ 

  bsp_board_init(BSP_INIT_LEDS);
  nrf_gpio_cfg_output(LED1);
  nrf_gpio_cfg_input(PIN_24, NRF_GPIO_PIN_PULLUP);
  nrf_gpio_pin_clear(LED1);

/**
service to others is the rent we pay
  for your room here on earth
**/


  // g_device.is_device_on = false;


  // use timer and toggle toggle state of LED
  // custom
  timers_init();




  // initialize modules
  user_uart_init(); 
  power_management_init();

  // bsp_hw_init configures DRDY as an input 
  bsp_hw_init();
  bsp_afe_init();

  // initialize BNO085 over I2C
  bno085_init();

  // set default mode to record data
  // g_device.mode = SYS_MODE_RECORD_DATA;
  setReports();

//home/ebenezer/Developer/aquatic-emg-sensor-main_x2/

  // custom
  application_timers_start();

  uint8_t imu_flags = 0;



  // enter MCU superloop
  for (;;) { 


    if(!nrf_gpio_pin_read(PIN_24)){
      nrf_delay_ms(100);
      if(!nrf_gpio_pin_read(PIN_24)){
        nrf_gpio_pin_toggle(LED1);
      }
    }



 
    if (bsp_afe_get_ecg(&emg_value) == BS_OK) {

      // get raw signal input

      float emg_value_raw = emg_value;
      //printf("\nfloat: %f\n", emg_value_raw);

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
      //printf("\nMA_SUM / MA_FILTER_SIZE: %f\n", MA_SUM / MA_FILTER_SIZE);

      // Assigning complex vector for FFT
      vector[emg_array_out_index++] = MA_SUM / MA_FILTER_SIZE;

      // reset array out index
      //if (emg_array_out_index >= EMG_SIGNAL_SIZE - 1)
      //        emg_array_out_index = 0;

      // a pointer to emg_array_out buffer
      float * c = emg_array_out;

      // print contents of emg_array_out ie MAV's
      while( *c != '\0'){
              emg_integrated += (*c < 0) ? (*c * -1) : *c;
              emg_ssi += *c * *c;
              c++;
      }

      emg_mean_absolute_value = emg_integrated / EMG_SIGNAL_SIZE;
      emg_variance = emg_ssi / (EMG_SIGNAL_SIZE - 1);
      emg_rms = sqrt(emg_ssi / EMG_SIGNAL_SIZE);



      //printf("\n\n");
      //printf("MA_SUM: %f\n", MA_SUM);
      //printf("\nemg_integrated: %f", emg_integrated);
      //printf("\nemg_ssi: %f", emg_ssi);
      //printf("\nemg_variance: %f", emg_variance);
      //printf("\nemg_rms: %f", emg_rms);
 

      MA_SUM = 0;
      emg_integrated = 0.0;
      emg_mean_absolute_value = 0.0;
      emg_ssi = 0.0;
      emg_rms = 0.0;
      emg_variance = 0.0;

      printf("array_out_index: %d\n", emg_array_out_index);

      // Calling frequency domain features
      if (emg_array_out_index == EMG_SIGNAL_SIZE)
      {
          emg_array_out_index_init = emg_array_out_index+1;
          fft(vector, emg_array_out_index_init);

          printf("in frequency domain:\n");

          for(size_t x = 0; x < emg_array_out_index_init; x++) {
                  printf("%f%+fi\n", creal(vector[x]), cimag(vector[x]));
          }

          meanfreq = mean(vector, emg_array_out_index_init, samplerate);
          medianfreq = median(vector, emg_array_out_index_init, samplerate);
          peakfrequ = peakfreq(vector, emg_array_out_index_init, samplerate);
          totalpower = powe(vector,emg_array_out_index_init);
          //printf("EMG VALUE RAW: %lf\n", emg_value_raw);
          printf("Mean Frequency: %lf\n", meanfreq);
          printf("Median Frequency: %lf\n", medianfreq);
          printf("Peak Frequency: %lf\n", peakfrequ);	
          printf("Total Power: %lf\n", totalpower);
          //printf("EMG Array Out Index %u\n",  emg_array_out_index);

      }

   // reset array out index
      if (emg_array_out_index >= EMG_SIGNAL_SIZE)
              emg_array_out_index = 0;
 





        //u8p_send(0, emg_value); 

        // sample IMU at 50 Hz = 1/20 sps 
        if (msec % 20 == 0) {
          // check if IMU was reset 
          if (bno085_was_reset()) setReports();

          // get data from IMU
          if (bno085_get_sensor_event(&sensor_value))
          {
            switch (sensor_value.sensorId)
            {
              case SH2_ACCELEROMETER:
                imu[0] = sensor_value.un.accelerometer.x * SCALE_FACTOR_ACC;
                imu[1] = sensor_value.un.accelerometer.y * SCALE_FACTOR_ACC;
                imu[2] = sensor_value.un.accelerometer.z * SCALE_FACTOR_ACC;
                imu_flags |= 0b001;
                break; 
              case SH2_GYROSCOPE_CALIBRATED:
                imu[3] = sensor_value.un.gyroscope.x * SCALE_FACTOR_GYR; 
                imu[4] = sensor_value.un.gyroscope.y * SCALE_FACTOR_GYR; 
                imu[5] = sensor_value.un.gyroscope.z * SCALE_FACTOR_GYR; 
                imu_flags |= 0b010;
                break;
              case SH2_MAGNETIC_FIELD_CALIBRATED:
                imu[6] = sensor_value.un.magneticField.x * SCALE_FACTOR_MAG;
                imu[7] = sensor_value.un.magneticField.y * SCALE_FACTOR_MAG;
                imu[8] = sensor_value.un.magneticField.z * SCALE_FACTOR_MAG;
                imu_flags |= 0b100;
                break; 
              default: break;
            } 

 

            //imu_show(imu);

          } // get data from IMU


        } // sample IMU at 50 Hz = 1/20 sps

        // increment ms counter 
        msec++;  
        // decrease semaphore 

       

    } // get emg_raw
    




   
  }
}


/**
 * @brief         Function for sending UART characters
 *
 * @param[in]     field   Field identifier 
 * @param[in]     value   16-bit value for transmission  
 *
 * @attention     None
 *
 * @return        None
 */
void u8p_send(uint8_t field, uint16_t value) {
    //printf("hello\n");
    // break 16-bit value into 6-6-4 groups 
    uint8_t bits[] = {
        (value >> 12) & 0xf,
        (value >> 6) & 0x3f,
        value & 0x3f
    };

    // break field ID into 2-3
    uint8_t fields[] = {
        (field >> 2) & 0x7,
        (field & 0x3) << 4
    };

    // identifier for a byte in the middle of sequence 
    static const uint8_t continue_byte = 0x80;

    static const uint8_t two_prefix = 0xc0;

    // identifier for EMG data header 
    static const uint8_t three_prefix = 0xe0;

    // identifier for other data header 
    static const uint8_t four_prefix = 0xf0;
    
    if (field) {
        // if field != 0, send 4-byte sequence, if not sending EMG
        TX_BYTE(four_prefix | fields[0]);
        TX_BYTE(continue_byte | fields[1] | bits[0]);
        TX_BYTE(continue_byte | bits[1]);
        TX_BYTE(continue_byte | bits[2]);
    } 
    else if ((value & ~0x7f) == 0) {
        // value of EMG is less than 127  
        TX_BYTE(value); 
    }
    else if ((value & ~0x7ff) == 0) {
        // value of EMG is between 128 and 2047 
        TX_BYTE(two_prefix | bits[1]);
        TX_BYTE(continue_byte | bits[2]);
    }
    else {
        // if field == 0, send 3-byte sequence, if sending EMG
        TX_BYTE(three_prefix | bits[0]);
        TX_BYTE(continue_byte | bits[1]);
        TX_BYTE(continue_byte | bits[2]);
    }
}

// Here is where you define the sensor outputs you want to receive
void setReports(void) {
  // NRF_LOG_INFO("Setting desired reports");
  static const sh2_SensorId_t sensors[] = {
    SH2_ACCELEROMETER,
    SH2_GYROSCOPE_CALIBRATED,
    SH2_MAGNETIC_FIELD_CALIBRATED
  };

  for (int i = 0; i < sizeof(sensors) / sizeof(sh2_SensorId_t); i++) {
    sh2_SensorId_t sensor = sensors[i];

    if (!bno085_enable_report(sensor, 5000)) {
      // NRF_LOG_INFO("Could not enable accelerometer");
      u8p_send(30, 0xE000 | sensor); 
    }
  }
}

/**
 * @brief         Function for assert macro callback.
 *
 * @param[in]     line_num     Line number of the failing ASSERT call.
 * @param[in]     p_file_name  File name of the failing ASSERT call.
 *
 * @attention     None
 *
 * @return        None
 */
void assert_nrf_callback(uint16_t line_num, const uint8_t *p_file_name)
{
  app_error_handler(DEAD_BEEF, line_num, p_file_name);
}

/**
 * @brief         Function for handling Queued Write Module errors.
 *
 * @param[in]     nrf_error   Error code containing information about what went wrong.
 *
 * @attention     None
 *
 * @return        None
 */
static void nrf_qwr_error_handler(uint32_t nrf_error)
{
  APP_ERROR_HANDLER(nrf_error);
}

/**
 * @brief         Function for putting the chip into sleep mode.
 *
 * @param[in]     None
 *
 * @attention     None
 *
 * @return        None
 */
static void sleep_mode_enter(void)
{
  uint32_t err_code = bsp_indication_set(BSP_INDICATE_IDLE);
  APP_ERROR_CHECK(err_code);

  // Prepare wakeup buttons.
  err_code = bsp_btn_ble_sleep_mode_prepare();
  APP_ERROR_CHECK(err_code);

  // Go to system-off mode (this function will not return; wakeup will cause a reset).
  err_code = sd_power_system_off();
  APP_ERROR_CHECK(err_code);
}

/**
 * @brief         Function for handling events from the BSP module.
 *
 * @param[in]     event   Event generated by button press.
 *
 * @attention     None
 *
 * @return        None
 */
void bsp_event_handler(bsp_event_t event)
{
  uint32_t err_code;
  switch (event)
  {
  case BSP_EVENT_SLEEP:
    sleep_mode_enter();
    break;

  case BSP_EVENT_DISCONNECT:
    err_code = sd_ble_gap_disconnect(m_conn_handle, BLE_HCI_REMOTE_USER_TERMINATED_CONNECTION);
    if (err_code != NRF_ERROR_INVALID_STATE)
    {
      APP_ERROR_CHECK(err_code);
    }
    break;

  case BSP_EVENT_WHITELIST_OFF:
    if (m_conn_handle == BLE_CONN_HANDLE_INVALID)
    {
      err_code = ble_advertising_restart_without_whitelist(&m_advertising);
      if (err_code != NRF_ERROR_INVALID_STATE)
      {
        APP_ERROR_CHECK(err_code);
      }
    }
    break;

  default:
    break;
  }
}

/**
 * @brief         Function for initializing the nrf log module.
 *
 * @param[in]     None
 *
 * @attention     None
 *
 * @return        None
 */
static void log_init(void)
{
  ret_code_t err_code = NRF_LOG_INIT(app_timer_cnt_get);
  APP_ERROR_CHECK(err_code);

  NRF_LOG_DEFAULT_BACKENDS_INIT();
}

/**
 * @brief         Function for initializing power management.
 *
 * @param[in]     None
 *
 * @attention     None
 *
 * @return        None
 */
static void power_management_init(void)
{
  ret_code_t err_code;
  err_code = nrf_pwr_mgmt_init();
  APP_ERROR_CHECK(err_code);
}

/**
 * @brief         Function for handling the idle state (main loop).
 *
 * @param[in]     None
 *
 * @attention     None
 *
 * @return        None
 */
static void idle_state_handle(void)
{
  if (NRF_LOG_PROCESS() == false)
  {
    nrf_pwr_mgmt_run();
  }
}

/**
 * @brief         Function for handling the Battery measurement timer timeout.
 *
 * @param[in]     p_context   Pointer to context
 *
 * @attention     None
 *
 * @return        None
 */
static void battery_level_meas_timeout_handler(void * p_context)
{
  UNUSED_PARAMETER(p_context);
  battery_level_update();
}

/**
 * @brief         Function for handling the Body temperature measurement timer timeout.
 *
 * @param[in]     p_context   Pointer to context
 *
 * @attention     None
 *
 * @return        None
 */
static void sensors_meas_timeout_handler(void * p_context)
{
  UNUSED_PARAMETER(p_context);
  // sensors_value_update();
}

/**
 * @brief         Function for handling the battery level update
 *
 * @param[in]     None
 *
 * @attention     None
 *
 * @return        None
 */
static void battery_level_update(void)
{
  uint8_t battery_level              = 0;
  static  uint8_t battery_cal_time   = 0;
  static  uint16_t sum_battery_level = 0;

  sys_bm_get_level_in_percent(&battery_level);
  sum_battery_level += battery_level;
  battery_cal_time ++;

  if (battery_cal_time >= 10)
  {
    battery_level = sum_battery_level / battery_cal_time;
    // NRF_LOG_INFO( "Battery avg : %d percent", battery_level);
    battery_cal_time  = 0;
    sum_battery_level = 0;

    // ble_bas_battery_level_update(&m_bas, battery_level, BLE_CONN_HANDLE_ALL);
  }
}


/** LED timeout handler **/
static void led_timeout_handler(void * p_contex)
{

  // every interval toggle state of led
  //nrf_gpio_pin_toggle(LED1);
  //printf("Hello");

}

/**@brief Function for starting applicatin timers

*/
static void timers_init(void)
{
  // initialize timer module
  //ret_code_t err_code;

  //err_code = app_timer_init();
  //APP_ERROR_CHECK(err_code);

  // Create timers.
  //err_code = app_timer_create(&m_led_timer_id, APP_TIMER_MODE_REPEATED, led_timeout_handler);
 // APP_ERROR_CHECK(err_code);
}

static void application_timers_start(void)
{
  //ret_code_t err_code;

  //err_code = app_timer_start(m_led_timer_id, LED_STATE_INTERVAL, NULL);
 // APP_ERROR_CHECK(err_code);
}







// start of frequency domain processing
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

















/* End of file -------------------------------------------------------- */
