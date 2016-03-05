#include "imu_sensor_fusion.h"
#include "imu_sensor.h"
#include "math.h"
#include "app.h"

#define GYROSCOPE_SENSITIVITY     4.375
#define ACCELERATOR_SENSITIVITY   0.061
#define GYRO_DEGREE_OFFSET        3.0
#define ACC_1_G                   1000.0

#define M_PI 3.1415926

#define dt 0.00125							// 1.25 ms sample rate!    

void complementary_filter(float acc_raw[3], float gyr_raw[3], float mag_raw[3], float *pitch, float *roll, float *yaw)
{
    float pitch_acc, roll_acc, yaw_mag;

    // Integrate the gyroscope data -> int(angularSpeed) = angle
    *pitch += gyr_raw[0] * dt; // Angle around the X-axis
    *roll -= gyr_raw[1] * dt;    // Angle around the Y-axis
    *yaw  += gyr_raw[2] * dt;
//    // Turning around the X axis results in a vector on the Y-axis
    pitch_acc = atan2f(acc_raw[0], acc_raw[2]) * 180 / M_PI;
    *pitch = *pitch * 0.98 + pitch_acc * 0.02;
    // Turning around the Y axis results in a vector on the X-axis
    roll_acc = atan2f(acc_raw[1], acc_raw[2]) * 180 / M_PI;
    *roll = *roll * 0.98 + roll_acc * 0.02;
    yaw_mag = atan2f(mag_raw[1], mag_raw[0]) * 180 / M_PI;
    *yaw = *yaw * 0.98 + yaw_mag * 0.02;
}

void imu_sensor_fusion_1(imu_sensor_data_t* sensor_raw, sensor_fusion_angle_t* sensor_angle,  imu_sensor_fusion_1_context_t* sensor_context)
{
    float k_acc, k_gyr, k_mag;
    float delta_gyr, delta_acc, delta_mag, delta_angle;
    float forceMagnitudeApprox, gyro_offset ;

    /*gyro factor*/
    forceMagnitudeApprox = sqrt(abs(sensor_raw->gyro[0]) * abs(sensor_raw->gyro[0]) +
                                abs(sensor_raw->gyro[1]) * abs(sensor_raw->gyro[1]) +
                                abs(sensor_raw->gyro[2]) * abs(sensor_raw->gyro[2]));

    gyro_offset = sqrt( (sensor_context->gyro_offset_x * sensor_context->gyro_offset_x) +
                        (sensor_context->gyro_offset_y * sensor_context->gyro_offset_y) +
                        (sensor_context->gyro_offset_z * sensor_context->gyro_offset_z));

    if(forceMagnitudeApprox > gyro_offset)
    {
        k_gyr = sensor_context->k_gyr_1;
    } else {
        k_gyr = sensor_context->k_gyr_2;
    }

    /*acc factor*/
    forceMagnitudeApprox = sqrt(abs(sensor_raw->acc[0]) * abs(sensor_raw->acc[0]) +
                                abs(sensor_raw->acc[1]) * abs(sensor_raw->acc[1]) +
                                abs(sensor_raw->acc[2]) * abs(sensor_raw->acc[2]));

    if((abs(forceMagnitudeApprox - ACC_1_G) / ACC_1_G) > 0.01 )
    {
        k_acc = sensor_context->k_acc_2;
    } else {
        k_acc = sensor_context->k_acc_1;
    }
    
    /*mag factor*/
    k_mag = sensor_context->k_mag_1;

    /*pitch*/
    sensor_raw->gyro[0] = sensor_raw->gyro[0] -  sensor_context->gyro_offset_x;
    delta_gyr =  sensor_raw->gyro[0] * dt;
    delta_acc = (atan2f(sensor_raw->acc[0], sensor_raw->acc[2]) * 180.0 / M_PI) - sensor_angle->pitch;
    delta_angle = k_gyr * delta_gyr + k_acc * delta_acc;
    sensor_context->gyro_offset_x += (   sensor_context->k_offset * (  ( sensor_raw->gyro[0] - (delta_angle / dt) ) -  sensor_context->gyro_offset_x  )   );
    sensor_angle->pitch += delta_angle;
    
    /*roll*/
    sensor_raw->gyro[1] = sensor_raw->gyro[1] - sensor_context->gyro_offset_y;
    delta_gyr =  sensor_raw->gyro[1] * dt;
    delta_acc = (atan2f(sensor_raw->acc[1], sensor_raw->acc[2]) * 180.0 / M_PI) - sensor_angle->roll;
    delta_angle = k_gyr * delta_gyr + k_acc * delta_acc;
    sensor_context->gyro_offset_y += (   sensor_context->k_offset * (  ( sensor_raw->gyro[1] - (delta_angle / dt) ) - sensor_context->gyro_offset_y  )   );
    sensor_angle->roll += delta_angle;

    /*yaw*/
    sensor_raw->gyro[2] = sensor_raw->gyro[2] - sensor_context->gyro_offset_z;
    delta_gyr =  sensor_raw->gyro[2] * dt;
    delta_mag = (atan2f(sensor_raw->mag[1], sensor_raw->mag[0]) * 180.0 / M_PI) - sensor_angle->yaw;
    delta_angle = k_gyr * delta_gyr + k_mag * delta_mag;
    sensor_context->gyro_offset_z += (   sensor_context->k_offset * (  ( sensor_raw->gyro[2] - (delta_angle / dt) ) - sensor_context->gyro_offset_z  )   );
    sensor_angle->yaw += delta_angle;
    
    
}



#define BETADEF		0.4f		// 2 * proportional gain

static float invSqrt(float x) {
	float halfx = 0.5f * x;
	float y = x;
	long i = *(long*)&y;
	i = 0x5f3759df - (i>>1);
	y = *(float*)&i;
	y = y * (1.5f - (halfx * y * y));
	return y;
}

///////madgwick algorithm
void MadgwickAHRSupdate(float* quat, float deltaT, float gx, float gy, float gz, float ax, float ay, float az, float mx, float my, float mz) {
	float q0=quat[0];
	float q1=quat[1];
	float q2=quat[2];
	float q3=quat[3];
	float recipNorm;
	float s0, s1, s2, s3;
	float qDot1, qDot2, qDot3, qDot4;
	float hx, hy;
	float _2q0mx, _2q0my, _2q0mz, _2q1mx, _2bx, _2bz, _4bx, _4bz, _2q0, _2q1, _2q2, _2q3, _2q0q2, _2q2q3, q0q0, q0q1, q0q2, q0q3, q1q1, q1q2, q1q3, q2q2, q2q3, q3q3;

	/*// Use IMU algorithm if magnetometer measurement invalid (avoids NaN in magnetometer normalisation)
	if((mx == 0.0f) && (my == 0.0f) && (mz == 0.0f)) {
		MadgwickAHRSupdateIMU(gx, gy, gz, ax, ay, az);
		return;
	}*/

	// Rate of change of quaternion from gyroscope
	qDot1 = 0.5f * (-q1 * gx - q2 * gy - q3 * gz);
	qDot2 = 0.5f * (q0 * gx + q2 * gz - q3 * gy);
	qDot3 = 0.5f * (q0 * gy - q1 * gz + q3 * gx);
	qDot4 = 0.5f * (q0 * gz + q1 * gy - q2 * gx);

	// Compute feedback only if accelerometer measurement valid (avoids NaN in accelerometer normalisation)
	if(!((ax == 0.0f) && (ay == 0.0f) && (az == 0.0f))) {

		// Normalise accelerometer measurement
		recipNorm = invSqrt(ax * ax + ay * ay + az * az);
		ax *= recipNorm;
		ay *= recipNorm;
		az *= recipNorm;   

		// Normalise magnetometer measurement
		recipNorm = invSqrt(mx * mx + my * my + mz * mz);
		mx *= recipNorm;
		my *= recipNorm;
		mz *= recipNorm;

		// Auxiliary variables to avoid repeated arithmetic
		_2q0mx = 2.0f * q0 * mx;
		_2q0my = 2.0f * q0 * my;
		_2q0mz = 2.0f * q0 * mz;
		_2q1mx = 2.0f * q1 * mx;
		_2q0 = 2.0f * q0;
		_2q1 = 2.0f * q1;
		_2q2 = 2.0f * q2;
		_2q3 = 2.0f * q3;
		_2q0q2 = 2.0f * q0 * q2;
		_2q2q3 = 2.0f * q2 * q3;
		q0q0 = q0 * q0;
		q0q1 = q0 * q1;
		q0q2 = q0 * q2;
		q0q3 = q0 * q3;
		q1q1 = q1 * q1;
		q1q2 = q1 * q2;
		q1q3 = q1 * q3;
		q2q2 = q2 * q2;
		q2q3 = q2 * q3;
		q3q3 = q3 * q3;

		// Reference direction of Earth's magnetic field
		hx = mx * q0q0 - _2q0my * q3 + _2q0mz * q2 + mx * q1q1 + _2q1 * my * q2 + _2q1 * mz * q3 - mx * q2q2 - mx * q3q3;
		hy = _2q0mx * q3 + my * q0q0 - _2q0mz * q1 + _2q1mx * q2 - my * q1q1 + my * q2q2 + _2q2 * mz * q3 - my * q3q3;
		_2bx = sqrt(hx * hx + hy * hy);
		_2bz = -_2q0mx * q2 + _2q0my * q1 + mz * q0q0 + _2q1mx * q3 - mz * q1q1 + _2q2 * my * q3 - mz * q2q2 + mz * q3q3;
		_4bx = 2.0f * _2bx;
		_4bz = 2.0f * _2bz;

		// Gradient decent algorithm corrective step
		s0 = -_2q2 * (2.0f * q1q3 - _2q0q2 - ax) + _2q1 * (2.0f * q0q1 + _2q2q3 - ay) - _2bz * q2 * (_2bx * (0.5f - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx) + (-_2bx * q3 + _2bz * q1) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my) + _2bx * q2 * (_2bx * (q0q2 + q1q3) + _2bz * (0.5f - q1q1 - q2q2) - mz);
		s1 = _2q3 * (2.0f * q1q3 - _2q0q2 - ax) + _2q0 * (2.0f * q0q1 + _2q2q3 - ay) - 4.0f * q1 * (1 - 2.0f * q1q1 - 2.0f * q2q2 - az) + _2bz * q3 * (_2bx * (0.5f - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx) + (_2bx * q2 + _2bz * q0) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my) + (_2bx * q3 - _4bz * q1) * (_2bx * (q0q2 + q1q3) + _2bz * (0.5f - q1q1 - q2q2) - mz);
		s2 = -_2q0 * (2.0f * q1q3 - _2q0q2 - ax) + _2q3 * (2.0f * q0q1 + _2q2q3 - ay) - 4.0f * q2 * (1 - 2.0f * q1q1 - 2.0f * q2q2 - az) + (-_4bx * q2 - _2bz * q0) * (_2bx * (0.5f - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx) + (_2bx * q1 + _2bz * q3) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my) + (_2bx * q0 - _4bz * q2) * (_2bx * (q0q2 + q1q3) + _2bz * (0.5f - q1q1 - q2q2) - mz);
		s3 = _2q1 * (2.0f * q1q3 - _2q0q2 - ax) + _2q2 * (2.0f * q0q1 + _2q2q3 - ay) + (-_4bx * q3 + _2bz * q1) * (_2bx * (0.5f - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx) + (-_2bx * q0 + _2bz * q2) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my) + _2bx * q1 * (_2bx * (q0q2 + q1q3) + _2bz * (0.5f - q1q1 - q2q2) - mz);
		recipNorm = invSqrt(s0 * s0 + s1 * s1 + s2 * s2 + s3 * s3); // normalise step magnitude
		s0 *= recipNorm;
		s1 *= recipNorm;
		s2 *= recipNorm;
		s3 *= recipNorm;

		// Apply feedback step
		qDot1 -= BETADEF * s0;
		qDot2 -= BETADEF * s1;
		qDot3 -= BETADEF * s2;
		qDot4 -= BETADEF * s3;
	}

	// Integrate rate of change of quaternion to yield quaternion
	q0 += qDot1 * (1.0f * deltaT);
	q1 += qDot2 * (1.0f * deltaT);
	q2 += qDot3 * (1.0f * deltaT);
	q3 += qDot4 * (1.0f * deltaT);

	// Normalise quaternion
	recipNorm = invSqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);
	q0 *= recipNorm;
	q1 *= recipNorm;
	q2 *= recipNorm;
	q3 *= recipNorm;
	quat[0]=q0;
	quat[1]=q1;
	quat[2]=q2;
	quat[3]=q3;
}

#define twoKpDef	(2.0f * 0.5f)	// 2 * proportional gain
#define twoKiDef	(2.0f * 0.0f)	// 2 * integral gain


static float twoKp = twoKpDef;											// 2 * proportional gain (Kp)
static float twoKi = twoKiDef;											// 2 * integral gain (Ki)
static float integralFBx = 0.0f,  integralFBy = 0.0f, integralFBz = 0.0f;	// integral error terms scaled by Ki

//====================================================================================================
// Functions

//---------------------------------------------------------------------------------------------------
// AHRS algorithm update

void MahonyAHRSupdate(float* quat, float deltaT, float gx, float gy, float gz, float ax, float ay, float az, float mx, float my, float mz) {
	float q0=quat[0];
	float q1=quat[1];
	float q2=quat[2];
	float q3=quat[3];
	float recipNorm;
  float q0q0, q0q1, q0q2, q0q3, q1q1, q1q2, q1q3, q2q2, q2q3, q3q3;  
	float hx, hy, bx, bz;
	float halfvx, halfvy, halfvz, halfwx, halfwy, halfwz;
	float halfex, halfey, halfez;
	float qa, qb, qc;

	// Compute feedback only if accelerometer measurement valid (avoids NaN in accelerometer normalisation)
	if(!((ax == 0.0f) && (ay == 0.0f) && (az == 0.0f))) {

		// Normalise accelerometer measurement
		recipNorm = invSqrt(ax * ax + ay * ay + az * az);
		ax *= recipNorm;
		ay *= recipNorm;
		az *= recipNorm;     

		// Normalise magnetometer measurement
		recipNorm = invSqrt(mx * mx + my * my + mz * mz);
		mx *= recipNorm;
		my *= recipNorm;
		mz *= recipNorm;   

        // Auxiliary variables to avoid repeated arithmetic
        q0q0 = q0 * q0;
        q0q1 = q0 * q1;
        q0q2 = q0 * q2;
        q0q3 = q0 * q3;
        q1q1 = q1 * q1;
        q1q2 = q1 * q2;
        q1q3 = q1 * q3;
        q2q2 = q2 * q2;
        q2q3 = q2 * q3;
        q3q3 = q3 * q3;   

        // Reference direction of Earth's magnetic field
        hx = 2.0f * (mx * (0.5f - q2q2 - q3q3) + my * (q1q2 - q0q3) + mz * (q1q3 + q0q2));
        hy = 2.0f * (mx * (q1q2 + q0q3) + my * (0.5f - q1q1 - q3q3) + mz * (q2q3 - q0q1));
        bx = sqrt(hx * hx + hy * hy);
        bz = 2.0f * (mx * (q1q3 - q0q2) + my * (q2q3 + q0q1) + mz * (0.5f - q1q1 - q2q2));

		// Estimated direction of gravity and magnetic field
		halfvx = q1q3 - q0q2;
		halfvy = q0q1 + q2q3;
		halfvz = q0q0 - 0.5f + q3q3;
        halfwx = bx * (0.5f - q2q2 - q3q3) + bz * (q1q3 - q0q2);
        halfwy = bx * (q1q2 - q0q3) + bz * (q0q1 + q2q3);
        halfwz = bx * (q0q2 + q1q3) + bz * (0.5f - q1q1 - q2q2);  
	
		// Error is sum of cross product between estimated direction and measured direction of field vectors
		halfex = (ay * halfvz - az * halfvy) + (my * halfwz - mz * halfwy);
		halfey = (az * halfvx - ax * halfvz) + (mz * halfwx - mx * halfwz);
		halfez = (ax * halfvy - ay * halfvx) + (mx * halfwy - my * halfwx);

		// Compute and apply integral feedback if enabled
		if(twoKi > 0.0f) {
			integralFBx += twoKi * halfex * (1.0f * deltaT);;	// integral error scaled by Ki
			integralFBy += twoKi * halfey * (1.0f * deltaT);
			integralFBz += twoKi * halfez * (1.0f * deltaT);
			gx += integralFBx;	// apply integral feedback
			gy += integralFBy;
			gz += integralFBz;
		}
		else {
			integralFBx = 0.0f;	// prevent integral windup
			integralFBy = 0.0f;
			integralFBz = 0.0f;
		}

		// Apply proportional feedback
		gx += twoKp * halfex;
		gy += twoKp * halfey;
		gz += twoKp * halfez;
	}
	
	// Integrate rate of change of quaternion
	gx *= (0.5f * (1.0f * deltaT));		// pre-multiply common factors
	gy *= (0.5f * (1.0f * deltaT));
	gz *= (0.5f * (1.0f * deltaT));
	qa = q0;
	qb = q1;
	qc = q2;
	q0 += (-qb * gx - qc * gy - q3 * gz);
	q1 += (qa * gx + qc * gz - q3 * gy);
	q2 += (qa * gy - qb * gz + q3 * gx);
	q3 += (qa * gz + qb * gy - qc * gx); 
	
	// Normalise quaternion
	recipNorm = invSqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);
	q0 *= recipNorm;
	q1 *= recipNorm;
	q2 *= recipNorm;
	q3 *= recipNorm;
	quat[0]=q0;
	quat[1]=q1;
	quat[2]=q2;
	quat[3]=q3;
}


