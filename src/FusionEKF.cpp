#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  /**
   * TODO: Finish initializing the FusionEKF.
   * TODO: Set the process and measurement noises
   */

  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  // Hj_ will be updated at every radar update

}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */

    // 2. set the process covariance matrix q

    ekf_.x_ = VectorXd(4);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // TODO: Convert radar from polar to cartesian coordinates 
      //         and initialize state.

        // input is (rho (angle), phi (hypotenuse), rho' (angle rate)
        float rho = measurement_pack.raw_measurements_[0],
                phi = measurement_pack.raw_measurements_[1], 
                rho_ = measurement_pack.raw_measurements_[2];
        if (rho > EIGEN_PI || rho < EIGEN_PI)
        {
            std::cout << "Rho is bigger/smaller than pi: " << rho << " which is " << rho/EIGEN_PI << " times of pi" << std::endl;
        }
        ekf_.x_ << cos(rho) * phi, sin(rho) * phi,
                    cos(rho) * rho_, sin(rho) * rho_;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
        ekf_.x_ << measurement_pack.raw_measurements_[0], 
                  measurement_pack.raw_measurements_[1], 
                  0, 
                  0;
    }

    // first measurement
    cout << "EKF: " << ekf_.x_ << endl;

    previous_timestamp_ = measurement_pack.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */

  /**
   * TODO: Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * TODO: Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
    float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
    previous_timestamp_ = measurement_pack.timestamp_;

    kf_.f_(0, 2) = dt;
    kf_.f_(1, 3) = dt;

    const float noise_ax = 9.0f;
    const float noise_ay = 9.0f;

    const float dt2 = dt*dt;
    const float dt3 = dt2 * dt / 2.0f;
    const float dt4 = dt3 * dt / 2.0f;

    // Q_: noise covariance matrix
    kf_.Q_ = MatrixXd(4,4);
    kf_.Q_ << dt4 * noise_ax, 0,              dt3 * noise_ax, 0,
              0             , dt4 * noise_ay, 0             , dt3 * noise_ay,
              dt3 * noise_ax, 0,              dt2 * noise_ax, 0,
              0             , dt3 * noise_ay, 0             , dt2 * noise_ay;

    ekf_.Predict();

  /**
   * Update
   */

  /**
   * TODO:
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      float rho = measurement_pack.raw_measurements_[0],
          phi = measurement_pack.raw_measurements_[1], 
          rho_ = measurement_pack.raw_measurements_[2];

	const auto Hj = tools.CalculateJacobian(ekf_.x_);

  } else {
      const auto& z = measurement_pack.raw_measurements_;
	  ekf_.H_ = H_laser;

      MatrixXd y = z - H_laser*ekf_.x_;
      MatrixXd S = H_laser*ekf_.P_*H_laser.transpose() + R_laser;
      MatrixXd K = ekf_.P_*H_laser.transpose()*S.inverse();

      x = x + K*y;
      P = (I - K*H_laser)*P;
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
