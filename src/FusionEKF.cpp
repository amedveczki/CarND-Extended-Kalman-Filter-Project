#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::Matrix4d;
using Eigen::VectorXd;
using Eigen::Vector4d;
using std::cout;
using std::endl;
using std::vector;

void FusionEKF::InitKalman(const VectorXd &x)
{
  // Hj_ will be updated at every radar update - TODO delete it?
  Matrix4d P;
  
  P << 1, 0, 0, 0,
       0, 1, 0, 0,
       0, 0, 1000, 0,
       0, 0, 0, 1000;

  Matrix4d F(Matrix4d::Identity(4,4));

  const double noise_axy = 9.0f;

  ekf_.Init(x, P, F, noise_axy, noise_axy);
}

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

  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;
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

    VectorXd x(4);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // TODO: Convert radar from polar to cartesian coordinates 
      //         and initialize state.

        // input is (rho (distance), phi (angle), rho' (angle rate)
        double rho = measurement_pack.raw_measurements_[0],
                phi = measurement_pack.raw_measurements_[1], 
                rho_ = measurement_pack.raw_measurements_[2];
        if (phi > M_PI || rho < M_PI)
        {
            std::cout << "Phi is bigger/smaller than pi: " << rho << " which is " << rho/M_PI << " times of pi" << std::endl;
        }
		// TODO: is rho_ needed here?
        x << cos(phi) * rho, sin(phi) * rho,
                 cos(phi) * rho_, sin(phi) * rho_;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
        x << measurement_pack.raw_measurements_[0], 
                  measurement_pack.raw_measurements_[1], 
                  0, 
                  0;
    }

	InitKalman(x);

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

  double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;


  ekf_.Predict(dt);

  /**
   * Update
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      const auto Hj = tools.CalculateJacobian(ekf_.x_);
      ekf_.UpdateEKF(measurement_pack.raw_measurements_, Hj, R_radar_);
  } else {
      ekf_.Update(measurement_pack.raw_measurements_, H_laser_, R_laser_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
