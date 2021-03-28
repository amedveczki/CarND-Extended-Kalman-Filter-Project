#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
			float noise_ax, float noise_ay)
    : noise_ax_(noise_ax)
    , noise_ay_(noise_ay)
{
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
}

void KalmanFilter::Predict(float dt) {
    F_(0, 2) = dt;
    F_(1, 3) = dt;

    const float dt2 = dt*dt;
    const float dt3 = dt2 * dt / 2.0f;
    const float dt4 = dt3 * dt / 2.0f;

    // Q: noise covariance matrix
    Matrix4d Q(kf_.Q_ << dt4 * noise_ax, 0,              dt3 * noise_ax, 0,
	      0             , dt4 * noise_ay, 0             , dt3 * noise_ay,
	      dt3 * noise_ax, 0,              dt2 * noise_ax, 0,
	      0             , dt3 * noise_ay, 0             , dt2 * noise_ay;
    x_ = F_ * x_;
    MatrixXd Ft = F_.transpose();
    P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z, const MatrixXd& H, const MatrixXd& R) {
  const MatrixXd y = z - H*x_;
  const MatrixXd S = H*P_*H.transpose() + R;
  const MatrixXd K = P_*H.transpose()*S.inverse();
  const MatrixXd I = MatrixXd::Identity(2, 2);


  x_ = x_ + K*y;
  P_ = (I - K*H_laser)*P;
}

void KalmanFilter::UpdateEKF(const VectorXd &z, const MatrixXd& H, const MatrixXd& R) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */

  const MatrixXd y = z - H*x_;
  const MatrixXd S = H*P_*H.transpose() + R;
  const MatrixXd K = P_*H.transpose()*S.inverse();
  const MatrixXd I = MatrixXd::Identity(2, 2);


  x_ = x_ + K*y;
  P_ = (I - K*H_laser)*P;
}
