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
			double noise_ax, double noise_ay)
    : noise_ax_(noise_ax)
    , noise_ay_(noise_ay)
{
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
}

void KalmanFilter::Predict(double dt) {
    F_(0, 2) = dt;
    F_(1, 3) = dt;

    const double dt2 = dt*dt;
    const double dt3 = dt2 * dt / 2.0f;
    const double dt4 = dt3 * dt / 2.0f;

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
  P_ = (I - K*H)*P;
}

void KalmanFilter::UpdateEKF(const VectorXd &z, const MatrixXd& H, const MatrixXd& R) {
  // We should convert our X state to (rho, phi, rho') (polar coordinates + speed)
  const double px = x_[0],
	      py = x_[1],
	      vx = x_[2],
	      vy = x_[3],

  const double rho = sqrt(px*px + py*py);
  const double phi = px == 0 ? 0 : atan2(py, px);
  const double rho_ = rho ? (px*vx = py*vy)/rho : 0;

  if (rho == 0)
      std::cout << " Warning - rho == 0, avoiding div zero!" << std::endl;

  const Vector3d z(rho, phi, rho_);

  MatrixXd y = z_ - H*x_;

  // Normalize
  if (y[1] < -EIGEN_PI)
      y[1] += 2*EIGEN_PI;

  if (y[1] > EIGEN_PI)
      y[1] -= 2*EIGEN_PI;

  const MatrixXd S = H*P_*H.transpose() + R;
  const MatrixXd K = P_*H.transpose()*S.inverse();
  const MatrixXd I = MatrixXd::Identity(2, 2);

  x_ = x_ + K*y;
  P_ = (I - K*H)*P;
}
