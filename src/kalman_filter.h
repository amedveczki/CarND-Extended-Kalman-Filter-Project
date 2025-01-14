#ifndef KALMAN_FILTER_H_
#define KALMAN_FILTER_H_

#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#pragma GCC diagnostic ignored "-Wint-in-bool-context"
#include "Eigen/Dense"
#pragma GCC diagnostic pop

class KalmanFilter {
 public:
  /**
   * Constructor
   */
  KalmanFilter();

  /**
   * Destructor
   */
  virtual ~KalmanFilter();

  /**
   * Init Initializes Kalman filter
   * @param x_in Initial state
   * @param P_in Initial state covariance
   * @param F_in Transition matrix (initial)
   * @param noise_ax Noise acceleration (x)
   * @param noise_ay Noise acceleration (y)
   */
  void Init(const Eigen::VectorXd &x_in, const Eigen::MatrixXd &P_in, const Eigen::MatrixXd &F_in,
			double noise_ax, double noise_ay);

  /**
   * Prediction Predicts the state and the state covariance
   * using the process model
   * @param dt Time between k and k+1 in s
   */
  void Predict(double dt);

  /**
   * Updates the state by using standard Kalman Filter equations
   * @param z The measurement at k+1
   * @param H The H measurement matrix for non-extended Kalman filtering
   * @param R The R measurement covariance matrix for non-extended Kalman filtering
   */
  void Update(const Eigen::VectorXd &z, const Eigen::MatrixXd& H, const Eigen::MatrixXd& R);

  /**
   * Updates the state by using Extended Kalman Filter equations
   * @param z The measurement at k+1
   * @param H The H measurement matrix for Extended Kalman filtering
   * @param R The R measurement covariance matrix for Extended Kalman filtering
   */
  void UpdateEKF(const Eigen::VectorXd &z, const Eigen::MatrixXd& H, const Eigen::MatrixXd& R);


  // state vector
  Eigen::VectorXd x_;

  // state covariance matrix
  Eigen::MatrixXd P_;

  // state transition matrix
  Eigen::MatrixXd F_;

  double noise_ax_;
  double noise_ay_;

/*
  // process covariance matrix
  Eigen::MatrixXd Q_;

  // measurement matrix
  Eigen::MatrixXd H_;

  // measurement covariance matrix
  Eigen::MatrixXd R_;
*/
};

#endif // KALMAN_FILTER_H_
