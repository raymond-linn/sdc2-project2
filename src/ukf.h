#ifndef UKF_H
#define UKF_H
#include "Eigen/Dense"
#include "measurement_package.h"
#include <vector>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  ///* time when the state is true, in us
  long previous_timestamp_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* Process Noise
  VectorXd process_noise_;

  ///* Radar covariance matrix
  MatrixXd R_radar_;

  ///* Lidar covariance matrix
  MatrixXd R_lidar_;

  ///* State dimension
  int n_x_;

  ///* Augmented state dimension
  int n_aug_;

  ///* Sigma dimension
  int n_sig_;

  ///* Radar measurement space dimension
  int n_z_radar_;

  ///* Laser measurement space dimension
  int n_z_laser_;

  ///* Sigma point spreading parameter
  double lambda_;

  ///* the current NIS for radar
  double NIS_radar_;

  ///* the current NIS for laser
  double NIS_laser_;

  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);

private:
  /**
   * Generate the Sigma points using ugmented coavariance
   * @param X
   * @param P
   * @param processNoise
   * @return Xsig_aug
   */
  MatrixXd GenerateSigmaPoints(const VectorXd& X, const MatrixXd& P, const VectorXd processNoise);

  /**
   * Predict the sigma points
   * @param Xsig_aug
   * @param delta_t
   * @return Xsig_pred
   */
  MatrixXd PredictSigmaPoints(const MatrixXd& Xsig_aug, double delta_t);

  /**
   * Predict mean
   * @param Xsig_pred
   * @param weights
   * @return x (predicted mean)
   */
  VectorXd PredictMean(const MatrixXd& Xsig_pred, const VectorXd& weights);

  /**
   * Predict Covariance
   * @param Xsig_pred
   * @param x
   * @param weights
   * @return P (predicted covariance)
   */
  MatrixXd PredictCovariance(const MatrixXd& Xsig_pred, const VectorXd& X, const VectorXd& weights);

  /**
   *
   * @param Xsig_pred
   * @return Zsig
   */
  MatrixXd TransformSigmaToLidarSpace(const MatrixXd& Xsig_pred);

  /**
   *
   * @param Xsig_pred
   * @return Zsig
   */
  MatrixXd TransformSigmaToRadarSpace(const MatrixXd& Xsig_pred);

  /**
   *
   * @param Xsig
   * @param x
   * @param Zsig
   * @param z
   * @param weights
   * @return Tc (cross correlation matrix)
   */
  MatrixXd CalculateCrossCorrelationMatrix(const MatrixXd& Xsig,
                                           const VectorXd& x,
                                           const MatrixXd& Zsig,
                                           const VectorXd& z,
                                           const VectorXd& weights);



};

#endif /* UKF_H */
