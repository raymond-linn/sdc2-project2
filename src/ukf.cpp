#include <iostream>
#include "ukf.h"

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ =  0.3; // 3.0;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.25; //M_PI/2;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15; //0.0225;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;//0.0225;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.15; //0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.0015; //0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.08; //0.3;

  // initialize the dimensions
  n_x_ = 5;
  n_aug_ = n_x_ + 2;
  n_sig_ = 2*n_aug_ + 1;
  n_z_radar_ = 3;
  n_z_laser_ = 2;

  lambda_ = 3 - n_aug_;

  is_initialized_ = false;

  // state vector
  x_ = VectorXd::Zero(n_x_);

  // Covariance matrix
  P_ = MatrixXd::Zero(n_x_, n_x_);
  P_ << 0.01, 0, 0, 0, 0,
          0, 0.01, 0 , 0, 0,
          0, 0, 0.01, 0, 0,
          0, 0, 0, 0.01, 0,
          0, 0, 0, 0, 0.01;

  // weights initialization
  weights_ = VectorXd(n_sig_);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int i=1; i < n_sig_; i++) {
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }

  // Radar Covariance Initialization
  R_radar_ = MatrixXd::Zero(n_z_radar_, n_z_radar_);
  R_radar_ << std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;

  // Lidar Covariance Initialization
  R_lidar_ = MatrixXd::Zero(n_z_laser_, n_z_laser_);
  R_lidar_ << std_laspx_*std_laspx_, 0,
              0, std_laspy_*std_laspy_;

  // initialize process noise
  process_noise_ = VectorXd::Zero(2);
  process_noise_ << pow(std_a_, 2), pow(std_yawdd_, 2);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  // check whether it is initialized
  if (!is_initialized_){

    // initialization
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
      double rho = meas_package.raw_measurements_[0];
      double phi = meas_package.raw_measurements_[1];
      double rho_dot = meas_package.raw_measurements_[2];
      double p_x, p_y;
      p_x = rho * cos(phi);
      p_y = rho * sin(phi);

      // check p_x and p_y are zeros then initialize them to be 1
      if (fabs(p_x) < 0.0001){
        p_x = 1;
      }
      if (fabs(p_y) < 0.0001){
        p_y = 1;
      }
      x_ << p_x, p_y, rho_dot, 0.0, 0.0;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER){
      double p_x, p_y;
      p_x = meas_package.raw_measurements_[0];
      p_y = meas_package.raw_measurements_[1];
      // check p_x and p_y are zeros then initialize them to be 1
      if (fabs(p_x) < 0.0001){
        p_x = 1;
      }
      if (fabs(p_y) < 0.0001){
        p_y = 1;
      }
      x_ << p_x, p_y, 0.0, 0.0, 0.0;
    }

    previous_timestamp_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }

  // if it is not first time
  double  dt = (meas_package.timestamp_ - previous_timestamp_)/1000000.0;
  previous_timestamp_ = meas_package.timestamp_;

  // Predict
  // dividing to smaller intervals to maintain stability
  while (dt > 0.1)
  {
    Prediction(0.05);
    dt -= 0.05;
  }

  Prediction(dt);

  // Update
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
    UpdateRadar(meas_package);
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    UpdateLidar(meas_package);
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

  // 1) Generate and Augment Sigma Points
  const MatrixXd Xsig = GenerateSigmaPoints(x_, P_, process_noise_);

  // 2) Predict Sigma Points
  Xsig_pred_ = PredictSigmaPoints(Xsig, delta_t);

  // 3) Predict Mean
  x_ = PredictMean(Xsig_pred_, weights_);

  // 4) Normalize
  while (x_(3)> M_PI) x_(3)-=2.*M_PI;
  while (x_(3)<-M_PI) x_(3)+=2.*M_PI;

  // 5) Predict Covariance
  P_ = PredictCovariance(Xsig_pred_, x_, weights_);
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

  const VectorXd z = meas_package.raw_measurements_;

  // 1) Transform to Lidar measurement space
  const MatrixXd Zsig = TransformSigmaToLidarSpace(Xsig_pred_);

  // 2) predict mean and covariance
  const VectorXd z_pred = PredictMean(Zsig, weights_);
  const MatrixXd S = PredictCovariance(Zsig, z_pred, weights_) + R_lidar_;

  // 3) calculate cross correlation  matrix
  MatrixXd Tc = CalculateCrossCorrelationMatrix(Xsig_pred_, x_, Zsig, z_pred, weights_);

  // 5) Kalman Gain
  MatrixXd K = Tc * S.inverse();
  // 6) residual
  VectorXd z_diff = z - z_pred;
  // 7) angle normalization
  //while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  //while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  // 8) update
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();

  // 9) Calculate NIS
  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {

  const VectorXd z = meas_package.raw_measurements_;

  // 1) Transform to Lidar measurement space
  const MatrixXd Zsig = TransformSigmaToRadarSpace(Xsig_pred_);

  // 2) predict mean and covariance
  const VectorXd z_pred = PredictMean(Zsig, weights_);
  const MatrixXd S = PredictCovariance(Zsig, z_pred, weights_) + R_radar_;

  // 3) calculate cross correlation  matrix
  MatrixXd Tc = CalculateCrossCorrelationMatrix(Xsig_pred_, x_, Zsig, z_pred, weights_);

  // 5) Kalman Gain
  const MatrixXd K  = Tc * S.inverse();

  // 6) residual
  VectorXd z_diff = z - z_pred;

  // 7) angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  // 8) update
  x_ = x_ + K * z_diff;

  // 9) Normalize
  while (x_(3)> M_PI) x_(3)-=2.*M_PI;
  while (x_(3)<-M_PI) x_(3)+=2.*M_PI;

  P_ = P_ - K * S * K.transpose();

  // 9) Calculate NIS
  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
}

// private function sections
/**
 * Generating augmented sigma points
 * @param X
 * @param P
 * @param processNoise
 * @return Xsig_aug
 */
MatrixXd UKF::GenerateSigmaPoints(const VectorXd& X, const MatrixXd& P, const VectorXd processNoise){

  // create augmented covariance Matrix
  MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);
  P_aug.topLeftCorner(5,5) = P;
  P_aug(5,5) = processNoise(0);
  P_aug(6,6) = processNoise(1);

  // create augmented state vector
  VectorXd x_aug = VectorXd::Zero(n_aug_);
  x_aug.head(n_x_) = X;

  // create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  // create sigma point Matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sig_);

  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }
  return Xsig_aug;
}

/**
 *
 * @param Xsig_aug
 * @param delta_t
 * @return Xsig_pred
 */
MatrixXd UKF::PredictSigmaPoints(const MatrixXd& Xsig_aug, double delta_t){

  //create matrix with predicted sigma points as columns
  MatrixXd Xsig_pred = MatrixXd(n_x_, n_sig_);

  //predict sigma points
  for (int i = 0; i< n_sig_; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
      px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
      py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
      px_p = p_x + v*delta_t*cos(yaw);
      py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred(0,i) = px_p;
    Xsig_pred(1,i) = py_p;
    Xsig_pred(2,i) = v_p;
    Xsig_pred(3,i) = yaw_p;
    Xsig_pred(4,i) = yawd_p;
  }
  return Xsig_pred;
}

/**
 *
 * @param Xsig_pred
 * @param weights
 * @return x (predicted mean)
 */
VectorXd UKF::PredictMean(const MatrixXd& Xsig_pred, const VectorXd& weights){

  VectorXd x = VectorXd::Zero(Xsig_pred.rows());
  for (int i = 0; i < n_sig_; i++){
    x = x + weights(i) * Xsig_pred.col(i);
  }
  return x;
}

/**
 *
 * @param Xsig_pred
 * @param x
 * @param weights
 * @return
 */
MatrixXd UKF::PredictCovariance(const MatrixXd& Xsig_pred, const VectorXd& X, const VectorXd& weights){

  MatrixXd P = MatrixXd::Zero(X.rows(), X.rows());
  for (int i = 0; i < n_sig_; i++){
    const VectorXd temp = Xsig_pred.col(i) - X;
    P = P + weights(i) * temp * temp.transpose();
  }
  return P;
}

/**
 *
 * @param Xsig_pred
 * @return Zsig
 */
MatrixXd UKF::TransformSigmaToLidarSpace(const MatrixXd& Xsig_pred){

  MatrixXd Zsig = MatrixXd(n_z_laser_, n_sig_);
  for (int i =0; i < n_sig_; i++){
    // extract values for better readibility
    double p_x = Xsig_pred(0,i);
    double p_y = Xsig_pred(1,i);

    // measurement model
    Zsig(0,i) = p_x;
    Zsig(1,i) = p_y;
  }
  return Zsig;
}

/**
 *
 * @param Xsig_pred
 * @return
 */
MatrixXd UKF::TransformSigmaToRadarSpace(const MatrixXd& Xsig_pred){

  MatrixXd Zsig = MatrixXd(n_z_radar_, n_sig_);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred(0,i);
    double p_y = Xsig_pred(1,i);
    double v  = Xsig_pred(2,i);
    double yaw = Xsig_pred(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    double rho = sqrt(p_x*p_x + p_y*p_y);
    double phi = atan2(p_y,p_x);
    double rho_dot = (p_x*v1 + p_y*v2 ) / rho;

    if (rho != rho){
      rho = 0;
    }

    if (phi != phi){
      phi = 0;
    }

    if (rho_dot != rho_dot) {
      rho_dot = 0;
    }
    Zsig(0,i) =  rho; //r
    Zsig(1,i) =  phi; //phi
    Zsig(2,i) =  rho_dot; //r_dot
  }
  return Zsig;
}

/**
 *
 * @param Xsig
 * @param x
 * @param Zsig
 * @param z
 * @param weights
 * @return Tc
 */
MatrixXd UKF::CalculateCrossCorrelationMatrix(const MatrixXd& Xsig,
                                         const VectorXd& x,
                                         const MatrixXd& Zsig,
                                         const VectorXd& z,
                                         const VectorXd& weights){

  MatrixXd Tc = MatrixXd::Zero(x.rows(), z.rows());
  for (int i = 0; i < n_sig_; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig.col(i) - x;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights(i) * x_diff * z_diff.transpose();
  }
  return Tc;
}