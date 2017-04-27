#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.4;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  is_initialized_ = false;

  //set state dimension
  n_x_ = 5;

  //augmented state dimension
  n_aug_ = 7;

  //define spreading parameter
  lambda_ = 3 - n_aug_;

  //predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  //sigma points weights
  double w0 = lambda_ / (lambda_ + n_aug_);
  double wi = 0.5 / (lambda_ + n_aug_);
  weights_ = VectorXd (2 * n_aug_ + 1);
  weights_.setConstant (wi);
  weights_ (0) = w0;

  P_.setIdentity ();

  lidar_measurement_function_ = [](const VectorXd& state_vector) -> VectorXd {
    VectorXd z (2);
    z (0) = state_vector (0);
    z (1) = state_vector (1);
    return z;
  };

  lidar_fix_measurement_ = [](const VectorXd& z) -> VectorXd {
    return z;
  };

  lidar_nis_saver_ = [this](double nis) {
    NIS_laser_ = nis;
  };

  radar_measurement_function_ = [](const VectorXd& state_vector) -> VectorXd {
    VectorXd z (3);
    float x = state_vector (0);
    float y = state_vector (1);
    float v = state_vector (2);
    float psi = state_vector (3);
    z (0) = std::sqrt (x * x + y * y);
    if (z(0) > 1e-4) {
      z (1) = atan2(y, x);
      z (2) =
          (x * cos (psi) * v +
           y * sin (psi) * v) /
          z (0);
    } else {
      z.setZero ();
    }
    return z;
  };

  radar_fix_measurement_ = [](const VectorXd& z) -> VectorXd {
    VectorXd zz (3);
    zz = z;
    zz (1) = fmod(zz(1) + M_PI, 2*M_PI) - M_PI;
    return zz;
  };

  radar_nis_saver_ = [this](double nis) {
    NIS_radar_ = nis;
  };

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  // code taken from EKF project
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    // first measurement

    // CTRV model, state is
    // x, y, v, psi, psi_dot
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      double rho = meas_package.raw_measurements_ [0];
      double phi = meas_package.raw_measurements_ [1];
      double rho_dot = meas_package.raw_measurements_ [2];
      double tg_phi = tan (phi);
      double x = rho / std::sqrt (tg_phi*tg_phi + 1);
      double y = rho / std::sqrt (tg_phi*tg_phi + 1) * tg_phi;
      x_ << x, y, rho_dot, phi, 0;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      double x = meas_package.raw_measurements_ [0];
      double y = meas_package.raw_measurements_ [1];
      double psi = atan2 (y, x);
      x_ << x, y, 0, psi, 0;
    }

    previous_timestamp_ = meas_package.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && !use_radar_) {
    return;
  }

  if (meas_package.sensor_type_ == MeasurementPackage::LASER && !use_laser_) {
    return;
  }

  //compute the time elapsed between the current and previous measurements
  float dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0; //dt - expressed in seconds
  previous_timestamp_ = meas_package.timestamp_;

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    Prediction(dt);
    UpdateRadar(meas_package);
  } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    Prediction(dt);
    UpdateLidar(meas_package);
  }

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //create augmented mean state
  x_aug.setZero ();
  x_aug.head (n_x_) = x_;

  //create augmented covariance matrix
  P_aug.setZero ();
  P_aug.block (0, 0, n_x_, n_x_) = P_;
  P_aug (n_aug_ - 2, n_aug_ - 2) = std_a_ * std_a_;
  P_aug (n_aug_ - 1, n_aug_ - 1) = std_yawdd_ * std_yawdd_;

  //create square root matrix
  MatrixXd A = P_aug.llt().matrixL();
  MatrixXd a = sqrt(lambda_ + n_aug_) * A;

  //create augmented sigma points
  Xsig_aug.colwise () = x_aug;
  Xsig_aug.block (0, 1, n_aug_, n_aug_) += a;
  Xsig_aug.block (0, 1 + n_aug_, n_aug_, n_aug_) -= a;

  //predict sigma points
  Xsig_pred_.setZero ();
  int sig_count = 2*n_aug_ + 1;
  for (int i=0; i<sig_count; i++) {
      VectorXd x_prev = Xsig_aug.block (0, i, n_aug_, 1);
      double x = x_prev (0);
      double y = x_prev (1);
      double v = x_prev (2);
      double psi = x_prev (3);
      double psi_dot = x_prev (4);
      double nu_a = x_prev (5);
      double nu_psi_dd = x_prev (6);
      double dt2 = delta_t*delta_t;
      VectorXd v_part (n_x_);
      v_part <<
        dt2 * cos (psi) * nu_a / 2,
        dt2 * sin (psi) * nu_a / 2,
        delta_t * nu_a,
        dt2 * nu_psi_dd / 2,
        delta_t * nu_psi_dd;
      //avoid division by zero
      //write predicted sigma points into right column
      if (fabs(psi_dot) < 1e-2) {
        Xsig_pred_.block (0, i, n_x_, 1) <<
            delta_t * v * cos (psi),
            delta_t * v * sin (psi),
            0,
            0,
            0;
      } else {
        Xsig_pred_.block (0, i, n_x_, 1) <<
            v / psi_dot * (sin (psi + psi_dot * delta_t) - sin (psi)),
            v / psi_dot * (- cos (psi + psi_dot * delta_t) + cos (psi)),
            0,
            delta_t * psi_dot,
            0;
      }
      Xsig_pred_.block (0, i, n_x_, 1) += x_prev.head (5) + v_part;
  }
}

/**
 * Updates the state and the state covariance matrix
 * @param {MeasurementPackage} meas_package
 * @param {int} n_z - measurement dimension
 * @param {MatrixXd} S_init - measurement covariance matrix of measurement noise
 * @param {mesaurement_function_t} function that maps from state vector to measurement vector
 * @param {fix_measurement_t} function for returning angle into -pi .. pi interval or another measurement vector correction during computation
 * @param {nis_saver_t} function that may save NIS value
 */
void UKF::Update(
    MeasurementPackage meas_package,
    int n_z,
    const MatrixXd& S_init,
    mesaurement_function_t measurement_function,
    fix_measurement_t fix_measurement,
    nis_saver_t nis_saver
) {

    //predict state
    VectorXd x (n_x_);
    x.setZero ();
    for (int i=0; i<2 * n_aug_ + 1; i++) {
        x += weights_ (i) * Xsig_pred_.block (0, i, n_x_, 1);
    }
    x(3) = fmod (x(3) + M_PI, 2*M_PI) - M_PI;
    //predict state covariance matrix
    MatrixXd P (n_x_, n_x_);
    P.setZero ();
    for (int i=0; i<2 * n_aug_ + 1; i++) {
        VectorXd Ai = Xsig_pred_.block (0, i, n_x_, 1) - x;
        Ai(3) = fmod (Ai(3) + M_PI, 2*M_PI) - M_PI;
        P += weights_ (i) * Ai * Ai.transpose ();
    }

    //create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

    //mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);

    //measurement covariance matrix S
    MatrixXd S = S_init;

    //transform sigma points into measurement space
    for (int i=0; i<2 * n_aug_ + 1; i++) {
      Zsig.col (i) = measurement_function (Xsig_pred_.col (i));
    }
    //calculate mean predicted measurement
    z_pred.setZero ();
    for (int i=0; i<2 * n_aug_ + 1; i++) {
        z_pred += weights_ (i) * Zsig.block (0, i, n_z, 1);
    }
    z_pred = fix_measurement (z_pred);

    for (int i=0; i<2 * n_aug_ + 1; i++) {
        VectorXd ai = Zsig.block (0, i, n_z, 1) - z_pred;
        ai = fix_measurement (ai);
        S += weights_ (i) * ai * ai.transpose ();
    }

    //create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z);

    //calculate cross correlation matrix
    Tc.setZero ();
    for (int i=0; i<2 * n_aug_ + 1; i++) {
        VectorXd xx = Xsig_pred_.block(0, i, n_x_, 1) - x;
        xx(3) = fmod (xx(3) + M_PI, 2*M_PI) - M_PI;

        VectorXd zz = Zsig.block (0, i, n_z, 1) - z_pred;
        zz = fix_measurement (zz);

        Tc += weights_ (i) * xx * zz.transpose ();
    }

    MatrixXd S_inv = S.inverse ();

    //calculate Kalman gain K;
    MatrixXd K = Tc * S.inverse ();

    //update state mean and covariance matrix
    VectorXd z_d = meas_package.raw_measurements_ - z_pred;
    z_d = fix_measurement (z_d);
    x_ = x + K * z_d;
    x_(3) = fmod (x_(3) + M_PI, 2*M_PI) - M_PI;
    P_ = P - K * S * K.transpose ();

    double nis = z_d.transpose () * S_inv * z_d;
    nis_saver (nis);
}


/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

  //calculate measurement covariance matrix S of measurement noise
  MatrixXd S (2, 2);
  S <<
      std_laspx_ * std_laspx_, 0,
      0, std_laspy_ * std_laspy_;

  return Update(
        meas_package,
        2,
        S,
        lidar_measurement_function_,
        lidar_fix_measurement_,
        lidar_nis_saver_);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {

  //calculate measurement covariance matrix S of measurement noise
  MatrixXd S (3, 3);
  S <<
      std_radr_ * std_radr_, 0, 0,
      0, std_radphi_ * std_radphi_, 0,
      0, 0, std_radrd_ * std_radrd_;

  return Update(
        meas_package,
        3,
        S,
        radar_measurement_function_,
        radar_fix_measurement_,
        radar_nis_saver_);
}
