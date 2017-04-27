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
  std_a_ = 0.2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.2;

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

//  std_radr_ = 0.3;
//  std_radphi_ = 0.0175;
//  std_radrd_ = 0.1;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  is_initialized_ = false;

  //set state dimension
  n_x_ = 5;

  n_aug_ = 7;

  //define spreading parameter
  lambda_ = 3 - n_aug_;

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  double w0 = lambda_ / (lambda_ + n_aug_);
  double wi = 0.5 / (lambda_ + n_aug_);
  weights_ = VectorXd (2 * n_aug_ + 1);
  weights_.setConstant (wi);
  weights_ (0) = w0;

cout << "weights: " << endl << weights_ << endl;

  P_.setIdentity ();
//  P_ (0, 0) = 1;
//  P_ (1, 1) = 1;
//  P_ (2, 2) = 1000;
//  P_ (3, 3) = 1;
//  P_ (4, 4) = 1000;
//  P_.setIdentity ();
//  P_ *= 1;

//  P_ <<     0.0043,   -0.0013,    0.0030,   -0.0022,   -0.0020,
//          -0.0013,    0.0077,    0.0011,    0.0071,    0.0060,
//           0.0030,    0.0011,    0.0054,    0.0007,    0.0008,
//          -0.0022,    0.0071,    0.0007,    0.0098,    0.0100,
//          -0.0020,    0.0060,    0.0008,    0.0100,    0.0123;

}

UKF::~UKF() {}

VectorXd ctrvToCartesian (VectorXd ctrv) {
  VectorXd c (4);
  c (0) = ctrv (0);
  c (1) = ctrv (1);
  c (2) = ctrv (2) * cos (ctrv (3));
  c (3) = ctrv (2) * sin (ctrv (3));
  return c;
}

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
      double x = rho / std::sqrt (tg_phi*tg_phi + 1) / tg_phi;
      double y = rho / std::sqrt (tg_phi*tg_phi + 1);
      x_ << x, y, rho_dot, phi, 0;
cout << "r init:" << x_ << endl;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      double x = meas_package.raw_measurements_ [0];
      double y = meas_package.raw_measurements_ [1];
      double psi = atan2 (y, x);
      x_ << x, y, 0, psi, 0;
cout << "l init:" << x_ << endl;
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

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  //compute the time elapsed between the current and previous measurements
  float dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0; //dt - expressed in seconds
  previous_timestamp_ = meas_package.timestamp_;

  /*****************************************************************************
   *  Update
   ****************************************************************************/

//  x_ <<   5.7441,
//         1.3800,
//         2.2049,
//         0.5015,
//         0.3528;

//  P_ <<     0.0043,   -0.0013,    0.0030,   -0.0022,   -0.0020,
//          -0.0013,    0.0077,    0.0011,    0.0071,    0.0060,
//           0.0030,    0.0011,    0.0054,    0.0007,    0.0008,
//          -0.0022,    0.0071,    0.0007,    0.0098,    0.0100,
//          -0.0020,    0.0060,    0.0008,    0.0100,    0.0123;
//  dt = 0.1;

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
//cout << "radar" << endl;
    Prediction(dt);
    UpdateRadar(meas_package);
  } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
//cout << "laser" << endl;
    Prediction(dt);
    UpdateLidar(meas_package);
  }

  // print the output
  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

//  //create sigma point matrix
//  MatrixXd Xsig = MatrixXd(n_x_, 2 * n_x_ + 1);

//  //calculate square root of P
//  MatrixXd A = P_.llt().matrixL();

//  //calculate sigma points ...
//  Xsig.colwise () = x_;
//  MatrixXd a = sqrt(lambda_ + n_x_) * A;
//  Xsig.block (0, 1, n_x_, n_x_) += a;
//  Xsig.block (0, 1 + n_x_, n_x_, n_x_) -= a;

cout << "dt: " << delta_t << endl;
  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //create augmented mean state
  x_aug.setZero ();
  x_aug.head (n_x_) = x_;
cout << "x_aug: " << endl << x_aug << endl;
  //create augmented covariance matrix
  P_aug.setZero ();
  P_aug.block (0, 0, n_x_, n_x_) = P_;
  P_aug (n_aug_ - 2, n_aug_ - 2) = std_a_ * std_a_;
  P_aug (n_aug_ - 1, n_aug_ - 1) = std_yawdd_ * std_yawdd_;
cout << "P_aug: " << endl << P_aug << endl;
  //create square root matrix
  MatrixXd A = P_aug.llt().matrixL();
cout << "A: " << endl << A << endl;
  MatrixXd a = sqrt(lambda_ + n_aug_) * A;
cout << "a: " << endl << a << endl;

//create augmented sigma points
  Xsig_aug.colwise () = x_aug;
  Xsig_aug.block (0, 1, n_aug_, n_aug_) += a;
  Xsig_aug.block (0, 1 + n_aug_, n_aug_, n_aug_) -= a;
cout << "Xsig_aug: " << endl << Xsig_aug << endl;

  //predict sigma points
  Xsig_pred_.setZero ();
  int sig_count = 2*n_aug_ + 1;
  for (int i=0; i<sig_count; i++) {
      VectorXd x_prev = Xsig_aug.block (0, i, n_aug_, 1);
      x_prev(3) = fmod (x_prev(3) + M_PI, 2*M_PI) - M_PI;
cout << "x_prev head (5): " << endl << x_prev.head (5) << endl;
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
      v_part(3) = fmod (v_part(3) + M_PI, 2*M_PI) - M_PI;
cout << "v part: " << endl << v_part << endl;
      //avoid division by zero
      //write predicted sigma points into right column
      if (fabs(psi_dot) < 1e-2) {
        Xsig_pred_.block (0, i, n_x_, 1) <<
            delta_t * v * cos (psi),
            delta_t * v * sin (psi),
            0,
            0,
            0;
cout << "process 1: " << endl << Xsig_pred_.block (0, i, n_x_, 1) << endl;
      } else {
        Xsig_pred_.block (0, i, n_x_, 1) <<
            v / psi_dot * (sin (psi + psi_dot * delta_t) - sin (psi)),
            v / psi_dot * (- cos (psi + psi_dot * delta_t) + cos (psi)),
            0,
            delta_t * psi_dot,
            0;
cout << "process 2: " << endl << Xsig_pred_.block (0, i, n_x_, 1) << endl;
      }
      Xsig_pred_.block (0, i, n_x_, 1) += x_prev.head (5) + v_part;
  }

cout << "Xsig_pred: " << endl << Xsig_pred_ << endl;

}

/**
 * Updates the state and the state covariance matrix
 * @param {MeasurementPackage} meas_package
 */
void UKF::Update(
    MeasurementPackage meas_package,
    int n_z,
    const MatrixXd& S_init,
    mesaurement_function_t measurement_function,
    fix_measurement_t fix_measurement
) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

    //predict state
    VectorXd x (n_x_);
    x.setZero ();
    for (int i=0; i<2 * n_aug_ + 1; i++) {
        x += weights_ (i) * Xsig_pred_.block (0, i, n_x_, 1);
    }
    //predict state covariance matrix
    MatrixXd P (n_x_, n_x_);
    P.setZero ();
    for (int i=0; i<2 * n_aug_ + 1; i++) {
        VectorXd Ai = Xsig_pred_.block (0, i, n_x_, 1) - x;
        Ai(3) = fmod (Ai(3) + M_PI, 2*M_PI) - M_PI;
        P += weights_ (i) * Ai * Ai.transpose ();
    }
cout << "predict x: " << endl << x << endl;
cout << "predict P: " << endl << P << endl;

    //create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

    //mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);

    //measurement covariance matrix S
    MatrixXd S = S_init;

    //transform sigma points into measurement space
    for (int i=0; i<2 * n_aug_ + 1; i++) {
      Zsig.col (i) = measurement_function (Xsig_pred_.col (i));
//      float x = Xsig_pred_ (0, i);
//      float y = Xsig_pred_ (1, i);
//      Zsig (0, i) = x;
//      Zsig (1, i) = y;
    }
    //calculate mean predicted measurement
cout << "update 1" << endl;
    z_pred.setZero ();
    for (int i=0; i<2 * n_aug_ + 1; i++) {
        z_pred += weights_ (i) * Zsig.block (0, i, n_z, 1);
    }
cout << "update 2" << endl;

cout << "update 3" << endl;
    for (int i=0; i<2 * n_aug_ + 1; i++) {
        VectorXd ai = Zsig.block (0, i, n_z, 1) - z_pred;
//        ai(1) = fmod (ai(1) + M_PI, 2*M_PI) - M_PI;
        ai = fix_measurement (ai);
cout << "ai: " << endl << ai << endl;
        S += weights_ (i) * ai * ai.transpose ();
    }
cout << "z_pred" << endl << z_pred << endl;
cout << "S" << endl << S << endl;

    //create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z);

    //calculate cross correlation matrix
    Tc.setZero ();
    for (int i=0; i<2 * n_aug_ + 1; i++) {
        VectorXd xx = Xsig_pred_.block(0, i, n_x_, 1) - x;
        xx(3) = fmod (xx(3) + M_PI, 2*M_PI) - M_PI;

        VectorXd zz = Zsig.block (0, i, n_z, 1) - z_pred;
//        zz(1) = fmod (zz(1) + M_PI, 2*M_PI) - M_PI;
        zz = fix_measurement (zz);

cout << "Tc i " << endl << weights_ (i) * xx * zz.transpose () << endl;
        Tc += weights_ (i) * xx * zz.transpose ();
    }
    //calculate Kalman gain K;
cout << "Tc" << endl << Tc << endl;
cout << "S inv" << endl << S.inverse () << endl;
    MatrixXd K = Tc * S.inverse ();
cout << "K" << endl << K << endl;
    //update state mean and covariance matrix
  cout << "measurement: " << endl;
  cout << meas_package.raw_measurements_ << endl;

    VectorXd z_d = meas_package.raw_measurements_ - z_pred;
//    z_d(1) = fmod (z_d(1) + M_PI, 2*M_PI) - M_PI;
    z_d = fix_measurement (z_d);
    x_ = x + K * z_d;
    x_(3) = fmod (x_(3) + M_PI, 2*M_PI) - M_PI;
  //  x_ += K * (z - z_pred);
    P_ = P - K * S * K.transpose ();
}


/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

  //calculate measurement covariance matrix S
  MatrixXd S (2, 2);
  S <<
      std_laspx_ * std_laspx_, 0,
      0, std_laspy_ * std_laspy_;

  return Update(
        meas_package,
        2,
        S,
        [](const VectorXd& state_vector) -> VectorXd {
          VectorXd z (2);
          z (0) = state_vector (0);
          z (1) = state_vector (1);
          return z;
        },
        [](const VectorXd& z) -> VectorXd {
          return z;
        }
        );


//  /**
//  TODO:

//  Complete this function! Use lidar data to update the belief about the object's
//  position. Modify the state vector, x_, and covariance, P_.

//  You'll also need to calculate the lidar NIS.
//  */

//    //predict state
//    VectorXd x (n_x_);
//    x.setZero ();
//    for (int i=0; i<2 * n_aug_ + 1; i++) {
//        x += weights_ (i) * Xsig_pred_.block <5, 1>(0, i);
//    }
//    //predict state covariance matrix
//    MatrixXd P (n_x_, n_x_);
//    P.setZero ();
//    for (int i=0; i<2 * n_aug_ + 1; i++) {
//        VectorXd Ai = Xsig_pred_.block <5, 1>(0, i) - x;
//        Ai(3) = fmod (Ai(3) + M_PI, 2*M_PI) - M_PI;
//        P += weights_ (i) * Ai * Ai.transpose ();
//    }
//cout << "predict x: " << endl << x << endl;
//cout << "predict P: " << endl << P << endl;

//    int n_z = 2;

//    //create matrix for sigma points in measurement space
//    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

//    //mean predicted measurement
//    VectorXd z_pred = VectorXd(n_z);

//    //measurement covariance matrix S
//    MatrixXd S = MatrixXd(n_z,n_z);

//    //transform sigma points into measurement space
//    for (int i=0; i<2 * n_aug_ + 1; i++) {
//      float x = Xsig_pred_ (0, i);
//      float y = Xsig_pred_ (1, i);
//      Zsig (0, i) = x;
//      Zsig (1, i) = y;
//    }
//    //calculate mean predicted measurement
//    z_pred.setZero ();
//    for (int i=0; i<2 * n_aug_ + 1; i++) {
//        z_pred += weights_ (i) * Zsig.block <2, 1> (0, i);
//    }
//    //calculate measurement covariance matrix S
//    S <<
//        std_laspx_ * std_laspx_, 0,
//        0, std_laspy_ * std_laspy_;

//    for (int i=0; i<2 * n_aug_ + 1; i++) {
//        VectorXd ai = Zsig.block <2, 1> (0, i) - z_pred;
//        ai(1) = fmod (ai(1) + M_PI, 2*M_PI) - M_PI;
//cout << "ai: " << endl << ai << endl;
//        S += weights_ (i) * ai * ai.transpose ();
//    }
//cout << "z_pred" << endl << z_pred << endl;
//cout << "S" << endl << S << endl;

//    //create matrix for cross correlation Tc
//    MatrixXd Tc = MatrixXd(n_x_, n_z);

//    //calculate cross correlation matrix
//    Tc.setZero ();
//    for (int i=0; i<2 * n_aug_ + 1; i++) {
//        VectorXd xx = Xsig_pred_.block(0, i, n_x_, 1) - x;
//        xx(3) = fmod (xx(3) + M_PI, 2*M_PI) - M_PI;

//        VectorXd zz = Zsig.block (0, i, n_z, 1) - z_pred;
//        zz(1) = fmod (zz(1) + M_PI, 2*M_PI) - M_PI;

//cout << "Tc i " << endl << weights_ (i) * xx * zz.transpose () << endl;
//        Tc += weights_ (i) * xx * zz.transpose ();
//    }
//    //calculate Kalman gain K;
//cout << "Tc" << endl << Tc << endl;
//cout << "S inv" << endl << S.inverse () << endl;
//    MatrixXd K = Tc * S.inverse ();
//cout << "K" << endl << K << endl;
//    //update state mean and covariance matrix
//  cout << "measurement: " << endl;
//  cout << meas_package.raw_measurements_ << endl;

//    VectorXd z_d = meas_package.raw_measurements_ - z_pred;
//    z_d(1) = fmod (z_d(1) + M_PI, 2*M_PI) - M_PI;
//    x_ = x + K * z_d;
//    x_(3) = fmod (x_(3) + M_PI, 2*M_PI) - M_PI;
//  //  x_ += K * (z - z_pred);
//    P_ = P - K * S * K.transpose ();
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {

  MatrixXd S (3, 3);
  S <<
      std_radr_ * std_radr_, 0, 0,
      0, std_radphi_ * std_radphi_, 0,
      0, 0, std_radrd_ * std_radrd_;

  return Update(
        meas_package,
        3,
        S,
        [](const VectorXd& state_vector) -> VectorXd {
cout << "--- m f 1" << endl;
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
cout << "--- m f 2" << endl;
          return z;
        },
        [](const VectorXd& z) -> VectorXd {
          VectorXd zz (3);
          zz = z;
          zz (1) = fmod(zz(1) + M_PI, 2*M_PI) - M_PI;
          return zz;
        }
        );

//  /**
//  TODO:

//  Complete this function! Use radar data to update the belief about the object's
//  position. Modify the state vector, x_, and covariance, P_.

//  You'll also need to calculate the radar NIS.
//  */

////  std_radr_ = 0.3;
////  std_radphi_ = 0.0175;
////  std_radrd_ = 0.1;

////  //create example matrix with predicted sigma points
////  Xsig_pred_ <<
////         5.9374,  6.0640,   5.925,  5.9436,  5.9266,  5.9374,  5.9389,  5.9374,  5.8106,  5.9457,  5.9310,  5.9465,  5.9374,  5.9359,  5.93744,
////           1.48,  1.4436,   1.660,  1.4934,  1.5036,    1.48,  1.4868,    1.48,  1.5271,  1.3104,  1.4787,  1.4674,    1.48,  1.4851,    1.486,
////          2.204,  2.2841,  2.2455,  2.2958,   2.204,   2.204,  2.2395,   2.204,  2.1256,  2.1642,  2.1139,   2.204,   2.204,  2.1702,   2.2049,
////         0.5367, 0.47338, 0.67809, 0.55455, 0.64364, 0.54337,  0.5367, 0.53851, 0.60017, 0.39546, 0.51900, 0.42991, 0.530188,  0.5367, 0.535048,
////          0.352, 0.29997, 0.46212, 0.37633,  0.4841, 0.41872,   0.352, 0.38744, 0.40562, 0.24347, 0.32926,  0.2214, 0.28687,   0.352, 0.318159;

//  VectorXd x (n_x_);
//  x.setZero ();
//  for (int i=0; i<2 * n_aug_ + 1; i++) {
//      x += weights_ (i) * Xsig_pred_.block <5, 1>(0, i);
//  }
//  //predict state covariance matrix
//  MatrixXd P (n_x_, n_x_);
//  P.setZero ();
//  for (int i=0; i<2 * n_aug_ + 1; i++) {
//      VectorXd Ai = Xsig_pred_.block <5, 1>(0, i) - x;
//      P += weights_ (i) * Ai * Ai.transpose ();
//  }
//cout << "predict x: " << endl << x << endl;
//cout << "predict P: " << endl << P << endl;

//  int n_z = 3;

//  //create matrix for sigma points in measurement space
//  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

//  //mean predicted measurement
//  VectorXd z_pred = VectorXd(n_z);

//  //measurement covariance matrix S
//  MatrixXd S = MatrixXd(n_z,n_z);

//  //transform sigma points into measurement space
//  for (int i=0; i<2 * n_aug_ + 1; i++) {
//    float x = Xsig_pred_ (0, i);
//    float y = Xsig_pred_ (1, i);
//    float v = Xsig_pred_ (2, i);
//    float psi = Xsig_pred_ (3, i);
//    Zsig (0, i) = std::sqrt (x * x + y * y);
//    if (Zsig (0, i) > 1e-4) {
//      Zsig (1, i) = atan2(y, x);
//      Zsig (2, i) =
//          (x * cos (psi) * v +
//           y * sin (psi) * v) /
//          Zsig (0, i);
//    } else {
//      Zsig.setZero ();
//    }

//  }
//  //calculate mean predicted measurement
//  z_pred.setZero ();
//  for (int i=0; i<2 * n_aug_ + 1; i++) {
//      z_pred += weights_ (i) * Zsig.block <3, 1> (0, i);
//  }
//  //calculate measurement covariance matrix S
//  S <<
//      std_radr_ * std_radr_, 0, 0,
//      0, std_radphi_ * std_radphi_, 0,
//      0, 0, std_radrd_ * std_radrd_;

//  for (int i=0; i<2 * n_aug_ + 1; i++) {
//      VectorXd ai = Zsig.block <3, 1> (0, i) - z_pred;
//cout << "ai: " << endl << ai << endl;
//      ai (1) = fmod(ai(1) + M_PI, 2*M_PI) - M_PI;
//      S += weights_ (i) * ai * ai.transpose ();
//  }
//cout << "z_pred" << endl << z_pred << endl;
//cout << "S" << endl << S << endl;

////x_ <<
////   5.93637,
////   1.49035,
////   2.20528,
////  0.536853,
////  0.353577;

////P_ <<
////0.0054342,  -0.002405,  0.0034157, -0.0034819, -0.00299378,
////-0.002405,    0.01084,   0.001492,  0.0098018,  0.00791091,
////0.0034157,   0.001492,  0.0058012, 0.00077863, 0.000792973,
////-0.0034819,  0.0098018, 0.00077863,   0.011923,   0.0112491,
////-0.0029937,  0.0079109, 0.00079297,   0.011249,   0.0126972;

////Xsig_pred_ <<
////       5.9374,  6.0640,   5.925,  5.9436,  5.9266,  5.9374,  5.9389,  5.9374,  5.8106,  5.9457,  5.9310,  5.9465,  5.9374,  5.9359,  5.93744,
////         1.48,  1.4436,   1.660,  1.4934,  1.5036,    1.48,  1.4868,    1.48,  1.5271,  1.3104,  1.4787,  1.4674,    1.48,  1.4851,    1.486,
////        2.204,  2.2841,  2.2455,  2.2958,   2.204,   2.204,  2.2395,   2.204,  2.1256,  2.1642,  2.1139,   2.204,   2.204,  2.1702,   2.2049,
////       0.5367, 0.47338, 0.67809, 0.55455, 0.64364, 0.54337,  0.5367, 0.53851, 0.60017, 0.39546, 0.51900, 0.42991, 0.530188,  0.5367, 0.535048,
////        0.352, 0.29997, 0.46212, 0.37633,  0.4841, 0.41872,   0.352, 0.38744, 0.40562, 0.24347, 0.32926,  0.2214, 0.28687,   0.352, 0.318159;

////Zsig <<
////    6.1190,  6.2334,  6.1531,  6.1283,  6.1143,  6.1190,  6.1221,  6.1190,  6.0079,  6.0883,  6.1125,  6.1248,  6.1190,  6.1188,  6.12057,
////   0.24428,  0.2337, 0.27316, 0.24616, 0.24846, 0.24428, 0.24530, 0.24428, 0.25700, 0.21692, 0.24433, 0.24193, 0.24428, 0.24515, 0.245239,
////    2.1104,  2.2188,  2.0639,   2.187,  2.0341,  2.1061,  2.1450,  2.1092,  2.0016,   2.129,  2.0346,  2.1651,  2.1145,  2.0786,  2.11295;

////z_pred <<
////    6.12155,
////   0.245993,
////    2.10313;

////S <<
////    0.0946171, -0.000139448,   0.00407016,
//// -0.000139448,  0.000617548, -0.000770652,
////   0.00407016, -0.000770652,    0.0180917;


//  //create matrix for cross correlation Tc
//  MatrixXd Tc = MatrixXd(n_x_, n_z);

//  //calculate cross correlation matrix
//  Tc.setZero ();
//  for (int i=0; i<2 * n_aug_ + 1; i++) {
//      VectorXd xx = Xsig_pred_.block(0, i, n_x_, 1) - x;
//      VectorXd zz = Zsig.block (0, i, n_z, 1) - z_pred;
//      Tc += weights_ (i) * xx * zz.transpose ();
//  }
//  //calculate Kalman gain K;
//  MatrixXd K = Tc * S.inverse ();
//  //update state mean and covariance matrix
//cout << "measurement: " << endl;
//cout << meas_package.raw_measurements_ << endl;

////VectorXd z = VectorXd(n_z);
////z <<
////    5.9214,   //rho in m
////    0.2187,   //phi in rad
////    2.0062;   //rho_dot in m/s
////cout << z << endl;

//  x_ = x + K * (meas_package.raw_measurements_ - z_pred);
////  x_ += K * (z - z_pred);
//  P_ = P - K * S * K.transpose ();
}
