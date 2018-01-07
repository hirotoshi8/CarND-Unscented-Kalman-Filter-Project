#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
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
  std_a_ = 3; //30

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 8;//30
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  
  this->is_initialized_ = false;

  this->n_x_ = 5;	///* State dimension
  this->n_aug_ = 7; ///* Augmented state dimension
  this->lambda_ = 3 - this->n_x_; ///* Sigma point spreading parameter
  this->n_sigma_ = 2 * n_aug_ + 1;

  /* State Vector */
  x_ = VectorXd::Zero(n_x_);
  
  /* State Covariance Matrix */
  P_ = MatrixXd::Identity(n_x_, n_x_);
  
  /* Sigma Point Matrix */
  Xsig_pred_ = MatrixXd::Zero(n_x_, n_sigma_);

  //Weights for Sigma point convert
  weights_ = VectorXd(n_sigma_);
  weights_(0) = lambda_ / (lambda_ + n_aug_);

  for (int i = 1; i < n_sigma_; i++) {
	  weights_(i) = 0.5 / (lambda_ + n_aug_);
  }

  //Process noise covariance matrix
  Q_ = MatrixXd::Zero(2, 2);
  Q_(0, 0) = std_a_*std_a_;
  Q_(1, 1) = std_yawdd_*std_yawdd_;

  //Measurement noise matrix R for rader
  R_rader_ = MatrixXd::Zero(3, 3);
  R_rader_(0, 0) = std_radr_*std_radr_;
  R_rader_(1, 1) = std_radphi_*std_radphi_;
  R_rader_(2, 2) = std_radrd_*std_radrd_;

  //Measurement noise matrix R for lidar
  R_lidar_ = MatrixXd::Zero(2, 2);
  R_lidar_(0, 0) = std_laspx_*std_laspx_;
  R_lidar_(1, 1) = std_laspy_*std_laspy_;

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
	std::cout << "********** Process Start ********** " << std::endl;
	/* calculate delta time [s] ??? */
	static double previous_timestamp_ = 0;

	if (!this->is_initialized_){
		std::cout << "********** UKF initialization ********** " << std::endl;
		
		if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
			/* Convert radar from polar to cartesian coordinates and initialize state */
			cout << "EKF with first Radar: " << endl;

			//Measurements
			double ro = meas_package.raw_measurements_[0];
			double theta = meas_package.raw_measurements_[1];
			double ro_dot = meas_package.raw_measurements_[2];

			double px = ro * cos(theta);
			double py = ro * sin(theta);

			//Initialize the state vector
			this->x_(0) = px;
			this->x_(1) = py;
		}
		else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
			/* Initialize state */
			cout << "EKF with first Laser: " << endl;

			//Measurements
			double px = meas_package.raw_measurements_[0];
			double py = meas_package.raw_measurements_[1];
			
			//Initialize the state vector
			this->x_(0) = px;
			this->x_(1) = py;
		}
		else {
			/* Nothing to do */
		}

		//complete initialization
		this->is_initialized_ = true;
		//update the time stamp
		previous_timestamp_ = meas_package.timestamp_;
		//debug
		//cout << "Complete the initialization"<<endl;
		return;
	}
	/* Calculate delta time [s]*/
	double delta_t = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;//dt - expressed in seconds
	previous_timestamp_ = meas_package.timestamp_;
	//cout << "delta_t " << delta_t << endl;

	/* Prediction */
	//debug
	//cout << "---------- Start Prediction ----------" << endl;
	this->Prediction(delta_t);

	//debug
	//cout << "---------- Start Update ----------" << endl;
	if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
		/* Update */
		if (this->use_radar_) {
			this->UpdateRadar(meas_package);
		}
	}
	else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
		/* Update */
		if (this->use_laser_) {
			this->UpdateLidar(meas_package);
		}
	}
	else {
		/* Nothing to do*/
	}
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

	//cout << "********** Make augumented Sigma Point ********** " << endl;
	
	/* Create  Augumented Sigma Point */
	//cout << "Create the augumented sigma point" << endl;
	VectorXd x_aug = VectorXd::Zero(n_aug_);
	MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);
	MatrixXd Xsig_aug = MatrixXd::Zero(n_aug_, n_sigma_);

	//create augmented mean state vector
	//cout << "Create the mean state vector" << endl;
	x_aug.head(5) = this->x_;
	x_aug(5) = 0; //nu_a
	x_aug(6) = 0; //nu_psi

	//create augmented covariance matrix
	//cout << "Create the covariance matrix" << endl;
	P_aug.topLeftCorner(n_x_, n_x_) = this->P_;
	P_aug.bottomRightCorner(2, 2) = this->Q_;

	//calculate squre root of P_aug (state covariance matrix)
	MatrixXd A = P_aug.llt().matrixL();

	//create augmented sigma points
	const double k = sqrt(lambda_ + n_x_);
	Xsig_aug.col(0) = x_aug;
	//cout << "Create augmented sigma points" << endl;
	for (int i = 0; i < n_aug_; i++) {
		Xsig_aug.col(i + 1) = x_aug + k*A.col(i);
		Xsig_aug.col(i + n_aug_ + 1) = x_aug - k*A.col(i);
	}


	//std::cout << "********** Predict Sigma Point ********** " << std::endl;
	/* Predict the Sigma Point */
	const double delta_t2 = delta_t*delta_t;

	for (int i = 0; i < n_sigma_; i++) {

		VectorXd x_mean = VectorXd(n_x_);
		VectorXd f = VectorXd(n_x_);
		VectorXd noise = VectorXd(n_x_);

		//State
		double px  = Xsig_aug(0, i);
		double py  = Xsig_aug(1, i);
		double v   = Xsig_aug(2, i);
		double psi = Xsig_aug(3, i);
		double psi_dot = Xsig_aug(4, i);
		//Noise
		double nu_a = Xsig_aug(5, i);
		double nu_psi_dd = Xsig_aug(6, i);

		//mean state
		x_mean(0) = px;
		x_mean(1) = py;
		x_mean(2) = v;
		x_mean(3) = psi;
		x_mean(4) = psi_dot;

		//state transfered
		if (fabs(psi_dot) > 0.001) {
			double k = v / psi_dot;
			f(0) = k*(sin(psi + psi_dot*delta_t) - sin(psi));
			f(1) = k*(-cos(psi + psi_dot*delta_t) + cos(psi));
		}
		else {
			f(0) = v*delta_t*cos(psi);
			f(1) = v*delta_t*sin(psi);
		}
		f(2) = 0;
		f(3) = psi_dot*delta_t;
		f(4) = 0;
		//noise transfered
		noise(0) = 0.5*delta_t2 * cos(psi) * nu_a;
		noise(1) = 0.5*delta_t2 * sin(psi) * nu_a;
		noise(2) = delta_t*nu_a;
		noise(3) = 0.5*delta_t2 * nu_psi_dd;
		noise(4) = delta_t * nu_psi_dd;

		this->Xsig_pred_.col(i) = x_mean + f + noise;
	}

	//std::cout << "********** Predict State Vector ********** " << std::endl;
	/* Predict the State Vector  */
	x_.fill(0);
	for (int i = 0; i < n_sigma_; i++) {
		x_ += weights_(i)*Xsig_pred_.col(i);
	}

	/* Predict the State Covariance Matrix */
	P_.fill(0);
	for (int i = 0; i < n_sigma_; i++) {
		VectorXd diff_X = Xsig_pred_.col(i) - x_;

		while (diff_X(3) > M_PI) diff_X(3) -= 2.0*M_PI;
		while (diff_X(3) <-M_PI) diff_X(3) += 2.0*M_PI;

		P_ += weights_(i)*diff_X*diff_X.transpose();
	}
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
	//Lidar dimension
	const int n_z_lidar = 5;

	//measurement
	VectorXd z = VectorXd::Zero(2);
	z(0) = meas_package.raw_measurements_(0);
	z(1) = meas_package.raw_measurements_(1);

	//Transfer matrix
	MatrixXd H_laser = MatrixXd(2,5);
	H_laser << 1, 0, 0, 0, 0,
		       0, 1, 0, 0, 0;

	MatrixXd Ht = H_laser.transpose();
	MatrixXd PHt = P_ * Ht; //[5x2] = [5x5]*[5x2]

	//Predict measure vector
	VectorXd z_pred = H_laser * x_; // [2x1] = [2x5] * [5x1]
	VectorXd y = z - z_pred;		// [2x1] = [2x1] - [2x1]
	MatrixXd S = H_laser* PHt + R_lidar_; // [2x2] = [2x5]*[5x2]*[2x2]
	MatrixXd Si = S.inverse();			  // [2x2]
	MatrixXd K = PHt * Si; //[5x2] = [5x2]*[2x2]

	//Update the posterior state
	x_ = x_ + (K * y); // [5x1] = [5x1] + [5x2]*[2x1]
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_laser) * P_; // [5x5] = [5x2]*[2x5] * [5x5]

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
	const int n_z = 3;
	//create matrix for sigma points in measurement space
	MatrixXd Zsig = MatrixXd(n_z, n_sigma_);

	//mean predicted measurement
	VectorXd z_pred = VectorXd(n_z);

	/* Predict Lidar measurement value */
	//std::cout << "********** Predict Rader Measurements ********** " << std::endl;
	//transform sigma points into measurement space
	for (int i = 0; i < n_sigma_; i++) {
		double px = Xsig_pred_(0, i);
		double py = Xsig_pred_(1, i);
		double v = Xsig_pred_(2, i);
		double psi = Xsig_pred_(3, i);
		double psi_dot = Xsig_pred_(4, i);

		double vx = v*cos(psi);
		double vy = v*sin(psi);

		double ro = sqrt(px*px + py*py);
		double theta = atan2(py, px); //-pi to +pi
		double theta_dot;
		//zero division check
		if (ro == 0) {
			theta_dot = 0;
		}
		else {
			theta_dot = (px*vx + py*vy) / ro;
		}

		Zsig(0, i) = ro;
		Zsig(1, i) = theta;
		Zsig(2, i) = theta_dot;
	}

	/*Predicted Mean state vector Z_hat */
	for (int i = 0; i < n_sigma_; i++) {
		z_pred += weights_(i) * Zsig.col(i);
	}

	/* Predicted innovation covariance matrix S */
	//measurement covariance matrix S
	MatrixXd S = MatrixXd(n_z, n_z);
	S.fill(0);
	for (int i = 0; i < n_sigma_; i++) {
		VectorXd diff_Z = Zsig.col(i) - z_pred;
		//angle normalization
		while (diff_Z(1)> M_PI) diff_Z(1) -= 2.*M_PI;
		while (diff_Z(1)<-M_PI) diff_Z(1) += 2.*M_PI;
		S = S + weights_(i) * (diff_Z * diff_Z.transpose());
	}
	//add noise to S covariance
	S = S + this->R_rader_;

	//std::cout << "********** Update state vector on Rader ********** " << std::endl;
	/* Update the state vector */
	//Measurement
	VectorXd z = VectorXd(n_z);
	z(0) = meas_package.raw_measurements_(0); //ro
	z(1) = meas_package.raw_measurements_(1); //theta
	z(2) = meas_package.raw_measurements_(2); //ro_dot

	//create matrix for cross correlation Tc
	MatrixXd Tc = MatrixXd(n_x_, n_z);

	//calculate cross correlation matrix
	Tc.fill(0);
	//std::cout << "calculate cross correlation matrix" << std::endl;
	for (int i = 0; i < n_sigma_; i++) {  //2n+1 simga points

		//residual
		VectorXd z_diff = Zsig.col(i) - z_pred;
		//angle normalization
		while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
		while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;

		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - this->x_;
		//angle normalization
		while (x_diff(3)> M_PI) x_diff(3) -= 2.*M_PI;
		while (x_diff(3)<-M_PI) x_diff(3) += 2.*M_PI;

		Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
	}

	//std::cout << "update kalman filter" << std::endl;
	//Kalman gain K;
	MatrixXd K = Tc * S.inverse();
	//residual
	VectorXd z_diff = z - z_pred;

	//angle normalization
	while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
	while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;

	//update state mean and covariance matrix
	this->x_ = x_ + K * z_diff;
	this->P_ = P_  -K * S * K.transpose();

	/* Calculate the NIS */


	/* debug */
	cout << "z_pred " <<endl << z_pred << endl << endl;
	cout << "S      " <<endl << S << endl << endl;
	cout << "Tc     " <<endl << Tc << endl << endl;
}
