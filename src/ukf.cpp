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
  use_laser_ = false;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3; //30

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;//30
  
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

  this->n_x_     = 5;	///* State dimension
  this->n_aug_   = 7; ///* Augmented state dimension
  this->lambda_  = 3 - this->n_aug_; ///* Sigma point spreading parameter
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
  Q_ << std_a_*std_a_,					   0,
					0, std_yawdd_*std_yawdd_;

  //Measurement noise matrix R for rader
  R_rader_ = MatrixXd::Zero(3, 3);
  R_rader_ << std_radr_*std_radr_,						 0,						0,
								0, std_radphi_*std_radphi_,                     0,
								0,						 0, std_radrd_*std_radrd_;
	  //Measurement noise matrix R for lidar
  R_lidar_ = MatrixXd::Zero(2, 2);
  R_lidar_ << std_laspx_*std_laspx_,                     0,
	                              0, std_laspy_*std_laspy_;

  //Transfer matrix for Lidar
  H_laser_ = MatrixXd(2, 5);
  H_laser_ << 1, 0, 0, 0, 0,
			  0, 1, 0, 0, 0;


}

UKF::~UKF() {}

void UKF::Normarize_angle(double &theta) {
	while (theta > M_PI) theta -= 2.0*M_PI;
	while (theta <-M_PI) theta += 2.0*M_PI;
	return;
}


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

	cout << "********** Process Start ********** " << std::endl;
	/* calculate delta time [s] ??? */
	static double previous_timestamp_ = 0;

	if (!this->is_initialized_){
		cout << "********** UKF initialization ********** " << std::endl;
		
		if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
			/* Convert radar from polar to cartesian coordinates and initialize state */
			cout << "UKF with first Radar: " << endl;

			//Measurements
			double ro     = meas_package.raw_measurements_[0];
			double theta  = meas_package.raw_measurements_[1];
			double ro_dot = meas_package.raw_measurements_[2];

			double px = ro * cos(theta);
			double py = ro * sin(theta);

			//double vx = ro_dot * cos(theta);
			//double vy = ro_dot * sin(theta);
			//double  v = sqrt(vx*vx + vy+vy);

			//Initialize the state vector
			this->x_(0) = px;
			this->x_(1) = py;

			//complete initialization
			this->is_initialized_ = true;
			//update the time stamp
			previous_timestamp_ = meas_package.timestamp_;
		}
		else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
			/* Initialize state */
			cout << "UKF with first Laser: " << endl;

			//Measurements
			double px = meas_package.raw_measurements_[0];
			double py = meas_package.raw_measurements_[1];			
			//double theta = atan2(py, px);

			//Initialize the state vector
			this->x_(0) = px;
			this->x_(1) = py;

			if (fabs(x_(0) < 0.001 && fabs(x_(1) < 0.001))) {
				this->x_(0) = px;
				this->x_(1) = py;
			}

			//complete initialization
			this->is_initialized_ = true;
			//update the time stamp
			previous_timestamp_ = meas_package.timestamp_;
		}
		else {
			/* Nothing to do */
		}
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
	if (this->use_radar_) {
		if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
			/* Update */
			this->UpdateRadar(meas_package);
		}
	}
	if (this->use_laser_) {
		if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
			/* Update */
			this->UpdateLidar(meas_package);
		}
	}
}


void UKF::make_sigma_point(MatrixXd *Xsig_aug) {
	/* Create  Augumented Sigma Point */
	//cout << "Create the augumented sigma point" << endl;
	VectorXd x_aug = VectorXd::Zero(n_aug_);
	x_aug.fill(0.);
	//create augmented mean state vector
	//cout << "Create the mean state vector" << endl;
	x_aug.head(5) = this->x_;
	//x_aug(5) = 0; //nu_a
	//x_aug(6) = 0; //nu_psi

	//create augmented covariance matrix
	//cout << "Create the covariance matrix" << endl;
	MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);
	P_aug.fill(0.);
	P_aug.topLeftCorner(n_x_, n_x_) = this->P_;
	P_aug(n_x_, n_x_)     = std_a_*std_a_;
	P_aug(n_x_+1, n_x_+1) = std_yawdd_*std_yawdd_;
	//P_aug.bottomRightCorner(2, 2) = this->Q_;

	//calculate squre root of P_aug (state covariance matrix)
	MatrixXd A = P_aug.llt().matrixL();

	//create augmented sigma points
	const double coeff = sqrt(lambda_ + n_aug_);
	
	//cout << "Create augmented sigma points" << endl;
	Xsig_aug->col(0) = x_aug;
	for (int i = 0; i < n_aug_; i++) {
		Xsig_aug->col(i + 1) = x_aug + coeff*A.col(i);
		Xsig_aug->col(i + n_aug_ + 1) = x_aug - coeff*A.col(i);
	}
	return;
}

void UKF::predict_sigma_point(const MatrixXd *Xsig_aug_in, double delta_t){
	/* Predict the Sigma Point */
	double delta_t2 = delta_t*delta_t;

	MatrixXd Xsig_aug = *Xsig_aug_in;

	for (int i = 0; i < n_sigma_; i++) {

		//State
		double px  = Xsig_aug(0, i);
		double py  = Xsig_aug(1, i);
		double v   = Xsig_aug(2, i);
		double yaw = Xsig_aug(3, i);
		double yaw_dot = Xsig_aug(4, i);
		//Noise
		double nu_a      = Xsig_aug(5, i);
		double nu_yaw_dd = Xsig_aug(6, i);
		//predicted state value
		double px_p, py_p;
		double coeff = v / yaw_dot;

		//state transfered
		if (fabs(yaw_dot) > 0.001) {
			px_p = px + coeff*( sin(yaw + yaw_dot*delta_t) - sin(yaw));
			py_p = py + coeff*(-cos(yaw + yaw_dot*delta_t) + cos(yaw));
		}
		else {
			px_p = px + v*delta_t*cos(yaw);
			py_p = py + v*delta_t*sin(yaw);
		}
		double v_p       = v;
		double yaw_p     = yaw + yaw_dot*delta_t;
		double yaw_dot_p = yaw_dot;

		//noise transfered
		px_p      += 0.5 * nu_a * cos(yaw) * delta_t2;
		py_p      += 0.5 * nu_a * sin(yaw) * delta_t2;
		v_p       += nu_a * delta_t;
		yaw_p     += 0.5 * nu_yaw_dd * delta_t2;
		yaw_dot_p += nu_yaw_dd * delta_t;

		this->Xsig_pred_(0, i) = px_p;
		this->Xsig_pred_(1, i) = py_p;
		this->Xsig_pred_(2, i) = v_p;
		this->Xsig_pred_(3, i) = yaw_p;
		this->Xsig_pred_(4, i) = yaw_dot_p;
	}
	return;
}

void UKF::predict_state_vector(void) {
	/* Predict the State Vector  */
	x_.fill(0);
	for (int i = 0; i < n_sigma_; i++) {
		x_ += weights_(i)*Xsig_pred_.col(i);
	}

	/* Predict the State Covariance Matrix */
	P_.fill(0);
	for (int i = 0; i < n_sigma_; i++) {
		VectorXd diff_X = Xsig_pred_.col(i) - x_;
		// Normarize angle[-pi to pi]
		Normarize_angle(diff_X(3));

		P_ += weights_(i) * diff_X * diff_X.transpose();
	}
	/* debug */
	cout << "Predict State Vector X: " << endl << x_ << endl << endl;
	cout << "Predict State Vector P: " << endl << P_ << endl <<endl;
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
	MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
	Xsig_aug.fill(0.);
	this->make_sigma_point(&Xsig_aug);

	//cout << "********** Predict Sigma Point ********** " << std::endl;
	this->predict_sigma_point(&Xsig_aug, delta_t);
	//std::cout << "********** Predict State Vector ********** " << std::endl;
	this->predict_state_vector();
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
	const int n_z_lidar = 2;

	VectorXd z_pred = H_laser_ * x_;            // [2x1] = [2x5] * [5x1]

	//Measurement
	VectorXd z = VectorXd::Zero(n_z_lidar);
	z << meas_package.raw_measurements_;

	const MatrixXd Ht = H_laser_.transpose();
	MatrixXd PHt = P_ * Ht; //[5x2] = [5x5]*[5x2]

	//Predict measure vector
	VectorXd z_diff = z - z_pred;		        // [2x1] = [2x1] - [2x1]
	MatrixXd S      = H_laser_* PHt + R_lidar_; // [2x2] = [2x5]*[5x2]*[2x2]
	MatrixXd Si     = S.inverse();			    // [2x2]
	MatrixXd K      = PHt * Si;                 //[5x2] = [5x2]*[2x2]

	//Update the posterior state
	x_ = x_ + (K * z_diff); // [5x1] = [5x1] + [5x2]*[2x1]
	
	const long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_laser_) * P_; // [5x5] = [5x2]*[2x5] * [5x5]

	//Calculate NIS for Lidar
	NIS_laser = z_diff.transpose()*Si*z_diff;
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

	/* Predict Lidar measurement value */
	//cout << "********** Predict Rader Measurements ********** " << std::endl;
	//transform sigma points into measurement space
	//create matrix for sigma points in measurement space
	MatrixXd Zsig = MatrixXd(n_z, n_sigma_);

	for (int i = 0; i < n_sigma_; i++) {
		double px  = Xsig_pred_(0, i);
		double py  = Xsig_pred_(1, i);
		double v   = Xsig_pred_(2, i);
		double psi = Xsig_pred_(3, i);

		double vx = v*cos(psi);
		double vy = v*sin(psi);

		Zsig(0, i) = sqrt(px*px + py*py);//ro
		Zsig(1, i) = atan2(py, px);      //psi[-pi to +pi]
		Zsig(2, i) = (px*vx + py*vy) / sqrt(px*px + py*py);//ro_dot
		//check zero divide
		if (Zsig(0, i) < 0.001){
			Zsig(2, i) = 0;
		}
	}


	/*Predicted Mean state vector Z_hat */
	//mean predicted measurement
	VectorXd z_pred = VectorXd(n_z);
	for (int i = 0; i < this->n_sigma_; i++) {
		z_pred = z_pred + weights_(i) * Zsig.col(i);
	}

	/* Predicted innovation covariance matrix S */
	//measurement covariance matrix S
	MatrixXd S = MatrixXd(n_z, n_z);
	S.fill(0);

	for (int i = 0; i < n_sigma_; i++) {
		VectorXd diff_Z = Zsig.col(i) - z_pred;
		// Normarize angle[-pi to pi]
		Normarize_angle(diff_Z(1));
		S = S + weights_(i) * (diff_Z * diff_Z.transpose());
	}
	//add noise to S covariance
	S = S + this->R_rader_;

	//cout << "********** Update state vector on Rader ********** " << std::endl;
	/* Update the state vector */
	//create matrix for cross correlation Tc
	MatrixXd Tc = MatrixXd(n_x_, n_z);
	//calculate cross correlation matrix
	Tc.fill(0);

	//std::cout << "calculate cross correlation matrix" << std::endl;
	for (int i = 0; i < n_sigma_; i++) {  //2n+1 simga points

		//residual
		VectorXd z_diff = Zsig.col(i) - z_pred;

		// Normarize angle[-pi to pi]
		Normarize_angle(z_diff(1));

		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - this->x_;

		// Normarize angle[-pi to pi]
		Normarize_angle(x_diff(3));

		Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
	}

	//std::cout << "update kalman filter" << std::endl;
	//Measurement
	VectorXd z = VectorXd(n_z);
	z << meas_package.raw_measurements_; //ro, theta, ro_dot

	MatrixXd Si = S.inverse();
	//Kalman gain K;
	MatrixXd K = Tc * Si;
	//residual
	VectorXd z_diff = z - z_pred;

	// Normarize angle[-pi to pi]
	Normarize_angle(z_diff(1));

	//update state mean and covariance matrix
	this->x_ = x_ + K * z_diff;
	this->P_ = P_ - K * S * K.transpose();

	//Calculate NIS for Rader
	NIS_radar = z_diff.transpose() * Si * z_diff;

	/* debug */
	cout << "z_pred_radar " <<endl << z_pred << endl << endl;
	cout << "S            " <<endl << S << endl << endl;
	cout << "Tc           " <<endl << Tc << endl << endl;
}
