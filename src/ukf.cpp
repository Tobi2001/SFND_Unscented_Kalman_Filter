#include <chrono>
#include <iostream>

#include "ukf.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF()
{
    // if this is false, laser measurements will be ignored (except during init)
    use_laser_ = true;

    // if this is false, radar measurements will be ignored (except during init)
    use_radar_ = true;

    // initial state vector
    x_ = VectorXd(5);

    // initial covariance matrix
    P_ = MatrixXd(5, 5);
    P_.setIdentity();


    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = 1.0;

    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = 0.6;

    /**
     * DO NOT MODIFY measurement noise values below.
     * These are provided by the sensor manufacturer.
     */

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

    /**
     * End DO NOT MODIFY section for measurement noise values
     */

    is_initialized_ = false;

    n_x_ = 5;

    n_aug_ = 7;

    n_z_radar_ = 3;

    n_z_lidar_ = 2;

    Xsig_pred_ = Eigen::MatrixXd(n_x_, 2 * n_aug_ + 1);
    Xsig_pred_.setZero();

    time_us_ = 0;

    lambda_ = 3 - n_x_;

    weights_ = VectorXd(2 * n_aug_ + 1);
    weights_(0) = lambda_ / (lambda_ + n_aug_);
    weights_.tail(2 * n_aug_).fill(0.5 / (n_aug_ + lambda_));


    R_radar_ = Eigen::MatrixXd(n_z_radar_, n_z_radar_);
    R_radar_ << std_radr_ * std_radr_, 0, 0,
        0, std_radphi_* std_radphi_, 0,
        0, 0, std_radrd_* std_radrd_;

    R_lidar_ = Eigen::MatrixXd(n_z_lidar_, n_z_lidar_);
    R_lidar_ << std_laspx_ * std_laspx_, 0,
        0, std_laspy_* std_laspy_;
}

UKF::~UKF() {}


void UKF::ProcessMeasurement(MeasurementPackage meas_package)
{
    if (!is_initialized_)
    {
        initState(meas_package);
        return;
    }

    double delta_t = microsecsToSecs(meas_package.timestamp_ - time_us_);
    if (meas_package.sensor_type_ == MeasurementPackage::SensorType::LASER ||
        meas_package.sensor_type_ == MeasurementPackage::SensorType::RADAR)
    {
        time_us_ = meas_package.timestamp_;
        Prediction(delta_t);

        if (meas_package.sensor_type_ == MeasurementPackage::SensorType::LASER)
        {
            UpdateLidar(meas_package);
        }
        else if (meas_package.sensor_type_ == MeasurementPackage::SensorType::RADAR)
        {
            UpdateRadar(meas_package);
        }
    }
    else
    {
        // given package is invalid
    }
}


void UKF::Prediction(double delta_t)
{
    predictSigmaPoints(generateSigmaPoints(), delta_t);
    predictState();
}


void UKF::UpdateLidar(MeasurementPackage meas_package)
{
    Eigen::MatrixXd Zsig = sigmaToLidarSpace();
    VectorXd z_pred = predictZ(Zsig, n_z_lidar_);
    MatrixXd S = genS(Zsig, z_pred, n_z_lidar_);
    S += R_lidar_;

    updateState(meas_package, Zsig, z_pred, S, n_z_lidar_);
}


void UKF::UpdateRadar(MeasurementPackage meas_package)
{
    Eigen::MatrixXd Zsig = sigmaToRadarSpace();
    VectorXd z_pred = predictZ(Zsig, n_z_radar_);
    MatrixXd S = genS(Zsig, z_pred, n_z_radar_, true);
    S += R_radar_;

    updateState(meas_package, Zsig, z_pred, S, n_z_radar_, true);
}


// Private methods

void UKF::initState(const MeasurementPackage& package)
{
    if (package.sensor_type_ == MeasurementPackage::SensorType::LASER)
    {
        x_ = stateFromLidar(package);
    }
    else if (package.sensor_type_ == MeasurementPackage::SensorType::RADAR)
    {
        x_ = stateFromRadar(package);
    }
    else
    {
        // unknown sensor type
        return;
    }

    time_us_ = package.timestamp_;
    is_initialized_ = true;
}

Eigen::VectorXd UKF::stateFromLidar(const MeasurementPackage& package) const
{
    Eigen::VectorXd state(5);
    state << package.raw_measurements_[0], package.raw_measurements_[1], 0, 0, 0;
    return state;
}

Eigen::VectorXd UKF::stateFromRadar(const MeasurementPackage& package) const
{
    Eigen::VectorXd state(3);
    double rho = package.raw_measurements_[0];
    double phi = package.raw_measurements_[1];
    double rho_dot = package.raw_measurements_[2];
    double vx = rho_dot * cos(phi);
    double vy = rho_dot * sin(phi);
    double v = sqrt(vx * vx + vy * vy);
    state << rho * cos(phi), rho* sin(phi), v, 0, 0;
    return state;
}

double UKF::microsecsToSecs(long long micros) const
{
    std::chrono::microseconds us(micros);
    auto secs = std::chrono::duration_cast<std::chrono::duration<double, std::chrono::seconds::period>>(us);
    return secs.count();
}

Eigen::MatrixXd UKF::generateSigmaPoints() const
{
    Eigen::VectorXd x_aug = VectorXd(n_aug_);
    x_aug.head(n_x_) = x_;
    x_aug.tail(n_aug_ - n_x_).setZero();
    Eigen::MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
    P_aug.setZero();
    P_aug.topLeftCorner(n_x_, n_x_) = P_;
    P_aug(n_x_, n_x_) = std_a_ * std_a_;
    P_aug(n_x_ + 1, n_x_ + 1) = std_yawdd_ * std_yawdd_;

    Eigen::MatrixXd result(n_aug_, 2 * n_aug_ + 1);
    Eigen::MatrixXd root = P_aug.llt().matrixL();

    result.col(0) = x_aug;
    double c = std::sqrt(lambda_ + n_aug_);
    for (int i = 0; i < 2 * n_aug_; ++i)
    {
        int sign = (i < n_aug_ ? 1 : -1);
        result.col(i + 1) = x_aug + sign * c * root.col(i % n_aug_);
    }
    return result;
}

void UKF::predictSigmaPoints(const Eigen::MatrixXd& sigmaPoints, const double delta_t)
{
    for (int i = 0; i < 2 * n_aug_ + 1; ++i)
    {
        Eigen::VectorXd upd(n_x_);
        upd << 0,
            0,
            0,
            sigmaPoints.col(i)(4) * delta_t,
            0;

        Eigen::VectorXd err(n_x_);
        err << 0.5 * delta_t * delta_t * cos(sigmaPoints.col(i)(3)) * sigmaPoints.col(i)(5),
            0.5 * delta_t * delta_t * sin(sigmaPoints.col(i)(3)) * sigmaPoints.col(i)(5),
            delta_t* sigmaPoints.col(i)(5),
            0.5 * delta_t * delta_t * sigmaPoints.col(i)(6),
            delta_t* sigmaPoints.col(i)(6);

        if (std::fabs(sigmaPoints.col(i)(4)) >= std::numeric_limits<double>::epsilon())
        {
            double v = sigmaPoints.col(i)(2) / sigmaPoints.col(i)(4);
            upd(0) = v * (sin(sigmaPoints.col(i)(3) + sigmaPoints.col(i)(4) * delta_t) - sin(sigmaPoints.col(i)(3)));
            upd(1) = v * (-cos(sigmaPoints.col(i)(3) + sigmaPoints.col(i)(4) * delta_t) + cos(sigmaPoints.col(i)(3)));
        }
        else
        {
            upd(0) = sigmaPoints.col(i)(2) * cos(sigmaPoints.col(i)(3)) * delta_t;
            upd(1) = sigmaPoints.col(i)(2) * sin(sigmaPoints.col(i)(3)) * delta_t;
        }
        Xsig_pred_.col(i) = sigmaPoints.col(i).head(n_x_) + upd + err;
    }
}

void UKF::normalizeAngle(Eigen::VectorXd& input, const int row) const
{
    while (input(row) > M_PI)
    {
        input(row) -= 2. * M_PI;
    }
    while (input(row) < -M_PI)
    {
        input(row) += 2. * M_PI;
    }
}

void UKF::predictState()
{
    x_.setZero();
    for (int i = 0; i < 2 * n_aug_ + 1; ++i)
    {
        x_ += weights_(i) * Xsig_pred_.col(i);
    }
    P_.setZero();
    for (int i = 0; i < 2 * n_aug_ + 1; ++i)
    {
        Eigen::VectorXd diff = Xsig_pred_.col(i) - x_;
        normalizeAngle(diff, 3);
        P_ += weights_(i) * (diff * diff.transpose());
    }
}

Eigen::MatrixXd UKF::sigmaToRadarSpace() const
{
    Eigen::MatrixXd Zsig = MatrixXd(n_z_radar_, 2 * n_aug_ + 1);
    for (int i = 0; i < 2 * n_aug_ + 1; ++i)
    {
        double px = Xsig_pred_.col(i)(0);
        double py = Xsig_pred_.col(i)(1);
        double v = Xsig_pred_.col(i)(2);
        double phi = Xsig_pred_.col(i)(3);
        double denom = std::sqrt(px * px + py * py);
        Zsig.col(i)(0) = denom;
        Zsig.col(i)(1) = 0.0;
        if (std::fabs(px) >= std::numeric_limits<double>::epsilon())
        {
            Zsig.col(i)(1) = std::atan2(py, px);
        }
        Zsig.col(i)(2) = 0.0;
        if (std::fabs(denom) >= std::numeric_limits<double>::epsilon())
        {
            Zsig.col(i)(2) = (px * cos(phi) * v + py * sin(phi) * v) / denom;
        }
    }
    return Zsig;
}

Eigen::MatrixXd UKF::sigmaToLidarSpace() const
{
    Eigen::MatrixXd Zsig = MatrixXd(n_z_lidar_, 2 * n_aug_ + 1);
    for (int i = 0; i < 2 * n_aug_ + 1; ++i)
    {
        Zsig.col(i) = Xsig_pred_.col(i).head(n_z_lidar_);
    }
    return Zsig;
}

Eigen::VectorXd UKF::predictZ(const Eigen::MatrixXd& Zsig, const int dim) const
{
    VectorXd z_pred = VectorXd(dim);
    z_pred.setZero();
    for (int i = 0; i < 2 * n_aug_ + 1; ++i)
    {
        z_pred += weights_(i) * Zsig.col(i);
    }
    return z_pred;
}

Eigen::MatrixXd UKF::genS(
    const Eigen::MatrixXd& Zsig, const Eigen::VectorXd& z_pred, const int dim, const bool norm_radar) const
{
    MatrixXd S = MatrixXd(dim, dim);
    S.fill(0);

    for (int i = 0; i < 2 * n_aug_ + 1; ++i)
    {
        VectorXd diff = Zsig.col(i) - z_pred;
        if (norm_radar)
        {
            normalizeAngle(diff, 1);
        }
        S += weights_(i) * diff * diff.transpose();
    }
    return S;
}

void UKF::updateState(
    const MeasurementPackage& package, const Eigen::MatrixXd& Zsig, const Eigen::VectorXd& z_pred,
    const Eigen::MatrixXd& S, const int dim, const bool norm_radar)
{
    MatrixXd Tc = MatrixXd(n_x_, dim);
    Tc.setZero();

    for (int i = 0; i < n_aug_ * 2 + 1; ++i)
    {
        VectorXd diffX = Xsig_pred_.col(i) - x_;
        normalizeAngle(diffX, 3);
        VectorXd diffZ = Zsig.col(i) - z_pred;
        if (norm_radar)
        {
            normalizeAngle(diffZ, 1);
        }
        Tc += weights_(i) * diffX * diffZ.transpose();
    }

    MatrixXd K = Tc * S.inverse();

    VectorXd z_diff = package.raw_measurements_ - z_pred;
    if (norm_radar)
    {
        normalizeAngle(z_diff, 1);
    }
    x_ = x_ + K * z_diff;
    P_ = P_ - K * S * K.transpose();
}
