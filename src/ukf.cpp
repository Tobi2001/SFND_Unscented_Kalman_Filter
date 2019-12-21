#include <chrono>

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

    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = 30;

    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = 30;

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

    /**
     * TODO: Complete the initialization. See ukf.h for other member properties.
     * Hint: one or more values initialized above might be wildly off...
     */
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package)
{
    if (!is_initialized_)
    {
        try
        {
            initState(meas_package);
        }
        catch (const std::runtime_error&)
        {
            throw;
        }
        return;
    }

    double delta_t = microsecsToSecs(time_us_ - meas_package.timestamp_);
    if (delta_t >= 0 &&
        (meas_package.sensor_type_ == MeasurementPackage::SensorType::LASER ||
        meas_package.sensor_type_ == MeasurementPackage::SensorType::RADAR))
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
    /**
     * TODO: Complete this function! Estimate the object's location.
     * Modify the state vector, x_. Predict sigma points, the state,
     * and the state covariance matrix.
     */
}

void UKF::UpdateLidar(MeasurementPackage meas_package)
{
    /**
     * TODO: Complete this function! Use lidar data to update the belief
     * about the object's position. Modify the state vector, x_, and
     * covariance, P_.
     * You can also calculate the lidar NIS, if desired.
     */
}

void UKF::UpdateRadar(MeasurementPackage meas_package)
{
    /**
     * TODO: Complete this function! Use radar data to update the belief
     * about the object's position. Modify the state vector, x_, and
     * covariance, P_.
     * You can also calculate the radar NIS, if desired.
     */
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
        throw std::runtime_error("Can't process measurement package: Sensor type is not supported");
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
