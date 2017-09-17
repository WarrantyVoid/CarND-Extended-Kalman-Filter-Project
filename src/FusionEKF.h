#ifndef FusionEKF_H_
#define FusionEKF_H_

#include "measurement_package.h"
#include "kalman_filter.h"

class FusionEKF
{
public:
  /**
  * Constructor.
  **/
  FusionEKF();

  /**
  * Destructor.
  **/
  virtual ~FusionEKF();

  /**
  * Run the whole flow of the Kalman Filter from here.
  **/
  void ProcessMeasurement(const MeasurementPackage &measurementPack);

  /**
  * Retrieves current estimate.
  **/
  TVector GetEstimate() const;

private:
  /** Kalman Filter update and prediction **/
  KalmanFilter mEkf;

  /** True, if the tracking toolbox was initialized, false otherwise (first measurement) **/
  bool mIsInitialized;

  /** Previous timestamp **/
  TTimeStamp mPreviousTimestamp;

  /** Sensor noise convariance matrix (laser) **/
  Eigen::MatrixXd mRLaser;

  /** Sensor noise convariance matrix (radar) **/
  Eigen::MatrixXd mRRadar;

  /** Measurement translation matrix (laser) **/
  Eigen::MatrixXd mHLaser;

  /** Measurement translation matrix (radar) **/
  Eigen::MatrixXd mHjRadar;
};

#endif /* FusionEKF_H_ */
