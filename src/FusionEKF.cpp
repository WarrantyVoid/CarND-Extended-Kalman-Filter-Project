#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF()
  : mEkf()
  , mIsInitialized(false)
  , mPreviousTimestamp(0ull)
  , mRLaser(2, 2)
  , mRRadar(3, 3)
  , mHLaser(2, 4)
  , mHjRadar(3, 4)
{
  //measurement covariance matrix - laser
  mRLaser << 0.0225, 0,
             0     , 0.0225;

  //measurement covariance matrix - radar
  mRRadar << 0.09, 0     , 0,
             0   , 0.0009, 0,
             0   , 0     , 0.09;

  //measurement translation matrix - laser
  mHLaser << 1, 0, 0, 0,
             0, 1, 0, 0;
}

FusionEKF::~FusionEKF()
{

}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurementPack)
{
  if (!mIsInitialized)
  {
    /*****************************************************************************
     *  Initialization
     ****************************************************************************/
    // Uncertainty and movement matrices are independent from sensor
    TMatrix P(4, 4);
    P << 1, 0, 0   , 0,
         0, 1, 0   , 0,
         0, 0, 1000, 0,
         0, 0, 0   , 1000;
    TMatrix F(4, 4);
    F << 1, 0, 1, 0,
         0, 1, 0, 1,
         0, 0, 1, 0,
         0, 0, 0, 1;
    TMatrix Q(4, 4);

    switch (measurementPack.sensorType)
    {
      case MeasurementPackage::RADAR:
      {
        // Convert x and y and v from polar
        TVector x(4);
        x(0) = measurementPack.rawMeasurements(0) * cos(measurementPack.rawMeasurements(1));
        x(1) = measurementPack.rawMeasurements(0) * sin(measurementPack.rawMeasurements(1));
        x(2) = measurementPack.rawMeasurements(2) * cos(measurementPack.rawMeasurements(1));
        x(3) = measurementPack.rawMeasurements(2) * sin(measurementPack.rawMeasurements(1));
        mEkf.Initialize(x, P, F);
      }
      break;
    case MeasurementPackage::LASER:
      {
        // Just use x and y directly
        TVector x(4);
        x(0) = measurementPack.rawMeasurements(0);
        x(1) = measurementPack.rawMeasurements(1);
        mEkf.Initialize(x, P, F);
      }
      break;
    }

    // done initializing
    mIsInitialized = true;
  }
  else
  {
    /*****************************************************************************
     *  Prediction
     ****************************************************************************/
    // predict state and new uncertainty given time delta and acceleration noise
    float dt = (measurementPack.timestamp - mPreviousTimestamp) / 1000000.0f;
    mEkf.Predict(dt, 9.0f, 9.0f);

    /*****************************************************************************
     *  Update
     ****************************************************************************/
    switch (measurementPack.sensorType)
    {
    case MeasurementPackage::RADAR:
      {
        // convert current prediction to polar and calculate measurement delta
        const TVector &x(mEkf.GetState());
        float roh = sqrt(x(0) * x(0) + x(1) * x(1));
        float phi = atan2(x(1), x(0));
        float rohdot = (x(0) * x(2) + x(1) * x(3)) / roh;
        TVector zd(3);
        zd << measurementPack.rawMeasurements(0) - roh,
              GetTools().CalculateAngleDelta(phi, measurementPack.rawMeasurements(1)),
              measurementPack.rawMeasurements(2) - rohdot;

        // calculate new approximation for measurement matrix
        mHjRadar = GetTools().CalculateRadarJacobian(x);

        // finally update
        mEkf.UpdateEKF(zd, mHjRadar, mRRadar);
      }
      break;
    case MeasurementPackage::LASER:
      // just update
      mEkf.Update(measurementPack.rawMeasurements, mHLaser, mRLaser);
      break;
    }
  }
  //std::cout << mEkf << std::endl;
  mPreviousTimestamp = measurementPack.timestamp;
}

TVector FusionEKF::GetEstimate() const
{
  return mEkf.GetState();
}

