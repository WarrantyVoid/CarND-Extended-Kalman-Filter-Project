#include "kalman_filter.h"

KalmanFilter::KalmanFilter()
: mX()
, mP()
, mF()
, mQ()
, mI()
{

}

KalmanFilter::~KalmanFilter()
{

}

void KalmanFilter::Initialize(TVector &x, TMatrix &P, TMatrix &F)
{
  mX = x;
  mP = P;
  mF = F;
  mI = TMatrix::Identity(mX.size(), mX.size());
}

void KalmanFilter::Predict(float dTime, float noise_ax, float noise_ay)
{
  // modify the F matrix so that the time is integrated
  mF(0, 2) = dTime;
  mF(1, 3) = dTime;

  // set the process covariance matrix Q based on time and noice
  float dt_2 = dTime * dTime;
  float dt_3 = dt_2 * dTime;
  float dt_4 = dt_3 * dTime;
  float n_ax = noise_ax;
  float n_ay = noise_ay;
  mQ = TMatrix(4, 4);
  mQ << dt_4/4*n_ax, 0          , dt_3/2*n_ax, 0,
        0          , dt_4/4*n_ay, 0          , dt_3/2*n_ay,
        dt_3/2*n_ax, 0          , dt_2*n_ax  , 0,
        0          , dt_3/2*n_ay, 0          , dt_2*n_ay;

  // predict
  mX = mF * mX;
  mP = mF * mP * mF.transpose() + mQ;
}

void KalmanFilter::UpdateKF(const TVector &z, const TMatrix &H, const TMatrix &R)
{
  Update(z - H * mX, H, R);
}

void KalmanFilter::UpdateEKF(const TVector &zd, const TMatrix &H, const TMatrix &R)
{
  Update(zd, H, R);
}

void KalmanFilter::Update(const TVector &y,const TMatrix &H, const TMatrix &R)
{
  TMatrix PHt = mP * H.transpose();
  TMatrix S = H * PHt + R;
  TMatrix K = PHt * S.inverse();

  mX = mX + (K * y);
  mP = (mI - K * H) * mP;
}
