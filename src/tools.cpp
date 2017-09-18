#include <iostream>
#include <limits>
#include <math.h>
#include "tools.h"

TVector Tools::CalculateRMSE(const TVectorList &estimations, const TVectorList &groundTruths) const
{
  TVector rmse(4);
  rmse << 0, 0, 0, 0;

  //check the validity of the inputs
  if(estimations.size() != groundTruths.size() || estimations.size() == 0)
  {
    std::cout << "Invalid estimation or ground_truth data" << std::endl;
    return rmse;
  }

  //accumulate squared residuals
  for(unsigned int i = 0; i < estimations.size(); ++i)
  {
    TVector residual = estimations[i] - groundTruths[i];
    residual = residual.array() * residual.array();
    rmse += residual;
  }

  //calculate the mean
  rmse = rmse / estimations.size();

  //calculate the squared root
  rmse = rmse.array().sqrt();

  //return the result
  return rmse;
}

float Tools::CalculateAngleDelta(float a1, float a2) const
{
  // normalize a1
  if (a1 > M_PI)
  {
    a1 -= M_PI * 2.0f;
  }
  else if (a1 < -M_PI)
  {
    a1 += M_PI * 2.0f;
  }

  // normalize a2
  if (a2 > M_PI)
  {
    a2 -= M_PI * 2.0f;
  }
  else if (a2 < -M_PI)
  {
    a2+= M_PI * 2.0f;
  }

  float aDelta = a2 - a1;

  // normalize delta
  if (aDelta > M_PI)
  {
    aDelta -= M_PI * 2.0f;
  }
  else if (aDelta < -M_PI)
  {
    aDelta += M_PI * 2.0f;
  }
  return aDelta;
}

TMatrix Tools::CalculateRadarJacobian(const TVector &xState) const
{
  TMatrix Hj(3, 4);

  //recover state parameters
  float px = xState(0);
  float py = xState(1);
  float vx = xState(2);
  float vy = xState(3);

  //pre-compute a set of terms to avoid repeated calculation
  float c1 = px*px+py*py;
  float c2 = sqrt(c1);
  float c3 = (c1*c2);
  float p1 = px*vy;
  float p2 = py*vx;

  //check division by zero
  if(isZero(c1))
  {
    std::cout << "CalculateJacobian () - Error - Division by Zero" << std::endl;
    Hj << 0, 0, 0, 0,
          0, 0, 0, 0;
    return Hj;
  }

  //compute the Jacobian matrix
  Hj << (px/c2)       , (py/c2)      , 0    , 0,
         -(py/c1)     , (px/c1)      , 0    , 0,
         py*(p2-p1)/c3, px*(p1-p2)/c3, px/c2, py/c2;

  return Hj;
}

bool Tools::isZero(float f) const
{
  return fabs(f) < std::numeric_limits<float>::epsilon();
}
