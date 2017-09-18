#ifndef KALMAN_FILTER_H_
#define KALMAN_FILTER_H_
#include "tools.h"

class KalmanFilter
{
public:
  /**
  * Constructor
  **/
  KalmanFilter();

  /**
  * Destructor
  **/
  virtual ~KalmanFilter();

  /**
  * Initializes Kalman filter
  * @param x_in Initial state
  * @param P_in Initial state covariance
  * @param F_in Transition matrix
  **/
  void Initialize(TVector &x, TMatrix &P, TMatrix &F);

  /**
  * Prediction Predicts the state and the state covariance
  * using the process model
  * @param dTime Time between k and k+1 in s
  * @param noise_ax Accelleration noise in x direction
  * @param noise_ay Accellaration noice in y direction
  **/
  void Predict(float dTime, float noise_ax, float noise_ay);

  /**
  * Updates the state by using standard Kalman Filter equations
  * @param z The measurement at k+1
  * @param H The measurement matrix
  * @param R The measurement convariance matrix
  **/
  void UpdateKF(const TVector &z, const TMatrix &H, const TMatrix &R);

  /**
  * Updates the state by using Extended Kalman Filter equations
  * @param dz The measurement delta between k and k+1
  * @param H The measurement matrix
  * @param R The measurement convariance matrix
  **/
  void UpdateEKF(const TVector &dz, const TMatrix &H, const TMatrix &R);

  /**
  * Retrieves currently predicted state x.
  **/
  inline const TVector &GetState() const;

protected:
  /**
  * Helper function to be used by update functions.
  **/
  void Update(const TVector &y, const TMatrix &H, const TMatrix &R);

  /**
  * Serialized filer state into output stream.
  **/
  friend std::ostream &operator<<(std::ostream &out, const KalmanFilter &kf)
  {
    out << "{ x=" << kf.mX << ", P=" << kf.mP << " }";
    return out;
  }

private:
  /** state vector **/
  TVector mX;

  /** state covariance matrix **/
  TMatrix mP;

  /** state transition matrix **/
  TMatrix mF;

  /** process covariance matrix **/
  TMatrix mQ;

  /** identity matrix **/
  TMatrix mI;
};

//=================================== Inlines ================================

const TVector &KalmanFilter::GetState() const
{
  return mX;
}

#endif /* KALMAN_FILTER_H_ */
