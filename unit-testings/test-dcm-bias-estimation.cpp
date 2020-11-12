#include <iostream>
#include <state-observation/dynamics-estimators/lipm-dcm-bias-estimator.hpp>
#include <state-observation/tools/probability-law-simulation.hpp>

using namespace stateObservation;

/// @brief runs the basic test
///
/// @return int : 0 if success, nonzero if fails
int testDcmBiasEstimator()
{

  double w0 = sqrt(cst::gravityConstant / 0.9);
  double dt = 0.005;

  double biasDriftPerSecondStd = 0.005;
  double zmpMeasurementErrorStd = 0.001;
  double dcmMeasurementErrorStd = 0.005;

  int signallength = int(120. / dt);

  ///////////////////////////////////////
  /// Build the ground truth signals
  ///////////////////////////////////////
  tools::ProbabilityLawSimulation ran;
  std::vector<double> dcm(signallength), bias(signallength), zmp(signallength);

  /// set the desired exponential convergence of the DCM
  double lambda = 2;

  /// initialize the dcm and bias to a random value
  dcm[0] = ran.getGaussianScalar(2, 0);
  double initBiasstd = 0.01;

  bias[0] = ran.getGaussianScalar(0.01, 0);
  double deviation = ran.getGaussianScalar(0.05, 0);
  zmp[0] = (lambda / w0 + 1) * dcm[0] + deviation;

  for(int i = 0; i < signallength - 1; ++i)
  {
    /// dcm dynamics
    dcm[i + 1] = dcm[i] + dt * w0 * (dcm[i] - zmp[i]);
    /// drift
    bias[i + 1] = bias[i] + ran.getGaussianScalar(biasDriftPerSecondStd * dt);
    /// set a noisy zmp to create  bounded drift of the DCM
    deviation += ran.getGaussianScalar(0.05, 0);
    zmp[i + 1] = (lambda / w0 + 1) * dcm[i] + ran.getGaussianScalar(0.05, 0) + deviation;
  }

  /////////////////////////////////
  /// Build the measurements
  /////////////////////////////////
  std::vector<double> dcm_m_unbiased(signallength), dcm_m(signallength), zmp_m(signallength);

  for(int i = 0; i < signallength; ++i)
  {
    dcm_m_unbiased[i] = dcm[i] + ran.getGaussianScalar(dcmMeasurementErrorStd);
    dcm_m[i] = dcm_m_unbiased[i] + bias[i];
    zmp_m[i] = zmp[i] + ran.getGaussianScalar(zmpMeasurementErrorStd);
  }

  /////////////////////////////////
  /// Run the estimator
  /////////////////////////////////
  std::vector<double> dcm_hat(signallength), bias_hat(signallength);

  double initbias = 0;

  LipmDcmBiasEstimator est(true, dcm_m[0], zmp_m[0], w0, dt, biasDriftPerSecondStd, zmpMeasurementErrorStd,
                           dcmMeasurementErrorStd, initbias, initBiasstd);

  IndexedVectorArray log;

  for(int i = 1; i < signallength; i++)
  {
    est.setInputs(dcm_m[i], zmp_m[i]);
    est.update();
    bias_hat[i] = est.getBias();
    dcm_hat[i] = est.getUnbiasedDCM();

    /// NaN detection
    if(dcm_hat[i] != dcm_hat[i])
    {
      return 1;
    }

    Vector6 log_k;

    log_k << i, bias_hat[i], bias[i], dcm_hat[i], dcm_m_unbiased[i], dcm[i];
    log.pushBack(log_k);
  }

  /// uncomment the following line to have a bit of a log
  // log.writeInFile("dcm.txt");

  /////////////////////////////////
  /// Check the rror
  /////////////////////////////////
  double error = 0;

  for(int i = signallength - 10; i < signallength; ++i)
  {
    error += fabs(dcm_hat[i] - dcm[i]) + fabs(bias_hat[i] - bias[i]);
  }

  std::cout << "Sum of error on the 10 last samples = " << error << std::endl;

  if(error > 0.02)
  {
    return 2;
  }
  else
  {
    return 0;
  }
}

int main()
{
  int exitCode;

  exitCode = testDcmBiasEstimator();

  if(exitCode != 0)
  {
    std::cout << "Failed, error code: " << exitCode << std::endl;
    return exitCode;
  }

  std::cout << "Succeeded" << std::endl;
  return exitCode;
}