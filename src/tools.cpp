#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using namespace std;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {


  // check the validity of the following inputs:
  if(estimations.size() == 0 || estimations.size() != ground_truth.size()){
    throw std::invalid_argument("CalculateRMSE() - estimation size and, or ground truth size are invalid");
  }

  VectorXd rmse(estimations[0].array().size());
  rmse.fill(0.0);

  //accumulate squared residuals
  for(int i=0; i < estimations.size(); ++i){
    // ... your code here
    VectorXd residual = estimations[i] - ground_truth[i];
    // coefficient-wise multiplication
    residual = residual.array()*residual.array();
    rmse += residual;
  }

  //calculate the mean
  rmse = rmse/estimations.size();

  //calculate the squared root
  rmse = rmse.array().sqrt();

  //return the result
  return rmse;
}

float Tools::CalculateNISConsistency(const std::vector<float> &nisValues, MeasurementPackage::SensorType  sensorType){
  // The laser and radar sensor devices Chi-Squared values from the Udacity class for 95 percentile
  // 2 degree of freedom (laser sensor) is: 5.991 and 3 degree of freedom (radar sensor) is: 7.815

  const float nisLaser = 5.991;
  const float nisRadar = 7.815;
  float nis;
  int nisCount = 0; // count for over 95 percentile

  if (sensorType == MeasurementPackage::RADAR){
    nis = nisRadar;
  }
  else if (sensorType == MeasurementPackage::LASER){
    nis = nisLaser;
  }

  for (int i = 0; i < nisValues.size(); i++){
    if (nisValues[i] > nis){
      nisCount++;
    }
  }

  return nisCount/(float)nisValues.size();
}
