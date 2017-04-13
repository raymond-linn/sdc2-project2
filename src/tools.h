//
// Created by Raymond Linn on 4/12/17.
//

#ifndef TOOLS_H
#define TOOLS_H

#include <vector>
#include "Eigen/Dense"
#include "measurement_package.h"

class Tools
{
public:
    /**
     * Constructor
     */
    Tools();

    /**
     * Destructor
     */
    virtual ~Tools();

    /**
     * A helper method to calculate RMSE
     * @param estimations
     * @param ground_truth
     * @return rmse
     */
    Eigen::VectorXd CalculateRMSE(const std::vector<Eigen::VectorXd>& estimations,
                                    const std::vector<Eigen::VectorXd>& ground_truth);

    /**
     *
     * @param nisVlaues
     * @param sensorType
     * @return nis
     */
    float CalculateNISConsistency(const std::vector<float>& nisVlaues, MeasurementPackage::SensorType sensorType);

};

#endif //TOOLS_H
