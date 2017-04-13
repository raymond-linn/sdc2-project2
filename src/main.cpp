
#include <iostream>
#include "Eigen/Dense"
#include <vector>
#include <iomanip>
#include "ukf.h"
#include "measurement_package.h"
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include "tools.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

void check_arguments(int argc, char* argv[]) {
  string usage_instructions = "Usage instructions: ";
  usage_instructions += argv[0];
  usage_instructions += " path/to/input.txt output.txt";

  bool has_valid_args = false;

  // make sure the user has provided input and output files
  if (argc == 1) {
    cerr << usage_instructions << endl;
  } else if (argc == 2) {
    cerr << "Please include an output file.\n" << usage_instructions << endl;
  } else if (argc == 3) {
    has_valid_args = true;
  } else if (argc > 3) {
    cerr << "Too many arguments.\n" << usage_instructions << endl;
  }

  if (!has_valid_args) {
    exit(EXIT_FAILURE);
  }
}

void check_files(ifstream& in_file, string& in_name,
                 ofstream& out_file, string& out_name) {
  if (!in_file.is_open()) {
    cerr << "Cannot open input file: " << in_name << endl;
    exit(EXIT_FAILURE);
  }

  if (!out_file.is_open()) {
    cerr << "Cannot open output file: " << out_name << endl;
    exit(EXIT_FAILURE);
  }
}

int main(int argc, char* argv[]) {

  check_arguments(argc, argv);

  string in_file_name_ = argv[1];
  ifstream in_file_(in_file_name_.c_str(), ifstream::in);

  string out_file_name_ = argv[2];
  ofstream out_file_(out_file_name_.c_str(), ofstream::out);

  check_files(in_file_, in_file_name_, out_file_, out_file_name_);

  /**********************************************
   *  Set Measurements                          *
   **********************************************/

  vector<MeasurementPackage> measurement_pack_list;
  // to Calculate RMSE
  vector<Eigen::VectorXd> estimations, ground_truths;
  // for ni values calculations
  vector<float> nisRadarValues, nisLidarValues;
  string line;

  // prep the measurement packages (each line represents a measurement at a
  // timestamp)
  while (getline(in_file_, line)) {
    string sensor_type;
    MeasurementPackage meas_package;
    istringstream iss(line);
    long timestamp;

    // reads first element from the current line
    iss >> sensor_type;

    if (sensor_type.compare("L") == 0) {
      // laser measurement

      // read measurements at this timestamp
      meas_package.sensor_type_ = MeasurementPackage::LASER;
      meas_package.raw_measurements_ = VectorXd(2);
      float px;
      float py;
      iss >> px;
      iss >> py;
      meas_package.raw_measurements_ << px, py;
      iss >> timestamp;
      meas_package.timestamp_ = timestamp;
      measurement_pack_list.push_back(meas_package);
    } else if (sensor_type.compare("R") == 0) {
      // radar measurement

      // read measurements at this timestamp
      meas_package.sensor_type_ = MeasurementPackage::RADAR;
      meas_package.raw_measurements_ = VectorXd(3);
      float ro;
      float theta;
      float ro_dot;
      iss >> ro;
      iss >> theta;
      iss >> ro_dot;
      meas_package.raw_measurements_ << ro, theta, ro_dot;
      iss >> timestamp;
      meas_package.timestamp_ = timestamp;
      measurement_pack_list.push_back(meas_package);
    }
    // populate the ground_truths vectors to calculate rmse
    float x_gt, y_gt, vx_gt, vy_gt;
    iss >> x_gt >> y_gt >> vx_gt >> vy_gt;
    VectorXd gt(4);
    gt << x_gt, y_gt, vx_gt, vy_gt;
    ground_truths.push_back(gt);
  }

  // Create a UKF instance
  UKF ukf;

  size_t number_of_measurements = measurement_pack_list.size();

  // start filtering from the second frame (the speed is unknown in the first
  // frame)
  for (size_t k = 0; k < number_of_measurements; ++k) {
    // set ground truth list to write out
    const auto& gt = ground_truths[k];
    // Call the UKF-based fusion
    ukf.ProcessMeasurement(measurement_pack_list[k]);

    // output the estimation
    const double x = ukf.x_(0);
    const double y = ukf.x_(1);
    const double v = ukf.x_(2);
    const double yaw = ukf.x_(3);
    const double yaw_d = ukf.x_(4);

    // write out to output file
    out_file_ << x << "\t" << y << "\t" << v << "\t" << yaw << "\t" << yaw_d << "\t";

    // populate estimations list to calculate rmse
    const double vx = v*cos(yaw);
    const double vy = v*sin(yaw);
    VectorXd est(4);
    est << x, y, vx, vy;
    estimations.push_back(est);

    // output the measurements
    if (measurement_pack_list[k].sensor_type_ == MeasurementPackage::LASER) {
      // output the estimation

      // p1 - meas
      out_file_ << measurement_pack_list[k].raw_measurements_(0) << "\t";

      // p2 - meas
      out_file_ << measurement_pack_list[k].raw_measurements_(1) << "\t";

      // write out nis values
      nisLidarValues.push_back(ukf.NIS_laser_);
    } else if (measurement_pack_list[k].sensor_type_ == MeasurementPackage::RADAR) {
      // output the estimation in the cartesian coordinates
      float ro = measurement_pack_list[k].raw_measurements_(0);
      float phi = measurement_pack_list[k].raw_measurements_(1);
      out_file_ << ro * cos(phi) << "\t"; // p1_meas
      out_file_ << ro * sin(phi) << "\t"; // p2_meas
      // write out nis values
      nisRadarValues.push_back(ukf.NIS_radar_);
    }

    // write out ground truth values
    out_file_ << gt(0) << "\t" << gt(1) << "\t" << gt(2) << "\t" << gt(3);
    out_file_ << "\n";
  }

  // calculate RMSE and NIS values
  Tools tools;
  cout << "Accuracy - RMSE:" << endl << tools.CalculateRMSE(estimations, ground_truths) << endl;
  cout << "NIS Radar:" << endl << setprecision(4) << setw(4)
       << tools.CalculateNISConsistency(nisRadarValues, MeasurementPackage::RADAR)*100 << "%" << endl;
  cout << "NIS Lidar:" << endl << setprecision(4) << setw(4)
       << tools.CalculateNISConsistency(nisLidarValues, MeasurementPackage::LASER)*100 << "%" << endl;

  // close files
  if (out_file_.is_open()) {
    out_file_.close();
  }

  if (in_file_.is_open()) {
    in_file_.close();
  }

  cout << "Done!" << endl;
  return 0;
}
