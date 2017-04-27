# Unscented Kalman Filter Project Starter Code
Self-Driving Car Engineer Nanodegree Program

![Unscented Kalman filter demo](https://github.com/parilo/CarND-Unscented-Kalman-Filter-Project/blob/master/demo.png "Unscented Kalman filter demo")

This Unscented Kalman filter is used for fusion of measurements from lidar and radar of one particular object. Files:

* ukf.cpp holds unscented Kalman filter code and parameters
* tools.cpp contains code for RMSE calculation over test data.
* data folder contains test measurements
* main.cpp reads measurements from specified testing data file and calculates RMSE (root mean squared error).
* UKF.ipynb may be used to visualize program output file with filtered trajectory and NIS (Normalized Innovation Squared) metric

## Dependencies

* cmake >= v3.5
* make >= v4.1
* gcc/g++ >= v5.4

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./UnscentedKF path/to/input.txt path/to/output.txt`. You can find
   some sample inputs in 'data/'.
    - eg. `./UnscentedKF ../data/obj_pose-laser-radar-synthetic-input.txt`
