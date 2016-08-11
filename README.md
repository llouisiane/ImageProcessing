# ImageProcessing

C++
Dependencies:
- [OpenCV](http://opencv.org)
- [Boost](http://www.boost.org/)


## segmentation
Detects bacteria on each picture (segmentation.cpp)
- use parameters from configuration file config.txt, command line argument: path of the configuration file
- input: raw data
- output: positions, angles (data_exp.txt) and sizes of rectangles (data_len.txt), pictures after preprocessing (e.g.: histogram equalization, binarization)

## tracking based on a kalman filter
Establishes trajectories over multiple frames (kalman.cpp)
- use parameters from configuration file config.txt, command line argument: path of the configuration file
- input: data_exp.txt, data_len.txt
- output:
updated positions and velocities estimate for each time step and trajectory index (updated_states.txt)
diagonal of the corresponding updated error covariance matrices (updated_error_covariance.txt)

## visualization
Enables to visualize the results (visualization.cpp)
- use parameters from configuration file config.txt, command line argument: path of the configuration file
- input: preprocessed pictures, updated_states.txt, updated_error_covariance.txt
- output: trajectories plotted on preprocessed images (pred_ORIGINAL_NAME)

