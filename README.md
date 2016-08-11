# ImageProcessing

Dependencies:

C++
- [OpenCV](http://opencv.org)
- [Boost](http://www.boost.org/)


## segmentation
- detects bacteria on each picture (segmentation.cpp)
- use parameters from configuration file config.txt, command line argument: path of the configuration file
- input: raw data
- output: positions, angles (data_exp.txt) and sizes of rectangles (data_len.txt), pictures after preprocessing (e.g.: histogram equalization, binarization)

## tracking based on a kalman filter
- establishes trajectories over multiple frames (kalman.cpp)
- use parameters from configuration file config.txt, command line argument: path of the configuration file
- input: data_exp.txt, data_len.txt
- output:
updated positions and velocities estimate for each time step and trajectory index (updated_states.txt)
diagonal of the corresponding updated error covariance matrices (updated_error_covariance.txt)

## visualization
- enables to visualize the results (visualization.cpp)
- use parameters from configuration file config.txt, command line argument: path of the configuration file
- input: preprocessed pictures, updated_states.txt, updated_error_covariance.txt
- output: trajectories plotted on preprocessed images (pred_ORIGINAL_NAME)

## step-by-step instructions for setting up OpenCV for imageprocessing
Linux:
- follow one of the guides, eg. [here](http://docs.opencv.org/3.1.0/d7/d9f/tutorial_linux_install.html)
- you have to include the shared library files (*.so) from the `bin` folder of the build to the PATH variable or copy them next to executable
- you have to change the "search directories" and "linker settings" in codeblocks (`Project -> Build options -> Search directories / Linker settings`) to accommodate for your file structure (and possibly add search directories and some `-lopencv_` type 'Other linker options')

