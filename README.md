# ImageProcessing

Dependencies:

C++
- [OpenCV](http://opencv.org)

Python
- [SciPy](http://www.scipy.org)
- [NumPy](http://www.numpy.org)
- [Matplotlib](http://matplotlib.org/)

## imageprocessing
detects bacteria, outputs positions, angles and sizes of rectangles

`ImageProcessing inpath filenametpl start stop X Y DX DY outpath`
<dl>
  <dt>inpath</dt>
  <dd>path with input images (including trailing slash)</dd>
  <dt>filenametpl</dt>
  <dd>filename template (printf syntax)</dd>
  <dt>start and stop</dt>
  <dd>ints used to format filenametpl</dd>
  <dt>X Y DX DY</dt>
  <dd>defines a subregion to read from the images</dd>
  <dt>outpath</dt>
  <dd>path for output (including trailing slash)</dd>
</dl>

More parameters in .cpp file (starting with "int width")  

[Raw data](https://www.youtube.com/watch?v=XVsikfZki0Q) in real time.

[![Segmented](http://img.youtube.com/vi/i4Po9AJZ46s/0.jpg)](https://www.youtube.com/watch?v=i4Po9AJZ46s)

## tracking
establishes frame-to-frame trajectories

Parameters in .cpp file (starting with "std::string dir")

[![Tracked](http://img.youtube.com/vi/vGFfmc9co-Y/0.jpg)](https://www.youtube.com/watch?v=vGFfmc9co-Y)

### ids_generate
assigns unique ids to particles

### spline_generate
smooths trajectories by splines (one parameter "smoothing")

## clustering
finds bacteria in the same cluster

## normalize
simple image normalization

## opticalflow
outputs optical flow grid and calculates correlation functions

## order of execution:
- imageprocessing onto raw images yields positions, orientations (data_exp.txt) and dimensions (data_len.txt)
- delete_duplicated_lines_data_exp yields data withouth duplicate output rectangles (due to unknown bug)
- tracking (yields [...]_tracked.txt)
- ids_generate (yields data_exp_tracked_ids.txt)
- splines_generate (yields [...]_splines.txt)
- final data files are data_exp_tracked_splines.txt, data_exp_tracked_ids.txt, data_len_tracked.txt (and data_vel_tracked_splines.txt)



## step-by-step instructions for setting up OpenCV for imageprocessing:
- download and install [Codeblocks](http://www.codeblocks.org) with mingw (includes gcc compiler)
- download and unpack [OpenCV](https://sourceforge.net/projects/opencvlibrary/files/) (tested 3.0)
- create symlink "opencv" to your OpenCV installation (super folder of "build", "sources" etc), e.g. under Windows using admin privileged command line and calling `mklink /J C:\[...\where your Codeblocks project is]\opencv C:\[...\where your OpenCV is]\opencv`
- download and install [CMake](https://cmake.org/)
- build opencv library with CMake: "[...]\opencv\sources" as source and "[...]\opencv\build\x86\mingw"(create folder first) as binaries destination, click "Configure", choose "MinGW Makefiles", after it finishes "Generate"
- add the `MinGW\bin` folder to your operating systems PATH variable
- in a command line, set your directory (using `cd`) to `[...]\opencv\build\x86`, run `mingw32-make`, and `mingw32-make install` after it finishes
- build and run imageprocessing.cbp
- alternatively [here](
https://mega.nz/#!3kdnHSxZ!lKF2blnbIrvAo7gzdzaepoIjdbDynHVd35vGuRkN9Ec) is our precompiled (for Windows 7 and newer) library