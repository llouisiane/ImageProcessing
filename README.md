# ImageProcessing

Dependencies:

C++
- [OpenCV](http://opencv.org)
- [Voro++](http://math.lbl.gov/voro++/)

Python
- [SciPy](http://www.scipy.org)
- [NumPy](http://www.numpy.org)
- [Matplotlib](http://matplotlib.org/)

## imageprocessing
detects bacteria, outputs positions, angles and sizes of rectangles
[Raw data](https://www.youtube.com/watch?v=XVsikfZki0Q) in real time.
[Segmented (after imageprocessing)](https://www.youtube.com/watch?v=i4Po9AJZ46s).

## tracking
establishes frame-to-frame [trajectories](https://www.youtube.com/watch?v=vGFfmc9co-Y)

### ids_generate
assigns unique ids to particles

### spline_generate
smooths trajectories by splines

## clustering
finds bacteria in the same cluster

## normalize
simple image normalization

## opticalflow
outputs optical flow grid and calculates correlation functions
