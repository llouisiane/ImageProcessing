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

[Raw data](https://www.youtube.com/watch?v=XVsikfZki0Q) in real time.

[![Segmented](http://img.youtube.com/vi/i4Po9AJZ46s/0.jpg)](https://www.youtube.com/watch?v=i4Po9AJZ46s)

## tracking
establishes frame-to-frame trajectories

[![Tracked](http://img.youtube.com/vi/vGFfmc9co-Y/0.jpg)](https://www.youtube.com/watch?v=vGFfmc9co-Y)

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
