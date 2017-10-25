# Depth from Defocus for Misaligned Camera Photos Subject to Parallax Effects
## Abstract
This project implements an attempt to solve the problem of determining depth from
defocus when given a focal stack (set of photos taken at different focuses) created from photos captured with
a hand-held camera with unknown calibration. Although the photos may initially suffer from significant parallax
due to camera/hand movement during capture, the image alignment steps explained in Section III of the report are shown to be
capable of compensating for such issues. Relative depth results are computed quickly, though suffer from some
shortcomings around depth edges as discussed in Section VII. The final computed depth map is also shown to
be sufficient for rendering visually pleasing synthetically refocused photos of the scene.

## Writeup
* [Final report](http://robertmahieu.com/docs/mahieu_ee367_final_report.pdf)
* [Poster](http://robertmahieu.com/docs/mahieu_ee367_final_poster.pdf)

## Instructions
To run the algorithm, simply run the script 'main.m' in MATLAB. This will produce results within a folder called 'results' which will be created in the execution directory. To edit the photo set used by the algorithm, change the code in the top section of 'main.m' to reflect the appropriate path and number of images.

## Dependencies
The code requires three external packages for performing:
* Affine image registration ([Image Alignment Toolbox v0.9.2](http://iatool.net/))
* Optical flow computation ([OpticalFlow by Ce Liu](https://people.csail.mit.edu/celiu/OpticalFlow/))
* Solving graph cuts ([GCoptimization v3.0](http://vision.csd.uwo.ca/code/))

These are all included in the 'dependencies' folder. The folders for each of these contain pre-built mex files that should run fine on Windows x64, but I've included the original package downloads just in case something goes wrong.
