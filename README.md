# Brightfield-biofilm-assay
Code for analyzing biofilms by bright-field microscopy

# Using the web app
## First time:
Install Julia from the following link: https://julialang.org/downloads/

## General Usage

`cd` into the project directory in terminal (e.g., `cd /Downloads/Brightfield-biofilm-assay/WebApp`) then run:

```bash
$> julia
`````

```julia
julia> using Pkg; Pkg.activate("."); Pkg.instantiate()
`````

Then run the app

```julia
julia> using GenieFramework
julia> Genie.loadapp() # load app
julia> up() # start server
`````

Open your browser and navigate to `http://localhost:8000/`. Sometimes `http://localhost:8000/?fake=true` is necessary.

# Inside the web app:
To upload, zip all images together with the *.zip extension. Select uploads by pressing "+". After upload, select which zipped file you would like to unzip for selecting Imin/Imax and for display. This is not strictly necessary if you don't have Imin/Imax, don't want to display raw images or test thresholding, and want to apply the same settings to all images. Select your Imin and Imax files from the options (if applicable). Optionally choose raw images to display (one image, unless
several images comprise a timelapse, in which case the files have to satisfy the naming convention "...int.tif"), where "int" is an integer. Change threshold using the slider. Default threshold value is 0.04, which we find works well in almost all cases. Choose images to test the threshold (optional). Same options as apply to the raw image display. If satisfied with the threshold, optionally select checkboxes. "Perform dust correction" is only applicable to timelapses for which
the first timepoint is effectively a dust-only image; that is, the cell density is sufficiently low that cells will not be considered dust. "Timelapse format" indicates that individual files are individual timepoints in a timelapse -- the criterion for which is described above. "Analyze all folders" means that you would like to use the same Imin/Imax (if relevant) and threshold to analyze all folders (zipped files) in batch mode. If this option is not selected, the zipped file that
was selected for settings calibration will be analyzed. Click "run analysis" to perform the analysis. Once this is finished, optionally display processed images/masks. At this stage, only single file selection is possible, because all timelapses have been aggregated into tif series. Finally, download processed images and numerical data (a zipped file containing these). Once done with the session, click the red button at the bottom, which will delete all uploads, display files,
and processed files from the WebApp directory.

# Known bugs:
The active tab for display "IMAGES" or "MASKS" will only update upon web page reload. Need to figure out how to auto-refresh.
