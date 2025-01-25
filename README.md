# Brightfield-biofilm-assay
Code for analyzing biofilms by bright-field microscopy

# Using the web app
## First time:
Install Julia from the following link: https://julialang.org/downloads/ then download this repository and unzip to a desired directory.

## General Usage

`cd` to your project directory in terminal (e.g., `cd /Downloads/Brightfield-biofilm-assay/WebApp`), then run:

```bash
$> julia
`````

```julia
julia> using Pkg; Pkg.activate("."); Pkg.instantiate()
`````
Note: Package installation may take several minutes. 

Then run the app

```julia
julia> using GenieFramework; Genie.loadapp(); up()
`````

To access the web app, open your browser and navigate to `http://localhost:8000/`. Sometimes `http://localhost:8000/?fake=true` is necessary if running the app through WSL.

# Inside the web app:
Data format: To prepare images for upload, compress tiff files to .zip format. **Do not zip a folder containing images, rather select all images and zip.** 

Upload data: Select uploads by pressing "+" within the "Upload image folders" box. 

Calibration: Within the "Select folder to calibrate settings" box, select which zipped file you would like to unzip for selecting Imin/Imax and for display. This is not strictly necessary if you don't have Imin/Imax, don't want to display raw images or test thresholding, and want to apply the same settings to all images. Select your Imin and Imax files from the options (if applicable). 

Display:  Optionally choose raw images to display (one image, unless
several images comprise a timelapse, in which case the files have to satisfy the naming convention "...int.tif"), where "int" is an integer. 

Threshold: To define biofilm vs background, change threshold using the slider. Default threshold value is 0.04. In general, the appropriate threshold will depend on microscope. Choose images to test the threshold (optional). Same options as apply to the raw image display. 

Optional checkboxes: 
"Perform dust correction": applicable to timelapses for which
the first timepoint is effectively a dust-only image; that is, the cell density is sufficiently low that cells will not be considered dust. 
"Timelapse format":  indicates that individual files are individual timepoints in a timelapse -- the criterion for which is described above. 
"Analyze all folders":  Use the same Imin/Imax (if relevant) and threshold to analyze all folders (zipped files) in batch mode. If this option is not selected, the zipped file that
was selected for settings calibration will be analyzed. 

Run analysis: Execute the analysis. Once this is finished, optionally display processed images/masks. At this stage, only single file selection is possible, because all timelapses have been aggregated into tif series. 

Download processed images and numerical data: Returns a zipped file with processed images and quantification.

Done session/cleanup: When an analysis session is complete, click the red button at the bottom, which will delete all uploads, display files, and processed files from the WebApp directory.

# Once done with the app:
Run:
```julia
julia> down(); exit()
`````

# Known bugs:
The active tab for display "IMAGES" or "MASKS" will only update upon web page reload. Need to figure out how to auto-refresh.
