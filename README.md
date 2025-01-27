# Brightfield-biofilm-assay
Code for analyzing biofilms by bright-field microscopy

# Using the web app
## First time:
Install Julia from the following link: https://julialang.org/downloads/ then download this repository and unzip to a desired directory.

## General Usage

`cd` to your project directory in terminal (e.g., `cd /Downloads/Brightfield-biofilm-assay/WebApp`). The path will depend on your operating system (the aforementioned command would work for Unix/Linux OS) and the method by which you downloaded this respository. For Windows OS, the command might look like `cd Downloads\Brightfield-biofilm\WebApp` -- that is, no slash before the parent directory and slash direction is reversed. If you downloaded the code directly (not via git), the path will look something like `Downloads\Brightfield-biofilm-assay-main\Brightfield-biofilm-assay-main\WebApp`. After changing directories, run:

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
**Data format:** To prepare images for upload, compress tiff files to .zip format. **Do not zip a folder containing images, rather select all images and zip.** If you have Imin and Imax images, these should be zipped together with your other images.

**Upload data:** Select uploads by pressing "+" within the "Upload image folders" box. 

**Calibration:** Within the "Select folder to calibrate settings" box, select which zipped file you would like to unzip for selecting Imin/Imax and for display. This is not strictly necessary if you don't have Imin/Imax, don't want to display raw images or test thresholding, and want to apply the same settings to all images. Select your Imin and Imax files from the options (if applicable). 

**Display:**  Optionally choose raw images to display: one image, unless
several images comprise a timelapse, in which case the files have to satisfy the naming convention "...int.tif"), where "int" is an integer. For example, you might have images with file names that look like "A1\_001.tif", "A1\_002.tif", etc. If the file name does not terminate in a number (which denotes the temporal order), or the portion upstream of the number is not the same for all images that comprise a timelapse, there is no way for the code to know that the images should be
analyzed as an ensemble. 

**Threshold:** To define biofilm vs background, change threshold using the slider. Default threshold value is 0.04. You can use the slider or arrow keys to adjust the value (the latter is easier for fine adjustments). In general, the appropriate threshold will depend on microscope. Choose images to test the threshold (optional). Same options as apply to the raw image display. The image testing might take a while depending on how many images you select. You will know that testing is done when images appear on the right-hand side of the GUI and when the progress indicator stops. 

**Optional checkboxes:** 
"Perform dust correction": applicable to timelapses for which
the first timepoint is effectively a dust-only image; that is, the cell density is sufficiently low that cells will not be considered dust. 
"Timelapse format":  indicates that individual files are individual timepoints in a timelapse -- the criterion for which is described above. 
"Analyze all folders":  Use the same Imin/Imax (if relevant) and threshold to analyze all folders (zipped files) in batch mode. If this option is not selected, the zipped file that
was selected for settings calibration will be analyzed. 

**Run analysis:** Execute the analysis. Once this is finished, optionally display processed images/masks. At this stage, only single file selection is possible, because all timelapses have been aggregated into tif series. This step might take a while, depending on how many images/datasets you are analyzing and the sizes of the constituent images. You can view a progress bar representing the percent of the analysis that has been completed in the terminal. Keep in mind, if you are analyzing
multiple sets of images, this progress bar will only indicate the progress on the image set it is currently processing.

**Download processed images and numerical data:** Returns a zipped file with processed images and quantification.

**Done session/cleanup:** When an analysis session is complete, click the red button at the bottom, which will delete all uploads, display files, and processed files from the WebApp directory.

# Once done with the app:
Run:
```julia
julia> down(); exit()
`````

# Known bugs:
The active tab for display "IMAGES" or "MASKS" will only update upon web page reload. Need to figure out how to auto-refresh.
