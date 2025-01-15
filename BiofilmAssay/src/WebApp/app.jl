##################################################
# TODO: 
# 1. Figure out why image updates only work when switching regimes (i.e. image -> gif or gif -> image, not image -> image)
# Reload fixes the issue, but how do I get it to do this automatically?
# It seems to work fine in the tab that is not currently active
# 2. Bug with clearing selected files 
# 3. Clean up redundant code in app.jl
# 4. Delete folders after use?
##################################################

using Main.Analysis
using Main.Displays
using GenieFramework
using FileIO
using TiffImages
@genietools

const FILE_PATH = joinpath("public", "uploads")
const DOWNLOADS_PATH = joinpath("public", "downloads")
const DISPLAYS_PATH = joinpath("public", "displays")
const DOWNLOADS_ZIP = joinpath("public", "data.zip")
mkpath(FILE_PATH)
mkpath(DOWNLOADS_PATH)
mkpath(DISPLAYS_PATH)
@app begin
    ##############################
    # Uploads 
    ##############################
    @out upfiles = readdir(FILE_PATH)
    @onchange fileuploads begin
        if ! isempty(fileuploads)
            @info "File was uploaded: " fileuploads
            filename = fileuploads["name"]

            try
                isdir(FILE_PATH) || mkpath(FILE_PATH)
                mv(fileuploads["path"], joinpath(FILE_PATH, filename), force=true)
            catch e
                @error "Error processing file: $e"
                notify(__model__,"Error processing file: $(fileuploads["name"])")
            end

            fileuploads = Dict{AbstractString,AbstractString}()
        end
        upfiles = readdir(FILE_PATH)
    end
    @event uploaded begin
        @info "uploaded"
        notify(__model__, "File was uploaded")
    end
    @event rejected begin
        @info "rejected"
        notify(__model__, "Please upload a valid file")
    end

    ##############################
    # Main
    ##############################
    @in selected_Imin = ""
    @in selected_Imax = ""
    @in dust_correction = false 
    @in batch_processing = false 
    @in fixed_thresh = 0.0400
    @in tab_selected = "Images"
    
    @in selected_raw_display_files = [""]
    @out upfiles = readdir(FILE_PATH)
    @in DisplayRawButtonProgress_process = false
    @out DisplayRawButtonProgress_processing = false
    @out image_display_path = ""
    @out mask_display_path = ""
    @onbutton DisplayRawButtonProgress_process begin
        DisplayRawButtonProgress_processing = true 
        DISPLAY_IMGPATH = ""
        DISPLAY_MASKPATH = ""
        image_display_path = "/$DISPLAY_IMGPATH#$(Base.time())" 
        mask_display_path = "/$DISPLAY_MASKPATH#$(Base.time())"
        filter!(e -> e != "", selected_raw_display_files)
        raw_image = load_raw_display(FILE_PATH, selected_raw_display_files)
        if length(selected_raw_display_files) > 0 
            if length(selected_raw_display_files) == 1 && length(size(raw_image)) == 2
                DISPLAY_IMGPATH = "displays/display.jpg"
                FileIO.save(joinpath(@__DIR__, "public", DISPLAY_IMGPATH), raw_image)
                image_display_path = "/$DISPLAY_IMGPATH#$(Base.time())" 
            else
                DISPLAY_IMGPATH = "displays/display.gif"
                FileIO.save(joinpath(@__DIR__, "public", DISPLAY_IMGPATH), raw_image, fps=10)
                image_display_path = "/$DISPLAY_IMGPATH#$(Base.time())" 
            end
        else
            image_display_path = ""
            mask_display_path = ""
            DISPLAY_IMGPATH = ""
            DISPLAY_MASKPATH = ""
            image_display_path = "/$DISPLAY_IMGPATH#$(Base.time())"
            mask_display_path = "/$DISPLAY_MASKPATH#$(Base.time())"
        end
        DisplayRawButtonProgress_processing = false
    end
    
    @in selected_test_files = [""]
    @out upfiles = readdir(FILE_PATH)
    @in TestButtonProgress_process = false
    @out TestButtonProgress_processing = false
    @onbutton TestButtonProgress_process begin
        TestButtonProgress_processing = true
        DISPLAY_IMGPATH = ""
        DISPLAY_MASKPATH = ""
        image_display_path = "/$DISPLAY_IMGPATH#$(Base.time())" 
        mask_display_path = "/$DISPLAY_MASKPATH#$(Base.time())"
        filter!(e -> e != "", selected_test_files)
        if length(selected_test_files) > 0 
            if length(selected_test_files) == 1
                @views dummy_img = TiffImages.load(joinpath(FILE_PATH, selected_test_files[1]))
                if length(size(dummy_img)) == 2
                    DISPLAY_IMGPATH = "displays/display.jpg"
                    DISPLAY_MASKPATH = "displays/mask_display.jpg"
                    @views normalized, overlay = image_test(FILE_PATH, selected_test_files[1], fixed_thresh, 101, 2)
                    FileIO.save(joinpath(@__DIR__, "public", DISPLAY_IMGPATH), normalized)
                    FileIO.save(joinpath(@__DIR__, "public", DISPLAY_MASKPATH), overlay)
                else
                    DISPLAY_IMGPATH = "displays/display.gif"
                    DISPLAY_MASKPATH = "displays/mask_display.gif"
                    normalized, overlay = timelapse_test(FILE_PATH, selected_test_files[1], fixed_thresh, 101, 2)
                    FileIO.save(joinpath(@__DIR__, "public", DISPLAY_IMGPATH), normalized, fps=10)
                    FileIO.save(joinpath(@__DIR__, "public", DISPLAY_MASKPATH), overlay, fps=10)
                end
            else
                DISPLAY_IMGPATH = "displays/display.gif"
                DISPLAY_MASKPATH = "displays/mask_display.gif"
                normalized, overlay = timelapse_test(FILE_PATH, selected_test_files, fixed_thresh, 101, 2)
                FileIO.save(joinpath(@__DIR__, "public", DISPLAY_IMGPATH), normalized, fps=10)
                FileIO.save(joinpath(@__DIR__, "public", DISPLAY_MASKPATH), overlay, fps=10)
            end
            image_display_path = "/$DISPLAY_IMGPATH#$(Base.time())" 
            mask_display_path = "/$DISPLAY_MASKPATH#$(Base.time())"
        else
            image_display_path = ""
            mask_display_path = ""
            DISPLAY_IMGPATH = ""
            DISPLAY_MASKPATH = ""
            image_display_path = "/$DISPLAY_IMGPATH#$(Base.time())"
            mask_display_path = "/$DISPLAY_MASKPATH#$(Base.time())"
        end
        TestButtonProgress_processing = false
    end
    
    @in AnalyzeButtonProgress_process = false
    @out AnalyzeButtonProgress_processing = false
    @out processed_images = readdir(DOWNLOADS_PATH)
    @onbutton AnalyzeButtonProgress_process begin
        AnalyzeButtonProgress_processing = true
        analysis_main(FILE_PATH, DOWNLOADS_PATH, upfiles, dust_correction, batch_processing, fixed_thresh, selected_Imin, selected_Imax)
        # Create a zipped version of the downloads folder
        run(`zip -r $DOWNLOADS_ZIP $DOWNLOADS_PATH`)
        processed_images = readdir(DOWNLOADS_PATH)
        AnalyzeButtonProgress_processing = false
    end

    @out processed_images = readdir(DOWNLOADS_PATH)
    @in selected_processed_image = ""
    @in DisplayProcessedButtonProgress_process = false
    @out DisplayProcessedButtonProgress_processing = false
    @onbutton DisplayProcessedButtonProgress_process begin
        DisplayProcessedButtonProgress_processing = true
        DISPLAY_IMGPATH = ""
        DISPLAY_MASKPATH = ""
        image_display_path = "/$DISPLAY_IMGPATH#$(Base.time())"
        mask_display_path = "/$DISPLAY_MASKPATH#$(Base.time())"
        filter!(e -> e != "", selected_raw_display_files)
        processed_image, overlay = load_processed_display(DOWNLOADS_PATH, selected_processed_image)
        if selected_processed_image != "" 
            if length(size(processed_image)) == 2
                DISPLAY_IMGPATH = "displays/display.jpg"
                DISPLAY_MASKPATH = "displays/mask_display.jpg"
                FileIO.save(joinpath(@__DIR__, "public", DISPLAY_IMGPATH), processed_image)
                FileIO.save(joinpath(@__DIR__, "public", DISPLAY_MASKPATH), overlay)
                image_display_path = "/$DISPLAY_IMGPATH#$(Base.time())" 
                mask_display_path = "/$DISPLAY_MASKPATH#$(Base.time())"
            else
                DISPLAY_IMGPATH = "displays/display.gif"
                DISPLAY_MASKPATH = "displays/mask_display.gif"
                FileIO.save(joinpath(@__DIR__, "public", DISPLAY_IMGPATH), processed_image, fps=10)
                FileIO.save(joinpath(@__DIR__, "public", DISPLAY_MASKPATH), overlay, fps=10)
                image_display_path = "/$DISPLAY_IMGPATH#$(Base.time())" 
                mask_display_path = "/$DISPLAY_MASKPATH#$(Base.time())"
            end
        else
            image_display_path = ""
            mask_display_path = ""
            DISPLAY_IMGPATH = ""
            DISPLAY_MASKPATH = ""
            image_display_path = "/$DISPLAY_IMGPATH#$(Base.time())"
            mask_display_path = "/$DISPLAY_MASKPATH#$(Base.time())"
        end
        selected_processed_image = ""
        DisplayProcessedButtonProgress_processing = false
    end
end

@page("/", "ui.jl")
