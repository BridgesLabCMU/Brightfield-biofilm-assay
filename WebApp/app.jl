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
using Main.utils
using GenieFramework
using FileIO
using TiffImages
@genietools

const FILE_PATH = joinpath("public", "uploads")
const DOWNLOADS_PATH = joinpath("public", "downloads")
const DISPLAYS_PATH = joinpath("public", "displays")
const DOWNLOADS_ZIP = joinpath("public", "downloads.zip")
if isdir(FILE_PATH)
    rm(FILE_PATH, recursive=true)
end
if isdir(DOWNLOADS_PATH)
    rm(DOWNLOADS_PATH, recursive=true)
end
if isdir(DISPLAYS_PATH)
    rm(DISPLAYS_PATH, recursive=true)
end
if isfile(DOWNLOADS_ZIP)
    rm(DOWNLOADS_ZIP)
end
mkpath(FILE_PATH)
mkpath(DOWNLOADS_PATH)
mkpath(DISPLAYS_PATH)
@app begin
    ##############################
    # Uploads 
    ##############################
    @out upfiles = readdir(FILE_PATH)
    @out zipped_files = [f for f in readdir(FILE_PATH) if occursin(".zip", f)]
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
        zipped_files = [f for f in readdir(FILE_PATH) if occursin(".zip", f)]
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
    @in analyze_folder = ""
    @in ConfirmButton_process = false
    @out ConfirmButton_processing = false
    @out analyze_folder_files = [""]
    @out folder = ""
    @onbutton ConfirmButton_process begin
        ConfirmButton_processing = true 
        folder = joinpath(FILE_PATH, splitext(analyze_folder)[1])
        if !isdir(folder)
            unzip_folder(joinpath(FILE_PATH, analyze_folder))
            mkpath(joinpath(DOWNLOADS_PATH, basename(folder)))
        end
        analyze_folder_files = readdir(joinpath(FILE_PATH, splitext(analyze_folder)[1]))
        ConfirmButton_processing = false 
    end

    @in selected_Imin = ""
    @in selected_Imax = ""
    @in dust_correction = false 
    @in batch_processing = false
    @in timelapse_flag = false 
    @in fixed_thresh = 0.0400
    @in tab_selected = "Images"
    
    @in selected_raw_display_files = [""]
    @in DisplayRawButtonProgress_process = false
    @out DisplayRawButtonProgress_processing = false
    @out image_display_path = ""
    @out mask_display_path = ""
    @onbutton DisplayRawButtonProgress_process begin
        DisplayRawButtonProgress_processing = true 
        folder = joinpath(FILE_PATH, splitext(analyze_folder)[1])
        DISPLAY_IMGPATH = ""
        DISPLAY_MASKPATH = ""
        image_display_path = "/$DISPLAY_IMGPATH#$(Base.time())" 
        mask_display_path = "/$DISPLAY_MASKPATH#$(Base.time())"
        filter!(e -> e != "", selected_raw_display_files)
        raw_image = load_raw_display(folder, selected_raw_display_files)
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
                @views dummy_img = TiffImages.load(joinpath(folder, selected_test_files[1]))
                if length(size(dummy_img)) == 2
                    DISPLAY_IMGPATH = "displays/display.jpg"
                    DISPLAY_MASKPATH = "displays/mask_display.jpg"
                    @views normalized, overlay = image_test(folder, selected_test_files[1], fixed_thresh, 101, 2)
                    FileIO.save(joinpath(@__DIR__, "public", DISPLAY_IMGPATH), normalized)
                    FileIO.save(joinpath(@__DIR__, "public", DISPLAY_MASKPATH), overlay)
                else
                    DISPLAY_IMGPATH = "displays/display.gif"
                    DISPLAY_MASKPATH = "displays/mask_display.gif"
                    normalized, overlay = timelapse_test(folder, selected_test_files[1], fixed_thresh, 101, 2)
                    FileIO.save(joinpath(@__DIR__, "public", DISPLAY_IMGPATH), normalized, fps=10)
                    FileIO.save(joinpath(@__DIR__, "public", DISPLAY_MASKPATH), overlay, fps=10)
                end
            else
                DISPLAY_IMGPATH = "displays/display.gif"
                DISPLAY_MASKPATH = "displays/mask_display.gif"
                normalized, overlay = timelapse_test(folder, selected_test_files, fixed_thresh, 101, 2)
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
        zipped_files = [f for f in readdir(FILE_PATH) if occursin("zip", f)]
        if batch_processing == true
            for zipped_file in zipped_files
                check_folder = joinpath(FILE_PATH, splitext(zipped_file)[1])
                if !isdir(check_folder)
                    unzip_folder(joinpath(FILE_PATH, zipped_file))
                    mkpath(joinpath(DOWNLOADS_PATH, basename(check_folder)))
                end
            end
            folders = [f for f in readdir(FILE_PATH, join=true) if !occursin("zip", f)]
            for f in folders
                analysis_main(f, joinpath(DOWNLOADS_PATH, basename(f)), readdir(f), dust_correction, timelapse_flag, fixed_thresh, selected_Imin, selected_Imax)
            end
        else
            if length(zipped_files) == 1  
                check_folder = joinpath(FILE_PATH, splitext(zipped_files[1])[1])
                if !isdir(check_folder)
                    unzip_folder(joinpath(FILE_PATH, zipped_file))
                    mkpath(joinpath(DOWNLOADS_PATH, basename(check_folder)))
                end
            end
            analysis_main(folder, joinpath(DOWNLOADS_PATH, basename(folder)), analyze_folder_files, dust_correction, timelapse_flag, fixed_thresh, selected_Imin, selected_Imax)
        end
        # Create a zipped version of the downloads folder
        write_zip(DOWNLOADS_ZIP, DOWNLOADS_PATH)
        processed_images = readdir(joinpath(DOWNLOADS_PATH, basename(folder)))
        AnalyzeButtonProgress_processing = false
    end

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
        processed_image, overlay = load_processed_display(joinpath(DOWNLOADS_PATH, basename(folder)), selected_processed_image)
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
    
    @in CleanupButtonProgress_process = false
    @out CleanupButtonProgress_processing = false
    @onbutton CleanupButtonProgress_process begin
        CleanupButtonProgress_processing = true 
        if isdir(FILE_PATH)
            rm(FILE_PATH, recursive=true)
            mkpath(FILE_PATH)
        end
        if isdir(DOWNLOADS_PATH)
            rm(DOWNLOADS_PATH, recursive=true)
            mkpath(DOWNLOADS_PATH)
        end
        if isdir(DISPLAYS_PATH)
            rm(DISPLAYS_PATH, recursive=true)
            mkpath(DISPLAYS_PATH)
        end
        if isfile(DOWNLOADS_ZIP)
            rm(DOWNLOADS_ZIP)
        end
        CleanupButtonProgress_processing = false 
        empty!(zipped_files)
        empty!(analyze_folder_files)
        analyze_folder = ""
        empty!(processed_images)
    end
end

@page("/", "ui.jl")
