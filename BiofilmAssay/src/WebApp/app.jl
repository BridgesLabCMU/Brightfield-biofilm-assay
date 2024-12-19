module App
using Main.Analysis.jl
using GenieFramework
@genietools

const FILE_PATH = joinpath("public", "uploads")
mkpath(FILE_PATH)

@app begin
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
    @in N = 0.03 
    @in ButtonProgress_process = false
    @in ButtonProgress_progress = 0.0
    @onbutton ButtonProgress_process begin
        for ButtonProgress_progress = 0:0.1:1
            @show ButtonProgress_progress
            sleep(0.5)
        end
        ButtonProgress_progress = 0.0
    end
    @in dust_correction = true
    @in bulk_processing = true
end

@page("/", "ui.jl")
end
