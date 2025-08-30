using Images
using CSV
using JSON
using TiffImages
using HistogramThresholding
using ImageView
using FLoops
using ImageTransformations
using CoordinateTransformations
using Interpolations
using StatsBase

function main()
    images_path = "/mnt/e/Spneumo_10-h_confocal_and_BF_growth-in-BioSpa/confocal"
    cell_threshold = 0.0003#0.0005
    xy_res = 0.325
    z_res = 0.5
    zstack_thresh = 0.02
    biomasses = []
    files = [f for f in readdir(images_path, join=true) if occursin("_noback.tif", f)]

    for noback_fn in files
        noback = TiffImages.load(noback_fn)

        _, _, _ = size(noback) 
        zstack = imfilter(
            dropdims(sum(noback, dims=3), dims=3),
            Kernel.gaussian(10)
        )

        otsu_thresh = find_threshold(zstack, Otsu())
        crop_mask = zstack .> max(zstack_thresh, otsu_thresh * 0.9)

        warped = imresize(noback, ratio=(1,1,z_res/xy_res))
        biomass = 0
        masks = similar(warped)

        for i in 1:size(warped, 3)
            @views slice_ = warped[:, :, i]
            otsu_thresh = find_threshold(slice_, Otsu())
            masks[:, :, i] = crop_mask .* (slice_ .> max(otsu_thresh, cell_threshold))
            @views biomass += sum(masks[:, :, i])
        end

        push!(biomasses, Float64(biomass))
        @show biomasses
        imshow(masks)

        # free memory
        noback = nothing
        denoised = nothing
        warped = nothing
        zstack = nothing
        crop_mask = nothing
    end

    data = Dict(zip(files, biomasses))
    CSV.write("../../Data/high_res_data_sp_new.csv", data)
end

main()
