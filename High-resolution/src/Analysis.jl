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
    images_path = "/mnt/e/Brightfield_paper/"
    numerical_data_path = "/mnt/e/Brightfield_paper/Data/"
    cell_threshold = 0.0005 #0.002
	xy_res = 0.325
    z_res = 0.5
    z_xy_ratio = z_res / xy_res 
    zstack_thresh = 0.09 #0.03
    biomasses = []
    files = [f for f in readdir(images_path, join=true) if (occursin("noback", f) && (occursin("D39", f) || occursin("SV36", f)))]
    for filename in files 
        image = TiffImages.load(filename) 
        height, width, slices = size(image)
        zstack = imfilter(dropdims(sum(image, dims=3), dims=3), Kernel.gaussian(10))
        otsu_thresh = find_threshold(zstack, Otsu())
        @show otsu_thresh
        @show filename
        crop_mask = zstack .> max(zstack_thresh, otsu_thresh*1.5)#max(zstack_thresh, otsu_thresh)
        warped = imresize(image, ratio=(1,1,z_res/xy_res)) 
        biomass = 0
        for i in 1:size(warped, 3)
            @views slice_ = warped[:, :, i]
            otsu_thresh = find_threshold(slice_, Otsu())
            @views biomass += sum(crop_mask .* (slice_ .> max(otsu_thresh, cell_threshold))) 
        end
        push!(biomasses, Float64(biomass))
        @show biomasses 
        image = nothing
        warped = nothing
        zstack = nothing
        crop_mask = nothing
    end
	data = Dict(zip(files, biomasses))
    CSV.write("../../Data/high_res_data_sp.csv", data)
end

main()
