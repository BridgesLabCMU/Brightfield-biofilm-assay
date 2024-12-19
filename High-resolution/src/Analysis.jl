using Images
using TiffImages: load
using HistogramThresholding
using ImageView
using FLoops
using ImageTransformations
using CoordinateTransformations
using Interpolations

function invtform(xy_res, z_res)
    M = [1 0 0 ; 0 1 0; 0 0 xy_res/z_res]
    LinearMap(M)
end 

function find_chunks(array, min_side, max_side)
    labeled_array = label_components(array)
    bounding_boxes = component_boxes(labeled_array)
    areas = component_lengths(labeled_array)
    biofilms = findall(x->x>1000, areas)
    bounding_boxes = bounding_boxes[biofilms]
    areas = areas[biofilms]
    popfirst!(bounding_boxes)
    popfirst!(areas)
    num_components = length(bounding_boxes)-1
    chunks = []
    for i in 1:num_components
        @views if areas[i] > 250000
            bounding_box = bounding_boxes[i]
            # Split the bounding box into equal-sized chunks
            # Delete the original from bounding boxes and add the chunks
        end
    end
    return chunks
end

function main()
    path = "/Volumes/T7 Shield/Brightfield_paper/"
    cell_threshold = 3.0e-5
	xy_res = 0.065 
    z_res = 0.3
    zstack_thresh = 0.4368
    biovolumes = []
    for filename in [f for f in readdir(path, join=true) if occursin(r"\.tif$", f)] 
        image = load(filename; lazyio=true) 
        zstack = imfilter(dropdims(sum(image, dims=3), dims=3), Kernel.gaussian(10))
        otsu_thresh = find_threshold(zstack, Otsu())
        if occursin("ara",filename)
            crop_mask = zstack .> max(zstack_thresh, otsu_thresh) 
        else
            @show otsu_thresh
            crop_mask = zstack .> otsu_thresh 
        end
        imshow(crop_mask)
        chunk_indices = find_chunks(crop_mask, 100, 500)
        biovolume = 0
        for chunk_index in chunk_indices
            chunk = image[chunk_index...,:]
            height, width, slices = size(chunk)
            warped_slices = round(Int, slices*z_res/xy_res)
            warped =  warpedview(image, invtform(xy_res, z_res), 
                      (1:height, 1:width, 1:warped_slices); 
                      method=Lanczos(4)) 
            @floop for i in 1:warped_slices 
                @show i
                @views slice_ = warped[:, :, i]
                otsu_thresh = find_threshold(slice_, Otsu())
                biovolume += sum(crop_mask .* (slice_ .> max(otsu_thresh/3, cell_threshold))) 
            end
        end
        push!(biovolumes, biovolume)
        @show biovolumes
    end
end

main()
