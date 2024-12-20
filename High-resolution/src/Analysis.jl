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

function split_indices(start_index, end_index) 
	total_distance = end_index - start_index  
	portions = ceil(total_distance / 500)  
	while total_distance / portions < 100 
		portions += 1 
	end  
	step = total_distance / portions 
	split_indices = [round(Int, start_index + i * step) for i in 1:portions-1] 
    prepend!(split_indices, [start_index])
    push!(split_indices, end_index)
	return split_indices 
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
            @views bounding_box = bounding_boxes[i]
            @views lower_row = first(bounding_box)[1] 
            @views upper_row = last(bounding_box)[1] 
            @views lower_column = first(bounding_box)[2] 
            @views upper_column = last(bounding_box)[2] 
            # Split the bounding box into equal-sized chunks
            row_splits = split_indices(lower_row, upper_row)
            column_splits = split_indices(lower_column, upper_column)
            # Delete the original from bounding boxes and add the chunks
            deleteat!(bounding_boxes, i)
            for j in 1:length(row_splits)-1
                for k in 1:length(column_splits)-1
                    push!(bounding_boxes, CartesianIndices((row_splits[j]:row_splits[j+1], column_splits[k]:column_splits[k+1])))
                end
            end
        end
    end
    return bounding_boxes 
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
        @show length(chunk_indices)
        for (j,chunk_index) in enumerate(chunk_indices)
            @show j
            chunk = image[chunk_index,:]
            @views crop_mask_chunk = crop_mask[chunk_index]
            height, width, slices = size(chunk)
            warped_slices = round(Int, slices*z_res/xy_res)
            warped = warpedview(chunk, invtform(xy_res, z_res), 
                      (1:height, 1:width, 1:warped_slices); 
                      method=Lanczos(4)) 
            for i in 1:warped_slices 
                @views slice_ = warped[:, :, i]
                otsu_thresh = find_threshold(slice_, Otsu())
                @views biovolume += sum(crop_mask_chunk .* (slice_ .> max(otsu_thresh/3, cell_threshold))) 
            end
            warped = nothing
            chunk = nothing
        end
        push!(biovolumes, biovolume)
        @show biovolumes
    end
end

main()
