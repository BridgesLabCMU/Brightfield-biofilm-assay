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

#=
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
=#

function main()
    images_path = "/mnt/f/Brightfield_paper/"
    numerical_data_path = "/mnt/f/Brightfield_paper/Data/"
    cell_threshold = 0.029
	xy_res = 0.065 
    z_res = 0.3
    z_xy_ratio = z_res / xy_res 
    zstack_thresh = 0.5
    thicknesses = []
    files = [f for f in readdir(images_path, join=true) if occursin(r"\.tif$", f) && !occursin("D39", f) && !occursin("36", f)]
    for filename in files 
        @show filename
        image = TiffImages.load(filename; lazyio=true) 
        height, width, slices = size(image)
        zstack = imfilter(dropdims(sum(image, dims=3), dims=3), Kernel.gaussian(10))
        otsu_thresh = find_threshold(zstack, Otsu())
        crop_mask = zstack .> max(zstack_thresh, otsu_thresh) 
        thickness = 0
        mask = Array{Gray{Bool}, 3}(undef, size(image))
        for i in 1:slices
            @views slice_ = image[:, :, i]
            otsu_thresh = find_threshold(slice_, Otsu())
            if i == 1 
                cell_threshold = (otsu_thresh/0.0477)*0.029 
            end
            mask[:,:,i] = crop_mask .* (slice_ .> max(otsu_thresh, cell_threshold))
            @views thickness += sum(mask[:,:,i]) 
        end
        #imshow(crop_mask) 
        thickness *= (z_xy_ratio*xy_res/(height*width))
        push!(thicknesses, thickness)
        @show thicknesses 
        imshow(mask)
    end
	data = Dict(zip(files, thicknesses))
    CSV.write("../../Data/high_res_data.csv", data)
end

main()
