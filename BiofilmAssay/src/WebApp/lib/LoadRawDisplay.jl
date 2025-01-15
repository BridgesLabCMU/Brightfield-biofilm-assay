module Displays 

export load_raw_display, load_processed_display, image_test, timelapse_test

using TiffImages
using Images
using NaturalSort
using IntegralArrays
using IntervalSets

function load_raw_display(file_path, filenames)
    ntimepoints = length(filenames)
    if ntimepoints == 0
        return nothing
    elseif ntimepoints == 1
        @views img = TiffImages.load(joinpath(file_path, filenames[1]))
    else
        @views dummy_img = TiffImages.load(joinpath(file_path, filenames[1]); mmap=true)
        height, width = size(dummy_img)
        img = Array{eltype(dummy_img), 3}(undef, height, width, ntimepoints)
        filenames = sort(filenames, lt=natural)
        for i in 1:ntimepoints
            @views img[:,:,i] = TiffImages.load(joinpath(file_path, filenames[i]))
        end
    end
    return img
end

function read_images!(ntimepoints, arr, file_path, files)
    @inbounds for t in 1:ntimepoints
        @views file = files[t]
        arr[:, :, t] = load(joinpath(file_path, file))
    end
    return nothing 
end

function mean_filter!(X, length_scale)
    iX = IntegralArray(X)
    @inbounds for i in CartesianIndex(1,1):CartesianIndex(size(X))
        x, y = i.I
        x_int = x±length_scale
        y_int = y±length_scale
        x_int = Interval(max(leftendpoint(x_int), 1),
                         min(rightendpoint(x_int), size(X)[1]))
        y_int = Interval(max(leftendpoint(y_int), 1),
                         min(rightendpoint(y_int), size(X)[2]))
        X[i] = iX[x_int, y_int]/(IntervalSets.width(x_int)*IntervalSets.width(y_int))
    end
    return nothing
end

function normalize_local_contrast_output(normalized, images, images_copy, blockDiameter, fpMean)
	length_scale = Int((blockDiameter-1)/2)
    if length(size(images)) == 3
        for t in 1:size(images, 3)
            img = images[:,:,t]
            img_copy = images_copy[:,:,t]
            mean_filter!(img_copy, length_scale)
            img = img - img_copy
            img .+= fpMean
            @. img[img < 0.0] = 0.0
            @. img[img > 1.0] = 1.0
            normalized[:,:,t] = img  
        end
    elseif length(size(images)) == 2
        img_copy = images_copy
        mean_filter!(img_copy, length_scale)
        images = images - img_copy
        images .+= fpMean
        @. images[images < 0.0] = 0.0
        @. images[images > 1.0] = 1.0
        normalized .= images
    else
        error("Dimension error in normalize_local_contrast_output")
    end
    return normalized 
end

function normalize_local_contrast(img, img_copy, blockDiameter)
	img = 1 .- img
	img_copy = 1 .- img_copy
	img_copy = imfilter(img_copy, Kernel.gaussian(blockDiameter))
	img = img - img_copy
    return img 
end

function compute_mask!(stack, masks, fixed_thresh, ntimepoints)
    @inbounds for t in 1:ntimepoints
        @views masks[:,:,t] = stack[:,:,t] .> fixed_thresh
    end
end

function load_processed_display(file_path, filename)
    if filename == "" 
        return nothing, nothing
    else
        if occursin("mask", filename)
            overlay = TiffImages.load(joinpath(file_path, filename))
            processed_filename = replace(filename, "mask" => "")
            processed_image = TiffImages.load(joinpath(file_path, processed_filename))
        else
            processed_image = TiffImages.load(joinpath(file_path, filename))
            overlay_filename = replace(filename, ".tif" => "mask.tif")
            overlay = TiffImages.load(joinpath(file_path, overlay_filename))
        end
    end
    return processed_image, overlay 
end

function mask_overlay(stack, masks, overlay)
    normalized = similar(stack)
	fpMax = maximum(stack)
	fpMin = minimum(stack)
	fpMean = (fpMax - fpMin) / 2.0 + fpMin
	normalized = normalize_local_contrast_output(normalized, stack, copy(stack), 101, fpMean)
	normalized = Gray{N0f8}.(normalized)
    @inbounds for i in CartesianIndices(normalized)
        gray_val = RGB{N0f8}(normalized[i], normalized[i], normalized[i])
        overlay[i] = masks[i] ? RGB{N0f8}(0,1,1) : gray_val
    end
    return normalized, overlay
end

function preprocess_noreg!(img_stack, normalized_stack, blockDiameter, sig)       
    @inbounds for t in 1:size(img_stack, 3)
        img = img_stack[:,:,t]
        img_copy = img_stack[:,:,t] 
        img_normalized = normalize_local_contrast(img, img_copy, 
                                    blockDiameter)
        normalized_stack[:,:,t] = imfilter(img_normalized, Kernel.gaussian(sig))
    end
end

function timelapse_test(file_path, filenames, fixed_thresh, blockDiameter, sig)
    if typeof(filenames) == String 
        @views stack = TiffImages.load(joinpath(file_path, filenames))
        height, width, ntimepoints = size(stack)
    else
        filenames = sort(filenames, lt=natural)
        @views trial_image = TiffImages.load(joinpath(file_path, filenames[1]); mmap=true)
        height, width = size(trial_image)
        ntimepoints = length(filenames)
        stack = Array{eltype(trial_image)}(undef, height, width, ntimepoints)
        read_images!(ntimepoints, stack, file_path, filenames)
    end
    stack = Float64.(stack)
    normalized_stack = similar(stack)
    preprocess_noreg!(stack, normalized_stack, blockDiameter, sig)
    masks = zeros(Bool, size(stack))
    compute_mask!(normalized_stack, masks, fixed_thresh, ntimepoints)
    overlay = zeros(RGB{N0f8}, size(stack)...)
    return mask_overlay(stack, masks, overlay)
end

function image_test(file_path, filename, fixed_thresh, blockDiameter, sig)
    image = TiffImages.load(joinpath(file_path, filename))
    image_copy = copy(image)
    image_normalized = normalize_local_contrast(image, image_copy, blockDiameter)
    image_normalized = imfilter(image_normalized, Kernel.gaussian(sig))
    mask = image_normalized .> fixed_thresh
    overlay = zeros(RGB{N0f8}, size(image)...)
    return mask_overlay(Float64.(image), mask, overlay)  
end

end # module
