module utils

export load_raw_display, load_processed_display, image_test, timelapse_test, analysis_main

using TiffImages
using NaturalSort
using IntegralArrays
using IntervalSets
using Images 
using StatsBase
using DataFrames
using CSV: write
using CoordinateTransformations
using AbstractFFTs
using Compat
using FFTW

round_odd(x) = div(x, 2) * 2 + 1
compmax(x) = length(x) > 1 ? maximum(x[1:end]) : 0

function load_raw_display(file_path, filenames)
    ntimepoints = length(filenames)
    if ntimepoints == 0
        return nothing
    elseif ntimepoints == 1
        @views img = TiffImages.load(joinpath(file_path, filenames[1]))
    else
        @views dummy_img = TiffImages.load(joinpath(file_path, filenames[1]); lazyio=true)
        height, width = size(dummy_img)
        img = Array{eltype(dummy_img), 3}(undef, height, width, ntimepoints)
        filenames = sort(filenames, lt=natural)
        for i in 1:ntimepoints
            @views img[:,:,i] = TiffImages.load(joinpath(file_path, filenames[i]))
        end
    end
    return img
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
    filenames = sort(filenames, lt=natural)
    @views trial_image = TiffImages.load(joinpath(file_path, filenames[1]); lazyio=true)
    height, width = size(trial_image)
    ntimepoints = length(filenames)
    stack = Array{eltype(trial_image)}(undef, height, width, ntimepoints)
    read_images!(ntimepoints, stack, file_path, filenames)
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

function phase_offset(source::AbstractArray, target::AbstractArray; kwargs...)
    plan = plan_fft(source)
    return phase_offset(plan, plan * source, plan * target; kwargs...)
end

function phase_offset(
    plan,
    source_freq::AbstractMatrix{<:Complex{T}},
    target_freq;
    upsample_factor = 1,
    normalize = false,
) where {T}
    image_product = @. source_freq * conj(target_freq)
    if normalize
        @. image_product /= max(abs(image_product), eps(T))
    end
    if isone(upsample_factor)
        cross_correlation = ifft!(image_product)
    else
        cross_correlation = plan \ image_product
    end
    maxima, maxidx = @compat findmax(abs, cross_correlation)
    shape = size(source_freq)
    midpoints = map(ax -> (first(ax) + last(ax)) / T(2), axes(source_freq))
    idxoffset = map(first, axes(cross_correlation))
    shift = @. T(ifelse(maxidx.I > midpoints, maxidx.I - shape, maxidx.I) - idxoffset)

    isone(upsample_factor) &&
        return (; shift, calculate_stats(maxima, source_freq, target_freq)...)

    shift = @. round(shift * upsample_factor) / T(upsample_factor)
    upsample_region_size = ceil(upsample_factor * T(1.5))
    dftshift = div(upsample_region_size, 2)
    sample_region_offset = @. dftshift - shift * upsample_factor
    cross_correlation = upsampled_dft(
        image_product,
        upsample_region_size,
        upsample_factor,
        sample_region_offset,
    )
    maxima, maxidx = @compat findmax(abs, cross_correlation)
    shift = @. shift + (maxidx.I - dftshift - idxoffset) / T(upsample_factor)

    stats = calculate_stats(maxima, source_freq, target_freq)
    return (; shift, stats...)
end

function upsampled_dft(
    data::AbstractMatrix{T},
    region_size,
    upsample_factor,
    offsets,
) where {T<:Complex}
    shiftrange = 1:region_size
    idxoffset = map(first, axes(data))
    sample_rate = inv(T(upsample_factor))
    freqs = fftfreq(size(data, 2), sample_rate)
    kernel = @. cis(-T(2π) * (shiftrange - offsets[2] - idxoffset[2]) * freqs')

    _data = kernel * data'

    freqs = fftfreq(size(data, 1), sample_rate)
    kernel = @. cis(T(2π) * (shiftrange - offsets[1] - idxoffset[1]) * freqs')
    _data = kernel * _data'
    return _data
end

function calculate_stats(crosscor_maxima, source_freq, target_freq)
    source_amp = mean(abs2, source_freq)
    target_amp = mean(abs2, target_freq)
    error = 1 - abs2(crosscor_maxima) / (source_amp * target_amp)
    phasediff = atan(imag(crosscor_maxima), real(crosscor_maxima))
    return (; error, phasediff)
end

function read_images!(file_path, ntimepoints, arr, files)
    @inbounds for t in 1:ntimepoints
        @views file = files[t]
        arr[:, :, t] = TiffImages.load(joinpath(file_path, file))
    end
    return nothing 
end

function compute_mask!(stack, masks, fixed_thresh, ntimepoints)
    @inbounds for t in 1:ntimepoints
        @views masks[:,:,t] = stack[:,:,t] .> fixed_thresh
    end
end

function dust_correct!(masks)
    rows, cols, nframes = size(masks)
    for col in 1:cols
        for row in 1:rows
            if masks[row, col, 1]  
                should_suppress = false
                for t in 2:nframes
                    if !masks[row, col, t]
                        should_suppress = true
                        break
                    end
                end
                if should_suppress
					masks[row, col, :] .= false
                end
            end
        end
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

function crop(img_stack)
    @views mask = .!any(isnan.(img_stack), dims=3)[:,:,1]
    @views mask_i = any(mask, dims=2)[:,1]
    @views mask_j = any(mask, dims=1)[1,:]
    i1 = findfirst(mask_i)
    i2 = findlast(mask_i)
    j1 = findfirst(mask_j)
    j2 = findlast(mask_j)
    cropped_stack = img_stack[i1:i2, j1:j2, :]
    return cropped_stack, (i1, i2, j1, j2)
end

function stack_preprocess(img_stack, normalized_stack, registered_stack, blockDiameter, nframes, mxshift, sig, Imin)       
    shifts = (0.0, 0.0) 
    @inbounds for t in 1:nframes
        img = img_stack[:,:,t]
        img_copy = img_stack[:,:,t] 
        img_normalized = normalize_local_contrast(img, img_copy, 
                                    blockDiameter)
        normalized_stack[:,:,t] = imfilter(img_normalized, Kernel.gaussian(sig))
        if t == 1
            registered_stack[:,:,t] = normalized_stack[:,:,t]
        else
            moving = normalized_stack[:,:,t]
            fixed = normalized_stack[:,:,t-1]
            shift, _, _ = phase_offset(fixed, moving, upsample_factor=1)
            if sqrt(shift[1]^2 + shift[2]^2) >= mxshift
                shift = Translation(shifts[1], shifts[2])
                registered_stack[:,:,t] = warp(moving, shift, axes(fixed))
                img_stack[:,:,t] = warp(img_stack[:,:,t], shift, axes(fixed))
            else
                shift = Tuple([-1*shift[1], -1*shift[2]])
                shift = shift .+ shifts
                shifts = shift
                shift = Translation(shift[1], shift[2])
                registered_stack[:,:,t] = warp(moving, shift, axes(fixed))
                img_stack[:,:,t] = warp(img_stack[:,:,t], shift, axes(fixed))
            end
        end
    end
    processed_stack, crop_indices = crop(registered_stack)
    row_min, row_max, col_min, col_max = crop_indices
    img_stack = img_stack[row_min:row_max, col_min:col_max, :]
    if Imin != nothing
        Imin = Imin[row_min:row_max, col_min:col_max]
    end
    return img_stack, processed_stack, Imin
end

function write_OD_images!(OD_images, output_dir, filename)
    if length(size(OD_images)) == 3 
        for t in 1:size(OD_images, 3)
            TiffImages.save("$output_dir/$filename"*"_OD"*"_t"*string(t)*".tif", 
                        OD_images[:,:,t])
        end
    elseif length(size(OD_images)) == 2
        TiffImages.save("$output_dir/$filename"*"_OD.tif", 
                    OD_images)
    end
end

function output_images!(stack, masks, overlay, OD_images, output_dir, filename)
    filename = split(filename, "/")[end]
    normalized = similar(stack)
	fpMax = maximum(stack)
	fpMin = minimum(stack)
	fpMean = (fpMax - fpMin) / 2.0 + fpMin
	normalized = normalize_local_contrast_output(normalized, stack, copy(stack), 101, fpMean)
	normalized = Gray{N0f8}.(normalized)
    TiffImages.save("$output_dir/$filename.tif", normalized)
    @inbounds for i in CartesianIndices(normalized)
        gray_val = RGB{N0f8}(normalized[i], normalized[i], normalized[i])
        overlay[i] = masks[i] ? RGB{N0f8}(0,1,1) : gray_val
    end
    TiffImages.save("$output_dir/$filename"*"mask.tif", overlay)
    if OD_images != nothing
        write_OD_images!(OD_images, output_dir, filename)
    end
end

function extract_base_and_ext(filename::String)
    # The batch criterion is that there is a number at the end of the filename
    matches = match(r"^(.*\D)\d+(\.[^\.]+)$", filename)
    if isnothing(matches)
        error("Filename $filename does not satisfy batch criteria")
    end
    base = matches.captures[1]
    ext = matches.captures[2]
    return base, ext
end

function timelapse_processing(images, blockDiameter, ntimepoints, shift_thresh, fixed_thresh, sig, dust_correction, output_dir, filename, Imin, Imax)
    images = Float64.(images)
    normalized_stack = similar(images)
    registered_stack = similar(images)
    images, output_stack, Imin = stack_preprocess(images, normalized_stack, registered_stack, blockDiameter, ntimepoints, shift_thresh, sig, Imin)
    masks = zeros(Bool, size(images))
    compute_mask!(output_stack, masks, fixed_thresh, ntimepoints)
    if dust_correction == true 
        dust_correct!(masks)
    end
    overlay = zeros(RGB{N0f8}, size(output_stack)...)
    biomasses = zeros(Float64, ntimepoints)
    OD_images = nothing
    if Imin != nothing
        OD_images = Array{Gray{Float32}, 3}(undef, size(images))
        for t in 1:ntimepoints
            if Imax != nothing
                OD_images[:,:,t] = @views (-1 .* log10.((images[:,:,t] .- Imin) ./ (Imax .- Imin)))
            else
                OD_images[:,:,t] = @views (-1 .* log10.((images[:,:,t] .- Imin) ./ (images[:,:,1] .- Imin)))
            end
            @inbounds biomasses[t] = @views mean(OD_images[:,:,t] .* masks[:,:,t])
        end
    else
        for t in 1:ntimepoints
            @inbounds biomasses[t] = @views mean((1 .- images[:,:,t]) .* masks[:,:,t])
        end
    end
    output_images!(images, masks, overlay, OD_images, output_dir, filename)
    return biomasses
end

function image_processing(image, blockDiameter, fixed_thresh, sig, output_dir, filename, Imin, Imax)
    ntimepoints = 1
    image = Float64.(image)
    img_copy = image 
    img_normalized = normalize_local_contrast(image, img_copy, blockDiameter)
    normalized_blurred = imfilter(img_normalized, Kernel.gaussian(sig))
    mask = normalized_blurred .> fixed_thresh
    overlay = zeros(RGB{N0f8}, size(image)...)
    OD_image = nothing
    biomass = nothing
    if Imin != nothing && Imax != nothing
        OD_image = Array{Gray{Float32}, 2}(undef, size(image))
        OD_image .= (-1 .* log10.((image .- Imin) ./ (Imax .- Imin)))
        biomass = mean(OD_image .* mask)
    else
        biomass = mean((1 .- image) .* mask)
    end
    output_images!(image, mask, overlay, OD_image, output_dir, filename)
    return biomass
end

function analysis_main(file_path, output_path, files, dust_correction, batch_processing, 
        fixed_thresh, Imin_path, Imax_path)
    sig = 2
    blockDiameter = 101 
    shift_thresh = 50

    Imin = nothing
    if Imin_path != ""
        Imin = TiffImages.load(joinpath(file_path, Imin_path))
    end
    
    Imax = nothing
    if Imax_path != ""
        Imax = TiffImages.load(joinpath(file_path, Imax_path))
    end

    BF_output_file = joinpath(output_path, "biomass.csv")
    analyzed = []
    biomass_data = []
    columns = []
    files = [f for f in files if f != Imin_path && f != Imax_path]
    for file in files
        if file ∉ analyzed
            test_image = TiffImages.load(joinpath(file_path, file); lazyio=true)
            img_dims = size(test_image)
            if length(filter(x -> x != 1, img_dims)) == 3
                height, width, ntimepoints = img_dims
                images = TiffImages.load(joinpath(file_path, file))
                target_base, target_ext = extract_base_and_ext(file)
                push!(columns, target_base)
                push!(biomass_data, timelapse_processing(images, 
                                                         blockDiameter,
                                                         ntimepoints,
                                                         shift_thresh,
                                                         fixed_thresh,
                                                         sig,
                                                         dust_correction,
                                                         output_path, target_base, Imin, Imax))
                push!(analyzed, file)
            elseif length(filter(x -> x != 1, img_dims)) == 2
                target_base, target_ext = extract_base_and_ext(file)
                matching_files = filter(file_name -> begin
                    base, ext = extract_base_and_ext(file_name)
                    base == target_base && ext == target_ext
                end, files)
                push!(columns, target_base)
                if length(matching_files) > 1 && batch_processing == true
                    timelapse_files = sort(matching_files, lt=natural)
                    ntimepoints = length(timelapse_files)
                    height, width = img_dims
                    images = Array{eltype(test_image), 3}(undef, height, width, ntimepoints)
                    read_images!(file_path, ntimepoints, images, timelapse_files)
                    push!(biomass_data, timelapse_processing(images, 
                                                             blockDiameter,
                                                             ntimepoints,
                                                             shift_thresh,
                                                             fixed_thresh,
                                                             sig,
                                                             dust_correction,
                                                             output_path, target_base, Imin, Imax))
                    for f in matching_files
                        push!(analyzed, f)
                    end
                else
                    height, width = img_dims
                    image = TiffImages.load(joinpath(file_path, file))
                    push!(biomass_data, image_processing(image, 
                                                         blockDiameter,
                                                         fixed_thresh, sig,
                                                         output_path, target_base, Imin, Imax))
                    push!(analyzed, file)
                end
            else
                error("Number of image dimensions must be either 2 or 3")
            end
        end
		max_length = maximum(length, biomass_data)
		padded = [vcat(l, fill(Missing, max_length - length(l))) for l in biomass_data]
		df = DataFrame(padded, Symbol.(columns))
        write(BF_output_file, df)
    end # loop over directories
end

end # module
