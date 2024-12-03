round_odd(x) = div(x, 2) * 2 + 1
compmax(x) = length(x) > 1 ? maximum(x[1:end]) : 0

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

function read_images!(ntimepoints, arr, files)
    @inbounds for t in 1:ntimepoints
        @views file = files[t]
        arr[:, :, t] = load(file)
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

function output_images!(stack, masks, overlay, dir, filename)
	flat_stack = vec(stack)
    img_min = quantile(flat_stack, 0.0035)
    img_max = quantile(flat_stack, 0.9965)
    adjust_histogram!(stack, LinearStretching(src_minval=img_min, src_maxval=img_max, dst_minval=0, dst_maxval=1))
	stack = Gray{N0f8}.(stack)
    save("$dir/Processed images/$filename.tif", stack)
    @inbounds for i in CartesianIndices(stack)
        gray_val = RGB{N0f8}(stack[i], stack[i], stack[i])
        overlay[i] = masks[i] ? RGB{N0f8}(0,1,1) : gray_val
    end
    save("$dir/Processed images/$filename"*"mask.tif", overlay)
end

function extract_base_and_ext(filename::String)
    matches = match(r"^(.*\D)\d+(\.[^\.]+)$", filename)
    if isnothing(matches)
        error("Filename $filename does not satisfy batch criteria")
    end
    base = matches.captures[1]
    ext = matches.captures[2]
    return base, ext
end

function timelapse_processing(images, blockDiameter, ntimepoints, shift_thresh, fixed_thresh, sig, dust_correction, dir, filename, Imin)
    images = Float64.(images)
    normalized_stack = similar(images)
    registered_stack = similar(images)
    images, output_stack, Imin = stack_preprocess(images, normalized_stack, registered_stack, blockDiameter, ntimepoints, shift_thresh, sig, Imin)
    masks = zeros(Bool, size(images))
    compute_mask!(output_stack, masks, fixed_thresh, ntimepoints)
    if dust_correction == "True"
        dust_correct!(masks)
    end
    overlay = zeros(RGB{N0f8}, size(output_stack)...)
    output_images!(images, masks, overlay, dir, filename)
    biomasses = zeros(Float64, ntimepoints)
    if Imin != nothing
        for t in 1:ntimepoints
            @inbounds biomasses[t] = @views mean((-1 .* log10.((images[:,:,t] .- Imin) ./ (images[:,:,1] .- Imin))) .* masks[:,:,t])
        end
    else
        for t in 1:ntimepoints
            @inbounds biomasses[t] = @views mean((1 .- images[:,:,t]) .* masks[:,:,t])
        end
    end
    return biomasses
end

function image_processing(image, blockDiameter, fixed_thresh, sig, dir, filename, Imin, Imax)
    ntimepoints = 1
    image = Float64.(image)
    img_copy = image 
    img_normalized = normalize_local_contrast(image, img_copy, blockDiameter)
    normalized_blurred = imfilter(img_normalized, Kernel.gaussian(sig))
    mask = normalized_blurred .> fixed_thresh
    overlay = zeros(RGB{N0f8}, size(image)...)
    output_images!(image, mask, overlay, dir, filename)
    if Imin != nothing && Imax != nothing
        biomass = mean((-1 .* log10.((image .- Imin) ./ (Imax .- Imin))) .* mask)
        return biomass
    else
        biomass = mean((1 .- image) .* mask)
        return biomass
    end
end

function analysis_main()
    config = parsefile("experiment_config.json")
    images_directories  = config["images_directory"]
    dust_correction = config["dust_correction"]
    batch = config["batch_processing"]
    fixed_thresh = config["fixed_thresh"] 
    Imin_path = config["Imin_path"]
    Imax_path = config["Imax_path"]
    sig = 2
    blockDiameter = 101 
    shift_thresh = 50

    Imin = nothing
    if Imin_path != ""
        Imin = load(Imin_path)
    end
    
    Imax = nothing
    if Imax_path != ""
        Imax = load(Imax_path)
    end

    @inbounds for k in eachindex(images_directories)
        dir = images_directories[k]
        if isdir("$dir/Processed images")
            rm("$dir/Processed images"; recursive = true)
        end
        if isdir("$dir/Numerical data")
            rm("$dir/Numerical data"; recursive = true)
        end
        mkdir("$dir/Processed images")
        mkdir("$dir/Numerical data")
        BF_output_file = "$dir/Numerical data/biomass.csv"

        analyzed = []
        biomass_data = []
        columns = []
		files = [f for f in readdir(dir, join=true) if occursin(r"\.tif$", f)]
        for file in files
            if file ∉ analyzed
                test_image = load(file; lazyio=true)
                img_dims = size(test_image)
                if length(filter(x -> x != 1, img_dims)) == 3
                    height, width, ntimepoints = img_dims
                    images = load(file)
                    target_base, target_ext = extract_base_and_ext(file)
                    push!(columns, target_base)
                    push!(biomass_data, timelapse_processing(images, 
                                                             blockDiameter,
                                                             ntimepoints,
                                                             shift_thresh,
                                                             fixed_thresh,
                                                             sig,
                                                             dust_correction,
                                                             dir, target_base, Imin))
                    push!(analyzed, file)
                elseif length(filter(x -> x != 1, img_dims)) == 2
                    target_base, target_ext = extract_base_and_ext(file)
                    matching_files = filter(file_name -> begin
                        base, ext = extract_base_and_ext(file_name)
                        base == target_base && ext == target_ext
                    end, files)
                    push!(columns, target_base)
                    if length(matching_files) > 1 && batch == "True"
                        timelapse_files = sort(matching_files, lt=natural)
                        ntimepoints = length(timelapse_files)
                        height, width = img_dims
                        images = Array{eltype(test_image), 3}(undef, height, width, ntimepoints)
                        read_images!(ntimepoints, images, timelapse_files)
                        push!(biomass_data, timelapse_processing(images, 
                                                                 blockDiameter,
                                                                 ntimepoints,
                                                                 shift_thresh,
                                                                 fixed_thresh,
                                                                 sig,
                                                                 dust_correction,
                                                                 dir, target_base, Imin))
                        for f in matching_files
                            push!(analyzed, f)
                        end
                    else
                        height, width = img_dims
                        image = load(file)
                        push!(biomass_data, image_processing(image, 
                                                             blockDiameter,
                                                             fixed_thresh, sig,
                                                             dir, target_base, Imin, Imax))
                        push!(analyzed, file)
                    end
                else
                    error("Number of image dimensions must be either 2 or 3")
                end
            end
        end
		max_length = maximum(length, biomass_data)
		padded = [vcat(l, fill(Missing, max_length - length(l))) for l in biomass_data]
		df = DataFrame(padded, Symbol.(columns))
        write(BF_output_file, df)
    end # loop over directories
end
