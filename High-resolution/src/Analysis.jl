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
    num_components = maximum(labeled_array)
    chunks = []
    for label in 1:num_components
        region_coords = findall(labeled_array .== label)
        rows = [x[1] for x in region_coords]
        cols = [x[2] for x in region_coords]
        row_min, row_max = minimum(rows), maximum(rows)
        col_min, col_max = minimum(cols), maximum(cols)
        row_len = row_max - row_min + 1
        col_len = col_max - col_min + 1
        if row_len > max_side || col_len > max_side
            for r_start in row_min:100:row_max
                for c_start in col_min:100:col_max
                    r_end = min(r_start + max_side - 1, row_max)
                    c_end = min(c_start + max_side - 1, col_max)
                    push!(chunks, (r_start, r_end, c_start, c_end))
                end
            end
        else
            r_start = max(row_min, 1)
            r_end = min(row_min + min_side - 1, size(array, 1))
            c_start = max(col_min, 1)
            c_end = min(col_min + min_side - 1, size(array, 2))
            push!(chunks, (r_start, r_end, c_start, c_end))
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
