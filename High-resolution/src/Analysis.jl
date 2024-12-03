using Images
using HistogramThresholding
using ImageView
using FLoops
using PythonCall
nd2 = pyimport("nd2")

function cell_mask!(image, mask, slices, cell_threshold)
    # Auto intensity threshold each z-slice independently
    @floop for i in 1:slices 
        @views slice_ = image[:, :, i]
        otsu_thresh = find_threshold(slice_, Otsu())
        if otsu_thresh/3 < cell_threshold 
            mask[:,:,i] = slice_ .> cell_threshold 
        else
            mask[:,:,i] = slice_ .> otsu_thresh/3
        end
    end
    imshow(mask)
    return nothing 
end

function crop_biofilms!(image, mask)
    # Modify the cell mask to isolate biofilm
    zstack = sum(image, dims=3)
    otsu_thresh = find_threshold(zstack, Otsu())
    crop_mask = zstack .> otsu_thresh 
    for i in 1:size(mask, 3)
        mask[:,:,i] .= mask[:,:,i] .& crop_mask
    end
end

function main()
    path = "/Volumes/T7 Shield/Brightfield_paper/"
    cell_threshold = 3.0e-5
    biovolumes = []
    for filename in readdir(path, join=true)
        image = pyconvert(Array, nd2.imread(filename))
        height, width, slices = size(image)
        mask = Array{Bool}(undef, height, width, slices)
        cell_mask!(image, mask, slices, cell_threshold)
        crop_biofilms!(image, mask)
        biovolume = sum(biofilm_image)
        push!(biovolumes, biovolume)
    end
end

main()
