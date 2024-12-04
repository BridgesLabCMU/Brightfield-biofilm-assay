function display_images!(stack, masks, overlay)
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
    imshow(overlay)
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

function timelapse_test!(test_images, fixed_thresh, blockDiameter, sig)
    trial_image = load(test_images[1]; lazyio=true)
    height, width = size(trial_image)
    ntimepoints = length(test_images)
    stack = Array{eltype(trial_image)}(undef, height, width, ntimepoints)
    read_images!(ntimepoints, stack, test_images)
    stack = Float64.(stack)
    normalized_stack = similar(stack)
    preprocess_noreg!(stack, normalized_stack, blockDiameter, sig)
    masks = zeros(Bool, size(stack))
    compute_mask!(normalized_stack, masks, fixed_thresh, ntimepoints)
    overlay = zeros(RGB{N0f8}, size(stack)...)
    display_images!(stack, masks, overlay)
    return nothing
end

function image_test!(image_path, fixed_thresh, blockDiameter, sig)
    image = load(image_path)
    image_copy = copy(image)
    image_normalized = normalize_local_contrast(image, image_copy, blockDiameter)
    image_normalized = imfilter(image_normalized, Kernel.gaussian(sig))
    mask = image_normalized .> fixed_thresh
    overlay = zeros(RGB{N0f8}, size(image)...)
    display_images!(Float64.(image), mask, overlay)
    return nothing
end

function GUI_main()
    win = GtkWindow("Experiment Configuration", 200, 300)
    vbox = GtkBox(:v)
    push!(win, vbox)
    
    directories = []
    test_images = []
    Imin = ""
    Imax = ""
    fixed_thresh = 0.03 

    select_dirs_button = GtkButton("Select experiment folders")
    push!(vbox, select_dirs_button)
    
    select_Imin_button = GtkButton("Select Imin file")
    push!(vbox, select_Imin_button)
    
    select_Imax_button = GtkButton("Select Imax file (for single image analysis)")
    push!(vbox, select_Imax_button)
    
    adj = GtkAdjustment(0.03, 0., 1., 0.001, 1.0, 0.0) 
    spin_button = GtkSpinButton(adj, 0.001, 3)
	Gtk4.value(spin_button, 0.03)
    push!(vbox, GtkLabel("Threshold:"))
    push!(vbox, spin_button)
    
    select_test_images_button = GtkButton("Select one tif or a tif series to test the threshold")
    push!(vbox, select_test_images_button)

    function select_directories(button)
        dlg = GtkFileChooserDialog(
            "Select Folders",
            win,
            Gtk4.FileChooserAction_SELECT_FOLDER,
                                (("_Cancel", Gtk4.ResponseType_CANCEL),
                                 ("_Open", Gtk4.ResponseType_ACCEPT)))
        
        dlgp = GtkFileChooser(dlg)
        Gtk4.G_.set_select_multiple(dlgp, true)
        signal_connect(dlg, "response") do widget, response_id
            if response_id == Gtk4.ResponseType_ACCEPT
                filename_list = Gtk4.G_.get_files(dlgp)
                sel = String[Gtk4.GLib.G_.get_path(Gtk4.GFile(f)) for f in Gtk4.GListModel(filename_list)]
                push!(directories, sel...)
                destroy(dlg)
            end
        end
        show(dlg)
    end

    signal_connect(select_dirs_button, "clicked") do widget
        select_directories(widget)
    end

    function select_Imin(button)
        dlg = GtkFileChooserDialog(
            "Select Imin file",
            win,
            Gtk4.FileChooserAction_OPEN,
                                (("_Cancel", Gtk4.ResponseType_CANCEL),
                                 ("_Open", Gtk4.ResponseType_ACCEPT)))
        
        dlgp = GtkFileChooser(dlg)
        Gtk4.G_.set_select_multiple(dlgp, false)
        signal_connect(dlg, "response") do widget, response_id
            if response_id == Gtk4.ResponseType_ACCEPT
                filename = Gtk4.G_.get_file(dlgp)
                sel = Gtk4.GLib.G_.get_path(Gtk4.GFile(filename))
                Imin = sel
                destroy(dlg)
            end
        end
        show(dlg)
    end

    signal_connect(select_Imin_button, "clicked") do widget
        select_Imin(widget)
    end

    function select_Imax(button)
        dlg = GtkFileChooserDialog(
            "Select Imax file",
            win,
            Gtk4.FileChooserAction_OPEN,
                                (("_Cancel", Gtk4.ResponseType_CANCEL),
                                 ("_Open", Gtk4.ResponseType_ACCEPT)))
        
        dlgp = GtkFileChooser(dlg)
        Gtk4.G_.set_select_multiple(dlgp, false)
        signal_connect(dlg, "response") do widget, response_id
            if response_id == Gtk4.ResponseType_ACCEPT
                filename = Gtk4.G_.get_file(dlgp)
                sel = Gtk4.GLib.G_.get_path(Gtk4.GFile(filename))
                Imax = sel
                destroy(dlg)
            end
        end
        show(dlg)
    end

    signal_connect(select_Imax_button, "clicked") do widget
        select_Imax(widget)
    end

    function select_test_images(button)
        dlg = GtkFileChooserDialog(
            "Select test image(s)",
            win,
            Gtk4.FileChooserAction_OPEN,
                                (("_Cancel", Gtk4.ResponseType_CANCEL),
                                 ("_Open", Gtk4.ResponseType_ACCEPT)))
        
        dlgp = GtkFileChooser(dlg)
        Gtk4.G_.set_select_multiple(dlgp, true)
        signal_connect(dlg, "response") do widget, response_id
            if response_id == Gtk4.ResponseType_ACCEPT
                filename_list = Gtk4.G_.get_files(dlgp)
                sel = String[Gtk4.GLib.G_.get_path(Gtk4.GFile(f)) for f in Gtk4.GListModel(filename_list)]
                push!(test_images, sel...)
                destroy(dlg)
            end
        end
        show(dlg)
    end
    
    signal_connect(select_test_images_button, "clicked") do widget
        select_test_images(widget)
    end
    
    test_button = GtkButton("Test")
    push!(vbox, test_button)

    function on_test(button)
        fixed_thresh = Gtk4.value(spin_button)
        ntimepoints = length(test_images)
        if ntimepoints == 1
            image_test!(test_images[1], fixed_thresh, 101, 2)
        else
            timelapse_test!(sort(test_images, lt=natural), fixed_thresh, 101, 2)
        end
        empty!(test_images)
    end
    
    signal_connect(test_button, "clicked") do widget
        on_test(widget)
    end

    dust_correction_checkbox = GtkCheckButton("Dust correction")
    batch_processing_checkbox = GtkCheckButton("Batch processing")
    push!(vbox, dust_correction_checkbox)
    push!(vbox, batch_processing_checkbox)

    done_button = GtkButton("Done")
    push!(vbox, done_button)

    function on_done(button)
        # Read value from spin button
        fixed_thresh = spin_button.value
        
        if dust_correction_checkbox.active == true
            dust_correction = "True" 
        else
            dust_correction = "False" 
        end
        if batch_processing_checkbox.active == true
            batch_processing = "True" 
        else
            batch_processing = "False" 
        end
        config = Dict(
            "images_directory" => directories,
            "Imin_path" => Imin,
            "Imax_path" => Imax,
            "fixed_thresh" => fixed_thresh,
            "dust_correction" => dust_correction,
            "batch_processing" => batch_processing 
        )
        
        open("experiment_config.json", "w") do f
            JSON.print(f, config)
        end
    end

    if !isinteractive()
        c = Condition()
        signal_connect(done_button, "clicked") do widget
            on_done(widget)
            notify(c)
        end
        @async Gtk4.GLib.glib_main()
        wait(c)
    else
        signal_connect(done_button, "clicked") do widget
            on_done(widget)
            close(win)
        end
    end

    show(win)
end
