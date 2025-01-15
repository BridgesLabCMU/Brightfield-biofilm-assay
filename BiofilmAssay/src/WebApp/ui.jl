row([
     cell(
          class="st-header q-pa-sm",
          style="max-width: 99%; width: 99%; margin: 0 auto;",
          [
    h1("Label-free biofilm analysis", class="my-4 text-center"),
   ])
])

row([
  cell(class="col-md-12 q-pa-sm", 
	  [
	  uploader( multiple = true,
		  accept = ".tif, .tiff, .ome.tif, .ome.tiff",
		  maxfilesize = 10000*10000*1000, # bytes
		  maxfiles = 38400, # 384-well plate, 100 timepoints
		  autoupload = true,
		  hideuploadbtn = true,
		  label = "Upload image files",
		  nothumbnails = true,
		  @on("rejected", :rejected),
		  @on("uploaded", :uploaded)
		  )
	  ])
])

row([
	cell(
          class="st-module col-6",
          [
           h6("Imin file"),
           Stipple.select(:selected_Imin; options=:upfiles, style="margin-bottom: 20px;"),
           h6("Imax file"),
           Stipple.select(:selected_Imax; options=:upfiles, style="margin-bottom: 40px;"),
            checkbox("Perform dust correction", :dust_correction),
            checkbox("Perform batch processing", :batch_processing, style="margin-bottom: 40px;"),

            h6("Choose file(s) to display (optional)"),
           Stipple.select(:selected_raw_display_files, options=:upfiles, clearable=true, hideselected = true, multiple = true,
                          counter=true),
            btn(
                "Display raw image(s)",
                @click(:DisplayRawButtonProgress_process),
                loading = :DisplayRawButtonProgress_process,
                percentage = :DisplayRawButtonProgress_progress,
                color = "blue", class="q-my-sm", style="margin-bottom: 40px;"
            ),

            h6("Set fixed threshold for masking"),
            slider(0.0000:0.0001:1.0000, :fixed_thresh),
            p("Fixed threshold = {{fixed_thresh}}", style="margin-bottom: 20px;"),
            h6("Choose file(s) to test threshold (optional)"),
           Stipple.select(:selected_test_files, options=:upfiles, clearable=true, hideselected = true, multiple = true,
                          counter=true),
            btn(
                "Test threshold",
                @click(:TestButtonProgress_process),
                loading = :TestButtonProgress_process,
                percentage = :TestButtonProgress_progress,
                color = "green", class="q-my-sm", style="margin-bottom: 40px;"
            ),
            btn(
                "Run analysis",
                @click(:AnalyzeButtonProgress_process),
                loading = :AnalyzeButtonProgress_process,
                percentage = :AnalyzeButtonProgress_progress,
                color = "orange", style="margin-bottom: 40px;"
            ),
            h6("Choose processed image(s) to display (optional)"),
           Stipple.select(:selected_processed_image, options=:processed_images, clearable=true, multiple = false),
            btn(
                "Display processed image(s)",
                @click(:DisplayProcessedButtonProgress_process),
                loading = :DisplayProcessedButtonProgress_process,
                percentage = :DisplayProcessedButtonProgress_progress,
                color = "blue", class="q-my-sm"
            ), 
            btn(
                "Download processed image(s) and numerical data",
                href = "/data.zip",
                color = "green", class="q-my-sm"
            )
          ]
         )
	cell(
          class="st-module col-6", [
			tabgroup(
				:tab_selected,
				inlinelabel = true,
				class = "bg-primary text-white shadow-2",
				[
					tab(name = "Images", icon = "images", label = "Images"),
					tab(name = "Masks", icon = "masks", label = "Masks"),
				],
			),
			tabpanels(
				:tab_selected,
				animated = true,
				var"transition-prev" = "scale",
				var"transition-next" = "scale",
				[
					tabpanel(name = "Images", [
                        imageview(
                            src = :image_display_path,
                            style = "height: 100%; (max-width: 750px)",
                        )
                        ]),
					tabpanel(name = "Masks", [
                        imageview(
                            src = :mask_display_path,
                            style = "height: 100%; (max-width: 750px)",
                        )
                    ]),
				],
			)]
          )
    ])
