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
		  accept = ".zip",
		  maxfilesize = 10000*10000*1000, # bytes
		  maxfiles = 38400, # 384-well plate, 100 timepoints
		  autoupload = true,
		  hideuploadbtn = true,
          label = "Upload image folders (zipped)", class = "tex_text",
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
           h6("Select folder (to set Imin/Imax, display images, and/or test threshold)"),
           Stipple.select(:analyze_folder; options=:zipped_files, style="margin-bottom: 20px;",
                         class="tex_text", clearable=true),
            btn(
                "Confirm",
                @click(:ConfirmButton_process),
                loading = :ConfirmButton_process,
                color = "blue", class="tex_text", style="margin-bottom: 20px;"
            )
	  ])
])

row([
	cell(
          class="st-module col-6",
          [
           h6("Imin file (optional)"),
           Stipple.select(:selected_Imin; options=:analyze_folder_files, style="margin-bottom: 20px;", class = "tex_text", clearable=true),
           h6("Imax file (optional)"),
           Stipple.select(:selected_Imax; options=:analyze_folder_files, style="margin-bottom: 40px;", class="tex_text", clearable=true),
            h6("Choose file(s) to display (optional)"),
           Stipple.select(:selected_raw_display_files, options=:analyze_folder_files, clearable=true, hideselected = true, multiple = true),
            btn(
                "Display raw image(s)",
                @click(:DisplayRawButtonProgress_process),
                loading = :DisplayRawButtonProgress_process,
                percentage = :DisplayRawButtonProgress_progress,
                color = "blue", class="tex_text", style="margin-bottom: 40px;"
            ),

            h6("Set fixed threshold for masking"),
            slider(0.0000:0.0001:1.0000, :fixed_thresh),
            p("Fixed threshold = {{fixed_thresh}}", style="margin-bottom: 20px;"),
            h6("Set block diameter for local contrast"),
            slider(100:1:500, :block_diameter),
            p("Block diameter = {{block_diameter}}", style="margin-bottom: 20px;"),
            h6("Choose file(s) to test threshold (optional)"),
           Stipple.select(:selected_test_files, options=:analyze_folder_files, clearable=true, hideselected = true, multiple = true),
            btn(
                "Test threshold",
                @click(:TestButtonProgress_process),
                loading = :TestButtonProgress_process,
                percentage = :TestButtonProgress_progress,
                color = "green", class="tex_text", style="margin-bottom: 40px;"
            ),
            checkbox("Perform dust correction", :dust_correction, class="tex_text"),
            checkbox("Timelapse format", :timelapse_flag, class="tex_text"),
            checkbox("Analyze all folders", :batch_processing, class="tex_text", style="margin-bottom: 40px;"),
            btn(
                "Run analysis",
                @click(:AnalyzeButtonProgress_process),
                loading = :AnalyzeButtonProgress_process,
                percentage = :AnalyzeButtonProgress_progress,
                color = "orange", class="tex_text", style="margin-bottom: 40px;"
            ),
            h6("Choose processed file to display (optional)"),
           Stipple.select(:selected_processed_image, options=:processed_images, clearable=true, multiple = false, style="margin-bottom: 20px"),
            btn(
                "Display processed image(s)",
                @click(:DisplayProcessedButtonProgress_process),
                loading = :DisplayProcessedButtonProgress_process,
                percentage = :DisplayProcessedButtonProgress_progress,
                color = "blue", class="tex_text", style="margin-bottom: 40px"
            ), 
            btn(
                "Download processed image(s) and numerical data",
                href = "/downloads.zip",
                color = "green", class="tex_text", style="margin-bottom: 40px"
            ), 
            btn(
                "Done session/cleanup",
                @click(:CleanupButtonProgress_process),
                loading = :CleanupButtonProgress_process,
                percentage = :CleanupButtonProgress_progress,
                color = "red", class="tex_text", style="margin-bottom: 40px"
            )
          ]
         )
	cell(
          class="st-module col-6", [
			tabgroup(
				:tab_selected,
				inlinelabel = true,
				class = "tex_text",
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
