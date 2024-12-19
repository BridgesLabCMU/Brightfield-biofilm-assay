# 1. Select directories where images are stored
# 2. Select Imin file and Imax file
# 3. Dust correction and batch processing checkboxes
# 4. fixed thresh slider
# 5. Test button that will run a test using the fixed thresh
heading("Label-free biofilm analysis")
row([
  cell(class="col-md-12", [
      uploader( multiple = true,
          accept = ".csv",
          maxfilesize = 1024*1024*1, # bytes
          maxfiles = 3,
          autoupload = true,
          hideuploadbtn = true,
          label = "Upload datasets",
          nothumbnails = true,
          style="max-width: 95%; width: 95%; margin: 0 auto;",

          @on("rejected", :rejected),
          @on("uploaded", :uploaded)
          )
      ])
 ])

row([     
     cell(class="st-module", [
            checkbox("Perform dust correction", :dust_correction),
          ])
     cell(class="st-module",[
            checkbox("Perform bulk processing", :bulk_processing)
          ]) 
])

row([
    cell(class="st-col col-3", [
        h1("Biofilm threshold"),
        slider(0:1, :N),
    ])
])

row([
    cell(class="st-col col-3", [
        btn(
            "Test threshold",
            @click(:ButtonProgress_process),
            loading = :ButtonProgress_process,
            percentage = :ButtonProgress_progress,
            color = "green",
        )
    ])
])
