module App
using Main.Analysis.jl
using GenieFramework, PlotlyBase
@genietools

@app begin
    #reactive code goes here
end

function ui()
    # 1. Select directories where images are stored
    # 2. Select Imin file and Imax file
    # 3. Dust correction and batch processing checkboxes
    # 4. fixed thresh slider
    # 5. Test button that will run a test using the fixed thresh
    row([
    cell(class="st-col col-3", [
        h1("A simple dashboard"),
        slider(1:1:1000, :N),
        p("The average of {{N}} random numbers is {{m}}", class="st-module"),
        ])
    ])
end

@page("/", ui)
end
