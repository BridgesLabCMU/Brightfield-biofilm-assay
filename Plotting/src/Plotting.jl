############
# To do: fix correlation coefficient calculation

using Makie 
using CairoMakie
using SwarmMakie
using LsqFit
using TiffImages
using CSV
using StatsBase
using DataFrames
using EasyFit
using Colors
using ColorSchemes: colorschemes
using PythonCall
CairoMakie.activate!(type="svg")

odr = pyimport("scipy.odr")

function extract_float(gray_string::String)
    m = match(r"Gray\{Float64\}\(([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\)", gray_string)
    if m !== nothing
        return parse(Float64, m.captures[1])
    else
        error("Invalid Gray{Float64} string format: $gray_string")
    end
end

function fig1B!(vols, vc_mass, pa_pf_mass, pf_wspf_mass, sp_mass)
    # Biovolume (x-axis) vs. brightfield biomass (y-axis)
    # Scatter + best-fit
    # Data are 1D dataframes where filenames/wells are the column names
    
    # Ensure vols are in the correct order
    substrings = ["D39", "36"]
	vols = filter(row -> !any(substr -> occursin(substr, row.first), substrings), vols)
	
	order_list = ["vpsL_1", "vpsL_2", "vpsL_3", 
                  "0ara_1", "0ara_2", "0ara_3", 
                  "0025ara_1", "0025ara_2", "0025ara_3",
                  "0035ara_1", "0035ara_2", "0035ara_3",
                  "005ara_1", "005ara_2", "005ara_3",
                  "01ara_1", "01ara_2", "01ara_3",
                  "Pf_nobiofilm_1", "Pf_nobiofilm_2", "Pf_nobiofilm_3",
                  "Pf_WT_1", "Pf_WT_2", "Pf_WT_3",
                  "Pf_wspF_1", "Pf_wspF_2", "Pf_wspF_3",
                  "Pa_pel_1", "Pa_pel_2", "Pa_pel_3",
                  "Pa_wspF_1", "Pa_wspF_2", "Pa_wspF_3",
                  ] 
	order_map = Dict(substring => idx for (idx, substring) in enumerate(order_list))
	vols.order = [order_map[substring] for row in eachrow(vols) for substring in order_list if occursin(substring, row.first)]
	vols = sort(vols, :order)
	select!(vols, Not(:order))

	function map_well_identifier(file_path)
		well_identifier = match(r"([A-Z])\d+", basename(file_path))[1][1]
		return Dict('A' => 1, 'B' => 2, 'C' => 3, 'D' => 4, 'E' => 5, 'F' => 6, 
                   'G' => 7, 'H' => 8)[well_identifier]
	end

    # Get averages and stds of vc_masses, pf_masses, sp_masses (correct order)
    pa_pf_mass = stack(pa_pf_mass, Not([]), variable_name = :FilePath, value_name = :Value)
    pa_pf_mass[!, :Well] = replace.(pa_pf_mass.FilePath, r".*/(.*?)_.*" => s"\1")
	pa_pf_mass[:, :condition] = map(map_well_identifier, pa_pf_mass.FilePath)
    grouped_values = [group.Value for group in groupby(pa_pf_mass, :condition)]
    pa_pf_averages = [mean(subarray) for subarray in grouped_values] 
    pa_pf_stds = [std(subarray) for subarray in grouped_values]
    pa_stds = pa_pf_stds[1:2]
    pa_averages = pa_pf_averages[1:2]
    pf_stds = [pa_pf_stds[4], pa_pf_stds[3], std(pf_wspf_mass[1,:])]
    pf_averages = [pa_pf_averages[4], pa_pf_averages[3], mean(pf_wspf_mass[1,:])]
    
    vc_mass = stack(vc_mass, Not([]), variable_name = :FilePath, value_name = :Value)
    vc_mass[!, :Well] = replace.(vc_mass.FilePath, r".*/(.*?)_.*" => s"\1")
	vc_mass[:, :condition] = map(map_well_identifier, vc_mass.FilePath)
    grouped_values = [group.Value for group in groupby(vc_mass, :condition)]
    vc_averages = [mean(subarray) for subarray in grouped_values] 
    vc_stds = [std(subarray) for subarray in grouped_values]

@pyexec """
def model(p, x):
    return p[0]*x 
""" => model 

    p0 = [(vc_averages[2]-vc_averages[1])/(vols[4,2]-vols[1,2])]
    lb = [0.0, 0.0]
    ub = [Inf, Inf]

    mass_avg = vcat(vc_averages, pf_averages, pa_averages)
    stds = vcat(vc_stds, pf_stds, pa_stds)
    col2 = vols.second
    vols_avg = [mean(col2[i:i+2]) for i in 1:3:length(col2) if i+2 <= length(col2)]
    vols_std = [std(col2[i:i+2]) for i in 1:3:length(col2) if i+2 <= length(col2)]
    @show fitlinear(mass_avg,vols_avg)

    linear = odr.Model(model)
    odrdata = odr.Data(mass_avg, vols_avg, wd=1/(stds/sqrt(9)), we=1/(vols_std/sqrt(3)))
    myodr = odr.ODR(odrdata, linear, beta0=p0)
    output = myodr.run()
    #fit = curve_fit(model, mass_avg, vols_avg, p0, lower=lb, upper=ub)
    #pstar = coef(fit)
    #@show fitlinear(mass_avg,vols_avg)
	#xbase = collect(range(minimum(vc_averages), maximum(vc_averages), 100))
    fig = Figure(size=(6*72,3*72))
    ax = Axis(fig[1, 1], xlabel="BF-biofilm biomass (a.u.)", ylabel="CF-biofilm biomass (μm³)")
    errorbars!(ax, vc_averages, vols_avg[1:6], vc_stds, color=Makie.wong_colors()[1], direction = :x)
    errorbars!(ax, pf_averages, vols_avg[7:9], pf_stds, color=Makie.wong_colors()[2], direction = :x)
    errorbars!(ax, pa_averages, vols_avg[10:11], pa_stds, color=Makie.wong_colors()[3], direction = :x)
    errorbars!(ax, vc_averages, vols_avg[1:6], vols_std[1:6], color=Makie.wong_colors()[1], direction = :y)
    errorbars!(ax, pf_averages, vols_avg[7:9], vols_std[7:9], color=Makie.wong_colors()[2], direction = :y)
    errorbars!(ax, pa_averages, vols_avg[10:11], vols_std[10:11], color=Makie.wong_colors()[3], direction = :y)
    scatter!(ax, vc_averages, vols_avg[1:6], color=Makie.wong_colors()[1], label=rich("V. cholerae"; font=:italic))
    scatter!(ax, pf_averages, vols_avg[7:9], color=Makie.wong_colors()[2], label=rich("P. fluorescens"; font=:italic))
    scatter!(ax, pa_averages, vols_avg[10:11], color=Makie.wong_colors()[3], label=rich("P. aeruginosa"; font=:italic))
    reg_params = pyconvert(Array, output.beta)
    reg_error = pyconvert(Array, output.sd_beta)
    avg = reg_params[1] .* mass_avg
    lines!(ax, mass_avg, avg, color="black")
    band!(ax, mass_avg, (reg_params[1]-reg_error[1]).*mass_avg, (reg_params[1]+reg_error[1]).*mass_avg, color=(:black, 0.2))
	ax.rightspinevisible = false
	ax.topspinevisible = false
    fig[1,2] = Legend(fig, ax, merge = true, unique = true, framevisible=false, labelsize=12, rowgap=0)
    save("figS1A.svg", fig)
end

function fig1A_OD_image(OD_image_path)
    image = Float64.(TiffImages.load(OD_image_path))
    fig = Figure(figure_padding = 0, fontsize=38)
    ax = Axis(fig[1, 1], aspect = DataAspect(), yticklabelsvisible = false, xticklabelsvisible = false, yticksvisible = false, xticksvisible = false)
    hm = heatmap!(ax, rotr90(image), colormap=:plasma)
    Colorbar(fig[end+1, :], hm, vertical = false, flipaxis = false, label = "OD", width = Relative(0.7))
    save("fig1A_OD_image.png", fig)
end

function fig1A_lineplot!(fig1A_data_path)
    data = DataFrame(CSV.File(fig1A_data_path))
    fig = Figure(size=(3*72, 2.7*72))
    ax = Axis(fig[1, 1], xlabel="Time (h)", ylabel="BF-biofilm biomass (a.u.)")
    columns = select(data, Cols.(contains.("/E"))) 
    avg = reduce(+, eachcol(columns)) ./ ncol(columns) 
    stdev = dropdims(std(Array(columns), dims=2), dims=2)  
    time = 0:0.5:nrow(columns)/2-0.5
    lines!(ax, time, avg, color=:black)
    band!(ax, time, avg-stdev, avg+stdev, color=(:black, 0.2))
    scatter!(ax, time, avg, color=:white, marker=:circle,  strokewidth=1)
	ax.rightspinevisible = false
	ax.topspinevisible = false
    save("fig1A_lineplot.svg", fig)
end

function EVOS_Nikon!(vols, EVOS, Nikon)
    substrings = ["D39", "36", "Pa", "Pf"]
	vols = filter(row -> !any(substr -> occursin(substr, row.first), substrings), vols)
	order_list = ["vpsL_1", "vpsL_2", "vpsL_3", 
                  "0ara_1", "0ara_2", "0ara_3", 
                  "0025ara_1", "0025ara_2", "0025ara_3",
                  "0035ara_1", "0035ara_2", "0035ara_3",
                  "005ara_1", "005ara_2", "005ara_3",
                  "01ara_1", "01ara_2", "01ara_3",
                  "Pf_nobiofilm_1", "Pf_nobiofilm_2", "Pf_nobiofilm_3",
                  "Pf_WT_1", "Pf_WT_2", "Pf_WT_3",
                  "Pf_wspF_1", "Pf_wspF_2", "Pf_wspF_3",
                  "Pa_pel_1", "Pa_pel_2", "Pa_pel_3",
                  "Pa_wspF_1", "Pa_wspF_2", "Pa_wspF_3",
                  ] 
	order_map = Dict(substring => idx for (idx, substring) in enumerate(order_list))
	vols.order = [order_map[substring] for row in eachrow(vols) for substring in order_list if occursin(substring, row.first)]
	vols = sort(vols, :order)
	select!(vols, Not(:order))
	
    order_list = ["vpsL_0001", "vpsL_0002", "vpsL_0003", 
                  "0ara_0001", "0ara_0002", "0ara_0003", 
                  "0025ara_0001", "0025ara_0002", "0025ara_0003",
                  "0035ara_0001", "0035ara_0002", "0035ara_0003",
                  "005ara_0001", "005ara_0002", "005ara_0003",
                  "01ara_0001", "01ara_0002", "01ara_0003",
                  ] 
	order_map = Dict(substring => idx for (idx, substring) in enumerate(order_list))

    EVOS = stack(EVOS, Not([]), variable_name = :FilePath, value_name = :Value)
    EVOS = EVOS[[!occursin("Im", str) for str in EVOS.FilePath], :]
	EVOS.order = [order_map[substring] for row in eachrow(EVOS) for substring in order_list if occursin(substring, row.FilePath)]
	EVOS = sort(EVOS, :order)
	select!(EVOS, Not(:order))
    col2 = EVOS.Value
    EVOS_averages = [mean(col2[i:i+2]) for i in 1:3:length(col2) if i+2 <= length(col2)]
    EVOS_stds = [std(col2[i:i+2]) for i in 1:3:length(col2) if i+2 <= length(col2)]

    p0 = [(EVOS_averages[2]-EVOS_averages[1])/(vols[4,2]-vols[1,2])]

    mass_avg = vcat(EVOS_averages)
    stds = vcat(EVOS_stds)
    vols = vols[[(!occursin("Pf", str) && !occursin("Pa",str)) for str in vols.first], :]
    col2 = vols.second
    vols_avg = [mean(col2[i:i+2]) for i in 1:3:length(col2) if i+2 <= length(col2)]
    vols_std = [std(col2[i:i+2]) for i in 1:3:length(col2) if i+2 <= length(col2)]

    @show fitlinear(mass_avg,vols_avg)
@pyexec """
def model(p, x):
    return p[0]*x 
""" => model 

    linear = odr.Model(model)
    odrdata = odr.Data(mass_avg, vols_avg, wd=1/(stds/sqrt(3)), we=1/(vols_std/sqrt(3)))
    myodr = odr.ODR(odrdata, linear, beta0=p0)
    output = myodr.run()

    fig = Figure(size=(6*72,3*72))
    ax = Axis(fig[1, 1], xlabel="BF-biofilm biomass (a.u.)", ylabel="CF-biofilm biomass (μm³)")
    errorbars!(ax, EVOS_averages, vols_avg[1:6], EVOS_stds, color=Makie.wong_colors()[1], direction = :x)
    errorbars!(ax, EVOS_averages, vols_avg[1:6], vols_std[1:6], color=Makie.wong_colors()[1], direction = :y)
    scatter!(ax, EVOS_averages, vols_avg[1:6], color=Makie.wong_colors()[1], label=rich("V. cholerae"; font=:italic))
	ax.rightspinevisible = false
	ax.topspinevisible = false
    fig[1,2] = Legend(fig, ax, merge = true, unique = true, framevisible=false, labelsize=12, rowgap=0)
    reg_params = pyconvert(Array, output.beta)
    reg_error = pyconvert(Array, output.sd_beta)
    avg = reg_params[1] .* mass_avg
    lines!(ax, mass_avg, avg, color="black")
    band!(ax, mass_avg, (reg_params[1]-reg_error[1]).*mass_avg, (reg_params[1]+reg_error[1]).*mass_avg, color=(:black, 0.2))
    save("figS1B.svg", fig)
end

function figS1C!(vc_mass, CV)
	order_list = [ 
                  "0ara_1", "0ara_2", "0ara_3", 
                  "0025ara_1", "0025ara_2", "0025ara_3",
                  "005ara_1", "005ara_2", "005ara_3",
                  "01ara_1", "01ara_2", "01ara_3",
                  ] 
	order_map = Dict(substring => idx for (idx, substring) in enumerate(order_list))
	CV.order = [order_map[substring] for row in eachrow(CV) for substring in order_list if occursin(substring, row.first)]
	CV = sort(CV, :order)
	select!(CV, Not(:order))
	
	function map_well_identifier(file_path)
		well_identifier = match(r"([A-Z])\d+", basename(file_path))[1][1]
		return Dict('A' => 1, 'B' => 2, 'C' => 3, 'D' => 4, 'E' => 5, 'F' => 6, 
                   'G' => 7, 'H' => 8)[well_identifier]
	end

    vc_mass = stack(vc_mass, Not([]), variable_name = :FilePath, value_name = :Value)
    vc_mass[!, :Well] = replace.(vc_mass.FilePath, r".*/(.*?)_.*" => s"\1")
	vc_mass[:, :condition] = map(map_well_identifier, vc_mass.FilePath)
    vc_mass = vc_mass[[(!occursin("B", str) && !occursin("E", str)) for str in vc_mass.Well],:]
    grouped_values = [group.Value for group in groupby(vc_mass, :condition)]
    vc_averages = [mean(subarray) for subarray in grouped_values] 
    vc_stds = [std(subarray) for subarray in grouped_values]

    p0 = [(vc_averages[2]-vc_averages[1])/(CV[4,2]-CV[1,2])]

    mass_avg = vcat(vc_averages)
    stds = vcat(vc_stds)
    col2 = CV.second
    CV_avg = [mean(col2[i:i+2]) for i in 1:3:length(col2) if i+2 <= length(col2)]
    CV_std = [std(col2[i:i+2]) for i in 1:3:length(col2) if i+2 <= length(col2)]
@pyexec """
def model(p, x):
    return p[0]*x 
""" => model 

    linear = odr.Model(model)
    odrdata = odr.Data(mass_avg, CV_avg, wd=1/(CV_avg/sqrt(9)), we=1/(CV_std/sqrt(3)))
    myodr = odr.ODR(odrdata, linear, beta0=p0)
    output = myodr.run()

    @show fitlinear(mass_avg,CV_avg)

    fig = Figure(size=(6*72,3*72))
    ax = Axis(fig[1, 1], xlabel="BF-biofilm biomass (a.u.)", ylabel="OD₅₉₀")
    errorbars!(ax, vc_averages, CV_avg, vc_stds, color=Makie.wong_colors()[1], direction = :x)
    errorbars!(ax, vc_averages, CV_avg, CV_std, color=Makie.wong_colors()[1], direction = :y)
    scatter!(ax, vc_averages, CV_avg, color=Makie.wong_colors()[1], label=rich("V. cholerae"; font=:italic))
	ax.rightspinevisible = false
	ax.topspinevisible = false
    reg_params = pyconvert(Array, output.beta)
    reg_error = pyconvert(Array, output.sd_beta)
    avg = reg_params[1] .* mass_avg
    lines!(ax, mass_avg, avg, color="black")
    band!(ax, mass_avg, (reg_params[1]-reg_error[1]).*mass_avg, (reg_params[1]+reg_error[1]).*mass_avg, color=(:black, 0.2))
    fig[1,2] = Legend(fig, ax, merge = true, unique = true, framevisible=false, labelsize=12, rowgap=0)
    save("figS1C.svg", fig)
end

function fig2_inoculum!(data)
    # 8 inocula, 9 replicates
    peaks = describe(data, :max)
    peaks = peaks[[(!occursin("A10", string(s)) && !occursin("A11", string(s)) && !occursin("A12", string(s))) for s in peaks.variable], :]
    peaks[!, :Well] = replace.(string.(peaks.variable), r".*/(.*?)_.*" => s"\1")

	function get_condition(well)
		letters = Dict('A' => 1, 'B' => 2, 'C' => 3, 'D' => 4,
					   'E' => 5, 'F' => 6, 'G' => 7, 'H' => 8)
		for (letter, num) in letters
			if occursin(letter, well)
				return num
			end
		end
		return missing 
	end

	# Add the new column "condition" to the DataFrame
	peaks.condition = [get_condition(well) for well in peaks.Well]

    grouped_values = [group.max for group in groupby(peaks, :condition)]
    data = vcat(grouped_values...)
    averages = [mean(subarray) for subarray in grouped_values] 
    maxes = [maximum(subarray) for subarray in grouped_values] 
    mins = [minimum(subarray) for subarray in grouped_values]

    category_num = Int.(1:8)
    category_num_swarm = Int.(repeat(1:8, inner = 9))
    fig = Figure(size=(4*72, 3.5*72))
    ax = Axis(fig[1, 1], xlabel="Inoculum (cells/mL)", ylabel="Peak BF-biofilm biomass (a.u.)")
    crossbar!(ax, category_num, averages, mins, maxes;
			  color=:grey, midlinecolor=:black)
    plt = beeswarm!(ax, category_num_swarm, Float64.(data), strokecolor=:black, color=:white, algorithm=UniformJitter(), strokewidth=1)
    ax.xticks=(1:8, [rich("10", superscript("7")), rich("10", superscript("6")), rich("10", superscript("5")), 
                     rich("10", superscript("4")), rich("10", superscript("3")), rich("10", superscript("2")), rich("10", superscript("1")), rich("10", superscript("0"))])
    ax.xticklabelrotation=45
    ax.xgridvisible = false
    ax.ygridvisible = false
	ax.rightspinevisible = false
	ax.topspinevisible = false
    save("fig2_inoculum.svg", fig)
end

function figS2_focal_plane!(data)
    # 7 focal planes, 3 replicates
    # c, d, e, f, g -> 5 plots
    rows = ["C", "D", "E"]
    fig = Figure(size=(10*72, 3*72))
    ga = fig[1, 1] = GridLayout()
    data_long = stack(data, Not([]), variable_name = :FilePath, value_name = :Value)
    data_long[!, :Well] = replace.(data_long.FilePath, r".*/(.*?)_.*" => s"\1")
    for (i,row) in enumerate(rows)
        ax = Axis(ga[1,i])
        data_sub = data_long[[occursin(row, str) for str in data_long.Well],:]
        data_sub[!, :XTick] = parse.(Int, replace.(data_sub.FilePath, r".*_(\d+)_.*" => s"\1"))
        grouped_values = [group.Value for group in groupby(data_sub, :XTick)]
		new_order = [4, 3, 2, 1, 5:length(grouped_values)...]
		grouped_values = grouped_values[new_order]
        focal_mean = mean(grouped_values[4])
        for j in 1:length(grouped_values)
            grouped_values[j] ./=  focal_mean
        end
        data = vcat(grouped_values...)
        averages = [mean(subarray) for subarray in grouped_values] 
        maxes = [maximum(subarray) for subarray in grouped_values] 
        mins = [minimum(subarray) for subarray in grouped_values]
        category_num = Int.(1:7)
        category_num_swarm = Int.(repeat(1:7, inner = 3))
        crossbar!(ax, category_num, averages, mins, maxes;
                  color=:white, midlinecolor=:black)
        plt = beeswarm!(ax, category_num_swarm, Float64.(data), strokecolor=:black, color=colorschemes[:Purples_3].colors[i], algorithm=UniformJitter(), strokewidth=1)
        ax.xticks=(1:7, ["-162", "-108", "-54", "0", "54", "108", "162"])
        ax.xticklabelrotation=45
        ax.xgridvisible = false
        ax.ygridvisible = false
        ax.rightspinevisible = false
        ax.topspinevisible = false
        if i == 2
            ax.xlabel="Relative focal plane (μm)"
            ax.title = "+0.025% Ara"
        end
        if i == 1
            ax.ylabel="BF-biofilm biomass (a.u.)"
            ax.title = "Untreated"
        end
        if i == 3
            ax.title = "+0.035% Ara"
        end
        #ylims!(ax, 0.0, 0.27)
    end
    linkaxes!(filter(x -> x isa Axis, fig.content)...)
    save("figS2_focal_plane.svg", fig)
end

function fig2_cps_nocps!(data)
    # 8 inocula, 9 replicates
    green = Colors.JULIA_LOGO_COLORS.green
    purple = Colors.JULIA_LOGO_COLORS.purple
    rename!(data, Dict(col => replace(col, r".*/(.*?)_.*" => s"\1") for col in names(data)))
    D39_WT = data[end, Cols("A1", "A2", "A3", "B1", "B2", "B3", "C1", "C2", "C3")]
    D39_dcps = data[end, Cols("A4", "A5", "A6", "B4", "B5", "B6", "C4", "C5", "C6")]
    SV36_WT = data[end, Cols("D1", "D2", "D3", "E1", "E2", "E3", "F1", "F2", "F3")]
    SV36_dcps = data[end, Cols("D4", "D5", "D6", "E4", "E5", "E6", "F4", "F5", "F6")]
    D39_data = [[values(D39_WT)...], [values(D39_dcps)...]]
    data = vcat(D39_data...)
    averages = [mean(subarray) for subarray in D39_data] 
    maxes = [maximum(subarray) for subarray in D39_data] 
    mins = [minimum(subarray) for subarray in D39_data]
    category_num = Int.(1:2)
    category_num_swarm = Int.(repeat(1:2, inner = 9))
    
    colormap1 = [purple, green]
    colormap2 = [purple, green]

    fig = Figure(size=(3*72, 3.5*72))
    ax = Axis(fig[1, 1], xlabel="", ylabel="BF-biofilm biomass (a.u.)")
    crossbar!(ax, category_num, averages, mins, maxes;
              color=:white, midlinecolor=[purple, green], colorrange=(1,2))
    plt = beeswarm!(ax, category_num_swarm, Float64.(data), strokecolor=:black, color=category_num_swarm, algorithm=UniformJitter(), strokewidth=1)
    plt.colormap[] = [purple, green]
    ax.xticks=(1:2, ["D39 wild-type", rich("D39 Δ", rich("cps"; font=:italic))])
    ax.xticklabelrotation=45
    ax.xgridvisible = false
    ax.ygridvisible = false
	ax.rightspinevisible = false
	ax.topspinevisible = false
    save("fig2_D39.svg", fig)
    
    SV36_data = [[values(SV36_WT)...], [values(SV36_dcps)...]]
    data = vcat(SV36_data...)
    averages = [mean(subarray) for subarray in SV36_data] 
    maxes = [maximum(subarray) for subarray in SV36_data] 
    mins = [minimum(subarray) for subarray in SV36_data]
    category_num = Int.(1:2)
    category_num_swarm = Int.(repeat(1:2, inner = 9))
    fig = Figure(size=(3*72, 3.5*72))
    ax = Axis(fig[1, 1], xlabel="", ylabel="BF-biofilm biomass (a.u.)")
    crossbar!(ax, category_num, averages, mins, maxes;
              color=:white, midlinecolor=[purple, green], colorrange=(1,2))
    plt = beeswarm!(ax, category_num_swarm, Float64.(data), strokecolor=:black, color=category_num_swarm, algorithm=UniformJitter(), strokewidth=1)
    plt.colormap[] = [purple, green] 
    ax.xticks=(1:2, ["SV36 wild-type", rich("SV36 ", rich("cps"; font=:italic), superscript("*"))])
    ax.xticklabelrotation=45
    ax.xgridvisible = false
    ax.ygridvisible = false
	ax.rightspinevisible = false
	ax.topspinevisible = false
    save("fig2_SV36.svg", fig)
end

function main()
    path_to_sans = "/root/.fonts/cmunss.ttf"
    path_to_sans_ital = "/root/.fonts/cmunsi.ttf"
    set_theme!(fonts = (; bold=path_to_sans, regular = path_to_sans, italic = path_to_sans_ital))
    data_folder = "../../Data/"
    #vols = DataFrame(CSV.File(joinpath(data_folder, "high_res_data.csv")))
    #vc_mass_path = "/mnt/e/Brightfield_paper/vc_fig1B/250120_4x_10x_plastic_Drawer1 final/4x/Numerical data/biomass.csv"
    #pa_pf_mass_path = "/mnt/e/Brightfield_paper/nonvc_fig1B/250120_4x_10x_plastic_2_Drawer2 final/4x/Numerical data/biomass.csv"
    #pf_wspf_path = "/mnt/e/Brightfield_paper/Pflu_wspF/4x/Numerical data/biomass.csv"
    #sp_mass_path = "/mnt/e/Brightfield_paper/sp_fig1B/4x/Numerical data/biomass.csv"
    #vc_mass = DataFrame(CSV.File(vc_mass_path))
    #pa_pf_mass = DataFrame(CSV.File(pa_pf_mass_path))
    #pf_wspf_mass = DataFrame(CSV.File(pf_wspf_path))
    #sp_mass = DataFrame(CSV.File(sp_mass_path))
    #fig1B!(vols, vc_mass, pa_pf_mass, pf_wspf_mass, nothing)

    #OD_image_path = "/mnt/f/Brightfield_paper/tests/OD_test.tif"
    #fig1A_OD_image(OD_image_path)  
    #fig1A_data_path = data_folder*"fig1A_lineplot.csv"
    #fig1A_lineplot!(fig1A_data_path)
    #focal_plane_path = "/mnt/e/Brightfield_paper/focal_plane/data/Numerical data/biomass.csv"
    #focal_plane_data = DataFrame(CSV.File(focal_plane_path))
    #figS2_focal_plane!(focal_plane_data)  
    #inoculum_path = "/mnt/e/Brightfield_paper/Inoculum/250205_4x_10x_plastic_PMR_Drawer1 05-Feb-2025 13-37-18/4x/Numerical data/biomass.csv"
    #inoculum_data = DataFrame(CSV.File(inoculum_path))
    #fig2_inoculum!(inoculum_data)
    
    #vols = DataFrame(CSV.File(joinpath(data_folder, "high_res_data.csv")))
    #vc_mass_path = "/mnt/e/Brightfield_paper/vc_fig1B/250120_4x_10x_plastic_Drawer1 final/10x/Numerical data/biomass.csv"
    #pa_pf_mass_path = "/mnt/e/Brightfield_paper/nonvc_fig1B/250120_4x_10x_plastic_2_Drawer2 final/10x/Numerical data/biomass.csv"
    #pf_wspf_path = "/mnt/e/Brightfield_paper/Pflu_wspF/10x/Numerical data/biomass.csv"
    #sp_mass_path = "/mnt/e/Brightfield_paper/sp_fig1B/10x/Numerical data/biomass.csv"
    #vc_mass = DataFrame(CSV.File(vc_mass_path))
    #pa_pf_mass = DataFrame(CSV.File(pa_pf_mass_path))
    #pf_wspf_mass = DataFrame(CSV.File(pf_wspf_path))
    #sp_mass = DataFrame(CSV.File(sp_mass_path))
    #fig1B!(vols, vc_mass, pa_pf_mass, pf_wspf_mass, nothing)
    
    #vols = DataFrame(CSV.File(joinpath(data_folder, "high_res_data.csv")))
    #EVOS_path = "/mnt/e/Brightfield_paper/EVOS/Numerical data/biomass.csv"
    #Nikon_path = "/mnt/e/Brightfield_paper/Nikon/Numerical data/biomass.csv"
    #EVOS = DataFrame(CSV.File(EVOS_path))
    #Nikon = DataFrame(CSV.File(Nikon_path))
    #EVOS_Nikon!(vols, EVOS, Nikon)
    
    #vc_mass_path = "/mnt/e/Brightfield_paper/vc_fig1B/250120_4x_10x_plastic_Drawer1 final/4x/Numerical data/biomass.csv"
    #vc_mass = DataFrame(CSV.File(vc_mass_path))
    #CV = DataFrame(CSV.File(joinpath(data_folder, "CV.csv")))
    #figS1C!(vc_mass, CV)
    
    #sp_brightfield_path = "/mnt/e/Brightfield_paper/250213_4x_10x_plastic_PMR_Drawer1 13-Feb-2025 12-24-15/4x/Numerical data/biomass.csv"
    #sp_brightfield = DataFrame(CSV.File(sp_brightfield_path))
    #fig2_cps_nocps!(sp_brightfield)   
end

main()
