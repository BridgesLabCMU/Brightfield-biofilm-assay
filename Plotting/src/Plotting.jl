using Makie 
using CairoMakie
using SwarmMakie
using LsqFit
using TiffImages
using CSV
using StatsBase
using DataFrames
using EasyFit
using ColorSchemes: colorschemes
CairoMakie.activate!(type="svg")

function extract_float(gray_string::String)
    m = match(r"Gray\{Float64\}\(([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\)", gray_string)
    if m !== nothing
        return parse(Float64, m.captures[1])
    else
        error("Invalid Gray{Float64} string format: $gray_string")
    end
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

function fig1B!(vols, vc_mass, pa_pf_mass,sp_mass)
    # Biovolume (x-axis) vs. brightfield biomass (y-axis)
    # Scatter + best-fit
    # Data are 1D dataframes where filenames/wells are the column names
    
    # Ensure vols are in the correct order
    #substrings = ["D39", "36", "Pa"]
	#vols = filter(row -> !any(substr -> occursin(substr, row.first), substrings), vols)
	
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
    pf_stds = [pa_pf_stds[4], pa_pf_stds[3]]
    pf_averages = [pa_pf_averages[4], pa_pf_averages[3]]
    
    vc_mass = stack(vc_mass, Not([]), variable_name = :FilePath, value_name = :Value)
    vc_mass[!, :Well] = replace.(vc_mass.FilePath, r".*/(.*?)_.*" => s"\1")
	vc_mass[:, :condition] = map(map_well_identifier, vc_mass.FilePath)
    grouped_values = [group.Value for group in groupby(vc_mass, :condition)]
    vc_averages = [mean(subarray) for subarray in grouped_values] 
    vc_stds = [std(subarray) for subarray in grouped_values]

    @. model(x, p) = p[1]*x + p[2] 

    p0 = [(vc_averages[2]-vc_averages[1])/(vols[4,2]-vols[1,2]), vc_averages[1]]
    lb = [0.0, 0.0]
    ub = [Inf, Inf]

    mass_avg = vcat(vc_averages, pf_averages, pa_averages)
    stds = vcat(vc_stds, pf_stds, pa_stds)
    vols = vols[[!occursin("Pf_wspF", str) for str in vols.first], :]
    col2 = vols.second
    vols_avg = [mean(col2[i:i+2]) for i in 1:3:length(col2) if i+2 <= length(col2)]
    vols_std = [std(col2[i:i+2]) for i in 1:3:length(col2) if i+2 <= length(col2)]
    fit = curve_fit(model, mass_avg, vols_avg, p0, lower=lb, upper=ub)
    pstar = coef(fit)

    @show fitlinear(mass_avg,vols_avg)

	xbase = collect(range(minimum(vc_averages), maximum(vc_averages), 100))

    fig = Figure(size=(6*72,3*72))
    ax = Axis(fig[1, 1], xlabel="BF-biofilm biomass (a.u.)", ylabel="CF-biofilm biomass (μm³)")
    errorbars!(ax, vc_averages, vols_avg[1:6], vc_stds, color=Makie.wong_colors()[1], direction = :x)
    errorbars!(ax, pf_averages, vols_avg[7:8], pf_stds, color=Makie.wong_colors()[2], direction = :x)
    errorbars!(ax, pa_averages, vols_avg[9:10], pa_stds, color=Makie.wong_colors()[3], direction = :x)
    errorbars!(ax, vc_averages, vols_avg[1:6], vols_std[1:6], color=Makie.wong_colors()[1], direction = :y)
    errorbars!(ax, pf_averages, vols_avg[7:8], vols_std[7:8], color=Makie.wong_colors()[2], direction = :y)
    errorbars!(ax, pa_averages, vols_avg[9:10], vols_std[9:10], color=Makie.wong_colors()[3], direction = :y)
    scatter!(ax, vc_averages, vols_avg[1:6], color=Makie.wong_colors()[1], label=rich("V. cholerae"; font=:italic))
    scatter!(ax, pf_averages, vols_avg[7:8], color=Makie.wong_colors()[2], label=rich("P. fluorescens"; font=:italic))
    scatter!(ax, pa_averages, vols_avg[9:10], color=Makie.wong_colors()[3], label=rich("P. aeruginosa"; font=:italic))
	lines!(ax, xbase, model.(xbase, (pstar,)), color="black")
	ax.rightspinevisible = false
	ax.topspinevisible = false
    fig[1,2] = Legend(fig, ax, merge = true, unique = true, framevisible=false, labelsize=12, rowgap=0)
    save("figS1A.svg", fig)
end

function EVOS_Nikon!(vols, EVOS, Nikon)
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
	
    order_list = ["vpsL_001", "vpsL_002", "vpsL_003", 
                  "0ara_001", "0ara_002", "0ara_003", 
                  "0025ara_001", "0025ara_002", "0025ara_003",
                  "0035ara_001", "0035ara_002", "0035ara_003",
                  "005ara_001", "005ara_002", "005ara_003",
                  "01ara_001", "01ara_002", "01ara_003",
                  ] 
	order_map = Dict(substring => idx for (idx, substring) in enumerate(order_list))

    Nikon = stack(Nikon, Not([]), variable_name = :FilePath, value_name = :Value)
    Nikon = Nikon[[!occursin("Im", str) for str in Nikon.FilePath], :]
	Nikon.order = [order_map[substring] for row in eachrow(Nikon) for substring in order_list if occursin(substring, row.FilePath)]
	Nikon = sort(Nikon, :order)
	select!(Nikon, Not(:order))
    col2 = Nikon.Value
    Nikon_averages = [mean(col2[i:i+2]) for i in 1:3:length(col2) if i+2 <= length(col2)]
    Nikon_stds = [std(col2[i:i+2]) for i in 1:3:length(col2) if i+2 <= length(col2)]

    @. model(x, p) = p[1]*x + p[2] 

    p0 = [(Nikon_averages[2]-Nikon_averages[1])/(vols[4,2]-vols[1,2]), Nikon_averages[1]]
    lb = [0.0, 0.0]
    ub = [Inf, Inf]

    mass_avg = vcat(Nikon_averages)
    stds = vcat(Nikon_stds)
    vols = vols[[(!occursin("Pf", str) && !occursin("Pa",str)) for str in vols.first], :]
    col2 = vols.second
    vols_avg = [mean(col2[i:i+2]) for i in 1:3:length(col2) if i+2 <= length(col2)]
    vols_std = [std(col2[i:i+2]) for i in 1:3:length(col2) if i+2 <= length(col2)]
    fit = curve_fit(model, mass_avg, vols_avg, p0, lower=lb, upper=ub)
    pstar = coef(fit)

    @show fitlinear(mass_avg,vols_avg)

	xbase = collect(range(minimum(Nikon_averages), maximum(Nikon_averages), 100))

    fig = Figure(size=(6*72,3*72))
    ax = Axis(fig[1, 1], xlabel="BF-biofilm biomass (a.u.)", ylabel="CF-biofilm biomass (μm³)")
    errorbars!(ax, Nikon_averages, vols_avg[1:6], Nikon_stds, color=Makie.wong_colors()[1], direction = :x)
    errorbars!(ax, Nikon_averages, vols_avg[1:6], vols_std[1:6], color=Makie.wong_colors()[1], direction = :y)
    scatter!(ax, Nikon_averages, vols_avg[1:6], color=Makie.wong_colors()[1], label=rich("V. cholerae"; font=:italic))
	lines!(ax, xbase, model.(xbase, (pstar,)), color="black")
	ax.rightspinevisible = false
	ax.topspinevisible = false
    fig[1,2] = Legend(fig, ax, merge = true, unique = true, framevisible=false, labelsize=12, rowgap=0)
    save("figS1B_2.svg", fig)
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

    @. model(x, p) = p[1]*x + p[2] 

    p0 = [(vc_averages[2]-vc_averages[1])/(CV[4,2]-CV[1,2]), vc_averages[1]]
    lb = [0.0, 0.0]
    ub = [Inf, Inf]

    mass_avg = vcat(vc_averages)
    stds = vcat(vc_stds)
    col2 = CV.second
    CV_avg = [mean(col2[i:i+2]) for i in 1:3:length(col2) if i+2 <= length(col2)]
    CV_std = [std(col2[i:i+2]) for i in 1:3:length(col2) if i+2 <= length(col2)]
    fit = curve_fit(model, mass_avg, CV_avg, p0, lower=lb, upper=ub)
    pstar = coef(fit)

    @show fitlinear(mass_avg,CV_avg)

	xbase = collect(range(minimum(vc_averages), maximum(vc_averages), 100))

    fig = Figure(size=(6*72,3*72))
    ax = Axis(fig[1, 1], xlabel="BF-biofilm biomass (a.u.)", ylabel="OD₅₉₀")
    errorbars!(ax, vc_averages, CV_avg, vc_stds, color=Makie.wong_colors()[1], direction = :x)
    errorbars!(ax, vc_averages, CV_avg, CV_std, color=Makie.wong_colors()[1], direction = :y)
    scatter!(ax, vc_averages, CV_avg, color=Makie.wong_colors()[1], label=rich("V. cholerae"; font=:italic))
	lines!(ax, xbase, model.(xbase, (pstar,)), color="black")
	ax.rightspinevisible = false
	ax.topspinevisible = false
    fig[1,2] = Legend(fig, ax, merge = true, unique = true, framevisible=false, labelsize=12, rowgap=0)
    save("figS1C.svg", fig)
end

function fig2_Vc!(data)
    # 6 inocula, 9 replicates
    peaks = describe(data, :max)
    peaks[!, :Well] = replace.(string.(peaks.variable), r".*/(.*?)_.*" => s"\1")

	function get_condition(well)
		if occursin(r"[17]", well)
			return 1
		elseif occursin(r"[28]", well)
			return 2
		elseif occursin(r"[39]", well)
			return 3
		elseif occursin(r"[410]", well)
			return 4
		elseif occursin(r"[511]", well)
			return 5
		elseif occursin(r"[612]", well)
			return 6
		else
			return missing
		end
	end

	# Add the new column "condition" to the DataFrame
	peaks.condition = [get_condition(well) for well in peaks.Well]


    grouped_values = [group.max for group in groupby(peaks, :condition)]
    data = vcat(grouped_values...)
    averages = [mean(subarray) for subarray in grouped_values] 
    maxes = [maximum(subarray) for subarray in grouped_values] 
    mins = [minimum(subarray) for subarray in grouped_values]

    category_num = Int.(1:6)
    category_num_swarm = Int.(repeat(1:6, inner = 9))
    fig = Figure(size=(3*72, 2.5*72))
    ax = Axis(fig[1, 1], xlabel="Inoculum (cells/mL)", ylabel="Peak biofilm OD (a.u.)")
    crossbar!(ax, category_num, averages, mins, maxes;
			  color=:white, midlinecolor=:black)
    plt = beeswarm!(ax, category_num_swarm, Float64.(data), strokecolor=:black, color=:white, algorithm=UniformJitter(), strokewidth=1)
    ax.xticks=(1:6, [rich("10", superscript("7")), rich("10", superscript("6")), rich("10", superscript("5")), 
                     rich("10", superscript("4")), rich("10", superscript("3")), rich("10", superscript("2"))])
    ax.xticklabelrotation=45
    ax.xgridvisible = false
    ax.ygridvisible = false
	ax.rightspinevisible = false
	ax.topspinevisible = false
    save("fig2_Vc.svg", fig)
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

function main()
    path_to_sans = "/root/.fonts/cmunss.ttf"
    path_to_sans_ital = "/root/.fonts/cmunsi.ttf"
    set_theme!(fonts = (; bold=path_to_sans, regular = path_to_sans, italic = path_to_sans_ital))
    data_folder = "../../Data/"
    #vols = DataFrame(CSV.File(joinpath(data_folder, "high_res_data.csv")))
    #vc_mass_path = "/mnt/e/Brightfield_paper/vc_fig1B/250120_4x_10x_plastic_Drawer1 final/4x/Numerical data/biomass.csv"
    #pa_pf_mass_path = "/mnt/e/Brightfield_paper/nonvc_fig1B/250120_4x_10x_plastic_2_Drawer2 final/4x/Numerical data/biomass.csv"
    #sp_mass_path = "/mnt/e/Brightfield_paper/sp_fig1B/4x/Numerical data/biomass.csv"
    #vc_mass = DataFrame(CSV.File(vc_mass_path))
    #pa_pf_mass = DataFrame(CSV.File(pa_pf_mass_path))
    #sp_mass = DataFrame(CSV.File(sp_mass_path))
    #fig1B!(vols, vc_mass, pa_pf_mass, nothing)

    #OD_image_path = "/mnt/f/Brightfield_paper/tests/OD_test.tif"
    #fig1A_OD_image(OD_image_path)  
    #fig1A_data_path = data_folder*"fig1A_lineplot.csv"
    #fig1A_lineplot!(fig1A_data_path)
    focal_plane_path = "/mnt/e/Brightfield_paper/focal_plane/data/Numerical data/biomass.csv"
    focal_plane_data = DataFrame(CSV.File(focal_plane_path))
    figS2_focal_plane!(focal_plane_data)  
    #inoculum_vc_path = "/mnt/f/Brightfield_paper/inoculum/250105_4x_bf_biospa_JP_Drawer3 05-Jan-2025 16-20-13/250105_Plate 1!PLATE_ID!_/Numerical data/biomass.csv"
    #inoculum_vc_data = DataFrame(CSV.File(inoculum_vc_path))
    #fig2_Vc!(inoculum_vc_data)
    
    #vols = DataFrame(CSV.File(joinpath(data_folder, "high_res_data.csv")))
    #vc_mass_path = "/mnt/e/Brightfield_paper/vc_fig1B/250120_4x_10x_plastic_Drawer1 final/10x/Numerical data/biomass.csv"
    #pa_pf_mass_path = "/mnt/e/Brightfield_paper/nonvc_fig1B/250120_4x_10x_plastic_2_Drawer2 final/10x/Numerical data/biomass.csv"
    #sp_mass_path = "/mnt/e/Brightfield_paper/sp_fig1B/4x/Numerical data/biomass.csv"
    #vc_mass = DataFrame(CSV.File(vc_mass_path))
    #pa_pf_mass = DataFrame(CSV.File(pa_pf_mass_path))
    #sp_mass = DataFrame(CSV.File(sp_mass_path))
    #fig1B!(vols, vc_mass, pa_pf_mass, nothing)
    
    #vols = DataFrame(CSV.File(joinpath(data_folder, "high_res_data.csv")))
    #EVOS_path = "/mnt/e/Brightfield_paper/EVOS/Numerical data/biomass.csv"
    #Nikon_path = "/mnt/e/Brightfield_paper/Nikon/Numerical data/biomass.csv"
    #sp_mass_path = "/mnt/e/Brightfield_paper/sp_fig1B/4x/Numerical data/biomass.csv"
    #EVOS = DataFrame(CSV.File(EVOS_path))
    #Nikon = DataFrame(CSV.File(Nikon_path))
    #sp_mass = DataFrame(CSV.File(sp_mass_path))
    #EVOS_Nikon!(vols, EVOS, Nikon)
    
    #vc_mass_path = "/mnt/e/Brightfield_paper/vc_fig1B/250120_4x_10x_plastic_Drawer1 final/4x/Numerical data/biomass.csv"
    #vc_mass = DataFrame(CSV.File(vc_mass_path))
    #CV = DataFrame(CSV.File(joinpath(data_folder, "CV.csv")))
    #figS1C!(vc_mass, CV)
    
end

main()
