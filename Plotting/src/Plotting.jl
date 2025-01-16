using Makie 
using CairoMakie
using SwarmMakie
using LsqFit
using TiffImages
using CSV
using StatsBase
using DataFrames
CairoMakie.activate!(type="svg")

function extract_float(gray_string::String)
    m = match(r"Gray\{Float64\}\(([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\)", gray_string)
    if m !== nothing
        return parse(Float64, m.captures[1])
    else
        error("Invalid Gray{Float64} string format: $gray_string")
    end
end


function fig1B!(vols, vc_mass, pf_mass, sp_mass)
    # Biovolume (x-axis) vs. brightfield biomass (y-axis)
    # Scatter + best-fit
    # Data are 1D dataframes where filenames/wells are the column names
    
    # Ensure vols are in the correct order
    substrings = ["D39", "36", "Pa"]
	vols = filter(row -> !any(substr -> occursin(substr, row.first), substrings), vols)
	
	order_list = ["vpsL", "0ara", "00125ara", "0025ara", "005ara", "01ara", "Pflu_no", "Pflu_WT",
                  "Pflu_wsp"] 
	order_map = Dict(substring => idx for (idx, substring) in enumerate(order_list))
	vols.order = [order_map[substring] for row in eachrow(vols) for substring in order_list if occursin(substring, row.first)]
	vols = sort(vols, :order)
	select!(vols, Not(:order))

	function map_well_identifier(file_path)
		well_identifier = match(r"([A-Z])\d+", basename(file_path))[1][1]
		return Dict('A' => 1, 'B' => 2, 'C' => 3, 'D' => 4, 'E' => 5, 'F' => 6, 
                   'G' => 7, 'H' => 8)[well_identifier]
	end

	function get_condition(well)
		if occursin("3", well)
			return 1
		elseif occursin("4", well)
			return 2
		elseif occursin("5", well)
			return 3
		elseif occursin("9", well)
			return 4 
		else
			return missing
		end
	end

    # Get averages and stds of vc_masses, pf_masses, sp_masses (correct order)
    pf_mass = stack(pf_mass, Not([]), variable_name = :FilePath, value_name = :Value)
    pf_mass[!, :Well] = replace.(pf_mass.FilePath, r".*/(.*?)_.*" => s"\1")
	pf_mass.condition = [get_condition(well) for well in pf_mass.Well]
    grouped_values = [group.Value for group in groupby(pf_mass, :condition)]
	grouped_values = [map(extract_float, inner_array) for inner_array in grouped_values] 
    vpsL = grouped_values[4]
    deleteat!(grouped_values, 4)
    data = vcat(grouped_values...)
    pf_averages = [mean(subarray) for subarray in grouped_values] 
    pf_stds = [std(subarray) for subarray in grouped_values]
    
    vc_mass = stack(vc_mass, Not([]), variable_name = :FilePath, value_name = :Value)
    vc_mass[!, :Well] = replace.(vc_mass.FilePath, r".*/(.*?)_.*" => s"\1")
	vc_mass[:, :condition] = map(map_well_identifier, vc_mass.FilePath)
    # Drop rows where condition = 6, 7, 8
    vc_mass = vc_mass[(vc_mass.condition .!= 6) .& (vc_mass.condition .!= 7) .& (vc_mass.condition .!= 8), :]
    grouped_values = [group.Value for group in groupby(vc_mass, :condition)]
    data = vcat(grouped_values..., vpsL)
    vc_averages = [mean(subarray) for subarray in grouped_values] 
    vc_stds = [std(subarray) for subarray in grouped_values]
    pushfirst!(vc_averages, mean(vpsL))
    pushfirst!(vc_stds, std(vpsL))

    @. model(x, p) = p[1]*x + p[2] 

    p0 = [(vc_averages[2]-vc_averages[1])/(vols[2,2]-vols[1,2]), vc_averages[1]]
    lb = [0.0, 0.0]
    ub = [Inf, Inf]

    mass_avg = vcat(vc_averages, pf_averages)
    fit = curve_fit(model, mass_avg, vols.second, p0, lower=lb, upper=ub)
    pstar = coef(fit)
	xbase = collect(range(minimum(vc_averages), maximum(vc_averages), 100))

    fig = Figure(size=(5*72,3*72))
    ax = Axis(fig[1, 1], xlabel="Biofilm OD (a.u.)", ylabel="Average thickness (µm)")
    errorbars!(ax, vc_averages, vols.second[1:6], vc_stds, color=Makie.wong_colors()[1], direction = :x)
    errorbars!(ax, pf_averages, vols.second[7:9], pf_stds, color=Makie.wong_colors()[2], direction = :x)
    scatter!(ax, vc_averages, vols.second[1:6], color=Makie.wong_colors()[1], label=rich("V. cholerae"; font=:italic))
    scatter!(ax, pf_averages, vols.second[7:9], color=Makie.wong_colors()[2], label=rich("P. fluorescens"; font=:italic))
	lines!(ax, xbase, model.(xbase, (pstar,)), color="black")
	ax.rightspinevisible = false
	ax.topspinevisible = false
    fig[1,2] = Legend(fig, ax, merge = true, unique = true, framevisible=false, labelsize=12, rowgap=0)
    save("fig1B.svg", fig)
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
    fig = Figure(size=(3*72, 2.5*72))
    ax = Axis(fig[1, 1], xlabel="Time (h)", ylabel="Biofilm OD (a.u.)")
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

function fig1C!(Vc_biomasses, Pa_biomasses, Sp_biomasses)
    # Peak biomass for biofilm former, non-former
    # 3 panels: V.c., P.a., S.p.
    # Jitter plot
    data = vcat(WT_growth, WT_dispersal)
    writedlm("$(plots_folder)/WT_radial_components.csv", data, ",")
    @show pvalue(UnequalVarianceTTest(Float64.(WT_growth), Float64.(WT_dispersal)))
    averages = [mean(data[1:5]), mean(data[6:10])]
    maxes = [maximum(data[1:5]), maximum(data[6:10])]
    mins = [minimum(data[1:5]), minimum(data[6:10])]
    category_num = Int.(1:2)
    category_num_swarm = Int.(repeat(1:2, inner = 5))
    fig = Figure()
    ax = Axis(fig[1, 1])
	crossbar!(ax, category_num, averages, mins, maxes;
			  color=:white, midlinecolor=:black)
    plt = beeswarm!(ax, category_num_swarm, Float64.(data), strokecolor=:black, color=:white, algorithm=UniformJitter(), strokewidth=1)
    ax.xticks=(1:2, unique(conditions))
    ax.xticklabelrotation=45
    ax.xlabel=""
    ax.ylabel="Radial displacement (µm)"
    ax.title=""
    ax.rightspinevisible = false
    ax.topspinevisible = false
    ax.xgridvisible = false
    ax.ygridvisible = false
    save("fig1C_Vc.svg", fig)
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
    fig = Figure(size=(3*72, 2.5*72))
    ax = Axis(fig[1, 1], xlabel="Relative focal plane (mm)", ylabel="Biofilm OD (a.u.)")
    data_long = stack(data, Not([]), variable_name = :FilePath, value_name = :Value)
    data_long[!, :Well] = replace.(data_long.FilePath, r".*/(.*?)_.*" => s"\1")
    data_long[!, :XTick] = parse.(Int, replace.(data_long.FilePath, r".*_(\d+)_.*" => s"\1"))
    grouped_values = [group.Value for group in groupby(data_long, :XTick)]
    data = vcat(grouped_values...)
    averages = [mean(subarray) for subarray in grouped_values] 
    maxes = [maximum(subarray) for subarray in grouped_values] 
    mins = [minimum(subarray) for subarray in grouped_values]
    category_num = Int.(1:7)
    category_num_swarm = Int.(repeat(1:7, inner = 3))
    crossbar!(ax, category_num, averages, mins, maxes;
			  color=:white, midlinecolor=:black)
    plt = beeswarm!(ax, category_num_swarm, Float64.(data), strokecolor=:black, color=:white, algorithm=UniformJitter(), strokewidth=1)
    ax.xticks=(1:7, ["-1.62", "-1.08", "-0.54", "0", "0.54", "1.08", "1.62"])
    ax.xticklabelrotation=45
    ax.xgridvisible = false
    ax.ygridvisible = false
	ax.rightspinevisible = false
	ax.topspinevisible = false
    ylims!(ax, 0.0, 0.27)
    save("figS2_focal_plane.svg", fig)
end

function main()
    path_to_sans = "/root/.fonts/cmunss.ttf"
    path_to_sans_ital = "/root/.fonts/cmunsi.ttf"
    set_theme!(fonts = (; bold=path_to_sans, regular = path_to_sans, italic = path_to_sans_ital))
    data_folder = "../../Data/"
    vols = DataFrame(CSV.File(joinpath(data_folder, "high_res_data.csv")))
    vc_mass_path = "/mnt/e/Brightfield_paper/vc_fig1B/4x/Numerical data/biomass.csv"
    pf_mass_path = "/mnt/e/Brightfield_paper/pf_fig1B/4x/Numerical data/biomass.csv"
    #sp_mass_path = "/mnt/e/Brightfield_paper/sp_fig1B/4x/Numerical data/biomass.csv"
    vc_mass = DataFrame(CSV.File(vc_mass_path))
    pf_mass = DataFrame(CSV.File(pf_mass_path))
    #sp_mass = DataFrame(CSV.File(sp_mass_path))
    fig1B!(vols, vc_mass, pf_mass, nothing)

    #OD_image_path = "/mnt/f/Brightfield_paper/tests/OD_test.tif"
    #fig1A_OD_image(OD_image_path)  
    #fig1A_data_path = data_folder*"fig1A_lineplot.csv"
    #fig1A_lineplot!(fig1A_data_path)
    #focal_plane_path = "/mnt/f/Brightfield_paper/focal_plane/data/Numerical data/biomass.csv"
    #focal_plane_data = DataFrame(CSV.File(focal_plane_path))
    #figS2_focal_plane!(focal_plane_data)  
    #inoculum_vc_path = "/mnt/f/Brightfield_paper/inoculum/250105_4x_bf_biospa_JP_Drawer3 05-Jan-2025 16-20-13/250105_Plate 1!PLATE_ID!_/Numerical data/biomass.csv"
    #inoculum_vc_data = DataFrame(CSV.File(inoculum_vc_path))
    #fig2_Vc!(inoculum_vc_data)
end

main()
