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
using ColorTypes
using ColorSchemes: colorschemes
using PythonCall
using Images: N0f8
CairoMakie.activate!(type="svg")

odr = pyimport("scipy.odr")

function gg_color_hue(n::Int)
	hues = range(15, stop=375, length=n+1)
	return [LCHuv(65, 100, h) for h in hues[1:end-1]]
end

function gg_color(i::Int, maxgroups::Int)
	hues = range(15, stop=375, length=maxgroups+1)
	return LCHuv(65, 100, hues[i])
end 

function extract_float(gray_string::String)
    m = match(r"Gray\{Float64\}\(([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\)", gray_string)
    if m !== nothing
        return parse(Float64, m.captures[1])
    else
        error("Invalid Gray{Float64} string format: $gray_string")
    end
end

function fig1B!(vols, vols_sp, vols_kleb, vc_mass, pa_pf_mass, pf_wspf_mass, sp_mass, kleb_mass)
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

    #order_list = ["D39_WT_1", "D39_WT_2", "D39_WT_3", "SV36_WT_1", "SV36_WT_2", "SV36_WT_3",
    #              "D39_dcps_1", "D39_dcps_2", "D39_dcps_3", "SV36_dcps_1", "SV36_dcps_2",
    #              "SV36_dcps_3"]
    order_list = ["D39_WT_nowash_10h_1", "D39_WT_nowash_10h_2", "D39_WT_nowash_10h_3",
                  "D39_Dcps_nowash_10h_1", "D39_Dcps_nowash_10h_2", "D39_Dcps_nowash_10h_3",
                  "SV36_WT_nowash_10h_1", "SV36_WT_nowash_10h_2", "SV36_WT_nowash_10h_3",
                  "SV36_cps-mut_nowash_10h_1", "SV36_cps-mut_nowash_10h_2", "SV36_cps-mut_nowash_10h_3"
                  ]
	order_map = Dict(substring => idx for (idx, substring) in enumerate(order_list))
	vols_sp.order = [order_map[substring] for row in eachrow(vols_sp) for substring in order_list if occursin(substring, row.first)]
	vols_sp = sort(vols_sp, :order)
	select!(vols_sp, Not(:order))
    vols = vcat(vols, vols_sp)
	
    order_list = ["Kleb_WT_1", "Kleb_WT_2", "Kleb_WT_3", "Kleb_mutant_1", "Kleb_mutant_2", 
                  "Kleb_mutant_3"]
    order_map = Dict(substring => idx for (idx, substring) in enumerate(order_list))
	vols_kleb.order = [order_map[substring] for row in eachrow(vols_kleb) for substring in order_list if occursin(substring, row.first)]
	vols_kleb = sort(vols_kleb, :order)
	select!(vols_kleb, Not(:order))
    vols = vcat(vols, vols_kleb)

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
    
	rename!(sp_mass, Symbol("Column5") => :Value)

	grp_order = ["D39 WT",        # 1
				 "D39 Dcps",      # 2
				 "SV36 WT",       # 3
				 "SV36 cps-mut"]  # 4  

	grouped_values = [
		collect(skipmissing(sp_mass[sp_mass.Strain .== g, :Value]))
		for g in grp_order
	]
	sp_averages = map(mean, grouped_values)
	sp_stds     = map(std,  grouped_values)

    rename!(kleb_mass, Dict(col => replace(col, r".*/(.*?)_.*" => s"\1") for col in names(kleb_mass)))
    kleb_WT = kleb_mass[1, Cols("A7", "A8", "A9", "B7", "B8", "B9", "C7", "C8", "C9")]
    kleb_mutant = kleb_mass[1, Cols("A10", "A11", "A12", "B10", "B11", "B12", "C10", "C11", "C12")]
    grouped_values = [[values(kleb_WT)...], [values(kleb_mutant)...]]
    kleb_averages = [mean(subarray) for subarray in grouped_values] 
    kleb_stds = [std(subarray) for subarray in grouped_values]

@pyexec """
def model(p, x):
    return p[0]*x 
""" => model 

    p0 = [(vc_averages[2]-vc_averages[1])/(vols[4,2]-vols[1,2])]
    lb = [0.0, 0.0]
    ub = [Inf, Inf]

    mass_avg = vcat(vc_averages, pf_averages, pa_averages, sp_averages, kleb_averages)
    stds = vcat(vc_stds, pf_stds, pa_stds, sp_stds, kleb_stds)
    col2 = vols.second
    vols_avg = [mean(col2[i:i+2]) for i in 1:3:length(col2) if i+2 <= length(col2)]
    vols_std = [std(col2[i:i+2]) for i in 1:3:length(col2) if i+2 <= length(col2)]
    @show fitlinear(mass_avg,vols_avg)

    linear = odr.Model(model)
    odrdata = odr.Data(mass_avg, vols_avg, wd=1/(stds/sqrt(9)), we=1/(vols_std/sqrt(3)))
    myodr = odr.ODR(odrdata, linear, beta0=p0)
    output = myodr.run()
    fig = Figure(size=(6*72,3*72))
    ax = Axis(fig[1, 1], xlabel="BF-biofilm biomass (a.u.)", ylabel="CF-biofilm biomass (μm³)")
    errorbars!(ax, vc_averages, vols_avg[1:6], vc_stds, color=Makie.wong_colors()[1], direction = :x)
    errorbars!(ax, pf_averages, vols_avg[7:9], pf_stds, color=Makie.wong_colors()[2], direction = :x)
    errorbars!(ax, pa_averages, vols_avg[10:11], pa_stds, color=Makie.wong_colors()[3], direction = :x)
    errorbars!(ax, sp_averages, vols_avg[12:15], sp_stds, color=Makie.wong_colors()[4], direction = :x)
    errorbars!(ax, [kleb_averages[2]], [vols_avg[17]], [kleb_stds[2]], color=Makie.wong_colors()[5], direction = :x)
    errorbars!(ax, vc_averages, vols_avg[1:6], vols_std[1:6], color=Makie.wong_colors()[1], direction = :y)
    errorbars!(ax, pf_averages, vols_avg[7:9], vols_std[7:9], color=Makie.wong_colors()[2], direction = :y)
    errorbars!(ax, pa_averages, vols_avg[10:11], vols_std[10:11], color=Makie.wong_colors()[3], direction = :y)
    errorbars!(ax, sp_averages, vols_avg[12:15], vols_std[12:15], color=Makie.wong_colors()[4], direction = :y)
    errorbars!(ax, [kleb_averages[2]], [vols_avg[17]], [vols_std[17]], color=Makie.wong_colors()[5], direction = :y)
    scatter!(ax, vc_averages, vols_avg[1:6], color=Makie.wong_colors()[1], label=rich("V. cholerae"; font=:italic))
    scatter!(ax, pf_averages, vols_avg[7:9], color=Makie.wong_colors()[2], label=rich("P. fluorescens"; font=:italic))
    scatter!(ax, pa_averages, vols_avg[10:11], color=Makie.wong_colors()[3], label=rich("P. aeruginosa"; font=:italic))
    scatter!(ax, sp_averages, vols_avg[12:15], color=Makie.wong_colors()[4], label=rich("S. pneumoniae"; font=:italic))
    scatter!(ax,[ kleb_averages[2]], [vols_avg[17]], color=Makie.wong_colors()[5], label=rich("K. pneumoniae"; font=:italic))
    reg_params = pyconvert(Array, output.beta)
    reg_error = pyconvert(Array, output.sd_beta)
    avg = reg_params[1] .* mass_avg
    lines!(ax, mass_avg, avg, color="black")
    band!(ax, mass_avg, (reg_params[1]-reg_error[1]).*mass_avg, (reg_params[1]+reg_error[1]).*mass_avg, color=(:black, 0.2))
	ax.rightspinevisible = false
	ax.topspinevisible = false
    ax.xgridvisible = false
    ax.ygridvisible = false
    fig[1,2] = Legend(fig, ax, merge = true, unique = true, framevisible=false, labelsize=12, rowgap=0)

    # 1) create the subset mask
	species_groups = vcat(
	  fill(1, length(vc_averages)),    # V. cholerae → group 1
	  fill(2, length(pf_averages)),    # P. fluorescens → group 2
	  fill(3, length(pa_averages)),    # P. aeruginosa → group 3
	  fill(4, length(sp_averages)),    # S. pneumoniae → group 4
	  fill(5, length(kleb_averages))   # K. pneumoniae → group 5
	)
	palette = Makie.wong_colors()[1:5]
    mask = mass_avg .<= 0.03

    # 2) extract subset vectors
    mass_sub   = mass_avg[mask]
    vols_sub   = vols_avg[mask]
    xerr_sub   = stds[mask]
    # make sure you still have vols_std defined as before:
    yerr_sub   = vols_std[mask]

    # 3) add an inset axis, 30% width & height of the main panel,
    #    anchored at the top-left of cell (1,1)
    inset = Axis(fig[1,1];
        width       = Relative(0.3),
        height      = Relative(0.3),
        halign = 0.1,
        valign = 0.9,
        xlabel      = "",              # no need to relabel
        ylabel      = "",
    )
    inset.xticklabelsvisible=false
    inset.yticklabelsvisible=false
    inset.xticksvisible = false
    inset.yticksvisible=false

	for g in 1:5
	  idx = findall(i -> mask[i] && species_groups[i] == g, eachindex(mask))
	  if !isempty(idx)
		# x‐errorbars
		errorbars!(inset,
		  mass_avg[idx], vols_avg[idx],
		  stds[idx],                  # x‐uncertainty
		  direction = :x,
		  color     = palette[g]
		)
		# y‐errorbars
		errorbars!(inset,
		  mass_avg[idx], vols_avg[idx],
		  vols_std[idx],              # y‐uncertainty
		  direction = :y,
		  color     = palette[g]
		)
		# points
		scatter!(inset,
		  mass_avg[idx], vols_avg[idx],
		  color     = palette[g],
		  markersize =10 
		)
	  end
	end
    lines!(inset, mass_sub, reg_params[1] .* mass_sub, color = :black)
    band!(inset,
                mass_sub,
                (reg_params[1] .- reg_error[1]) .* mass_sub,
                (reg_params[1] .+ reg_error[1]) .* mass_sub,
                color = (:black, 0.2)
    )
	inset.rightspinevisible = false
	inset.topspinevisible = false
    inset.xgridvisible = false
    inset.ygridvisible = false
    save("fig1B.svg", fig)
end

function figS1A!(vols, vols_sp, vols_kleb, vc_mass, pa_pf_mass, pf_wspf_mass, sp_mass, kleb_mass)
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

    #order_list = ["D39_WT_1", "D39_WT_2", "D39_WT_3", "SV36_WT_1", "SV36_WT_2", "SV36_WT_3",
    #              "D39_dcps_1", "D39_dcps_2", "D39_dcps_3", "SV36_dcps_1", "SV36_dcps_2",
    #              "SV36_dcps_3"]
    order_list = ["D39_WT_nowash_10h_1", "D39_WT_nowash_10h_2", "D39_WT_nowash_10h_3",
                  "D39_Dcps_nowash_10h_1", "D39_Dcps_nowash_10h_2", "D39_Dcps_nowash_10h_3",
                  "SV36_WT_nowash_10h_1", "SV36_WT_nowash_10h_2", "SV36_WT_nowash_10h_3",
                  "SV36_cps-mut_nowash_10h_1", "SV36_cps-mut_nowash_10h_2", "SV36_cps-mut_nowash_10h_3"
                  ]
	order_map = Dict(substring => idx for (idx, substring) in enumerate(order_list))
	vols_sp.order = [order_map[substring] for row in eachrow(vols_sp) for substring in order_list if occursin(substring, row.first)]
	vols_sp = sort(vols_sp, :order)
	select!(vols_sp, Not(:order))
    vols = vcat(vols, vols_sp)
	
    order_list = ["Kleb_WT_1", "Kleb_WT_2", "Kleb_WT_3", "Kleb_mutant_1", "Kleb_mutant_2", 
                  "Kleb_mutant_3"]
    order_map = Dict(substring => idx for (idx, substring) in enumerate(order_list))
	vols_kleb.order = [order_map[substring] for row in eachrow(vols_kleb) for substring in order_list if occursin(substring, row.first)]
	vols_kleb = sort(vols_kleb, :order)
	select!(vols_kleb, Not(:order))
    vols = vcat(vols, vols_kleb)

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
    
    rename!(sp_mass, Dict(col => replace(col, r".*/(.*?)_.*" => s"\1") for col in names(sp_mass)))
    D39_WT = sp_mass[1, Cols("A1", "A2", "A3", "B1", "B2", "B3", "C1", "C2", "C3")]
    D39_dcps = sp_mass[1, Cols("A4", "A5", "A6", "B4", "B5", "B6", "C4", "C5", "C6")]
    SV36_WT = sp_mass[1, Cols("D1", "D2", "D3", "E1", "E2", "E3", "F1", "F2", "F3")]
    SV36_dcps = sp_mass[1, Cols("D4", "D5", "D6", "E4", "E5", "E6", "F4", "F5", "F6")]
    grouped_values = [[values(D39_WT)...], [values(D39_dcps)...],[values(SV36_WT)...],  
                      [values(SV36_dcps)...]]	
	sp_averages = map(mean, grouped_values)
	sp_stds     = map(std,  grouped_values)

    rename!(kleb_mass, Dict(col => replace(col, r".*/(.*?)_.*" => s"\1") for col in names(kleb_mass)))
    kleb_WT = kleb_mass[1, Cols("A7", "A8", "A9", "B7", "B8", "B9", "C7", "C8", "C9")]
    kleb_mutant = kleb_mass[1, Cols("A10", "A11", "A12", "B10", "B11", "B12", "C10", "C11", "C12")]
    grouped_values = [[values(kleb_WT)...], [values(kleb_mutant)...]]
    kleb_averages = [mean(subarray) for subarray in grouped_values] 
    kleb_stds = [std(subarray) for subarray in grouped_values]

@pyexec """
def model(p, x):
    return p[0]*x 
""" => model 

    p0 = [(vc_averages[2]-vc_averages[1])/(vols[4,2]-vols[1,2])]
    lb = [0.0, 0.0]
    ub = [Inf, Inf]

    mass_avg = vcat(vc_averages, pf_averages, pa_averages, sp_averages, kleb_averages)
    stds = vcat(vc_stds, pf_stds, pa_stds, sp_stds, kleb_stds)
    col2 = vols.second
    vols_avg = [mean(col2[i:i+2]) for i in 1:3:length(col2) if i+2 <= length(col2)]
    vols_std = [std(col2[i:i+2]) for i in 1:3:length(col2) if i+2 <= length(col2)]
    @show fitlinear(mass_avg,vols_avg)

    linear = odr.Model(model)
    odrdata = odr.Data(mass_avg, vols_avg, wd=1/(stds/sqrt(9)), we=1/(vols_std/sqrt(3)))
    myodr = odr.ODR(odrdata, linear, beta0=p0)
    output = myodr.run()
    fig = Figure(size=(6*72,3*72))
    ax = Axis(fig[1, 1], xlabel="BF-biofilm biomass (a.u.)", ylabel="CF-biofilm biomass (μm³)")
    errorbars!(ax, vc_averages, vols_avg[1:6], vc_stds, color=Makie.wong_colors()[1], direction = :x)
    errorbars!(ax, pf_averages, vols_avg[7:9], pf_stds, color=Makie.wong_colors()[2], direction = :x)
    errorbars!(ax, pa_averages, vols_avg[10:11], pa_stds, color=Makie.wong_colors()[3], direction = :x)
    errorbars!(ax, sp_averages, vols_avg[12:15], sp_stds, color=Makie.wong_colors()[4], direction = :x)
    errorbars!(ax, [kleb_averages[2]], [vols_avg[17]], [kleb_stds[2]], color=Makie.wong_colors()[5], direction = :x)
    errorbars!(ax, vc_averages, vols_avg[1:6], vols_std[1:6], color=Makie.wong_colors()[1], direction = :y)
    errorbars!(ax, pf_averages, vols_avg[7:9], vols_std[7:9], color=Makie.wong_colors()[2], direction = :y)
    errorbars!(ax, pa_averages, vols_avg[10:11], vols_std[10:11], color=Makie.wong_colors()[3], direction = :y)
    errorbars!(ax, sp_averages, vols_avg[12:15], vols_std[12:15], color=Makie.wong_colors()[4], direction = :y)
    errorbars!(ax, [kleb_averages[2]], [vols_avg[17]], [vols_std[17]], color=Makie.wong_colors()[5], direction = :y)
    scatter!(ax, vc_averages, vols_avg[1:6], color=Makie.wong_colors()[1], label=rich("V. cholerae"; font=:italic))
    scatter!(ax, pf_averages, vols_avg[7:9], color=Makie.wong_colors()[2], label=rich("P. fluorescens"; font=:italic))
    scatter!(ax, pa_averages, vols_avg[10:11], color=Makie.wong_colors()[3], label=rich("P. aeruginosa"; font=:italic))
    scatter!(ax, sp_averages, vols_avg[12:15], color=Makie.wong_colors()[4], label=rich("S. pneumoniae"; font=:italic))
    scatter!(ax,[ kleb_averages[2]], [vols_avg[17]], color=Makie.wong_colors()[5], label=rich("K. pneumoniae"; font=:italic))
    reg_params = pyconvert(Array, output.beta)
    reg_error = pyconvert(Array, output.sd_beta)
    avg = reg_params[1] .* mass_avg
    lines!(ax, mass_avg, avg, color="black")
    band!(ax, mass_avg, (reg_params[1]-reg_error[1]).*mass_avg, (reg_params[1]+reg_error[1]).*mass_avg, color=(:black, 0.2))
	ax.rightspinevisible = false
	ax.topspinevisible = false
    ax.xgridvisible = false
    ax.ygridvisible = false
    fig[1,2] = Legend(fig, ax, merge = true, unique = true, framevisible=false, labelsize=12, rowgap=0)

    # 1) create the subset mask
	species_groups = vcat(
	  fill(1, length(vc_averages)),    # V. cholerae → group 1
	  fill(2, length(pf_averages)),    # P. fluorescens → group 2
	  fill(3, length(pa_averages)),    # P. aeruginosa → group 3
	  fill(4, length(sp_averages)),    # S. pneumoniae → group 4
	  fill(5, length(kleb_averages))   # K. pneumoniae → group 5
	)
	palette = Makie.wong_colors()[1:5]
    mask = mass_avg .<= 0.03

    # 2) extract subset vectors
    mass_sub   = mass_avg[mask]
    vols_sub   = vols_avg[mask]
    xerr_sub   = stds[mask]
    # make sure you still have vols_std defined as before:
    yerr_sub   = vols_std[mask]

    # 3) add an inset axis, 30% width & height of the main panel,
    #    anchored at the top-left of cell (1,1)
    inset = Axis(fig[1,1];
        width       = Relative(0.3),
        height      = Relative(0.3),
        halign = 0.1,
        valign = 0.9,
        xlabel      = "",              # no need to relabel
        ylabel      = "",
    )
    inset.xticklabelsvisible=false
    inset.yticklabelsvisible=false
    inset.xticksvisible = false
    inset.yticksvisible=false

	for g in 1:5
	  idx = findall(i -> mask[i] && species_groups[i] == g, eachindex(mask))
	  if !isempty(idx)
		# x‐errorbars
		errorbars!(inset,
		  mass_avg[idx], vols_avg[idx],
		  stds[idx],                  # x‐uncertainty
		  direction = :x,
		  color     = palette[g]
		)
		# y‐errorbars
		errorbars!(inset,
		  mass_avg[idx], vols_avg[idx],
		  vols_std[idx],              # y‐uncertainty
		  direction = :y,
		  color     = palette[g]
		)
		# points
		scatter!(inset,
		  mass_avg[idx], vols_avg[idx],
		  color     = palette[g],
		  markersize =10 
		)
	  end
	end
    lines!(inset, mass_sub, reg_params[1] .* mass_sub, color = :black)
    band!(inset,
                mass_sub,
                (reg_params[1] .- reg_error[1]) .* mass_sub,
                (reg_params[1] .+ reg_error[1]) .* mass_sub,
                color = (:black, 0.2)
    )
	inset.rightspinevisible = false
	inset.topspinevisible = false
    inset.xgridvisible = false
    inset.ygridvisible = false
    save("figS1A.svg", fig)
end

function fig1A_OD_image(OD_image_path, mask_path)
    image = Float64.(TiffImages.load(OD_image_path))
    fig = Figure(figure_padding = 0, fontsize=38)
    ax = Axis(fig[1, 1], aspect = DataAspect(), yticklabelsvisible = false, xticklabelsvisible = false, yticksvisible = false, xticksvisible = false)
    hm = heatmap!(ax, rotr90(image), colormap=:plasma)
    Colorbar(fig[end+1, :], hm, vertical = false, flipaxis = false, label = "OD", width = Relative(0.7))
    save("fig1A_OD_image.png", fig)

    mask = TiffImages.load(mask_path)
	CYAN = RGB{N0f8}(0, 1, 1)
	mask = mask .== CYAN
    fig = Figure(figure_padding = 0, fontsize=38)
    ax = Axis(fig[1, 1], aspect = DataAspect(), yticklabelsvisible = false, xticklabelsvisible = false, yticksvisible = false, xticksvisible = false)
    image .*= mask
    hm = heatmap!(ax, rotr90(image), colormap=:plasma)
    save("fig1A_OD_image_masked.png", fig)
end

function fig1A_lineplot!(fig1A_data_path)
    data = DataFrame(CSV.File(fig1A_data_path))
    fig = Figure(size=(6.5*72, 2.7*72))
    ax = Axis(fig[1, 1], xlabel="Time (h)", ylabel="BF-biofilm biomass (a.u.)")
    vpvc_columns = select(data, Cols.(contains.("/E"))) 
    vpvc_avg = reduce(+, eachcol(vpvc_columns)) ./ ncol(vpvc_columns) 
    vpvc_stdev = dropdims(std(Array(vpvc_columns), dims=2), dims=2)  
    time = 0:0.5:nrow(vpvc_columns)/2-0.5
    lines!(ax, time, vpvc_avg, color=:black)
    band!(ax, time, vpvc_avg-vpvc_stdev, vpvc_avg+vpvc_stdev, color=(:black, 0.2))
    scatter!(ax, time, vpvc_avg, color=:white, marker=:circle,  strokewidth=1, 
             label = rich(rich("Pbad-vpvC"; font=:italic), superscript("W240R"), "+0.1% Ara"))
    vpsl_columns = select(data, Cols.(contains.("/H"))) 
    vpsl_avg = reduce(+, eachcol(vpsl_columns)) ./ ncol(vpsl_columns) 
    vpsl_stdev = dropdims(std(Array(vpsl_columns), dims=2), dims=2)  
    lines!(ax, time, vpsl_avg, color=:black)
    band!(ax, time, vpsl_avg-vpsl_stdev, vpsl_avg+vpsl_stdev, color=(:black, 0.2))
    scatter!(ax, time, vpsl_avg, color=:white, marker=:utriangle,  strokewidth=1, 
             label = rich("Δ", rich("vpsL"; font=:italic)))
    fig[1,2] = Legend(fig, ax, framevisible=false, rowgap=0)
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

    fig = Figure(size=(3.7*72,3*72))
    ax = CairoMakie.Axis(fig[1, 1], xlabel="BF-biofilm biomass (a.u.)", ylabel="CF-biofilm biomass (μm³)")
    errorbars!(ax, EVOS_averages, vols_avg[1:6], EVOS_stds, color=Makie.wong_colors()[1], direction = :x)
    errorbars!(ax, EVOS_averages, vols_avg[1:6], vols_std[1:6], color=Makie.wong_colors()[1], direction = :y)
    scatter!(ax, EVOS_averages, vols_avg[1:6], color=Makie.wong_colors()[1], label=rich("V. cholerae"; font=:italic))
	ax.rightspinevisible = false
	ax.topspinevisible = false
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

    fig = Figure(size=(3*72,3*72))
    ax = CairoMakie.Axis(fig[1, 1], xlabel="BF-biofilm biomass (a.u.)", ylabel="OD₅₉₀")
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
    rows = ["C", "D", "G"]
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
            ax.title = "+0.1% Ara"
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

function detect_outliers_iqr(data)
    q1 = quantile(data, 0.25)
    q3 = quantile(data, 0.75)
    iqr = q3 - q1

    lower_bound = q1 - 1.5 * iqr
    upper_bound = q3 + 1.5 * iqr

    outliers = filter(x -> x < lower_bound || x > upper_bound, data)
    return outliers, lower_bound, upper_bound
end

function fig3_stats!(screen_plate_paths)
	wells_to_remove = Dict(
		1 => ["A1"], 2 => ["A1"], 3 => [], 4 => ["C6"], 5 => ["A1"], 
		6 => ["A1", "B6", "H1"], 7 => ["A1", "G6"], 8 => ["A1"], 9 => ["A1"], 
		10 => ["A1", "A9"], 11 => ["A1", "C5", "D11"], 12 => ["A1", "C3", "C7"], 
		13 => [], 14 => [], 15 => ["A1"], 16 => ["A1", "B6", "H2"], 
		17 => ["A1", "F10"], 18 => ["A1"], 19 => ["C9"], 20 => ["B1"], 
		21 => ["A1", "C6"], 22 => [], 23 => ["E5"], 24 => ["A1", "F12"], 
		25 => ["A1"]
	)
    all_biofilm = []
    all_planktonic = []
    for (i,plate_path) in enumerate(screen_plate_paths)
        biofilm = DataFrame(CSV.File(plate_path*"/Numerical data/biomass.csv"))
        planktonic = DataFrame(CSV.File(plate_path*"/Numerical data/planktonic.csv"))
		rename!(biofilm, Dict(col => replace(col, r".*/(.*?)_.*" => s"\1") for col in names(biofilm)))
        rename!(planktonic, Dict(col => replace(col, r".*/(.*?)_.*" => s"\1") for col in names(planktonic)))
        planktonic = planktonic[:,Not(wells_to_remove[i])]
        biofilm = biofilm[:,Not(wells_to_remove[i])]
        append!(all_planktonic, values(planktonic[1,:]))
        append!(all_biofilm, values(biofilm[1,:]))
    end
    sort!(all_biofilm, rev=true)
    sort!(all_planktonic, rev=true)
    biofilm_o, biofilm_lb, biofilm_ub = detect_outliers_iqr(all_biofilm)
    planktonic_o, planktonic_lb, planktonic_ub = detect_outliers_iqr(all_planktonic)
	fig = Figure(size=(4*72, 3.5*72))
	ax = Axis(fig[1, 1], xlabel="Transposon rank", ylabel="BF-Biofilm biomass (a.u.)")
	scatter!(ax, 1:length(all_biofilm), all_biofilm, color=Makie.wong_colors()[5], markersize=2)
    hlines!(ax, biofilm_lb, color=:grey, linestyle=:dash)
    hlines!(ax, biofilm_ub, color=:grey, linestyle=:dash)
    ax.xgridvisible = false
    ax.ygridvisible = false
	ax.rightspinevisible = false
	ax.topspinevisible = false
    save("fig3_biofilm.svg", fig)
	
    fig = Figure(size=(4*72, 3.5*72))
    ax = Axis(fig[1, 1], xlabel="BF-Biofilm biomass (a.u.)", ylabel="Frequency")
    hist!(ax, all_biofilm, color=(Makie.wong_colors()[5], 0.5), bins=30)
    vlines!(ax, biofilm_lb, color=:grey, linestyle=:dash)
    vlines!(ax, biofilm_ub, color=:grey, linestyle=:dash)
    ax.xgridvisible = false
    ax.ygridvisible = false
	ax.rightspinevisible = false
	ax.topspinevisible = false
    save("fig3_biofilm_hist.svg", fig)

    fig = Figure(size=(4*72, 3.5*72))
	ax = Axis(fig[1, 1], xlabel="Transposon rank", ylabel="Planktonic biomass (a.u.)")
	scatter!(ax, 1:length(all_planktonic), all_planktonic, color=Makie.wong_colors()[6], markersize=2)
    hlines!(ax, planktonic_lb, color=:grey, linestyle=:dash)
    hlines!(ax, planktonic_ub, color=:grey, linestyle=:dash)
    ax.xgridvisible = false
    ax.ygridvisible = false
	ax.rightspinevisible = false
	ax.topspinevisible = false
    save("fig3_planktonic.svg", fig)
    
    fig = Figure(size=(4*72, 3.5*72))
    ax = Axis(fig[1, 1], xlabel="Planktonic biomass (a.u.)", ylabel="Frequency")
    hist!(ax, all_planktonic, color=(Makie.wong_colors()[6], 0.5), bins=30)
    vlines!(ax, planktonic_lb, color=:grey, linestyle=:dash)
    vlines!(ax, planktonic_ub, color=:grey, linestyle=:dash)
    ax.xgridvisible = false
    ax.ygridvisible = false
	ax.rightspinevisible = false
	ax.topspinevisible = false
    save("fig3_planktonic_hist.svg", fig)
end

function wash_nowash!(wash, nowash)
    # Rename columns for easier handling
    rename!(wash, Dict(col => replace(col, r".*/(.*?)_.*" => s"\1") for col in names(wash)))
    rename!(nowash, Dict(col => replace(col, r".*/(.*?)_.*" => s"\1") for col in names(nowash)))

    # Define conditions and tick labels; only using "WT", "ara0025", and "ara01"
    conditions = ["WT", "ara0025", "ara01"]
    tick_labels = Dict("WT" => "Untreated", "ara0025" => "+0.025% Ara", "ara01" => "+0.1% Ara")
    
    # Define column mapping for the remaining conditions.
    col_mapping = Dict(
        "WT" => "C", "ara0025" => "D", "ara01" => "G"
    )

    # Initialize arrays for the beeswarm (jitter) data and summary statistics.
    data_x = Float64[]          # x (group) positions for each measurement
    data_y = Float64[]          # measured values for each point
    treatment_flag = Int[]      # treatment flag: 1 for wash, 2 for nowash

    # Group summary arrays for crossbar! (one entry per treatment group)
    group_summary_x = Float64[]
    group_means = Float64[]
    group_mins = Float64[]
    group_maxes = Float64[]
    group_treatment = Int[]

    # Prepare x-axis tick positions and labels (one tick per condition)
    condition_ticks = Float64[]
    condition_labels_vec = String[]

    # Assign a unique group ID for each treatment within a condition.
    group_id = 1
    for cond in conditions
        # Extract values using the appropriate column mapping and suffixes.
        wash_vals = values(wash[end, Cols(col_mapping[cond] .* ["2", "3", "4", "5", "6", "7", "8", "9", "10"])])
        nowash_vals = values(nowash[end, Cols(col_mapping[cond] .* ["2", "3", "4", "5", "6", "7", "8", "9", "10"])])
    
        # --- For the wash treatment ---
        n_wash = length(wash_vals)
        append!(data_y, wash_vals)
        append!(data_x, fill(group_id, n_wash))
        append!(treatment_flag, fill(1, n_wash))
    
        push!(group_summary_x, group_id)
        push!(group_means, mean(wash_vals))
        push!(group_mins, minimum(wash_vals))
        push!(group_maxes, maximum(wash_vals))
        push!(group_treatment, 1)
    
        group_id += 1  # Increment group id for nowash treatment
    
        # --- For the nowash treatment ---
        n_nowash = length(nowash_vals)
        append!(data_y, nowash_vals)
        append!(data_x, fill(group_id, n_nowash))
        append!(treatment_flag, fill(2, n_nowash))
    
        push!(group_summary_x, group_id)
        push!(group_means, mean(nowash_vals))
        push!(group_mins, minimum(nowash_vals))
        push!(group_maxes, maximum(nowash_vals))
        push!(group_treatment, 2)
    
        # Place tick at midpoint between the two treatment groups for this condition
        push!(condition_ticks, group_id - 0.5)
        push!(condition_labels_vec, tick_labels[cond])
    
        group_id += 1  # Increment for the next condition
    end

    # Define colors for the treatments – using the first two colors from Makie’s Paired_3 colormap
    wash_color = Makie.to_colormap(:Paired_3)[1]
    nowash_color = Makie.to_colormap(:Paired_3)[2]

    # Create a figure with two slots: main axis at fig[1,1] and legend axis at fig[1,2].
    fig = Figure(size = (6*72, 3.2*72))
    ax = Axis(fig[1, 1],
        xlabel = "",
        ylabel = "BF-biofilm biomass (a.u.)",
        xticks = (condition_ticks, condition_labels_vec)
    )

    # Plot crossbars (with midlines showing the means) for each treatment group.
    midlinecolors = [ (t == 1 ? wash_color : nowash_color) for t in group_treatment ]
    crossbar!(ax, group_summary_x, group_means, group_mins, group_maxes;
              color = :white, midlinecolor = midlinecolors)

    # Plot the beeswarm (jitter) points.
    plt = beeswarm!(ax, data_x, data_y,
        color = treatment_flag,
        algorithm = UniformJitter(),
        strokecolor = :black,
        strokewidth = 1
    )
    plt.colormap[] = [wash_color, nowash_color]

    # Add dummy scatter plots (with no data) to generate legend entries.
    d_wash = scatter!(ax, Float64[], Float64[]; color=wash_color, markersize=8)
    d_nowash = scatter!(ax, Float64[], Float64[]; color=nowash_color, markersize=8)

    # Create the legend in the second axis slot, fig[1,2]
    lg = Legend(fig[1, 2], [d_wash, d_nowash], ["With washing", "Without washing"])

    # General styling for the main axis.
    ax.xticklabelrotation = 45
    ax.xgridvisible = false
    ax.ygridvisible = false
    ax.rightspinevisible = false
    ax.topspinevisible = false

    save("wash_nowash_plot.svg", fig)
end


function sample_timelapses!(vc_data, sp_data)
    fig = Figure(size=(5*72, 2.7*72))
    ax = Axis(fig[1, 1], xlabel="Time (h)", ylabel="BF-biofilm biomass (a.u.)")
    rename!(sp_data, Dict(col => replace(col, r".*/(.*?)_.*" => s"\1") for col in names(sp_data)))
    D39_WT = sp_data[:, Cols("A1", "A2", "A3", "B1", "B2", "B3", "C1", "C2", "C3")]
    D39_dcps = sp_data[:, Cols("A4", "A5", "A6", "B4", "B5", "B6", "C4", "C5", "C6")]
    SV36_WT = sp_data[:, Cols("D1", "D2", "D3", "E1", "E2", "E3", "F1", "F2", "F3")]
    SV36_dcps = sp_data[:, Cols("D4", "D5", "D6", "E4", "E5", "E6", "F4", "F5", "F6")]

	# S. pneumoniae 
    time = 0:0.5:nrow(D39_WT)/2-0.5
    avg = reduce(+, eachcol(D39_WT)) ./ ncol(D39_WT) 
    stdev = dropdims(std(Array(D39_WT), dims=2), dims=2)  
    lines!(ax, time, avg, color=gg_color(1,4), label="D39 wild-type")
    band!(ax, time, avg-stdev, avg+stdev, color=(gg_color(1,4), 0.3))
    scatter!(ax, time, avg, color=gg_color(1,4), marker=:circle,  strokewidth=1)
    
    avg = reduce(+, eachcol(D39_dcps)) ./ ncol(D39_dcps) 
    stdev = dropdims(std(Array(D39_dcps), dims=2), dims=2)  
    lines!(ax, time, avg, color=gg_color(2,4), label=rich("D39 Δ", rich("cps"; font=:italic)))
    band!(ax, time, avg-stdev, avg+stdev, color=(gg_color(2,4), 0.3))
    scatter!(ax, time, avg, color=gg_color(2,4), marker=:circle,  strokewidth=1)
	
    avg = reduce(+, eachcol(SV36_WT)) ./ ncol(SV36_WT) 
    stdev = dropdims(std(Array(SV36_WT), dims=2), dims=2)  
    lines!(ax, time, avg, color=gg_color(3,4), label="SV36 wild-type")
    band!(ax, time, avg-stdev, avg+stdev, color=(gg_color(3,4), 0.3))
    scatter!(ax, time, avg, color=gg_color(3,4), marker=:circle,  strokewidth=1)
	
    avg = reduce(+, eachcol(SV36_dcps)) ./ ncol(SV36_dcps) 
    stdev = dropdims(std(Array(SV36_dcps), dims=2), dims=2)  
    lines!(ax, time, avg, color=gg_color(4,4), label=rich("SV36 ", rich("cps"; font=:italic), superscript("*")))
    band!(ax, time, avg-stdev, avg+stdev, color=(gg_color(4,4), 0.3))
    scatter!(ax, time, avg, color=gg_color(4,4), marker=:circle,  strokewidth=1)

    fig[1,2] = Legend(fig, ax, merge = true, unique = true, framevisible=false, labelsize=12, rowgap=0)

	ax.rightspinevisible = false
	ax.topspinevisible = false
    save("fig1C.svg", fig)
end

function main()
    path_to_sans = "/root/.fonts/cmunss.ttf"
    path_to_sans_ital = "/root/.fonts/cmunsi.ttf"
    set_theme!(fonts = (; bold=path_to_sans, regular = path_to_sans, italic = path_to_sans_ital))
    data_folder = "../../Data/"
    #correction = 0.325^3
    #vols = DataFrame(CSV.File(joinpath(data_folder, "high_res_data.csv")))
    #vols = mapcols(col -> eltype(col) <: Number ? col .* correction : col, vols)
    #vols_sp = DataFrame(CSV.File(joinpath(data_folder, "high_res_data_sp_new.csv")))
    #vols_sp = mapcols(col -> eltype(col) <: Number ? col .* correction : col, vols_sp)
    #vols_kl = DataFrame(CSV.File(joinpath(data_folder, "high_res_data_kleb.csv")))
    #vols_kl = mapcols(col -> eltype(col) <: Number ? col .* correction : col, vols_kl)
    #vc_mass_path = joinpath(data_folder, "vc_fig2C/4x/biomass.csv")
    #pa_pf_mass_path = joinpath(data_folder, "pa_pf_fig2C/4x/biomass.csv")
    #pf_wspf_path = joinpath(data_folder, "pf_wspF_fig2C/4x/biomass.csv")
    #sp_mass_path = joinpath(data_folder, "sp_fig2C/4x/BF_biomass_at_10-h.csv")
    #kleb_mass_path = joinpath(data_folder, "kp_fig2C/4x/biomass.csv")
    #vc_mass = DataFrame(CSV.File(vc_mass_path))
    #pa_pf_mass = DataFrame(CSV.File(pa_pf_mass_path))
    #pf_wspf_mass = DataFrame(CSV.File(pf_wspf_path))
    #sp_mass = DataFrame(CSV.File(sp_mass_path))
    #kleb_mass = last(DataFrame(CSV.File(kleb_mass_path)), 1)
    #fig1B!(vols, vols_sp, vols_kl, vc_mass, pa_pf_mass, pf_wspf_mass, sp_mass, kleb_mass)
    
	#vc_mass_path = "/mnt/e/Brightfield_paper/vc_fig1B/250119_4x_10x_plastic_Drawer1 19-Jan-2025 17-39-44/4x/Numerical data/biomass.csv"
    #sp_mass_path = "/mnt/e/Brightfield_paper/sp_fig1B/4x/Numerical data/biomass.csv"
    #vc_mass = DataFrame(CSV.File(vc_mass_path))
    #sp_mass = DataFrame(CSV.File(sp_mass_path))
	#sample_timelapses!(vc_mass, sp_mass)

    #OD_image_path = "/mnt/e/Brightfield_paper/tests/OD_test.tif"
    #mask_path = "/mnt/e/Brightfield_paper/tests/example_mask.tif"
    #fig1A_OD_image(OD_image_path, mask_path)  
    #fig1A_data_path = data_folder*"fig1A_lineplot.csv"
    #fig1A_lineplot!(fig1A_data_path)
    #focal_plane_path = "/mnt/e/Brightfield_paper/focal_plane/data/Numerical data/biomass.csv"
    #focal_plane_data = DataFrame(CSV.File(focal_plane_path))
    #figS2_focal_plane!(focal_plane_data)  
    #inoculum_path = "/mnt/e/Brightfield_paper/Inoculum/250205_4x_10x_plastic_PMR_Drawer1 05-Feb-2025 13-37-18/4x/Numerical data/biomass.csv"
    #inoculum_data = DataFrame(CSV.File(inoculum_path))
    #fig2_inoculum!(inoculum_data)
    
    correction = 0.325^3
    vols = DataFrame(CSV.File(joinpath(data_folder, "high_res_data.csv")))
    vols = mapcols(col -> eltype(col) <: Number ? col .* correction : col, vols)
    vols_sp = DataFrame(CSV.File(joinpath(data_folder, "high_res_data_sp_new.csv")))
    vols_sp = mapcols(col -> eltype(col) <: Number ? col .* correction : col, vols_sp)
    vols_kl = DataFrame(CSV.File(joinpath(data_folder, "high_res_data_kleb.csv")))
    vols_kl = mapcols(col -> eltype(col) <: Number ? col .* correction : col, vols_kl)
    vc_mass_path = joinpath(data_folder, "vc_fig2C/10x/biomass.csv")
    pa_pf_mass_path = joinpath(data_folder, "pa_pf_fig2C/10x/biomass.csv")
    pf_wspf_path = joinpath(data_folder, "pf_wspF_fig2C/10x/biomass.csv")
    sp_mass_path = joinpath(data_folder, "sp_fig2C/10x/biomass.csv")
    kleb_mass_path = joinpath(data_folder, "kp_fig2C/10x/biomass.csv")
    vc_mass = DataFrame(CSV.File(vc_mass_path))
    pa_pf_mass = DataFrame(CSV.File(pa_pf_mass_path))
    pf_wspf_mass = DataFrame(CSV.File(pf_wspf_path))
    sp_mass = DataFrame(CSV.File(sp_mass_path))[21:21,:]
    kleb_mass = last(DataFrame(CSV.File(kleb_mass_path)), 1)
    figS1A!(vols, vols_sp, vols_kl, vc_mass, pa_pf_mass, pf_wspf_mass, sp_mass, kleb_mass)
    
    #vols = DataFrame(CSV.File(joinpath(data_folder, "high_res_data.csv"))) .* 0.325^3
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
    

    #screen_path = "/mnt/e/Brightfield_paper/Spneumo_screen"
    #screen_plates = [f for f in readdir(screen_path, join=true)]
    #fig3_stats!(screen_plates)
    
    #wash_data_path = "/mnt/e/Brightfield_paper/vc_fig1B/250120_4x_10x_plastic_Drawer1 final/4x/Numerical data/biomass.csv" 
    #nowash_data_path = "/mnt/e/Brightfield_paper/vc_fig1B/250119_4x_10x_plastic_Drawer1 19-Jan-2025 17-39-44/4x/Numerical data/biomass.csv"
    #wash_data = DataFrame(CSV.File(wash_data_path))
    #nowash_data = DataFrame(CSV.File(nowash_data_path))
    #wash_nowash!(wash_data, nowash_data)
end

main()
