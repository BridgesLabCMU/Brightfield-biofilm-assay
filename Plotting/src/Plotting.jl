using Makie 
using CairoMakie
using LsqFit
using TiffImages
using CSV
using StatsBase
using DataFrames
CairoMakie.activate!(type="svg")
path_to_sans = "/mnt/c/Downloads/cmunss.ttf"
set_theme!(fonts = (; bold=path_to_sans, regular = path_to_sans))

function fig1B!(Vc_biovolumes, Vc_biomasses, Pa_biovolumes, Pa_biomasses, Sp_biovolumes, Sp_biomasses)
    # Biovolume (x-axis) vs. brightfield biomass (y-axis)
    # 3 panels: V.c., P.a., S.p.
    # Scatter + best-fit
    Vc_mass_averages = [mean(Vc_biomasses[1:5]), mean(Vc_biomasses[6:10]), mean(Vc_biomasses[11:15]), 
                        mean(Vc_biomasses[16:20]), mean(Vc_biomasses[21:25])]
    Pa_mass_averages = 
    Sp_mass_averages =
    Vc_mass_stds = [std(Vc_biomasses[1:5]), std(Vc_biomasses[6:10]), std(Vc_biomasses[11:15]), 
                    std(Vc_biomasses[16:20]), std(Vc_biomasses[21:25])]
    Pa_mass_stds =
    Sp_mass_stds =
    @. model(x, p) = p[1]*x + p[2] 

    p0 = [(Vc_biomasses[2]-Vc_biomasses[1])/(Vc_biovolumes[2]-Vc_biovolumes[1]), Vc_biomasses[1]]
    lb = [0.0, 0.0]
    ub = [Inf, Inf]

    fit = curve_fit(model, x, y, p0, lower=lb, upper=ub)
    pstar = coef(fit)
	xbase = collect(range(minimum(x), maximum(x), 100))

    fig = Figure(size=(300,300))
    ax = Axis(fig[1, 1], xlabel="Total biofilm biovolume (µm³)", ylabel="Biofilm biomass (a.u.)")
    scatter!(ax, Vc_biovolumes, Vc_biomasses, color=:black)
	lines!(ax, xbase, model.(xbase, (pstar,)), color="black")
    save("fig1B_Vc.svg", fig)
end

function fig1A_OD_image(OD_image_path)
    image = Float64.(TiffImages.load(OD_image_path))
    fig = Figure()
    ax = Axis(fig[1, 1])
    hm = heatmap!(ax, image, colormap=:plasma)
	Colorbar(fig[:, end+1], hm)
    save("fig1A_OD_image.svg", fig)
end

function fig1A_lineplot!(fig1A_data_path)
    data = DataFrame(CSV.File(fig1A_data_path))
    fig = Figure(size=(3.5*72, 3*72))
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
    fig
    #save("fig1A_lineplot.svg", fig)
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

function main()
    data_folder = "../../Data/"
    OD_image_path = "/mnt/f/Brightfield_paper/tests/OD_test.tif"
    fig1A_data_path = data_folder*"fig1A_lineplot.csv"
    fig1A_lineplot!(fig1A_data_path)
end

main()
