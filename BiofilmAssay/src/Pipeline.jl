module Pipeline

export pipeline

using Gtk4
using JSON
using NaturalSort
using TiffImages: load, save
using IntegralArrays
using IntervalSets
using Images 
using ImageView
using StatsBase
using DataFrames
using CSV: write
using CoordinateTransformations
using AbstractFFTs
using Compat
using FFTW

include("GUI.jl")
include("Analysis.jl")

function pipeline()
    #GUI_main()
    analysis_main()
end

end # module
