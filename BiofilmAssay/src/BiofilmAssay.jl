include("Pipeline.jl")
using .Pipeline

function julia_main()#::Cint
	try
        Pipeline.pipeline()
        return 0 # if things finished successfully
	catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
	end
end
julia_main()
#end 