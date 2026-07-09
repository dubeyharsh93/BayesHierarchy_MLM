"""
Main module for `MetabolomicsWorkbenchAPI.jl` -- a Julia wrapper to create HTTP requests and interact 
with datasets through the REST API of https://www.metabolomicsworkbench.org/.
"""



module PANSTEATITISstudy

    using MatrixLM, MetabolomicsWorkbenchAPI, DataFrames
    using RecipesBase
    using GLM, StatsBase, Distributions

    include("./utils.jl")    
    export tstat2pval, pval2qval
        
    include("./mLinearModel.jl")    
    export getCoefs, getPermPvals, funStandardize, funStandardize!

    include("./myPlots.jl")    
    export MLMheatmap, mlmheatmap, ConfidencePlot, confidenceplot

end