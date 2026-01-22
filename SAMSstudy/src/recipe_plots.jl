# using Plots

# create a type MLMheatmap for dispatch.
# An object of type `MLMheatmap` has the field args which is the tuple of arguments 
# the plot function is invoked with, which can be either `mlmheatmap(x,...)` 
# or plot(x, seriestype = :mlmheatmap).
# @userplot MLMheatmap

mutable struct MLMheatmap{AbstractType}
    args::Any                                      
end
    
mlmheatmap(args...; kw...) = RecipesBase.plot(MLMheatmap{typeof(args[1])}(args); kw...)


@recipe function f(x::MLMheatmap; annotationargs = ())
    # get the input arguments stored in x.args, here we expect only one  
    mB = x.args[1] 
    mBval = x.args[2] 
    
    # return error message if the input arugment is different from `AbstractMatrix`
    typeof(mB) <: AbstractMatrix || error("Pass a Matrix as the arg to heatmap")

    # get size of input matrix
    rows, cols = size(mB)
    
    # turn off the background grid
    grid := false                      
    
    # remove x-axis and y-axis
    xaxis := false
    yaxis := false
    
#     myXlabel = x.args[1]
#     myYlabel = x.args[2]
    
    ##########################
    # SERIES SHOWING HEATMAP #
    ##########################
    
    @series begin                      
        seriestype := :heatmap   
        # myXlabel
        # myYlabel
        mB
    end

    #################################
    # SERIES SHOWING VERTICAL LINES #
    #################################
    
    # each line is added as its own series
    for i in 0:cols         
        @series begin
            seriestype := :path
            primary := false          # to avoid showing the lines in a legend
            linecolor --> :black
            [i, i] .+ 0.5, [0, rows] .+ 0.5  # start[x,y], end[x,y] values of lines
        end
    end
    
    ###################################
    # SERIES SHOWING HORIZONTAL LINES #
    ###################################    

    for i in 0:rows
        @series begin
            seriestype := :path
            primary := false
            linecolor --> :black
            [0, cols] .+ 0.5, [i,i] .+ 0.5
        end
    end

    ############################
    # SERIES ANNOTATING VALUES #
    ############################
    
    
    # @series begin
    #     seriestype := :scatter
    #     markeralpha --> 0.0
    #     markersize --> 1
    #     seriescolor := RGBA(0,0,0,0.0)        # do
    #     series_annotations := text.(vec(round.(mBval,digits=2)), annotationargs...)
    #     # series_annotations := text.(1:(cols*rows), annotationargs...)
    #     primary := false
    #     repeat(1:cols, inner = rows), repeat(1:rows, outer = cols)
    # end
end

# ################## 
# CONFIDENCE PLOT #
# ##################

mutable struct ConfidencePlot{AbstractType}
    args::Any                                      
end

confidenceplot(args...; kw...) = RecipesBase.plot(ConfidencePlot{typeof(args[1])}(args); kw...)

@recipe function f(h::ConfidencePlot) 
    # check types of the input arguments
    if length(h.args) != 3 || !(typeof(h.args[1]) <: AbstractVector) ||
        !(typeof(h.args[2]) <: AbstractVector) || !(typeof(h.args[3]) <: AbstractVector)
        error("Confidence Plots should be given three vectors.  Got: $(typeof(h.args))")
    end
    # Note: if confidence or interval not symmetric, then ε should be a vector of tuple.
    # collect(zip(vCI1, vCI2));
    
    x, y, ε = h.args
    x_vline = [0]
    v_significant = (x+ε).*(x-ε) .> 0

    # set a default value for an attribute with `-->`
    # xlabel --> "Effect size"
    # ylabel --> "Variables"
    markershape --> :rect
         
    # set up the subplots
    legend := false
    link := :both
    # framestyle := [:none :axes :none]
    # yaxis := false 
    y_foreground_color_axis := :white
    y_foreground_color_border := :white
    ylims := (-0.5,length(y)+0.25)
    grid := false


    # main confidence plot
    @series begin
        seriestype := :scatter
        xerror := ε
        linecolor := :black#nothing
        # get the seriescolor passed by the user
        c = get(plotattributes, :seriescolor, :auto)
        # highlight big errors, otherwise use the user-defined color
        # markercolor --> ifelse.(v_significant, :blue, c)
        zcolor --> v_significant
        color -->  cgrad([:lightblue, :blue])
        
        x, y
    end
    
    # vertical red line at x = 0
    @series begin
        seriestype := :vline
        linecolor := :red
        primary := false
        alpha := 0.5
        x_vline
    end
end
