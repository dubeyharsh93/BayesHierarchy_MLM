#=
Synopsosis: `demog.jl`.

Functions list:


- getDemographicTable:
    Returns a dataframe containing the demography summary table.
    
=#

"""
**getDemographicTable** -*Function*

    getDemographicTable (dfIndividuals::DataFrame) => DataFrame

Returns a dataframe containing the demography summary table.

"""
function getDemographicTable(dfIndividuals::DataFrame)

    # Group by sex
    gdf = groupby(dfIndividuals, :Sex);

    # Get mean values
    mymean(X) = mean(skipmissing(X))
    df1a = DataFrames.combine(gdf, [:Age, :BMI, :SmokingPackYears, :PercentEmphysema] .=> mymean)
    df1a[:,2:end] = round.((df1a[:,2:end]); digits = 1)
    rename!(df1a, Dict(:Age_mymean => "Age", :BMI_mymean => "BMI",
                       :SmokingPackYears_mymean => "SmokingPackYears", :PercentEmphysema_mymean => "PercentEmphysema"));

    # Get standard deviation values
    mystd(X) = std(skipmissing(X))
    df1b = DataFrames.combine(gdf, [:Age, :BMI, :SmokingPackYears, :PercentEmphysema] .=> mystd)
    df1b[:,2:end] = round.((df1b[:,2:end]); digits = 1)
    rename!(df1b, Dict(:Age_mystd => "Age", :BMI_mystd => "BMI",
                       :SmokingPackYears_mystd => "SmokingPackYears", :PercentEmphysema_mystd => "PercentEmphysema"));

    # Join mean and standard deviation values
    dfDem1 = string.(df1a[:,2:end]).*repeat(["("], size(df1a,1),size(df1a,2)-1).* 
             string.(df1b[:,2:end]).*repeat([")"], size(df1a,1),size(df1a,2)-1);
    insertcols!(dfDem1, 1, :Sex => df1a.Sex, :Participants => DataFrames.combine(gdf, nrow)[:,2])

    # Get sum values
    df2a = DataFrames.combine(gdf, [:NHW, :CurrentSmoker, :COPD] .=> sum)

    # Get percentage values
    df2b = round.((df2a[:,2:end]./ dfDem1.Participants).*100, digits = 1)

    # Join sum and percentage values
    dfDem2 = string.(df2a[:,2:end]).*repeat(["("], size(df2a,1),size(df2a,2)-1).* 
             string.(df2b[:,1:end]).*repeat([")"], size(df2a,1),size(df2a,2)-1)
    insertcols!(dfDem2, 1, :Sex => df2a.Sex)
    rename!(dfDem2, Dict(:NHW_sum => "NHW", :CurrentSmoker_sum => "CurrentSmoker",
                       :COPD_sum => "COPD"))
    # Join demographics dataframes
    dfDem = leftjoin(dfDem1, dfDem2, on = :Sex )

    # Pivot table
    dfDem = permutedims(dfDem, 1, "Variable") 

	####################
	# Compute P-values #
	####################

	pval = zeros(8) *1.0
	# grp = groupby(dfAgren, :flc);
	pval[2] = UnequalVarianceTTest(gdf[1].Age[:], gdf[2].Age[:]) |> x -> pvalue(x)
	pval[3] = UnequalVarianceTTest(gdf[1].BMI[:], gdf[2].BMI[:]) |> x -> pvalue(x)
	pval[4] = UnequalVarianceTTest(gdf[1].SmokingPackYears[:], gdf[2].SmokingPackYears[:]) |> x -> pvalue(x)
	pval[5] = UnequalVarianceTTest(gdf[1].PercentEmphysema[:], gdf[2].PercentEmphysema[:]) |> x -> pvalue(x)
	
	pval[6] = freqtable(dfIndividuals, :NHW, :Sex) |> 
		x -> ChisqTest(x) |>
		x -> pvalue(x)
	pval[7] = freqtable(dfIndividuals, :CurrentSmoker, :Sex) |> 
		x -> ChisqTest(x) |>
		x -> pvalue(x)
	pval[8] = freqtable(dfIndividuals, :COPD, :Sex) |> 
		x -> ChisqTest(x) |>
		x -> pvalue(x)

	pval = pval |>
		x -> round.(x, digits = 9) |>
		x -> string.(x)
	pval[1] = " "

	dfDem.Pvalue = pval;
	
    return dfDem
end