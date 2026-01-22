#=
Synopsosis: `utils_copd_spiro.jl`.

Type:
- Cohort: structure containing response, independent and references data

Functions list:
- get_data:
    Extracts data from the specified cohort.

- get_demog_table:
    Returns demographic table in Markdown format.

-get_superclass_table:
    Returns super class table in Markdown format.   
    
- get_subclass_table:
    Returns sub class table in Markdown format.

- getCoefs_copd_spiro:
    Returns coefficients, confidence of interval, T-stats.
    
=#

# using DataFrames, CSV, StatsBase
########
# TYPE #
########

struct CohortData
    # Individuals #
    dfInd::DataFrame

    # Metabolites responses #
    dfMet::DataFrame
    
    # Outcome matrix #
    mY::Matrix{Float64}
     
    # Number of columns in the outcome matrix #
    m::Int64

    # Metabolites annotations #
    dfRef::DataFrame   
end


struct MLMCoeffs
    # Estimated coefficients #
    coef::Matrix{Float64}

    # Confidence of Interval #
    ci::Matrix{Float64}
    
    # T statistics of coefficients #
    tstats::Matrix{Float64}
     
    # Variance of estimated coefficients #
    var::Matrix{Float64}
end

"""
*get_data** -*Function*

    get_data(mycohort::String) => CohortData

Extracts data from the specified cohort and returns a `CohortData` object.

***Arguments***


- `mycohort` name of the directory containing the data of the cohort.

"""
function get_data(mycohort::String)
    ###############
    # Individuals #
    ###############

    fileIndividuals = joinpath(@__DIR__,"..","data","processed",mycohort,
                               mycohort*"_ClinicalCovariates.csv");
    dfInd = CSV.read(fileIndividuals, DataFrame);
    # keep complete cases
    dfInd = dfInd[findall(completecases(dfInd)), :];

    # change Union{Missing, Float64} to Float64
    dfInd.SmokingPackYears = float.(dfInd.SmokingPackYears);
    dfInd.PercentEmphysema = float.(dfInd.PercentEmphysema);
    # names of the categorical variables
    # vIndivCatNames = ["Sex", "NHW", "CurrentSmoker", "COPD", "FinalGold"]

    #########################
    # Metabolites responses #
    #########################

    # Get metabolite file
    fileMeta = joinpath(@__DIR__,"..","data","processed",mycohort,"inl2_Meta.csv");
    dfMet = CSV.read(fileMeta, DataFrame);

    funStandardize!(dfMet, tocenter = true)

    # Order compID columns
    dfMet = dfMet[!, vcat([:SampleID],Symbol.(sort(names(dfMet)[2:end])))]
    # get the outcome matrix 
    mY = Matrix(dfMet[:,2:end]);
    # get number of columns in the outcome matrix
    m = size(mY, 2);

    ###########################
    # Metabolites annotations #
    ###########################	

    # Get reference metabolite file
    fileRef = joinpath(@__DIR__,"..","data","processed",mycohort,"refMeta.csv");
    dfRef = CSV.read(fileRef, DataFrame);

    # Sort according to compID 
    sort!(dfRef, [:CompID]);

    # Keep only references accoring metabolites in the response dataframe
    # - Get CompID from the outcome dataframe
    dfMetID = DataFrame(CompID = names(dfMet)[2:end]);
    # Use join function to filter the reference dataframe
    dfRef = rightjoin(dfRef,dfMetID, on = :CompID);

    return CohortData(dfInd, dfMet, mY, m, dfRef)
end


"""
*get_demog_table** -*Function*

get_demog_table(df::DataFrame) => Markdown.MD

Returns demographic table in Markdown format.

***Arguments***

- `df` dataframe containing the individuals covariates informations.

"""
function get_demog_table(df::DataFrame)
    # get demographic dataframe
    dfDemographic = getDemographicTable(df)	

	return latexify(dfDemographic; env=:mdtable, latex=false)
end


"""
*get_superclass_table** -*Function*

get_superclass_table(df::DataFrame) => Markdown.MD

Returns super class table in Markdown format.

***Arguments***

- `df` dataframe containing the metabolites references attributes informations.

"""
function get_superclass_table(df::DataFrame)
    # freq table
	vFreq = freqtable(df.SuperPathway)
	dfFreq = DataFrame(Class = names(vFreq)[1], Count = vFreq)
	# join dataframes
	dfFreq =leftjoin(rename(sort(unique(df[:,[:SuperPathway, :SuperClassID]])),
					 [:SuperPathway, :SuperClassID].=> [:Class, :ID]), 
                     dfFreq,
					 on = :Class)
    
	return latexify(dfFreq; env=:mdtable, latex=false)
end
   


"""
*get_subclass_table** -*Function*

get_subclass_table(df::DataFrame) => Markdown.MD

Returns sub class table in Markdown format.

***Arguments***

- `df` dataframe containing the metabolites references attributes informations.

"""
function get_subclass_table(df::DataFrame)
	# freq table
	vFreq = freqtable(df.SubPathway)
	dfFreq = DataFrame(SubClass = names(vFreq)[1], Count = vFreq)
	# join dataframes
	dfFreq =leftjoin(rename(unique(df[:,[:SubPathway, :SubClassID]]),
					  		[:SubPathway, :SubClassID].=> [:SubClass, :ID]),
					    dfFreq,
                        on = :SubClass)
	return dfFreq = sort(dfFreq, [:ID])

end


"""
**getCoefs_copd_spiro** -*Function*.

    getCoefs_copd_spiro(Y, X, Z) => MLMCoeffs

Returns coefficients, confidence of interval, T-stats.

"""
function getcoeffs_mlm(Y, X, Z; hasXIntercept = true, hasZIntercept = true)
    coef, ci, tstats, var = getCoefs(Y, X, Z; 
        hasXIntercept = hasXIntercept, hasZIntercept = hasZIntercept);
    
    return MLMCoeffs(coef, ci, tstats, var)    
end


"""
**getCoefs_copd_spiro** -*Function*.

    getCoefs_copd_spiro(Y, X, Z) => MLMCoeffs

Returns coefficients, confidence of interval, T-stats.

"""
function compare_mlm(c::MLMCoeffs, s::MLMCoeffs)
    
    
    
    
    coef, ci, tstats, var = getCoefs(Y, X, Z);
    
    return MLMCoeffs(coef, ci, tstats, var)    
end

"""
**get_modules_coefs** -*Function*.

    get_modules_coefs(chrt::CohortData, X) => MLMCoeffs

Returns coefficients, confidence of interval, T-stats for Z based on modules.

"""
function get_modules_coefs(chrt::CohortData, X)
	######################
	# Metabolite Classes #
	######################
	# nameModules = ["blue", "red", "turquoise", "brown", "yellow","green",
	# 			"magenta", "black", "greenyellow", "purple", "pink"]
	nameModules =  ["pink", "purple", "greenyellow", "black", "magenta", "green",
		"yellow", "brown", "turquoise", "red", "blue"]

	# remove cols with NA class Metabolites in mY
	idx_noNA = findall(chrt.dfRef.SuperPathway .!= "NA")
	mY_noNA = chrt.mY[:, idx_noNA];
	dfRef_noNA = filter(row-> row.SuperPathway != "NA", chrt.dfRef)
	dfZmodule = select(dfRef_noNA, [:CompID, :SubPathway, :SuperPathway])
	
	# dfZmodule = select(dfRef, [:CompID, :SubPathway, :SuperPathway])
	
	# blue	
	idx_blue = findall(occursin.(r"Dicarboxylate|Monohydroxy|Long Chain|Medium Chain|Acyl Carnitines|Endocannabinoid",
            dfZmodule.SubPathway))
	idx_blue = vcat(idx_blue, findall(occursin.(r"Nucleotide", dfZmodule.SuperPathway)))
	
	# red
	idx_red = findall(occursin.(r"Ceramide|Sphingomyelin", dfZmodule.SubPathway))
	
	# turquoise
	idx_turquoise = findall(
        occursin.(r"Tryptophan|Glutamate|Histidine|Glycine|Serine|Threonine|Methionine|Cysteine|SAM|Taurine|Polyamine|Urea|Arginine|Proline|TCA Cycle|Leucine|Isoleucine|Valine",
        dfZmodule.SubPathway))
	idx_turquoise = vcat(idx_turquoise, findall(occursin.(r"Xenobiotics", dfZmodule.SuperPathway)))
	# NOTES: Branched Chain Amino Acids: |Leucine|Isoleucine|Valine
		
	# brown
	idx_brown = findall(occursin.(r"Gamma-glutamyl|Glutamate|Urea|Lysine|Methionine|Phenylalanine|Bile Acid|Acyl Choline|Lysophospholipid|Leucine|Isoleucine|Valine", dfZmodule.SubPathway))
	# NOTES: Branched Chain Amino Acids: |Leucine|Isoleucine|Valine
		
	# yellow
	idx_yellow = findall(occursin.(r"Benzoate|Xanthine Metabolism|Food", dfZmodule.SubPathway))
	
	# green
	idx_green = findall(occursin.(r"Lysophospholipid|Phosphatidylcholine|Phosphatidylinositol|Plasmalogen", dfZmodule.SubPathway))
	
	# magenta
	idx_magenta = findall(occursin.(r"Androgenic|Pregnenolone|Corticosteroid|Progestin", dfZmodule.SubPathway))
	
	# black
	idx_black = findall(occursin.(r"Diacylglycerol|Phosphatidylethanolamine|Acyl Carnitine", dfZmodule.SubPathway))
	
	# greenyellow
	idx_greenyellow = findall(occursin.(r"Cofactors and Vitamins", dfZmodule.SuperPathway))
	
	# purple
	idx_purple = findall(occursin.(r"Acetylated Peptides|Benzoate|Secondary Bile Acid", dfZmodule.SubPathway))
	
	# pink
	idx_pink = findall(occursin.(r"Chemical|Dipeptide|Hemoglobin", dfZmodule.SubPathway))
	
	# idx_modules = (idx_blue, idx_red, idx_turquoise, idx_brown, idx_yellow, idx_green,
	# 				idx_magenta, idx_black, idx_greenyellow, idx_purple, idx_pink);
	idx_modules = (idx_pink, idx_purple, idx_greenyellow, idx_black, idx_magenta, idx_green, idx_yellow, idx_brown, idx_turquoise, idx_red, idx_blue);
	
	
	# Generate Z modules matrix
	# mZmdl = zeros(Int,size(dfRef, 1), length(idx_modules))
	mZmdl = zeros(Int,size(dfRef_noNA, 1), length(idx_modules))
	
	for i in 1:length(idx_modules)
		mZmdl[idx_modules[i] ,i] .= 1
	end

    return getcoeffs_mlm(mY_noNA, X, mZmdl), mZmdl	
end