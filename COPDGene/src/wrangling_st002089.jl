#=
Author: Gregory Farage
Date: 2022-03-21
Synopsosis: `preprocessing.jl` contains function and to wrangle and preprocess metabolomics data.
Data set information:
- METABOLOMICS WORKBENCH Suni16_20220208_123148 DATATRACK_ID:3071 STUDY_ID:ST002089 ANALYSIS_ID:AN003410 PROJECT_ID:PR000907
- METABOLOMICS WORKBENCH Suni16_20220208_123148 DATATRACK_ID:3071 STUDY_ID:ST002089 ANALYSIS_ID:AN003411 PROJECT_ID:PR000907
- METABOLOMICS WORKBENCH Suni16_20220208_123148 DATATRACK_ID:3071 STUDY_ID:ST002089 ANALYSIS_ID:AN003412 PROJECT_ID:PR000907
- METABOLOMICS WORKBENCH Suni16_20220208_123148 DATATRACK_ID:3071 STUDY_ID:ST002089 ANALYSIS_ID:AN003413 PROJECT_ID:PR000907

Functions list:

- getMetaRef: 
    Returns a data frame containing reference ID for the metabolomics.

- getMetaboData:
    Returns a dataframe containing metabolites values for each individual sample,
    where the rows are metabolites and the columns are individual samples.

- getTransposedMetaboData:
    Returns the transposed dataframe containing metabolites values for each individual
    sample, where the columns are metabolites and the rows are individual samples.

=#

#############
# Functions #
#############

"""
**getMetaboRef** -*Function*.

getMetaboRef(fileName::String) => DataFrame

Returns a data frame containing reference ID for the metabolomics.

"""
function getMetaboRef(fileName::String)
    
    # localize targeted dataset in the text file
    vIOFile = readlines(fileName);
    startLine = findall(occursin.(r"METABOLITES_START", vIOFile))[1]+1;
    endLine =  findall(occursin.(r"METABOLITES_END", vIOFile))[1]-1;
    
    # read column names
    dfHeaderRef = CSV.read(fileName, DataFrame; header=false, skipto = startLine, limit = (1), delim = '	');
    
    # read data set
    dfMetaboRef = CSV.read(fileName, DataFrame; 
        header=false, skipto = (startLine+1), limit = (endLine-startLine), delim = '	');
    
    rename!(dfMetaboRef, Symbol.(collect(dfHeaderRef[1, :])));

    return dfMetaboRef
end


"""
**getMetaboData** -*Function*.

getMetaboData(fileName::String) => DataFrame

Returns a dataframe containing metabolites values for each individual sample, where the rows are metabolites and the columns are individual samples.

"""
function getMetaboData(fileName::String)
    
    # localize targeted dataset in the text file
    vIOFile = readlines(fileName);
    startLine = findall(occursin.(r"MS_METABOLITE_DATA_START", vIOFile))[1]+1;
    endLine =  findall(occursin.(r"MS_METABOLITE_DATA_END", vIOFile))[1]-1;
    
    # read column names
    dfMetaboHeader = CSV.read(fileName, DataFrame; header=false, skipto = startLine, limit = (1), delim = '	');
    vHeader = String.(collect(dfMetaboHeader[1,:]));
    vHeader[1] = "metabolite_name";
    
    # read data set 
    dfMetaboData = CSV.read(fileName, DataFrame; 
        header=false, missingstring=["-", "NA"],
        skipto = (startLine+2), limit = (endLine-startLine-1), 
        delim = '	');
    
    rename!(dfMetaboData, Symbol.(vHeader));

    return dfMetaboData
end


"""
**getTransposedMetaboData** -*Function*.

getTransposedMetaboData(fileName::String; colNamesRef::String = "CHEM_ID") => DataFrame

Returns the transposed dataframe containing metabolites values for each individual sample, where the columns are metabolites and the rows are individual samples.

- colNamesRef indicates what columns names are used to reference the metabolite names. By default, Chem ID is used ("CHEM-ID"); possible other id formats are "LIB_ID", "COMP_ID" and "CHRO_LIB_ENTRY_ID".

"""
function getTransposedMetaboData(fileName::String; colNamesRef::String = "CHEM_ID")
    dfMetaboRef = getMetaboRef(fileName)
    dfMetaboData = getMetaboData(fileName)
    
    df = dfMetaboRef |>
        x -> select(x, [1, 2]) |>
        x -> leftjoin(x, dfMetaboData, on = :metabolite_name) |>
        x -> select(x, Not([1]))
    df[!, 1] = string.(df[!, 1])
    df = permutedims(df, 1)

    rename!(df, Dict(Symbol.(colNamesRef) => "Samples_ID"))

    return df
end

"""
**splitVecString** -*Function*.

splitVecString(vec::Vector{String}, sepMain::String = "|", sepSub::String = ":") => Matrix{String}

Splits a vector of string into a matrix of string according to the main separator character, `speMain`, and the sublevel separator charater, `sepSub`.

"""
function splitVecString(vec::Vector{String}, sepMain::String = "|", sepSub::String = ":")
    n = length(vec)
    numCol = sum(occursin.(sepMain, split(vec[1])))+1
    matString = repeat([""], n, numCol)

    for i in 1:n
        tmp = replace(vec[i], sepMain => sepSub)
        matString[i, :] = String.(reshape(strip.(split(tmp, sepSub)),2, numCol))[2, :]
    end
    
    return matString
end


"""
**getCovariatesNames** -*Function*.

getCovariatesNames(strNames::String, sepMain::String = "|", sepSub::String = ":") => Vector{Symbol}

Extracts covariates names according to the main separator character, `speMain`, and the sublevel separator charater, `sepSub`.

"""
function getCovariatesNames(strNames::String, sepMain::String = "|", sepSub::String = ":")
    
    numCol = sum(occursin.(sepMain, split(strNames)))+1
    covarNames = strNames |>
                    x ->  replace(x, sepMain => sepSub) |>
                    x ->  replace(x, "." => "_") |>
                    x -> Symbol.(reshape(strip.(split(x, sepSub)),2, numCol))[1, :]
    
    return covarNames
end


"""
**getIndivData** -*Function*.

getIndivData(fileName::String) => DataFrame

Returns a dataframe containing individuals covariates.

"""
function getIndivData(fileName::String)
    
    # localize targeted dataset in the text file
    vIOFile = readlines(fileName);
    startLine = findall(occursin.(r"#SUBJECT_SAMPLE_FACTORS:", vIOFile))[1]+1;
    endLine =  findall(occursin.(r"#COLLECTION", vIOFile))[1]-1;
    
    # read data set    
    dfIndividuals = CSV.read(fileName, DataFrame; 
        header=false, missingstring=["-", "NA"],
        skipto = (startLine), limit = (endLine-startLine+1), 
        delim = '	');
    select!(dfIndividuals, Not([1]));
    
    # get names of columns to remove
    rmNames = Symbol.(names(dfIndividuals)[end-1:end])
    
    # part 1: corvariates containing string type values 
    matIndivCovariates = splitVecString(String.(dfIndividuals[:,rmNames[1]]))
    matIndivCovariates = replace(matIndivCovariates, "NA"=> missing)
    newColNames = getCovariatesNames(String(dfIndividuals[1,rmNames[1]]))
    for i in 1:length(newColNames)
       dfIndividuals[:, newColNames[i]] = matIndivCovariates[:, i];
    end
    
    # part 2: covariates containing numerical values
    matIndivCovariates = splitVecString(String.(dfIndividuals[:,rmNames[2]]), ";", "=")
    matIndivCovariates = replace(matIndivCovariates, "NA"=> missing, "-"=> missing)
    matIndivCovariates = passmissing(parse).(Float64, matIndivCovariates)
    newColNames = getCovariatesNames(String(dfIndividuals[1,rmNames[2]]), ";", "=")
    for i in 1:length(newColNames)
       dfIndividuals[:, newColNames[i]] =   matIndivCovariates[:, i];
    end   
            
    # clean columns name and rename appropriately
    select!(dfIndividuals, Not(rmNames))
    rename!(dfIndividuals, Dict(:Column2 => :subject, :Column3 => :sample));
    
    return dfIndividuals
end