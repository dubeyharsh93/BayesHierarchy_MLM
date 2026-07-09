#=
Author: Gregory Farage
Date: 2022-09-01
Synopsosis: `wrangling_util.jl` contains function and to wrangle metabolomics data.

Functions list:

- print_variables_missing:
    Print variables' name containing missing values and shows how many missing data.

- print_df_size
    Print dataframe size.

- summary_variables_missing
    Returns a dataframe that summarizes missing data in the given dataframe.

- getMetaRef: 
    Returns a data frame containing reference ID for the metabolomics.

- getMetaboData:
    Returns a dataframe containing metabolites values for each individual sample,
    where the rows are metabolites and the columns are individual samples.

- getTransposedMetaboData:
    Returns the transposed dataframe containing metabolites values for each individual
    sample, where the columns are metabolites and the rows are individual samples.

=#

# ############
# Functions #
# ############

"""
    print_variables_missing(df::DataFrame) => String

Print variables' name containing missing values and shows how many missing data.

"""
function print_variables_missing(df::DataFrame)
    
    vMissing = map(eachcol(df)) do col
                   sum(ismissing.(col))
               end
    
    idxColMiss = findall(vMissing .!= 0)
    
    for i in idxColMiss
        println("$(names(df)[i]) contains $(vMissing[i]) missing values.")
    end
    
end


"""
    print_df_size(df::DataFrame) => String

Print dataframe size.

"""
function print_df_size(df::DataFrame)
    
    println("The dataframe contains $(size(df, 1)) rows and $(size(df, 2)) columns")
    
end


"""
    summary_variables_missing(df::DataFrame) => DataFrame

Returns a dataframe that summarizes missing data in the given dataframe.

"""
function summary_variables_missing(df::DataFrame)
    
    vMissing = map(eachcol(df)) do col
                   sum(ismissing.(col))
               end
    
    if sum(vMissing) != 0
        dfSummary = DataFrame(Variables = names(df), 
                          Missing_Count = vMissing, 
                          Missing_Percentage = round.(100 .*vMissing./size(df,1), digits=2))
        return dfSummary
    else
        println("No missing data.")
    end
end









"""
**getMetaboRef** -*Function*.

getMetaboRef(strFile::String, isTextFile::Bool = false) => DataFrame

Returns a data frame containing reference ID for the metabolomics.

Arguments
- `strFile` contains the file name or the content of the file as a string.
- `isTextFile` indicates if `strFile` contains file name path (`true`) or the content of the file (`false`). 

"""
function getMetaboRef(filename::String, isTextFile::Bool = false)
    
    # localize targeted dataset in the text file
    if isTextFile
        vIOFile = readlines(fileName);
    else
        vIOFile = split(fileName, '\n');
    end
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

getMetaboData(fileName::String; addString2colNames::String = "") => DataFrame

Returns a dataframe containing metabolites values for each individual sample, where the rows are metabolites and the columns are individual samples.

Arguments

- `addString2colNames` contains a string (default is "") to prefix to potential numerical column names.

"""
function getMetaboData(fileName::String; addString2colNames::String = "")
    
    # localize targeted dataset in the text file
    vIOFile = readlines(fileName);
    startLine = findall(occursin.(r"MS_METABOLITE_DATA_START", vIOFile))[1]+1;
    endLine =  findall(occursin.(r"MS_METABOLITE_DATA_END", vIOFile))[1]-1;
    
    # read column names
    dfMetaboHeader = CSV.read(fileName, DataFrame; header=false, skipto = startLine, limit = (1), delim = '	');
    vHeader = string.(collect(dfMetaboHeader[1,:]));
    vHeader[1] = "metabolite_name";
    
    if !isempty(addString2colNames)
        vHeader[2:end] .= addString2colNames.*vHeader[2:end]   
    end
    
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

Splits a vector of string into a matrix of string according to the main separator character, `sepMain`, and the sublevel separator charater, `sepSub`.

---

splitVecString(s::String, sepMain::String = "|", sepSub::String = ":") => Matrix{String}

Splits a string into a matrix of string according to the main separator character, `sepMain`, and the sublevel separator charater, `sepSub`.

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

function splitVecString(s::String, sepMain::String = "|", sepSub::String = ":")
    
    numCol = sum(occursin.(sepMain, split(s)))+1
    matString = repeat([""], 1, numCol)

    tmp = replace(s, sepMain => sepSub)
    matString[1, :] = String.(reshape(strip.(split(tmp, sepSub)),2, numCol))[2, :]
    
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
                    x ->  reshape(strip.(split(x, sepSub)),2, numCol)[1, :] |>
                    x ->  Symbol.(replace.(x, " " => "_", ":" => "_", "-" => "_"))
    return covarNames
end


"""
**getIndivData** -*Function*.

getIndivData(fileName::String) => DataFrame

Returns a dataframe containing individuals covariates.

"""
function getIndivData(fileName::String)#, isTextFile::Bool = false)
    
    # # localize targeted dataset in the text file
    # if isTextFile
    #     vIOFile = readlines(fileName);
    # else
    #     vIOFile = split(fileName, '\n');
    # end
           vIOFile = split(fileName, '\n');
    # localize targeted dataset in the text file
    startLine = findall(occursin.(r"#SUBJECT_SAMPLE_FACTORS:", vIOFile))[1]+1;
    endLine =  findall(occursin.(r"#COLLECTION", vIOFile))[1]-1;
    
    # read data set    
    dfIndividuals = CSV.read(map(IOBuffer, string.(vIOFile[(startLine):(endLine)])), DataFrame; 
        header=false, missingstring=["NA"], delim = '	');
    # dfIndividuals = CSV.read(fileName, DataFrame; 
    #     header=false, missingstring=["-", "NA"],
    #     skipto = (startLine), limit = (endLine-startLine+1), 
    #     delim = '	');
    
    
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
    matIndivCovariates = replace(matIndivCovariates, "NA"=> missing, "-"=> "0")
    # matIndivCovariates = passmissing(parse).(Float64, matIndivCovariates)
    newColNames = getCovariatesNames(String(dfIndividuals[1,rmNames[2]]), ";", "=")
    for i in 1:length(newColNames)
       dfIndividuals[:, newColNames[i]] =   matIndivCovariates[:, i];
    end   
            
    # clean columns name and rename appropriately
    select!(dfIndividuals, Not(rmNames))
    rename!(dfIndividuals, Dict(:Column2 => :subject, :Column3 => :sample));
    
    return dfIndividuals
end



function new_getIndivData(fileText::String)#, isTextFile::Bool = false)
    
    # # localize targeted dataset in the text file
    # if isTextFile
    #     vIOFile = readlines(fileName);
    # else
    #     vIOFile = split(fileName, '\n');
    # end
    vIOFile = split(fileText, '\n');
    # localize targeted dataset in the text file
    startLine = findall(occursin.(r"#SUBJECT_SAMPLE_FACTORS:", vIOFile))[1]+1;
    endLine =  findall(occursin.(r"#COLLECTION", vIOFile))[1]-1;
    
    # read data set    
    dfIndividuals = CSV.read(map(IOBuffer, string.(vIOFile[(startLine):(endLine)])),
                                DataFrame; header=false, delim = '	');
                                # missingstring=["-", "NA"], delim = '	');
    # dfIndividuals = CSV.read(fileName, DataFrame; 
    #     header=false, missingstring=["-", "NA"],
    #     skipto = (startLine), limit = (endLine-startLine+1), 
    #     delim = '	');
    
    # Remove first column 
    select!(dfIndividuals, Not([1]));
    
    # get names of columns containing categorical and continuous values
    rmNames = Symbol.(names(dfIndividuals)[end-1:end])
    
    
    
    # part 1: corvariates containing string type values 
    matIndivCovariates = splitVecString(String.(dfIndividuals[:,rmNames[1]]))
    matIndivCovariates = replace(matIndivCovariates, "NA"=> missing)
    newColNames = getCovariatesNames(String(dfIndividuals[1,rmNames[1]]))
    for i in 1:length(newColNames)
       dfIndividuals[:, newColNames[i]] = matIndivCovariates[:, i];
    end

#     # part 2: covariates containing numerical values
#     matIndivCovariates = splitVecString(String.(dfIndividuals[:,rmNames[2]]), ";", "=")
#     matIndivCovariates = replace(matIndivCovariates, "NA"=> missing, "-"=> missing)
#     # matIndivCovariates = passmissing(parse).(Float64, matIndivCovariates)
#     newColNames = getCovariatesNames(String(dfIndividuals[1,rmNames[2]]), ";", "=")
#     for i in 1:length(newColNames)
#        dfIndividuals[:, newColNames[i]] =   matIndivCovariates[:, i];
#     end   

#     # clean columns name and rename appropriately
#     select!(dfIndividuals, Not(rmNames))
#     rename!(dfIndividuals, Dict(:Column2 => :subject, :Column3 => :sample));

    return dfIndividuals
end

function json_new_getIndivData(fileText::String)#, isTextFile::Bool = false)
    
    jsonFile = JSON3.read(fileText)
    
    dfIndividuals = DataFrame(jsonFile[:SUBJECT_SAMPLE_FACTORS])
    
    DataFrame.(dfSample.Factors);
    dfIndividuals =  hcat(select(dfIndividuals, ["Sample ID"]),
                            jsonbuild_df_names(DataFrame.(dfSample.Factors)),
                            jsonbuild_df_names(DataFrame.(dfSample."Additional sample data")))
    
    # Remove first column 
    # select!(dfIndividuals, Not([1]));
    
    # get names of columns containing categorical and continuous values
    # rmNames = Symbol.(names(dfIndividuals)[end-1:end])
    
    
    
    # part 1: corvariates containing string type values 
    # matIndivCovariates = splitVecString(String.(dfIndividuals[:,rmNames[1]]))
#     matIndivCovariates = replace(matIndivCovariates, "NA"=> missing)
#     newColNames = getCovariatesNames(String(dfIndividuals[1,rmNames[1]]))
#     for i in 1:length(newColNames)
#        dfIndividuals[:, newColNames[i]] = matIndivCovariates[:, i];
#     end

#     # part 2: covariates containing numerical values
#     matIndivCovariates = splitVecString(String.(dfIndividuals[:,rmNames[2]]), ";", "=")
#     matIndivCovariates = replace(matIndivCovariates, "NA"=> missing, "-"=> missing)
#     # matIndivCovariates = passmissing(parse).(Float64, matIndivCovariates)
#     newColNames = getCovariatesNames(String(dfIndividuals[1,rmNames[2]]), ";", "=")
#     for i in 1:length(newColNames)
#        dfIndividuals[:, newColNames[i]] =   matIndivCovariates[:, i];
#     end   

#     # clean columns name and rename appropriately
#     select!(dfIndividuals, Not(rmNames))
#     rename!(dfIndividuals, Dict(:Column2 => :subject, :Column3 => :sample));

    return dfIndividuals
end





# """
# **get_variable_names** -*Function*.

# get_variable_names(strNames::String, sepMain::String = "|", sepSub::String = ":") => Vector{Symbol}

# Return variables name as a `Symbol` vector. Extracts variable names according to the main separator character, `speMain`, and the sublevel separator charater, `sepSub`. If the following characters {' ', ':', '-'} are included in the varaiable names, they will be replaced by the character '_'.   

# ---
# get_variable_names(strNames::Vectro{String}, sepMain::String = "|", sepSub::String = ":") => Vector{Symbol}

# Return the union of the unique variables name from a vector of `String`. Extracts variable names according to the main separator character, `speMain`, and the sublevel separator charater, `sepSub`. If the following characters {' ', ':', '-'} are included in the varaiable names, they will be replaced by the character '_'.   


# """
# function get_variables_names(strNames::String, sepMain::String = "|", sepSub::String = ":")
    
#     numCol = sum(occursin.(sepMain, split(strNames)))+1
#     varNames = strNames |>
#                     x ->  replace(x, sepMain => sepSub) |>
#                     x ->  replace(x, "." => "_") |>
#                     x ->  reshape(strip.(split(x, sepSub)),2, numCol)[1, :] |>
#                     x ->  Symbol.(replace.(x, " " => "_", ":" => "_", "-" => "_"))
#     return varNames
# end

# function get_variables_names(vStrNames::Vector{String}, sepMain::String = "|", sepSub::String = ":")
    
#     varNames = []
#     for i  in 1:length(vStrNames)
#         varNames = unique(vcat(varNames, get_variables_names(vStrNames[i])))
#     end
    
#     return varNames
# end

# function jsonget_variables_names(vDF::Vector{DataFrame}, sepMain::String = "|", sepSub::String = ":")

#     varNames = unique(reduce(vcat, names.(vDF)))
        
#     return varNames
# end

# function jsonbuild_df_names(vDF::Vector{DataFrame}, sepMain::String = "|", sepSub::String = ":")
    
#     refVarNames = jsonget_variables_names(vDF, sepMain, sepSub)
        
#     for i in 1:length(vDF)
#         vMissNames = setdiff(refVarNames, names(vDF[i]))
#         vDF[i] = hcat(vDF[i], DataFrame(repeat(["NA"], 1,length(vMissNames)), vMissNames))
#     end
    
#     df = reduce(vcat,vDF)
    
#     return df
# end
