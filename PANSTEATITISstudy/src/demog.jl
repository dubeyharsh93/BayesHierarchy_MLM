#=
Synopsosis: `demog.jl`.

Functions list:


- getDemographicST001052:
    Returns a dataframe containing the demography summary table.
    
=#
using MetabolomicsWorkbenchAPI


function xtrctHistScore(vecH) 
    idxnotmissing = findall(.!(ismissing.(vecH)));
    vScore = Vector{Union{Missing, Int}}(undef, length(vecH));
    for i in idxnotmissing
        if isdigit(vecH[i][1])
           vScore[i] = parse(Int, vecH[i][1]) 
        end
    end
    
    return vScore
end

"""
**getDemographicST001052** -*Function*

getDemographicST001052 () => DataFrame

Returns a dataframe containing the demography summary table.

"""
function getDemographicST001052(str_study::String = "ST001052")

    dfIndividuals =  fetch_samples(str_study);

    # dfIndividuals = ifelse.(dfIndividuals .== "-", "0", dfIndividuals)
    dfIndividuals = ifelse.(dfIndividuals .== "-", missing, dfIndividuals);
    
    select!(dfIndividuals, Symbol.(["Sample ID",
                                "Group",
                                "Gender",
                                "Annuli",
                                "Age",
                                "WEIGHT (KG)",
                                "LENGTH (CM)",
                                "Histology Adipose",
    ]));
    
      
    rename!(dfIndividuals, Dict(
        :Group => "Status",
        :Gender => "Sex",
        Symbol("Histology Adipose") => "Histological Score",
    ));
    
    
    dfIndividuals.:"Histological Score" .= xtrctHistScore(dfIndividuals.:"Histological Score");
    
    # filter complete cases
    idxComplete = findall(completecases(dfIndividuals))
    dfIndividuals = dfIndividuals[idxComplete, :]
    
    idxDiseased = findall(occursin.("D", dfIndividuals.Status)) ;
    idxHealthy = findall(occursin.("H", dfIndividuals.Status));
    
    dfIndividuals.Status[idxDiseased] .= "Diseased";
    dfIndividuals.Status[idxHealthy] .= "Healthy";
    
    idxMale = findall(occursin.("M", dfIndividuals.Sex)) ;
    idxFemale = findall(occursin.("F", dfIndividuals.Sex));
    
    dfIndividuals.Sex[idxMale] .= "Male";
    dfIndividuals.Sex[idxFemale] .= "Female";
    
    
    # insertcols!(dfIndividuals, 3, :GroupStatus => occursin.("D", dfIndividuals.Status));
    
    # Categorical Variables
    catVar = ["Status", "Sex"]
    vDFcat = Vector{DataFrame}(undef, length(catVar));
    
    for i in 1:length(catVar)
        dfCat = combine(groupby(dfIndividuals, [Symbol(catVar[i])]), nrow => :Count)
        sort!(dfCat, [Symbol(catVar[i])]) 
        dfCat = DataFrame(vcat([names(dfCat)[1] " "], Matrix(dfCat)), ["Clinical Features", "Count/ mean(SD)"])
        vDFcat[i] = dfCat 
    end
    
    # Continuous Variables
    contVar = ["Annuli", "Age", "WEIGHT (KG)", "LENGTH (CM)", "Histological Score"]
    vDFcont = Vector{DataFrame}(undef, length(contVar));
    
    for i in 1:length(contVar)
        # calculate mean
        vVar = dfIndividuals[:,Symbol(contVar[i])];
        idxNotMiss = findall(.!ismissing.(vVar))
        vVar = parse.(Float32, string.(vVar[idxNotMiss]))
        myMean = vVar |> mean |> x->round(x, digits = 2)
        # calculate SD
        myStd = vVar |> std |> x->round(x, digits = 2)
    
        dfCont  = DataFrame([contVar[i] string(myMean,"(", myStd, ")")], ["Clinical Features", "Count/ mean(SD)"])
        vDFcont[i] = dfCont 
    end
    vDF = vcat(vDFcat, vDFcont)
    
    dfDem = reduce(vcat, vDF)
    
    return dfDem
end



"""
**getDemographicST001059** -*Function*

getDemographicST001059 () => DataFrame

Returns a dataframe containing the demography summary table.

"""
function getDemographicST001059(str_study::String = "ST001059")

    dfIndividuals =  fetch_samples(str_study);

    # dfIndividuals = ifelse.(dfIndividuals .== "-", "0", dfIndividuals)
    dfIndividuals = ifelse.(dfIndividuals .== "-", missing, dfIndividuals);
    
    select!(dfIndividuals, Symbol.(["Sample ID",
                                "Group",
                                "Gender",
                                "Annuli",
                                "Age",
                                "WEIGHT (KG)",
                                "LENGTH (CM)",
                                "VET SCORE (Adipose)",
    ]));
    
    # filter complete cases
    idxComplete = findall(completecases(dfIndividuals))
    dfIndividuals = dfIndividuals[idxComplete, :]
    
    idxDiseased = findall(occursin.("D", dfIndividuals.Group)) ;
    idxHealthy = findall(occursin.("H", dfIndividuals.Group));
    
    dfIndividuals.Group[idxDiseased] .= "Diseased";
    dfIndividuals.Group[idxHealthy] .= "Healthy";
    
    idxMale = findall(occursin.("M", dfIndividuals.Gender)) ;
    idxFemale = findall(occursin.("F", dfIndividuals.Gender));
    
    dfIndividuals.Gender[idxMale] .= "Male";
    dfIndividuals.Gender[idxFemale] .= "Female";
    
    
    # insertcols!(dfIndividuals, 3, :GroupStatus => occursin.("D", dfIndividuals.Group));
    
    # Categorical Variables
    catVar = ["Group", "Gender"]
    vDFcat = Vector{DataFrame}(undef, length(catVar));
    
    for i in 1:length(catVar)
        dfCat = combine(groupby(dfIndividuals, [Symbol(catVar[i])]), nrow => :Count)
        sort!(dfCat, [Symbol(catVar[i])]) 
        dfCat = DataFrame(vcat([names(dfCat)[1] " "], Matrix(dfCat)), ["Clinical Features", "Count/ mean(SD)"])
        vDFcat[i] = dfCat 
    end
    
    # Continuous Variables
    contVar = ["Annuli", "Age", "WEIGHT (KG)", "LENGTH (CM)", "VET SCORE (Adipose)"]
    vDFcont = Vector{DataFrame}(undef, length(contVar));
    
    for i in 1:length(contVar)
        # calculate mean
        vVar = dfIndividuals[:,Symbol(contVar[i])];
        idxNotMiss = findall(.!ismissing.(vVar))
        vVar = parse.(Float32, string.(vVar[idxNotMiss]))
        myMean = vVar |> mean |> x->round(x, digits = 2)
        # calculate SD
        myStd = vVar |> std |> x->round(x, digits = 2)
    
        dfCont  = DataFrame([contVar[i] string(myMean,"(", myStd, ")")], ["Clinical Features", "Count/ mean(SD)"])
        vDFcont[i] = dfCont 
    end
    vDF = vcat(vDFcat, vDFcont)
    
    dfDem = reduce(vcat, vDF)
    
    return dfDem
end