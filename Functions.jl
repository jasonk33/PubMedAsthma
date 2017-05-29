using BioMedQuery
using BioMedQuery.UMLS
using BioMedQuery.Processes
using BioMedQuery.Entrez
using BioMedQuery.Entrez.DB
using MySQL
using DataFrames
using StatsBase
using PlotlyJS
using RCall
using AssociationRules



function split_rule!(dat)
    n = size(dat, 1)
    dat[:lhs] = Array{String,1}(n)
    dat[:rhs] = Array{String,1}(n)
    for i = 1:n
        dat[i, :lhs], dat[i, :rhs] = split(dat[i, :rules], " => ")
    end
end

function apriori2(dat::DataFrame, supp = 0.2, conf = 0.01, minlen = 1, maxlen = 10, minlift = 1.2)
    @rput supp  # add user input to r environment
    @rput conf  # add user input to r environment
    @rput minlen  # add user input to r environment
    @rput maxlen  # add user input to r environment
    @rput minlift  # add user input to r environment
    @rput dat  # add user input to r environment
    R"library('arules')"    # load r library to create association rules
    R"dat2 <- as(as(dat, 'matrix'),'itemMatrix')"   # get data into proper format for apriori call
    R"rules1 <- apriori(dat2, parameter = list(supp = supp, conf = conf, minlen = minlen, maxlen = maxlen), control = list(verbose = FALSE))"   # use apriori algorithm to get association rules for data
    R"rules1 <- if (length(rules1) == 0) data.frame() else rules1"  # check if any rules were created
    R"rules1 <- character_columns(as(rules1, \"data.frame\"))"  # formatting
    R"rules_sub <- subset(rules1, subset = lift > minlift)" # filter by minimum lift
    rules_df = @rget rules_sub; # get dataframe from R
    R"rm(dat, dat2, rules1, rules_sub, supp, conf, minlen, maxlen, minlift)"    # clean up R environment
    split_rule!(rules_df);  # formatting
    rules_df = rules_df[:, [:lhs, :rhs, :support, :confidence, :lift]]  # create columns
    for i in 1:length(rules_df[1])
        lhs = rules_df[i,1][2:end-1]
        rhs = rules_df[i,2][2:end-1]
        rules_df[i,1] = lhs
        rules_df[i,2] = rhs
    end # remove brakcets around lhs and rhs strings in dataframe
    rules_df[:chi_squared] = length(dat[1]).*rules_df[:,:support].*(rules_df[:,:lift]-1).^2.*rules_df[:,:support].*rules_df[:,:confidence]./(rules_df[:,:confidence]-rules_df[:,:support])./(rules_df[:,:lift]-rules_df[:,:confidence]) # add column with chi squared statistic for each rule
    sort!(rules_df, cols = :chi_squared, rev=true)  # sort rules in descending order of chi squared statistic
    rules_df # return the association rules datarame
end

function arules_viz(itemset, method="grouped", num_rules=50)
    @rput itemset  # add user input to r environment
    @rput num_rules  # add user input to r environment
    @rput method  # add user input to r environment
    R"library('arules')"    # load r library to create association rules
    R"library('arulesViz')" # load r library to create association rules visualizations
    R"data <- as(as(itemset, 'matrix'),'itemMatrix')"   # get data into proper format
    R"rules <- apriori(data, parameter = list(supp = .01, conf = .01, minlen = 2, maxlen = 3), control = list(verbose = FALSE))"  # use apriori algorithm to get association rules for data
    plot = R"capture.output(plot(rules, method=method, control=list(k=num_rules)))" # create plot (either grouped or graph) for association rules
    R"rm(itemset, data, rules, num_rules, method)"   # clean up R environment
    return plot # return the visualization
end

function get_mesh_semantics_filtered(db)
    mesh_descriptor = mysql_execute(db, "SELECT * FROM mesh_descriptor;") #get mesh descriptor data from MySQL
    mesh_heading = mysql_execute(db, "SELECT * FROM mesh_heading;") #get header data from MySQL
    mesh2umls = mysql_execute(db, "SELECT * FROM mesh2umls;") #get data umls data from MySQL
    rename!(mesh_descriptor, [:id, :name], [:did, :mesh_descriptor]) #change columns name for join
    data = sort(join(mesh_descriptor, mesh_heading, on = :did)[:,[:mesh_descriptor, :pmid]], cols=:pmid) #data of mesh terms for each article
    arts = unique(data[:pmid]) #article ID's
    rename!(mesh2umls, :mesh, :mesh_descriptor) #change column name for join
    data_semantic = join(data, mesh2umls, on = :mesh_descriptor) #data of semantic types for each article
    semantic_counts = sort(collect(zip(values(countmap(data_semantic[:umls])),keys(countmap(data_semantic[:umls])))),rev=true) #counts for semantic types
    semantics=DataFrame(Any,0,2)
    for i in 1:length(semantic_counts)
        semantic = [semantic_counts[i][1],semantic_counts[i][2]]
        push!(semantics, semantic)
    end #get semantic counts into usable form
    frequency = DataArray(Float64, length(semantics[2]))
    for i in 1:length(semantics[2])
        freq = length(unique(data_semantic[data_semantic[:umls] .== semantics[i,2],2]))/length(arts)
        frequency[i] = freq
    end #calculate frequency for each semantic type
    semantics[:freq]=frequency #add frequencies to data
    sort!(semantics, cols = :freq, rev = true) #sort by frequency
    rename!(semantics, [:x1, :x2], [:count, :semantic_type]) #rename columns
    mesh_counts = sort(collect(zip(values(countmap(data[:mesh_descriptor])),keys(countmap(data[:mesh_descriptor])))),rev=true) #counts for mesh descriptors
    mesh_descrips=DataFrame(Any,0,2)
    for i in 1:length(mesh_counts)
        mesh_descrip = [mesh_counts[i][1],mesh_counts[i][2]]
        push!(mesh_descrips, mesh_descrip)
    end #get mesh counts into usable form
    frequency = DataArray(Float64, length(mesh_descrips[2]))
    for i in 1:length(mesh_descrips[2])
        freq = length(unique(data[data[:mesh_descriptor] .== mesh_descrips[i,2],2]))/length(arts)
        frequency[i] = freq
    end #calculate frequency for each mesh term
    mesh_descrips[:freq]=frequency #add frequencies to data
    sort!(mesh_descrips, cols = :freq, rev = true) #sort by frequency
    rename!(mesh_descrips, [:x1, :x2], [:count, :mesh_descriptor]) #rename columns
    umls_filtered = mesh2umls[(mesh2umls[:umls] .== "Disease or Syndrome") | (mesh2umls[:umls] .== "Mental or Behavioral Dysfunction") | (mesh2umls[:umls] .== "Neoplastic Process"),:] #filter by semantic type
    data_filtered = join(data, umls_filtered, on = :mesh_descriptor) #new data after filtering
    arts_filtered = unique(data_filtered[:pmid]) #article ID's after filtering
    mesh_counts_filtered = sort(collect(zip(values(countmap(data_filtered[:mesh_descriptor])),keys(countmap(data_filtered[:mesh_descriptor])))),rev=true) #counts for mesh descriptors after filtering
    mesh_descrips_filtered=DataFrame(Any,0,2)
    for i in 1:length(mesh_counts_filtered)
        mesh_descrip_filtered = [mesh_counts_filtered[i][1],mesh_counts_filtered[i][2]]
        push!(mesh_descrips_filtered, mesh_descrip_filtered)
    end #get mesh counts into usable form
    frequency = DataArray(Float64, length(mesh_descrips_filtered[2]))
    for i in 1:length(mesh_descrips_filtered[2])
        freq = length(unique(data_filtered[data_filtered[:mesh_descriptor] .== mesh_descrips_filtered[i,2],2]))/length(arts_filtered)
        frequency[i] = freq
    end #calculate frequency for each mesh term
    mesh_descrips_filtered[:freq]=frequency #add frequencies to data
    sort!(mesh_descrips_filtered, cols = :freq, rev = true) #sort by frequency
    rename!(mesh_descrips_filtered, [:x1, :x2], [:count, :mesh_descriptor]) #rename columns
    return mesh_descrips,semantics, mesh_descrips_filtered
end

function occurances_to_itemsets(des_ind_dict, disease_occurances)
    name_dict = sort(collect(des_ind_dict), by=x->x[2])
    col_names = DataArray(String, length(des_ind_dict))
    for i in 1:length(des_ind_dict)
        name = name_dict[i][1]
        col_names[i] = name
    end
    itemsets = DataFrame(Matrix(convert(Array{Int64}, disease_occurances')))
    names!(itemsets, [symbol(col_names[i]) for i in 1:length(col_names)])
    return itemsets
end

function semantic_fold(semantic_ped, semantic_adult, term::String, verbose=true)
    try
        fold = (semantic_ped[semantic_ped[2].==term,:freq]/semantic_adult[semantic_adult[2].==term,:freq])[1]
    catch
        error("The term $term is not in the semantic data sets")
    end
    if fold < 1
        fold = inv(fold)
        if verbose == true
            println("The term $term is associated more with adult asthma")
        end
    else
        if verbose == true
            println("The term $term is associated more with pediatric asthma")
        end
    end
    return round(fold,2)
end

function mesh_fold(mesh_ped, mesh_adult)

    mesh_set = join(mesh_ped,mesh_adult,on=:mesh_descriptor)
    mesh_set_filt = mesh_set[mesh_set[:freq] + mesh_set[:freq_1] .>.005,:]

    mesh_set_filt[:fold_ped]=Array(Float64, length(mesh_set_filt[1]))
    for i in mesh_set_filt[2]
        j = mesh_set_filt[mesh_set_filt[2].==i,1]
        mesh_set_filt[mesh_set_filt[2].==i,:fold_ped] = (mesh_set_filt[mesh_set_filt[2].==i,:freq]/mesh_set_filt[mesh_set_filt[2].==i,:freq_1])[1]
    end

    mesh_set_filt[:fold_adult]=Array(Float64, length(mesh_set_filt[1]))
    for i in mesh_set_filt[2]
        j = mesh_set_filt[mesh_set_filt[2].==i,1]
        mesh_set_filt[mesh_set_filt[2].==i,:fold_adult] = (mesh_set_filt[mesh_set_filt[2].==i,:freq_1]/mesh_set_filt[mesh_set_filt[2].==i,:freq])[1]
    end

    folds_mesh = sort(mesh_set_filt[(mesh_set_filt[:fold_ped].>5) | (mesh_set_filt[:fold_adult].>5),[2,6,7]],cols=2,rev=true)

    return folds_mesh[folds_mesh[2].>1,1:2],sort(folds_mesh[folds_mesh[3].>1,[1,3]],cols=2,rev=true)
end
