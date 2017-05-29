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
    #Uses the arules package in R to call the apriori algorithm to generate association rules for data
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
    #Uses the arulesViz package in R to create group and graph based visualization of assocation rules
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
    #utility function that takes in a mysql database and return a dataframes of mesh descriptors before and after filtering, and umls semantic types, as well as repective counts and frequencies
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
    #Changes the sparse occurance matrix into a dataframe ready to be used in the apriori call
    name_dict = sort(collect(des_ind_dict), by=x->x[2]) #sort the indices from the dictionary
    col_names = DataArray(String, length(des_ind_dict)) #initialize empty vector to be used for column names
        name = name_dict[i][1]
        col_names[i] = name
    end #match the indices with the column names
    itemsets = DataFrame(Matrix(convert(Array{Int64}, disease_occurances')))    #convert the sparse occurance matrix into a dataframe
    names!(itemsets, [symbol(col_names[i]) for i in 1:length(col_names)])   #add the appropriate column names to the dataframe
    return itemsets #return the object as a dataframe
end

function semantic_fold(semantic_ped, semantic_adult, term::String, verbose=true)
    #Calculates the fold difference for specific umls semantic type between two sets of data
    try
        fold = (semantic_ped[semantic_ped[2].==term,:freq]/semantic_adult[semantic_adult[2].==term,:freq])[1]   #calculate fold difference for term
    catch
        error("The term $term is not in the semantic data sets")    #check if term entered is in both sets of semantic terms
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
    end #check to see which type of asthma the term is more associated with
    return round(fold,2)    #round final fold difference
end

function mesh_fold(mesh_ped, mesh_adult)
    #Gets dataframe of fold differences between semantic terms in two data sets (filtering applied within function)
    mesh_set = join(mesh_ped,mesh_adult,on=:mesh_descriptor)    #join both mesh descriptor dataframes
    mesh_set_filt = mesh_set[mesh_set[:freq] + mesh_set[:freq_1] .>.005,:]  #filter for only terms where total frequency is above 0.5%

    mesh_set_filt[:fold_ped]=Array(Float64, length(mesh_set_filt[1]))   #add column to store fold differences for pediatric terms
    for i in mesh_set_filt[2]
        mesh_set_filt[mesh_set_filt[2].==i,:fold_ped] = (mesh_set_filt[mesh_set_filt[2].==i,:freq]/mesh_set_filt[mesh_set_filt[2].==i,:freq_1])[1]
    end #get fold differences

    mesh_set_filt[:fold_adult]=Array(Float64, length(mesh_set_filt[1]))   #add column to store fold differences for adult terms
    for i in mesh_set_filt[2]
        mesh_set_filt[mesh_set_filt[2].==i,:fold_adult] = (mesh_set_filt[mesh_set_filt[2].==i,:freq_1]/mesh_set_filt[mesh_set_filt[2].==i,:freq])[1]
    end #get fold differences

    folds_mesh = sort(mesh_set_filt[(mesh_set_filt[:fold_ped].>5) | (mesh_set_filt[:fold_adult].>5),[2,6,7]],cols=2,rev=true)   #filter for fold differences above 5

    return folds_mesh[folds_mesh[2].>1,1:2],sort(folds_mesh[folds_mesh[3].>1,[1,3]],cols=2,rev=true)    #seperate fold differences for pediatric and adult terms
end

function filter_mesh_by_concepts(db, umls_concepts...)
    #Updated on bcbi BioMedQuery github
    if length(umls_concepts) == 1
        uc = string("'", replace(umls_concepts, "'", "''") , "'")
        query  = mysql_execute(db, "SELECT mesh FROM mesh2umls
        WHERE umls LIKE $uc ")
    else
        query_1 = string(" '", umls_concepts[1], "'")
        queries = DataArray(String, length(umls_concepts)-1)
        for i in 2:length(umls_concepts)
            query = string(", '", umls_concepts[i], "'")
            queries[i-1] = query
        end
        query_2 = join(queries)
        query_joined = string("SELECT mesh FROM mesh2umls WHERE umls IN (", query_1, query_2, " )")
        query  = mysql_execute(db, query_joined)
    end
    #return data array
    return get_value(query.columns[1])
end

function lhs_rhs_vals(association_rules)
    #Grouped Matrix plot data in Julia (outdated)
    lhs_counts = collect(zip(values(countmap(association_rules[1])),keys(countmap(association_rules[1]))))
    lhs_before=DataFrame(term=[])
    for i in 1:length(lhs_counts)
        left = [lhs_counts[i][2]]
        push!(lhs_before, left)
    end

    rhs_counts = collect(zip(values(countmap(association_rules[2])),keys(countmap(association_rules[2]))))
    rhs_before=DataFrame(term=[])
    for i in 1:length(rhs_counts)
        right = [rhs_counts[i][2]]
        push!(rhs_before, right)
    end

    lhs = DataArray(String,length(lhs_counts))
    for i in 1:length(lhs_counts)
        name = string(lhs_counts[i][2], " : ", lhs_counts[i][1])
        lhs[i]=name
    end

    rhs = DataArray(String,length(rhs_counts))
    for i in 1:length(rhs_counts)
        name = string(rhs_counts[i][2], " : ", rhs_counts[i][1])
        rhs[i] = name
    end

    x_vals = DataArray(String, length(association_rules[1]))
    for i in 1:length(association_rules[1])
        x_val = lhs[association_rules[i,1].==lhs_before[1],1][1]
        x_vals[i] = x_val
    end

    y_vals = DataArray(String, length(association_rules[1]))
    for i in 1:length(association_rules[1])
        y_val = rhs[association_rules[i,2].==rhs_before[1],1][1]
        y_vals[i] = y_val
    end

    return x_vals, y_vals
end

function change_bubbles(x_vals, y_vals)
    #Change where bubble lie in grouped matrix visualization of association rules
    dat2 = DataFrame(x=x_vals, y=y_vals)
    counts = []
    for i in 1:length(dat2[1])
        yes = isequal(mode(x_vals), dat2[i,1]) | isequal(mode(y_vals), dat2[i,2])
        if yes == false
            count = i
            push!(counts, count)
        end
    end
    cc = convert(DataArray{Int}, counts)
    datt = dat2[cc,:]
    counts = []
    for i in 1:length(dat2[1])
        yes = isequal(mode(x_vals), dat2[i,1]) | isequal(mode(y_vals), dat2[i,2])
        if yes == true
            count = i
            push!(counts, count)
        end
    end
    cc = convert(DataArray{Int}, counts)
    dattt = append!(datt,dat2[cc,:])
    return dattt[1],dattt[2]
end
