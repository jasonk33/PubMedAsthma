using BioMedQuery
pass =


credentials = Credentials("jasonk33", pass)
tgt = get_tgt(credentials)
term = "obesity"
query = Dict("string"=>term, "searchType"=>"exact" )
all_results= search_umls(tgt, query)


email= "jason_katz@brown.edu"
search_term="(obesity[MeSH Major Topic]) AND (\"2010\"[Date - Publication] : \"2012\"[Date - Publication])"
max_articles = 20
overwrite=true
verbose = false

results_dir = "./results"
if !isdir(results_dir)
    mkdir(results_dir)
end

    host="Jasons-MacBook-Air-2.local" #If want to hide - use enviroment variables instead
    mysql_usr="jasonk33"
    #mysql_pswd=""
    dbname="pubmed_obesity_2010_2012"
    config = Dict(:host=>host,
                     :dbname=>dbname,
                     :username=>mysql_usr,
                     :pswd=>mysql_pswd,
                     :overwrite=>overwrite)
    save_func = save_efetch_mysql

@time begin
    db = pubmed_search_and_save(email, search_term, max_articles,
    save_func, config, verbose)
end

using MySQL
con = mysql_connect(host, mysql_usr, mysql_pswd, dbname)
using DataFrames
using DataStructures
mysql_execute(con, "SHOW TABLES;")

article = mysql_execute(con, "SELECT * FROM article;")

author = mysql_execute(con, "SELECT * FROM author;")

author2article = mysql_execute(con, "SELECT * FROM author2article;")

mesh_descriptor = mysql_execute(con, "SELECT * FROM mesh_descriptor;")

mesh_heading = mysql_execute(con, "SELECT * FROM mesh_heading;")

mesh_qualifier = mysql_execute(con, "SELECT * FROM mesh_qualifier;")

rename!(mesh_descriptor, :id, :did)
rename!(mesh_descriptor, :name, :mesh_descriptor)

data = sort(join(mesh_descriptor, mesh_heading, on = :did)[:,[:mesh_descriptor, :pmid]], cols=:pmid)

occur = sort(collect(zip(values(countmap(data[:mesh_descriptor])),keys(countmap(data[:mesh_descriptor])))),rev=true)

using BioMedQuery.Processes
using BioMedQuery.UMLS
append = false
@time begin
    map_mesh_to_umls_async!(db, credentials; append_results=append)
end

using BioMedQuery.DBUtils

db_query(con, "SELECT * FROM author;")


mysql_execute(con, "SHOW TABLES;")




mesh2umls = mysql_execute(con, "SELECT * FROM mesh2umls;")





max_descrips = sum(mesh_heading[:pmid] .== mode(mesh_heading[:pmid]))
test = DataFrame(AbstractString,0,max_descrips)
for i in 1:length(arts)
    transaction = convert(Array{AbstractString}, unique(data[1][data[:pmid].==arts[i],1]))
    dat = DataArray(AbstractString, max_descrips)
    dat[1:length(transaction)] = transaction
    push!(test, dat)
end
for i in 1:length(test)
    test[ isna(test[:,i]), i] = " "
end


apriori(test[1:5,1:16], .4, .5)




groceries = ["milk", "bread", "eggs", "apples", "oranges", "beer"]
transacts = [sample(groceries, 4, replace = false) for x in 1:1000]
transacts_df = DataFrame(transactions(transacts))
apriori(transacts_df[1:10,:], .4, .5)





######################################functions###########################

using BioMedQuery.Processes
using BioMedQuery.Entrez
using MySQL
using DataFrames
using StatsBase
using RCall

email= "jason_katz@brown.edu"
host="Jasons-MacBook-Air-2.local"
mysql_usr="jasonk33"

function search_and_related(MeSH_Major_Topic, MeSH_Term_1="none", MeSH_Term_2="none" ; num_articles=20, year_start=2000, year_end=2017)
    #using BioMedQuery.Processes
    #using BioMedQuery.Entrez
    #using MySQL
    #using DataFrames
    #using StatsBase
    if year_start > year_end
      error("Starting year must be before ending year")
    end
    if (MeSH_Term_1 == "none") & (MeSH_Term_2 != "none")
      error("Enter MeSH Term 1 before MeSH Term 2")
    end
    if MeSH_Term_1 == "none"
      search_term="($MeSH_Major_Topic[MeSH Major Topic]) AND (\"$year_start\"[Date - Publication] : \"$year_end\"[Date - Publication])"
    elseif MeSH_Term_2 == "none"
      search_term="(($MeSH_Major_Topic[MeSH Major Topic]) AND $MeSH_Term_1[MeSH Terms]) AND (\"$year_start\"[Date - Publication] : \"$year_end\"[Date - Publication])"
    else
      search_term="((($MeSH_Major_Topic[MeSH Major Topic]) AND $MeSH_Term_1[MeSH Terms]) AND $MeSH_Term_2[MeSH Terms]) AND (\"$year_start\"[Date - Publication] : \"$year_end\"[Date - Publication])"
    end
    overwrite=true
    verbose = false
    global email #Configure Locally Before
    global host #Configure Locally Before
    global mysql_usr #Configure Locally Before
    mysql_pswd = ENV["MYSQL_PSWD"] #Configure Locally Before
    dbname="pubmed_$MeSH_Major_Topic"
    config = Dict(:host=>host,
                     :dbname=>dbname,
                     :username=>mysql_usr,
                     :pswd=>mysql_pswd,
                     :overwrite=>overwrite)
    db=nothing
    try
        @time begin
        db = pubmed_search_and_save(email, search_term, num_articles, save_efetch_mysql, config, verbose)
        end
    catch
        error("Not enough articles matching criteria, either reduce the number of articles, or change some of the criteria")
    end
    #con = mysql_connect(host, mysql_usr, mysql_pswd, dbname)
    mesh_descriptor = mysql_execute(db, "SELECT * FROM mesh_descriptor;")
    mesh_heading = mysql_execute(db, "SELECT * FROM mesh_heading;")
    rename!(mesh_descriptor, :id, :did)
    rename!(mesh_descriptor, :name, :mesh_descriptor)
    data = sort(join(mesh_descriptor, mesh_heading, on = :did)[:,[:mesh_descriptor, :pmid]], cols=:pmid)
    occur = sort(collect(zip(values(countmap(data[:mesh_descriptor])),keys(countmap(data[:mesh_descriptor])))),rev=true)
    #mysql_disconnect(con)
    return [data, occur]
end

data = search_and_related("alcohol", num_articles=100)

function arules(data ; support_min=.01, confidence_min=.6, lift_min=1.2, right_side = "none", left_side = "none", left_side_2 = "none")
    #Get data using search_and_related function
    #using DataFrames
    #using RCall (requires R to be installed)
    if (left_side == "none") & (left_side_2 != "none")
      error("Must input first left side term before inputing a second left side term")
    end
    arts = unique(data[1][:pmid])
    transaction = convert(Array{AbstractString}, unique(data[1][data[1][:pmid].==arts[1],1]))
    descriptors = convert(Array{AbstractString}, unique(data[1][:mesh_descriptor]))
    itemsets = DataFrame(Any, 0,length(descriptors))
    names!(itemsets, [symbol(descriptors[i]) for i in 1:length(descriptors)])
    for j = 1:length(arts)
        contains=Vector(length(descriptors))
    for i = 1:length(descriptors)
        contain = Int(in(descriptors[i],convert(Array{AbstractString}, unique(data[1][data[1][:pmid].==arts[j],1]))))
        contains[i] = contain
    end
      push!(itemsets, contains)
    end
    writetable("output.csv",itemsets)
    R"library('arules')"
    path = homedir()
    @rput path
    R"data_R = read.csv(paste(path, '/output.csv', sep = ''))"
    R"iMat = as(as(data_R, 'matrix'),'itemMatrix')"
    @rput support_min
    @rput confidence_min
    @rput lift_min
    @rput right_side
    @rput left_side
    @rput left_side_2
    R"rules <- apriori(iMat, parameter = list(support = support_min, confidence = confidence_min))"
    if (right_side == "none") & (left_side == "none") & (left_side_2 == "none")
      R"rules_sub = subset(rules, subset = lift > lift_min)"
    elseif (right_side != "none") & (left_side != "none") & (left_side_2 != "none")
      R"rules_sub = subset(rules, subset = rhs %in% right_side & lhs %in% left_side & lhs %in% left_side_2 & lift > lift_min)"
    elseif (right_side != "none") & (left_side != "none") & (left_side_2 == "none")
      R"rules_sub = subset(rules, subset = rhs %in% right_side & lhs %in% left_side & lift > lift_min)"
    elseif (right_side !== "none") & (left_side == "none") & (left_side_2 == "none")
      R"rules_sub = subset(rules, subset = rhs %in% right_side & lift > lift_min)"
    elseif (right_side == "none") & (left_side !== "none") & (left_side_2 !== "none")
      R"rules_sub = subset(rules, subset = lhs %in% left_side & lhs %in% left_side_2 & lift > lift_min)"
    elseif (right_side == "none") & (left_side !== "none") & (left_side_2 == "none")
      R"rules_sub = subset(rules, subset = lhs %in% left_side & lift > lift_min)"
    else
      error("Something is wrong with function")
    end
    R"num_rules = length(rules_sub)"
    @rget num_rules
    if num_rules == 0
      error("No rules matching criteria, reduce minimum criteria in function call")
    end
    R"inspect(rules_sub)"
end

arules(data, support_min=.1, confidence_min=.25, lift_min=.6, left_side="middle.aged", left_side_2="adult", right_side="male")


###############################################################################



#using AssociationRules
using StatsBase
# for sample() function

# simulate transactions
groceries = ["milk", "bread", "eggs", "apples", "oranges", "beer"]
transactions = [sample(groceries, 4, replace = false) for x in 1:1000]
transactions_df = DataFrame(transactions')







function filter_mesh_by_concepts(db, umls_semantic_type_1, umls_semantic_type_2 = "none", umls_semantic_type_3 = "none", umls_semantic_type_4 = "none", umls_semantic_type_5 = "none", umls_semantic_type_6 = "none", umls_semantic_type_7 = "none", umls_semantic_type_8 = "none", umls_semantic_type_9= "none", umls_semantic_type_10 = "none")

    mesh_descriptor=[]
    mesh_heading=[]
    mesh2umls=[]
    try
        mesh_descriptor = mysql_execute(db, "SELECT * FROM mesh_descriptor;") #get mesh descriptor data from MySQL
    catch
        error("mesh_descriptor table does not exitst in database, use pubmed_search_and_save() function to create your database")
    end
    try
        mesh_heading = mysql_execute(db, "SELECT * FROM mesh_heading;") #get header data from MySQL
    catch
        error("mesh_heading table does not exist in database, use pubmed_search_and_save() function to create your database")
    end
    try
        mesh2umls = mysql_execute(db, "SELECT * FROM mesh2umls;") #get data umls data from MySQL
    catch
        error("You must first add umls terms to database use map_mesh_to_umls_async!() function")
    end
    rename!(mesh_descriptor, [:id, :name], [:did, :mesh_descriptor]) #change columns name for join
    data = sort(join(mesh_descriptor, mesh_heading, on = :did)[:,[:mesh_descriptor, :pmid]], cols=:pmid) #data of mesh terms for each article
    rename!(mesh2umls, :mesh, :mesh_descriptor) #change column name for join
    if umls_semantic_type_2 == "none"
        umls_filtered = mesh2umls[mesh2umls[:umls] .== "$umls_semantic_type_1",:] #filter by semantic type
    elseif umls_semantic_type_3 == "none"
        umls_filtered = mesh2umls[(mesh2umls[:umls] .== "$umls_semantic_type_1") | (mesh2umls[:umls] .== "$umls_semantic_type_2"),:] #filter by semantic type
    elseif umls_semantic_type_4 == "none"
        umls_filtered = mesh2umls[(mesh2umls[:umls] .== "$umls_semantic_type_1") | (mesh2umls[:umls] .== "$umls_semantic_type_2") | (mesh2umls[:umls] .== "$umls_semantic_type_3"),:] #filter by semantic type
    elseif umls_semantic_type_5 == "none"
        umls_filtered = mesh2umls[(mesh2umls[:umls] .== "$umls_semantic_type_1") | (mesh2umls[:umls] .== "$umls_semantic_type_2") | (mesh2umls[:umls] .== "$umls_semantic_type_3") | (mesh2umls[:umls] .== "$umls_semantic_type_4"),:] #filter by semantic type
    else
        umls_filtered = mesh2umls[(mesh2umls[:umls] .== "$umls_semantic_type_1") | (mesh2umls[:umls] .== "$umls_semantic_type_2") | (mesh2umls[:umls] .== "$umls_semantic_type_3") | (mesh2umls[:umls] .== "$umls_semantic_type_4") | (mesh2umls[:umls] .== "$umls_semantic_type_5"),:] #filter by semantic type
    end
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
    return mesh_descrips_filtered

end



itemsets_filtered = DataFrame(Any, 0,length(mesh_descrips_filtered[2])) #initialize empty dataframe
names!(itemsets_filtered, [symbol(mesh_descrips_filtered[i,2]) for i in 1:length(mesh_descrips_filtered[2])]) #add column names
@time begin
    for j = 1:length(arts_filtered)
        contains_filtered=Vector(length(mesh_descrips_filtered[2]))
        for i = 1:length(mesh_descrips_filtered[2])
            contain_filtered = Int(in(mesh_descrips_filtered[i,2],convert(Array{String}, unique(data_filtered[data_filtered[:pmid].==arts_filtered[j],1])))) #get 1 or 0 for each item in a transaction
            contains_filtered[i] = contain_filtered #add each item's value (1 or 0) to the transaction vector
        end
        push!(itemsets_filtered, contains_filtered) #add all transaction vectors to dataframe
    end
end


x_vals,y_vals=lhs_rhs_vals(association_rules);

plot(scatter(x=x_vals, y=y_vals, mode="markers", marker_size=500.*association_rules[:support], marker_color=association_rules[:chi_squared]),Layout(height=600, width=600,margin=[1,1,1,1],title="Grouped Matrix for 43 Rules - Adult Asthma",yaxis_title="RHS",xaxis_title="LHS"))



plot(scatter(x=association_rules_ped[1], y=association_rules_ped[2], mode="markers", marker_size=750.*association_rules_ped[:support], marker_color=association_rules_ped[:chi_squared]),Layout(height=600, width=600,margin=[1,1,1,1],title="Grouped Matrix for $(length(association_rules_ped[1])) Rules - Pediatric Asthma",yaxis_title="RHS",xaxis_title="LHS"))
plot(scatter(x=association_rules_adult[1], y=association_rules_adult[2], mode="markers", marker_size=500.*association_rules_adult[:support], marker_color=association_rules_adult[:chi_squared]),Layout(height=600, width=600,margin=[1,1,1,1],title="Grouped Matrix for $(length(association_rules_adult[1])) Rules - Adult Asthma",yaxis_title="RHS",xaxis_title="LHS"))




########adultasthmascript##########
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
using Rsvg


search_term = "asthma[mh] AND adult[mh] NOT (infant[mh] OR child[mh] OR adolescent[mh]) AND (1800[Date - Publication] : 2/13/2014[Date - Publication])" #search term for adult asthma
email= "jason_katz@brown.edu"
host="Jasons-MacBook-Air-2.local"
mysql_usr="jasonk33"
credentials = Credentials("jasonk33", ENV["UMLS_PSWD"])
overwrite=true
verbose = false
append = false
max_articles = 25000
dbname="pubmed_asthma_adult"
config = Dict(:host=>host,:dbname=>dbname,:username=>mysql_usr,:pswd=>ENV["MYSQL_PSWD"],:overwrite=>overwrite)

@time begin
    pubmed_search_and_save(email, search_term, max_articles, save_efetch_mysql, config, verbose) #get info from articles
end
@time begin
    map_mesh_to_umls_async!(db, credentials; append_results=append) #get semantic types for articles
end

db = mysql_connect(host, mysql_usr, ENV["MYSQL_PSWD"], dbname)

mesh_descrips,semantics,mesh_descrips_filtered=get_mesh_semantics_filtered(db)

des_ind_dict, disease_occurances = umls_semantic_occurrences_2(db, "Disease or Syndrome", "Mental or Behavioral Dysfunction", "Neoplastic Process")

itemsets = occurances_to_itemsets(des_ind_dict, disease_occurances)

association_rules = apriori2(itemsets, .01, .01, 2, 3, 0)
############################################################Results############################################################
unique(mysql_execute(db, "SELECT * FROM mesh_heading")[:pmid]) #Number of articles matching search criteria
length(mesh_descrips[:,2]) #Number of unique mesh descriptors
length(semantics[:,3]) #Number of unique semantic types
plot(bar(x=mesh_descrips[1:25,:mesh_descriptor], y=mesh_descrips[1:25,:freq])) #Top 25 mesh descriptors
plot(bar(x=semantics[1:10,:semantic_type], y=semantics[1:10,:freq])) #Top 10 semantic types
length(mesh_descrips_filtered[:,:mesh_descriptor]) #Number of unique mesh descriptors after filtering
plot(bar(x=mesh_descrips_filtered[1:25,:mesh_descriptor], y=mesh_descrips_filtered[1:25,:freq])) #Top 25 mesh descriptors after filtering
############################################################Results############################################################




q=plot(scatter(x=association_rules[1], y=association_rules[2], mode="markers", marker_size=500.*association_rules[:support], marker_color=association_rules[:chi_squared]),Layout(height=600, width=600,margin=[1,1,1,1],title="Grouped Matrix for 43 Rules - Adult Asthma",yaxis_title="RHS",xaxis_title="LHS"))







#####################################################################################

semantic_ped = readtable("/Users/JasonKatz/Desktop/Code/PubMedAsthma/Data/Pediatric/Semantics.csv")
semantic_adult = readtable("/Users/JasonKatz/Desktop/Code/PubMedAsthma/Data/Adult/Semantics.csv")
mesh_ped = readtable("/Users/JasonKatz/Desktop/Code/PubMedAsthma/Data/Pediatric/MeshDescripsFiltered.csv")
mesh_adult = readtable("/Users/JasonKatz/Desktop/Code/PubMedAsthma/Data/Adult/MeshDescripsFiltered.csv")


semantic_fold(semantic_ped, semantic_adult, "Age Group")
semantic_fold(semantic_ped, semantic_adult, "Population Group")
semantic_fold(semantic_ped, semantic_adult, "Human")
semantic_fold(semantic_ped, semantic_adult, "Disease or Syndrome")
semantic_fold(semantic_ped, semantic_adult, "Organism Attribute")


semantic_fold(semantic_ped, semantic_adult, "Family Group")
semantic_fold(semantic_ped, semantic_adult, "Environmental Effect of Humans")
semantic_fold(semantic_ped, semantic_adult, "Conceptual Entity")
semantic_fold(semantic_ped, semantic_adult, "Organization")
semantic_fold(semantic_ped, semantic_adult, "Regulation or Law")


semantic_fold(semantic_ped, semantic_adult, "Anatomical Abnormality")
semantic_fold(semantic_ped, semantic_adult, "Organophosphorus Compound")
semantic_fold(semantic_ped, semantic_adult, "Body Substance")
semantic_fold(semantic_ped, semantic_adult, "Cell Function")
semantic_fold(semantic_ped, semantic_adult, "Neoplastic Process")

mesh_fold(mesh_ped,mesh_adult)
#############################




function lhs_rhs_vals(association_rules)
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


function filter_mesh_by_concepts(db, umls_concepts...)
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
