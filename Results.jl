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
using Rsvg
#using DataTables

db_ped = mysql_connect("Jasons-MacBook-Air-2.local", "jasonk33", ENV["MYSQL_PSWD"], "pubmed_asthma_pediatric")  #connect to mysql database for pediatric asthma
db_adult = mysql_connect("Jasons-MacBook-Air-2.local", "jasonk33", ENV["MYSQL_PSWD"], "pubmed_asthma_adult")  #connect to mysql database for pediatric asthma

mesh_descrips_ped,semantics_ped,mesh_descrips_filtered_ped=get_mesh_semantics_filtered(db_ped); #get mesh descriptors before and after filtering, and umls semantic types, as well as repective counts and frequencies
mesh_descrips_adult,semantics_adult,mesh_descrips_filtered_adult=get_mesh_semantics_filtered(db_adult); #get mesh descriptors before and after filtering, and umls semantic types, as well as repective counts and frequencies
des_ind_dict_ped, disease_occurances_ped = umls_semantic_occurrences(db_ped, "Disease or Syndrome", "Mental or Behavioral Dysfunction", "Neoplastic Process");  #get sparse occurance matrix after filtering by umls semantic type
des_ind_dict_adult, disease_occurances_adult = umls_semantic_occurrences(db_adult, "Disease or Syndrome", "Mental or Behavioral Dysfunction", "Neoplastic Process");  #get sparse occurance matrix after filtering by umls semantic type
itemsets_ped = occurances_to_itemsets(des_ind_dict_ped, disease_occurances_ped) #convert sparse occurance matrix into dataframe
itemsets_adult = occurances_to_itemsets(des_ind_dict_adult, disease_occurances_adult) #convert sparse occurance matrix into dataframe
association_rules_ped = apriori2(itemsets_ped, .01, .01, 2, 3, 0)   #get association rules for pediatric asthma filtered mesh descriptors
association_rules_adult = apriori2(itemsets_adult, .01, .01, 2, 3, 0)   #get association rules for adult asthma filtered mesh descriptors
plot(bar(x=mesh_descrips_ped[1:25,:mesh_descriptor], y=mesh_descrips_ped[1:25,:freq]))  #plot top 25 mesh descriptors for pediatric asthma
plot(bar(x=mesh_descrips_adult[1:25,:mesh_descriptor], y=mesh_descrips_adult[1:25,:freq]))  #plot top 25 mesh descriptors for adult asthma
plot(bar(x=semantics_ped[1:10,:semantic_type], y=semantics_ped[1:10,:freq]))    #plot top 10 umls semantic types for pediatric asthma
plot(bar(x=semantics_adult[1:10,:semantic_type], y=semantics_adult[1:10,:freq]))    #plot top 10 umls semantic types for adult asthma
plot(bar(x=mesh_descrips_filtered_ped[1:25,:mesh_descriptor], y=mesh_descrips_filtered_ped[1:25,:freq]))  #plot top 25 mesh descriptors after filtering for pediatric asthma
plot(bar(x=mesh_descrips_filtered_adult[1:25,:mesh_descriptor], y=mesh_descrips_filtered_adult[1:25,:freq]))  #plot top 25 mesh descriptors after filtering for adult asthma


folds=[]
for i in join(semantics_ped, semantics_adult, on = :semantic_type)[2]
    fold = semantic_fold(semantics_ped, semantics_adult,i,false)
    push!(folds,fold)
end
round(mean(folds.<2),2) #check what percentage of semantic terms are within two folds of each other for pediatric versus adult asthma data

isequal(sort(semantics_ped[1:5,2]),sort(semantics_adult[1:5,2]))    #check if the top 5 semantic types for pediatric and adult asthma are the same
semantics_ped[1:5,2]    #top 5 umls semantic types (same for both pediatric and adult asthma)


for i in ["Family Group","Environmental Effect of Humans","Conceptual Entity","Organization","Regulation or Law"]
    println(semantic_fold(semantics_ped, semantics_adult,i)," is the fold difference")
end #get fold differences for semantic types mentioned in the article

for i in ["Anatomical Abnormality","Organophosphorus Compound","Body Substance","Cell Function","Neoplastic Process"]
    try println(semantic_fold(semantics_ped, semantics_adult,i)," is the fold difference")
    catch
        println("The term $i is not in the semantic data sets")
    end
end #get fold differences for semantic types mentioned in the article

mesh_descrips_filtered_ped[mesh_descrips_filtered_ped[:freq] .> .02,:]  #get mesh descriptors (after filtering) that have frequencies above 2% for pediatric asthma
mesh_descrips_filtered_adult[mesh_descrips_filtered_adult[:freq] .> .02,:]  #get mesh descriptors (after filtering) that have frequencies above 2% for adult asthma


mesh_fold(mesh_descrips_filtered_ped,mesh_descrips_filtered_adult)  #get semantic terms with the biggest fold differences

arules_viz(itemsets_ped)    #grouped matrix visualization of association rules for pediatric asthma
arules_viz(itemsets_adult)    #grouped matrix visualization of association rules for adult asthma

arules_viz(itemsets_ped, "graph")    #graph based visualization of association rules for pediatric asthma
arules_viz(itemsets_adult, "graph")    #graph based visualization of association rules for adult asthma


mysql_disconnect(db_ped)    #disconnect from mysql
mysql_disconnect(db_adult)    #disconnect from mysql
