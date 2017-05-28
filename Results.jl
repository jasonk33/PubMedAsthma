using BioMedQuery
using BioMedQuery.UMLS
using BioMedQuery.Processes
using BioMedQuery.Entrez
using BioMedQuery.Entrez.DB
using MySQL
using DataFrames
using DataTables
using StatsBase
using PlotlyJS
using RCall
using AssociationRules
using Rsvg

db_ped = mysql_connect("Jasons-MacBook-Air-2.local", "jasonk33", ENV["MYSQL_PSWD"], "pubmed_asthma_pediatric")
db_adult = mysql_connect("Jasons-MacBook-Air-2.local", "jasonk33", ENV["MYSQL_PSWD"], "pubmed_asthma_adult")

mesh_descrips_ped,semantics_ped,mesh_descrips_filtered_ped=get_mesh_semantics_filtered(db_ped);
mesh_descrips_adult,semantics_adult,mesh_descrips_filtered_adult=get_mesh_semantics_filtered(db_adult);
des_ind_dict_ped, disease_occurances_ped = umls_semantic_occurrences(db_ped, "Disease or Syndrome", "Mental or Behavioral Dysfunction", "Neoplastic Process");
des_ind_dict_adult, disease_occurances_adult = umls_semantic_occurrences(db_adult, "Disease or Syndrome", "Mental or Behavioral Dysfunction", "Neoplastic Process");
itemsets_ped = occurances_to_itemsets(des_ind_dict_ped, disease_occurances_ped)
itemsets_adult = occurances_to_itemsets(des_ind_dict_adult, disease_occurances_adult)
association_rules_ped = apriori2(itemsets_ped, .01, .01, 2, 3, 0)
association_rules_adult = apriori2(itemsets_adult, .01, .01, 2, 3, 0)
plot(bar(x=mesh_descrips_ped[1:25,:mesh_descriptor], y=mesh_descrips_ped[1:25,:freq]))
plot(bar(x=mesh_descrips_adult[1:25,:mesh_descriptor], y=mesh_descrips_adult[1:25,:freq]))
plot(bar(x=semantics_ped[1:10,:semantic_type], y=semantics_ped[1:10,:freq]))
plot(bar(x=semantics_adult[1:10,:semantic_type], y=semantics_adult[1:10,:freq]))
plot(bar(x=mesh_descrips_filtered_ped[1:25,:mesh_descriptor], y=mesh_descrips_filtered_ped[1:25,:freq]))
plot(bar(x=mesh_descrips_filtered_adult[1:25,:mesh_descriptor], y=mesh_descrips_filtered_adult[1:25,:freq]))


folds=[]
for i in join(semantics_ped, semantics_adult, on = :semantic_type)[2]
    fold = semantic_fold(semantics_ped, semantics_adult,i,false)
    push!(folds,fold)
end
round(mean(folds.<2),2)

isequal(sort(semantics_ped[1:5,2]),sort(semantics_adult[1:5,2]))
semantics_ped[1:5,2]


for i in ["Family Group","Environmental Effect of Humans","Conceptual Entity","Organization","Regulation or Law"]
    println(semantic_fold(semantics_ped, semantics_adult,i)," is the fold difference")
end

for i in ["Anatomical Abnormality","Organophosphorus Compound","Body Substance","Cell Function","Neoplastic Process"]
    try println(semantic_fold(semantics_ped, semantics_adult,i)," is the fold difference")
    catch
        println("The term $i is not in the semantic data sets")
    end
end

mesh_descrips_filtered_ped[mesh_descrips_filtered_ped[:freq] .> .02,:]
mesh_descrips_filtered_adult[mesh_descrips_filtered_adult[:freq] .> .02,:]


mesh_fold(mesh_descrips_filtered_ped,mesh_descrips_filtered_adult)



arules_viz(itemsets_ped)
arules_viz(itemsets_adult)

arules_viz(itemsets_ped, "graph")
arules_viz(itemsets_adult, "graph")


mysql_disconnect(db_ped)
mysql_disconnect(db_adult)
