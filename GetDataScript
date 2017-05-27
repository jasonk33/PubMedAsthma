using BioMedQuery.UMLS
using BioMedQuery.Processes


email= "jason_katz@brown.edu"
host="Jasons-MacBook-Air-2.local"
mysql_usr="jasonk33"
credentials = Credentials("jasonk33", ENV["UMLS_PSWD"])
overwrite=true
verbose = false
append = false
max_articles = 25000

search_term_adult = "asthma[mh] AND adult[mh] NOT (infant[mh] OR child[mh] OR adolescent[mh]) AND (1800[Date - Publication] : 2/13/2014[Date - Publication])"
dbname_adult="pubmed_asthma_adult"
config = Dict(:host=>host,:dbname=>dbname_adult,:username=>mysql_usr,:pswd=>ENV["MYSQL_PSWD"],:overwrite=>overwrite)

@time begin
    db_adult = pubmed_search_and_save(email, search_term_adult, max_articles, save_efetch_mysql, config, verbose)
end
@time begin
    map_mesh_to_umls_async!(db_adult, credentials; append_results=append)
end

search_term_pediatric = "asthma[mh] AND (infant[mh] OR child[mh] OR adolescent[mh]) NOT adult[mh] AND (1800[Date - Publication] : 2/13/2014[Date - Publication])" #search term for pediatric asthma
dbname_pediatric="pubmed_asthma_pediatric"
config = Dict(:host=>host,:dbname=>dbname_pediatric,:username=>mysql_usr,:pswd=>ENV["MYSQL_PSWD"],:overwrite=>overwrite)

@time begin
    db_ped = pubmed_search_and_save(email, search_term_pediatric, max_articles, save_efetch_mysql, config, verbose) #get info from articles
end
@time begin
    map_mesh_to_umls_async!(db_ped, credentials; append_results=append)
end

mysql_disconnect(db_ped)
mysql_disconnect(db_adult)
