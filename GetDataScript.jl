using BioMedQuery.UMLS
using BioMedQuery.Processes


email= "jason_katz@brown.edu"   #account for querying UMLS
host="Jasons-MacBook-Air-2.local"   #host for mysql
mysql_usr="jasonk33"    #mysql username
credentials = Credentials("jasonk33", ENV["UMLS_PSWD"]) #credentials for mapping umls semantic types to mesh terms
verbose = false #function input
max_articles = 25000    #max articles to be returned in search

search_term_adult = "asthma[mh] AND adult[mh] NOT (infant[mh] OR child[mh] OR adolescent[mh]) AND (1800[Date - Publication] : 2/13/2014[Date - Publication])"   #search term for adult asthma
dbname_adult="pubmed_asthma_adult"  #database name to store search results
config = Dict(:host=>host,:dbname=>dbname_adult,:username=>mysql_usr,:pswd=>ENV["MYSQL_PSWD"],:overwrite=>true) #formatting

@time begin
    db_adult = pubmed_search_and_save(email, search_term_adult, max_articles, save_efetch_mysql, config, verbose)
end #search for adult asthma articles and save in mysql database
@time begin
    map_mesh_to_umls_async!(db_adult, credentials; append_results=false)
end #map umls semantic types to mesh descriptors

search_term_pediatric = "asthma[mh] AND (infant[mh] OR child[mh] OR adolescent[mh]) NOT adult[mh] AND (1800[Date - Publication] : 2/13/2014[Date - Publication])" #search term for pediatric asthma
dbname_pediatric="pubmed_asthma_pediatric"  #database name to store search results
config = Dict(:host=>host,:dbname=>dbname_pediatric,:username=>mysql_usr,:pswd=>ENV["MYSQL_PSWD"],:overwrite=>overwrite)    #formatting

@time begin
    db_ped = pubmed_search_and_save(email, search_term_pediatric, max_articles, save_efetch_mysql, config, verbose) #get info from articles
end #search for pediatric asthma articles and save in mysql database
@time begin
    map_mesh_to_umls_async!(db_ped, credentials; append_results=append)
end #map umls semantic types to mesh descriptors

mysql_disconnect(db_ped)    #disconnect from mysql
mysql_disconnect(db_adult)    #disconnect from mysql
