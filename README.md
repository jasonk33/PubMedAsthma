# Mining and Visualizing MeSH-based Associations in PubMed

The goal of this project was to replicate the results obtained in the PubMedMiner article, which can be found here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4419975/

Two searches were conducted using the PubMed database, one for articles pertaining to pediatric asthma, and a second for adult asthma. Data relating to mesh terms used in the articles, as well as the UMLS semantic types associated with each term were extracted and stored in mysql databases. The mesh term were then filtered by the semantic types: (1) “Disease or Syndrome,” (2) “Mental or Behavioral Dysfunction,” and (3) “Neoplastic Process,”. Comparisons were then made between mesh terms and semantic types found in the different sets of articles. Next, association rules for the filtered mesh terms were generated using the apriori algorithm. Both grouped-matrix and graph based visualization were produced to further analyze the association rules. 

Here are the observed differences between the results in the article, and our results:

  Pediatric Asthma Search:
  
    Ours: 24,380 articles, 7,385 unique MeSH descriptors, 117 semantic types, 1,016 MeSH descriptors after filtering, 32 association rules
    
    Theirs: 23,448 articles, 6,524 unique MeSH descriptors, 123 semantic types, 669 MeSH descriptors after filtering, 35 association rules
    
  Adult Asthma Search:
  
    Ours: 23,513 articles, 8,892 unique MeSH descriptors, 121 semantic types, 1,213 MeSH descriptors after filtering, 43 association rules
    
    Theirs: 22,839 articles, 8,688 unique MeSH descriptors, 127 semantic types, 948 MeSH descriptors after filtering, 34 association rules
