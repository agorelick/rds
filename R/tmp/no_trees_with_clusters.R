no_trees_with_clusters <-
function(k, clusters){
    
#    Calculates the number of trees where m1, m2, ..., mC, metastases form a cluster (have one root)
#    k: number of remaining cancer samples of the same phylogeny
#    clusters: list of number of metastases samples per C clusters
#    returns: number of trees where the given clusters have independent origin
    
    n = k + sum(clusters)
    no_trees = 1
    for (c in clusters){
        
        # this is an addition to the python script to allow for RDS caculcations for over 450 metastases samples     
        c <- gmp::as.bigz(c)
        
        no_trees <-  no_trees * (factorial(2 * c - 3) / (2^(c-2) * factorial(c - 2)))
        
    }
    no_trees <-  no_trees * (factorial((2*k) + (2 * length(clusters) - 3)) 
                / (2^(k + (length(clusters) - 2)) * factorial(k + (length(clusters) - 2))))  
    
    no_trees    
}
