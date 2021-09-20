cluster_probability <-
function(k, m) {
    
#   Calculates the probability that a tree where all m metastases form a cluster (have one root) exists by chance
#   k: number of remaining cancer samples of the same phylogeny
#   m: number of metastases samples
#   returns: probabilility of cluster by chance
    
    n = k + m
    
    no_trees = calculate_no_common_trees(k, m)
    
    no_trees / calculate_no_trees(n)
}
