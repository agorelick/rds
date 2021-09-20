rds <-
function(k,m,l) {
    #Calculates the root diversity scores as the probability that a tree 
    #where l metastases form a cluster (have one root) exists by chance 
    #as defined by Equation (2)
    #:param k: number of primary tumor region samples
    #:param m: total number of metastases samples
    #:param l: number of metastases that actually cluster together
    #:return probabilility that a cluster of size l arises by chance

    if(l==1) {
        return(1)
    } else if(l < 1) {
        stop('Size of largest cluster l can not be smaller than 1')
    } else {
        n <- k + m
        no_trees <- as.numeric(no_trees_min_cluster(k, m, l))
    
        # uses the result from equation (S1) and implements equation (S2)
        no_trees / as.numeric(calculate_no_trees(n))
    }
}
