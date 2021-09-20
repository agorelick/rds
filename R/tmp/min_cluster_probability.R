min_cluster_probability <-
function(k, m, l) {
    
#    Calculates the probability that a tree where all m metastases form a cluster (have one root) exists by chance
#    k: number of primary tumor region samples
#    m: total number of metastases samples
#    l: number of metastases that actually cluster together
#    returns: probabilility that a cluster of size l arises by chance
    
    
    if (m == 1) {
        NA }
    else if (l == 1){
         1
        }  else {
    # this is an addition to the python script to allow for RDS caculcations for over 10,000 samples      
    k <- gmp::as.bigz(k)     
        
    n = k + m
      
    no_trees = no_trees_min_cluster(k, m, l)
        
    if (is.infinite(no_trees) | is.infinite(calculate_no_trees(n))) { 
    stop("The number of trees calculated are too large for R to handle. Reduce input KLM or change script.")
    }   
        
    # if the number of trees gets too large, R can't display it and will difine it as infinite.
    # to adapt the script to work with R, the RDS  
    # uses the result from equation (S1) and implements equation (S2)
        as.numeric(no_trees / calculate_no_trees(n)) 
        }
 }
