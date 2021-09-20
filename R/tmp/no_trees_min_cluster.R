no_trees_min_cluster <-
function(k, m, l){
    
#    Calculates the number of trees where m metastases form a cluster of size at least l
#    Implements Formula (S1) from SI...
#    k: number of primary tumor region samples
#    m: number of metastases samples
#    l: minimal cluster size
#    returns: number of distinct tree topologies where l out of m metastases form a cluster
    
    
    n = k + m
    no_trees = 0
    
    # this IF statement is an addition to the python script. If floor(m / l) == 0 python will not calculate the 
    # number of trees but rather set it to zero. That behavior is reproduced here.
    if ((floor(m / l))!=0) {
    for (no_clusters in 1:(floor(m / l))) {
        no_trees  <-  no_trees + ((-1)^(no_clusters-1) * multibinomial(m, l, no_clusters) / factorial(no_clusters) 
                     * no_trees_with_clusters(n - (no_clusters * l), rep.int(l, times = no_clusters)))
        }
        }
    no_trees

}
