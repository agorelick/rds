##' calculate_no_trees
##'
##' Calculate the total number of distinct trees given by Edwards and Cavalli-Sforza, 1963
##'
##' @param n number of samples
##' @return number of distinct trees
##' @export
calculate_no_trees <- function(n) {
    factorial((2*n)-3)/((2^(n-2))*factorial(n-2))
}


##' calculate_no_common_trees
##' 
##' Calculate the number of monophyletic trees given by Theorem 3 in the SI of Reiter et al. 2020 Nat Gen.
##' @param k number of remaining samples
##' @param m number of metastases samples
##' @return number of monophyletic trees
##' @export
calculate_no_common_trees <- function(k, m) {
    ((factorial(((2*k)-1)))/((2^(k-1))*(factorial(k-1))))*((factorial((2*m)-3))/((2^(m-2))*factorial(m-2))) 
}


##' cluster_probability
##' 
##' Calculate the probability that a tree where all m metastases form a cluster (have one root) exists by chance
##' @param k number of remaining cancer samples of the same phylogeny
##' @param m number of metastases samples
##' @return probabilility of cluster by chance
##' @export 
cluster_probability <- function(k, m) {
    n = k + m
    no_trees = calculate_no_common_trees(k, m)
    no_trees / calculate_no_trees(n)
}


##' multibinomial
##' 
##' Calculate the product of the binomials
##' @param m number of remaining cancer samples
##' @param l size of observed metastases cluster
##' @param cluster_numbers vector of possible coexisting clusters of size l for a given phylogeny
##' @return product of binomials
##' @export
multibinomial <- function(m, l, cluster_numbers) {
    mb = 1
    for (c in 0:(cluster_numbers-1)) {
        mb <- mb*(choose(m - c * l, l))
    }
    mb
}


##' no_trees_with_clusters
##' 
##' Calculates the number of trees where m1, m2, ..., mC, metastases form a cluster (have one root)
##' @param k number of remaining cancer samples of the same phylogeny
##' @param clusters vector of number of metastases samples per C clusters
##' @return number of trees where the given clusters have independent origin
##' @export
no_trees_with_clusters <- function(k, clusters){
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


##' no_trees_min_cluster
##' 
##' Calculates the number of trees where m metastases form a cluster of size at least l. Implements Formula (S1) from SI.
##' @param k number of primary tumor region samples
##' @param m number of metastases samples
##' @param l minimal cluster size
##' @return number of distinct tree topologies where l out of m metastases form a cluster
##' @export
no_trees_min_cluster <- function(k, m, l){
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


##' min_cluster_probability
##' 
##' Calculates the probability that a tree where all m metastases form a cluster (have one root) exists by chance
##' @param k number of primary tumor region samples
##' @param m total number of metastases samples
##' @param l number of metastases that actually cluster together
##' @return probabilility that a cluster of size l arises by chance
##' @export
min_cluster_probability <- function(k, m, l) {
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


##' rds
##' 
##' Calculates the root diversity scores as the probability that a tree where l metastases form a cluster (have one root) exists by chance as defined by Equation (2).
##' @param k number of primary tumor region samples
##' @param m total number of metastases samples
##' @param l number of metastases that actually cluster together
##' @return probabilility that a cluster of size l arises by chance
##' @export
rds <- function(k,m,l) {
    RDS_vector <- Vectorize(min_cluster_probability)
    out <- RDS_vector(k, m, l)
    out
}



