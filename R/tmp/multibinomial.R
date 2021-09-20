multibinomial <-
function(m, l, cluster_numbers) {
   
#     Calculate the product of the binomials
#     m: number of remaining cancer samples
#     l: size of observed metastases cluster
#     cluster_numbers: list of possible coexisting clusters of size l for a given phylogeny
#     returns: product of binomials
    
    mb = 1
    for (c in 0:(cluster_numbers-1)) {
        mb <- mb*(choose(m - c * l, l))
    }
    mb
}
