calculate_no_trees <-
function(n) {
    
#    Function calculates the total number of distinct trees given by Edwards and Cavalli-Sforza, 1963
#    n: number of samples
#    returns: number of distinct trees
        
    factorial((2*n)-3)/((2^(n-2))*factorial(n-2))
}
