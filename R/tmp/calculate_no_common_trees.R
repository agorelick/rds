calculate_no_common_trees <-
function(k, m) {
    
#    Function calculates the number of monophyletic trees given by Theorem 3 in the SI of Reiter et al. 2020 Nat Gen.
#    k: number of remaining samples
#    m: number of metastases samples
#    returns: number of distinct trees
   
    ((factorial(((2*k)-1)))/((2^(k-1))*(factorial(k-1))))*((factorial((2*m)-3))/((2^(m-2))*factorial(m-2))) 
}
