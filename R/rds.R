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
        if(class(k)=="bigz") c <- gmp::as.bigz(c)
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
min_cluster_probability <- function(k, m, l, use_gmp) {
    if (l > m) {
        stop('L is larger than M. Make sure KLM values are accurate to continue.')
    }
    if (m == 1) {
        NA 
    } else if (l == 1){
        1
    }  else {
        # this is an addition to the python script to allow for RDS caculcations for over 10,000 samples
        if(use_gmp) {
            require(gmp)
            k <- as.bigz(k)
        }
        n = k + m
        no_trees = no_trees_min_cluster(k, m, l)
        if (is.infinite(no_trees) | is.infinite(calculate_no_trees(n))) {
            stop("The number of trees calculated are too large for R to handle. Reduce input KLM or use rds(...,use_gmp=T)")
        }
        # if the number of trees gets too large, R can't display it and will difine it as infinite.
        # to adapt the script to work with R, the RDS
        # uses the result from equation (S1) and implements equation (S2)
        as.numeric(no_trees / calculate_no_trees(n))
    }
}


##' rds
##' 
##' Main function to run RDS calculations on a cancer phylogeny (phylo object). Must use Naxerova Lab Normal/Primary/Metastasis tip label conventions.
##'
##' @param phy A phylo-class cancer phylogeny
##' @param drop_na Omit rows in the output data.frame with NA RDS values (detault=T)
##' @param use_gmp Use gmp R package to calculate RDS for large trees. Required for trees with more around 30 tips (default=T)
##' @return data.frame where rows are metastasis sample-types and columns are k, l, m, and RDS values.
##' @export
rds <- function(phy, drop_na=T, use_gmp=T){
    require(stringr)

    ## root the tree at the normal
    if(!is.rooted(phy)) {
        node.number <- grep('^N',phy$tip.label)
        phy <- phytools::reroot(phy, node.number=node.number)
    }

    # get sample type of each sample
    st <- sapply(phy$tip.label,function(x) unlist(strsplit(x,"[0-9]+"))[1])
    st[grep('PT[a-z]',st)] <- 'PT'
    st[which(st=="PerOv")] <- "Per"

    # select tumor samples
    sel_t  <- str_detect(st, "^N", negate=TRUE)          
    # tumor samples
    tt <- st[sel_t]

    # get number of samples for each tumor type
    df_tt <- data.frame(table(tt))

    # get metastasis types in the tree
    sel_m <- df_tt[,1]
    tumor_types <- as.character(df_tt[sel_m,1])

    ntt <- integer()
    for (tumor in df_tt[,1]) {
        ntt[tumor] <- df_tt[which(df_tt[,1]==tumor),2]
    }

    k <- integer()
    l <- integer()
    m <- integer()
    rootID <- length(phy$tip.label) + 1

    # skip klm calculation if only one sample type is present             
    if (length(tumor_types)<2) {return()}  

    for (met in tumor_types) {

        m[met] <- ntt[met]
        k[met] <- sum(ntt) - m[met]
        l[met] <- 1

        for (i in rootID:max(phy$edge[,1])) {
            clade <- extract.clade(phy, i)
            labels <- clade$tip.label

            #change to original KLM script, added map function because the original version would count 
            # samples with a different suffixes as different samples (LN1a wouldn't count as LN)
            labels_st <- sapply(labels,function(x) strsplit(x,"[0-9]+")[[1]][1])
            labels_st[grep('PT[a-z]',labels_st)] <- 'PT'
            labels_st[which(labels_st=="PerOv")] <- "Per"

            clade_size <- length(labels)
            # get the size of the largest clade formed by a given met

            if (length(which(labels_st==met))==clade_size) {
                if (clade_size > l[met]) {
                    l[met] = clade_size
                }
            }
        }
    }

    df <- data.frame(k=k,l=l,m=m)
    df <- cbind(type=rownames(df), df)

    for(i in 1:nrow(df)) {
        df$RDS[i] <- calculate_rds(k=df$k[i],m=df$m[i],l=df$l[i], use_gmp=use_gmp)
    }

    if(drop_na) df <- df[!is.na(df$RDS),]
    rownames(df) <- NULL
    return(df)
}  



##' calculate_rds
##' 
##' Calculates the root diversity scores as the probability that a tree where l metastases form a cluster (have one root) exists by chance as defined by Equation (2).
##' @param k number of primary tumor region samples
##' @param m total number of metastases samples
##' @param l number of metastases that actually cluster together
##' @param use_gmp Use gmp R package to calculate RDS for large trees. Required for trees with more around 30 tips (default=F)
##' @return probabilility that a cluster of size l arises by chance
##' @export
calculate_rds <- function(k,m,l, use_gmp) {
    RDS_vector <- Vectorize(min_cluster_probability)
    out <- RDS_vector(k, m, l, use_gmp=use_gmp)
    out
}



