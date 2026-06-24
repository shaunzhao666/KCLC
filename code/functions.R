library(mvtnorm)
library(MASS)
library(cluster)
library(MultiPhen)
### functions related to simulation ###
lambdas_setup <- function(beta, scene, K){
  # set up the 6 scenes lambdas in paper
  # beta: control effect size
  # scene: {1, 2, 3, 4, 5, 6}
  # K: total number of phenotypes
  if(!(scene %in% c(1, 2, 3, 4, 5, 6))){
    stop("the scene number should be 1, 2, 3, 4, 5 or 6")
  }else{
    M <- ifelse(scene%%2==1, 50, 100)
    K_star <- K/M
    if (scene <= 2) {
      lambda <- c(replicate(M-1, rep(0, K_star), simplify = FALSE) , list(beta*c(1:K_star)))
    }else if (scene <= 4 & scene > 2){
      lambda <-c(replicate(M-2, rep(0, K_star), simplify = FALSE), 
                 list(((2*beta)/(K_star+1))*c(1:K_star)), 
                 list(2*beta*c(rep(1, K_star/2), rep(0, K_star/2))))
    }else{
      lambda <- c(replicate(M-3, rep(0, K_star), simplify = FALSE), 
                  list(c(c(1:(K_star/2)), c((K_star/2):1))*(beta/(1+K_star/2))), 
                  list(((2*beta)/(K_star+1))*c(1:K_star)), list(2*beta*c(rep(1, K_star/2), rep(0, K_star/2))))
    }
    return(lambda)
  }
}

simulation <- function(n, MAF, lambda, rho_a, rho_w, c = sqrt(0.5), seed = NA){
  # simulation according to the description in paper. All categories have same number of phenotypes.
  # n: number of subjects;
  # cat: {(11, ..., 1k*), ..., (M1, ..., Mk*)} show categories
  # MAF: minor allele frequency
  # lambda: {(lambda11...lambda1k*), ...(lambdamM1, ....lambamMk*)} effect size
  # rho_a: control correlation across phenotypic categories
  # rho_w: controls the correlation among phenotypes within each category
  # c: constant
  # seed: optional
  counter <- 1
  cat <- lapply(lambda, function(x) {
    n <- length(x)
    vals <- counter:(counter + n - 1)
    counter <<- counter + n
    return(vals)
  })
  
  M <- length(cat) # M: number of categories;
  K_star <- lengths(cat)[1] # number of phenotypes in each category
  K <- sum(lengths(cat)) # K: total number of phenotypes
  
  # generate x (genotypes)
  if(!is.na(seed)){
    set.seed(seed)
  }
  x <- matrix(sample(0:2, size = n, replace = TRUE, prob = c((1-MAF)^2, 2*MAF*(1-MAF), MAF^2)), nrow = n, byrow = TRUE )#R^n
  
  # generate A
  A <- matrix(1, nrow = M, ncol = M)
  Sigma_a <- (1 - rho_a)*diag(M) + rho_a * A # M*M matrix 
  if(!is.na(seed)){
    set.seed(seed)
  }
  A <- rmvnorm(n, mean=rep(0, M), sigma=Sigma_a) # n*M matrix
  
  # generate Sigma_w 
  Sigma_w <- matrix(NA, nrow=K_star, ncol=K_star)
  for(i in 1:K_star){
    for(j in 1:K_star){
      Sigma_w[i, j] <- rho_w^abs(i - j)
    }
  }
  
  # initialize y
  y <- matrix(NA, nrow=n, ncol=K)
  if(!is.na(seed)){
    set.seed(seed)
  }
  subseeds <- sample(0:10000, size = M)
  for(m in 1:M){
    # generate E
    set.seed(subseeds[m])
    E <- rmvnorm(n, mean = rep(0, K_star), sigma = Sigma_w) # n*K_star matrix
    ym <- x %*% matrix(lambda[[m]], nrow=1, byrow = TRUE) + c*matrix(A[, m], nrow = n, byrow = TRUE)%*%matrix(1, nrow=1, ncol=K_star) + sqrt(1-c^2)*E
    y[, ((m - 1) * K_star + 1):(m * K_star)] <- ym
  }
  return(list(y = y, x = x, cat = cat))
}



### Kmeans ###
kmeans_Silhouette <- function(y, cat){
  # using the Silhouette score  to find the optimum k (find the elbow) 
  M <- length(cat)
  K <- ncol(y)
  n <- nrow(y)
  nested_list <- list() # a list of M groups of clusters
  for(m in 1:M){
    y.sub <- y[, cat[[m]]]
    
    # find the optimal k
    maxK <- ncol(y.sub)-1
    ss.avg <- c()
    foroptimalclusters <- list()
    for(k in 2:maxK){
      clusters <- kmeans(t(y.sub), centers = k)
      foroptimalclusters[[k]] <- clusters
      ss <- silhouette(clusters$cluster, dist = dist(t(y.sub)))
      ss.avg <- c(ss.avg, mean(ss[,3]))
      
    }
    optimal_k <- which.max(ss.avg) + 1
    optimalcluster <- foroptimalclusters[[optimal_k]]
    clusters <- lapply(1:optimal_k, function(x) as.list(which(optimalcluster$cluster == x)))
    
    nested_list[[m]] <- clusters
  }
  return(nested_list)
}


### HCM ###
dissimilarity_average <- function(x, y) {
  dissimilarity <- c()
  for (i in 1:ncol(x)){
    for (j in 1:ncol(y)){
      dissimilarity <- c(dissimilarity, 1-cor(x[, i], y[, j])) # correlation
    }
  }
  return(mean(dissimilarity))
}

compute_distance_matrix <- function(data, clusters){
  # clusters <- list([1], [[2], [3]], [[4], [5], [6]]) like this
  n <- length(clusters)
  dist_matrix <- matrix(0, n, n)
  for( i in 1:(n-1)){
    for( j in (i+1):n){
      dist_matrix[i, j] <- dissimilarity_average(data[,unlist(clusters[[i]]), drop = FALSE], data[, unlist(clusters[[j]]), drop = FALSE])
      dist_matrix[j, i] <- dist_matrix[i, j]
    }
  }
  return(dist_matrix)
}

HCM <- function(y, cat){
  M <- length(cat)
  K <- ncol(y)
  n <- nrow(y)
  nested_list <- list() # a list of M groups of clusters
  for(m in 1:M){
    y.sub <- y[, cat[[m]]]
    clusters <- lapply(1:ncol(y.sub), function(x) list(x))# initialize
    h <- c()
    while(length(clusters) > 1){
      dissimilarity <- compute_distance_matrix(y.sub, clusters)
      diag(dissimilarity) <- NA
      hb <- min(dissimilarity, na.rm = TRUE)
      h <- c(h, hb)
      indices <- which(dissimilarity == hb, arr.ind = TRUE)
      new_cluster <- c(clusters[[indices[1, 1]]], clusters[[indices[1, 2]]])
      clusters <- c(clusters[-c(indices[1, 1], indices[1, 2])], list(new_cluster))
      
    }
    diff_h <- h[-1] - h[1:(length(h) - 1)]
    k_hat <- which.max(diff_h)  # when to stop
    clusters <-lapply(1:ncol(y.sub), function(x) list(x))# initialize
    h_hat <- h[k_hat]
    hb <- 0
    while(hb < h_hat){
      dissmilarity <- compute_distance_matrix(y.sub, clusters)
      diag(dissmilarity) <- NA
      hb <- min(dissmilarity, na.rm = TRUE)
      indices <- which(dissmilarity == hb, arr.ind = TRUE)
      new_cluster <- c(clusters[[indices[1, 1]]], clusters[[indices[1, 2]]])
      clusters <- c(clusters[-c(indices[1, 1], indices[1, 2])], list(new_cluster))
    }
    Lm <- length(clusters)
    nested_list[[m]] <- clusters
  }
  return(nested_list)
}


## CLC ###
CLC <- function(y, x, clusters, cat){
  # clusters: M groups of clusters
  M <- length(cat)
  K <- ncol(y)
  n <- nrow(y)
  p <- c()
  for(m in 1:M){
    y.sub <- y[, cat[[m]]]
    m.clusters <- clusters[[m]]
    K_star <- length(cat[[m]])
    # T calculation
    x_bar <- mean(x)
    U.m <- t(y.sub) %*% matrix(x-x_bar, nrow = n, byrow = TRUE)  # a K_star*1 matrix
    var_x <- var(x)
    var_y <- sapply(1:ncol(y.sub), function(i) var(y.sub[, i]))
    V.m <- var_y*var_x[1, 1]*n 
    T.m <- matrix(sapply(1:K_star, function(i) U.m[i, 1]/sqrt(V.m[i])), nrow = K_star, byrow = TRUE) # a K_star*1 matrix
    
    # generate the B based on clusters
    B.m <- matrix(0, nrow = K_star, ncol = length(m.clusters))
    for(i in 1:K_star){
      one.loc <- which(sapply(m.clusters, function(x)  i %in% x))
      B.m[i, one.loc] <- 1
    }
    # similarity matrix of ym (correlation)
    Sigma.m <- cor(y.sub)
    W.m <- t(B.m) %*% ginv(Sigma.m) # a Lm*K_star matrix, sigma.m general inverse
    
    # T^Lm_CLC calculation
    T_CLC.m <- t(W.m %*% T.m) %*% ginv(W.m %*% Sigma.m %*% t(W.m)) %*% (W.m %*% T.m)
    
    p.m <- 1-pchisq(T_CLC.m, df = length(m.clusters))
    
    p <- c(p, p.m)
  }
  return(p)
}

### MANOVA ###
MANOVA.test <- function(y, x, cat){
  K <- ncol(y)
  n <- nrow(y)
  M <- length(cat)
  p <- c()
  for(m in 1:M){
    y.sub <- y[, cat[[m]]]
    manova.m <- summary(manova(y.sub~x))
    p <- c(p, manova.m$stats[1, "Pr(>F)"])
  }
  return(p)
}

### multiPhen ###
multiphen.test <- function(y, x, cat){
  K <- ncol(y)
  n <- nrow(y)
  M <- length(cat)
  p <- c()
  row.names(x) <- c(1:(dim(x)[1]))
  colnames(x) <- c(1:(dim(x)[2]))
  opts = mPhen.options(c("regression","pheno.input"))
  
  for(m in 1:M){
    y.sub <- y[, cat[[m]]]
    row.names(y.sub) <- c(1:(dim(y.sub)[1]))
    colnames(y.sub) <- c(1:(dim(y.sub)[2]))
    mphen.model <- mPhen(x , y.sub, opts = opts)
    p <- c(p, mphen.model$Results[1, 1, , ]["JointModel", "pvalue"])
  }
  return(p)
}

### TATES ###
rho_r_approximation <- function(r){
  return(-0.0008-0.0023*r+0.6226*r^2+0.0149*r^3+0.1095*r^4-0.0219*r^5+0.2179*r^6)
} # approximation

mej_cal <- function(y, j){
  if(j == 1){
    rho <- c(1)
  }else{
    r <- cor(y[, 1:j])
    non.diag.mask <- lower.tri(r) | upper.tri(r)
    rho <- r
    rho[non.diag.mask] <- rho_r_approximation(r[non.diag.mask])
  }
  lambda <- eigen(rho)$values
  return(j - sum(lambda>1))
}

TATES.test <- function(y, x, cat){
  K <- ncol(y)
  n <- nrow(y)
  M <- length(cat)
  p <- c()
  for(m in 1:M){
    y.sub <- y[, cat[[m]]]
    K_star <- length(cat[[m]])
    # calculate p by GWAS
    p.values <- sapply(1:K_star, function(i) {
      model <- lm(y.sub[, i] ~ x)
      summary_model <- summary(model)
      summary_model$coefficients["x", "Pr(>|t|)"]
    })
    # sort p and rearrange phenotype by sorted p
    sorted.indices <- order(p.values)
    p.ascend <- p.values[sorted.indices]
    sorted.y.sub <- y.sub[, sorted.indices]
    
    me <- mej_cal(sorted.y.sub, K_star)
    PT.Candidate <- lapply(1:K_star, function(j){
      me*p.ascend[j]/mej_cal(sorted.y.sub, j)
    })
    p <- c(p, min(unlist(PT.Candidate)))
  }
  return(p)
}


### FC ###
T_estimation <- function(P, alpha){
  alpha_m0 <- alpha/length(P)
  # print(length(P))
  p_desc <- sort(P, decreasing = TRUE)
  t_choose_bound<- c(p_desc, 0)
  j <- 1
  candidate <- 0
  I_M_vec <- c(rev(seq_along(P)), 1)
  tests <- (t_choose_bound <= I_M_vec*alpha_m0)
  max_value <- which(tests == TRUE)[1]
  t <- I_M_vec[max_value] * alpha_m0
  return(t)
}




