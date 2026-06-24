source("./code/functions.R")
library(optparse)
### input parameters ###
option_list <- list(
  make_option(c("-e", "--seed"), 
              type = "integer", 
              default = NA,  # Default values
              help = "seed", 
              metavar = "integer"),
  make_option(c("-s", "--scene"), 
              type = "integer", 
              default = 1,  # Default values
              help = "Scene number (1-6)", 
              metavar = "integer"),
  make_option(c("-r", "--replicatenum"), 
              type = "numeric", 
              default = 200,  
              help = "number of replicates", 
              metavar = "numeric"),
  make_option(c("-K", "--Knum"), 
              type = "numeric", 
              default = 1000,  
              help = "total phenotypes", 
              metavar = "numeric"),
  make_option(c("-n", "--nsample"), 
              type = "numeric", 
              default = 2000,  
              help = "sample size", 
              metavar = "numeric"),
  make_option(c("-b", "--beta"), 
              type = "numeric", 
              default = 0.005,  
              help = "beta", 
              metavar = "numeric"),
  make_option(c("-f", "--MAF"), 
              type = "numeric", 
              default = 0.3,  
              help = "minor allele frequency", 
              metavar = "numeric"),
  make_option(c("-a", "--rhoa"), 
              type = "numeric", 
              default = 0.2,  
              help = "rho_a:control correlation across phenotypic categories", 
              metavar = "numeric"),
  make_option(c("-w", "--rhow"), 
              type = "numeric", 
              default = 0.3,  
              help = "rho_w:controls the correlation among phenotypes within each category", 
              metavar = "numeric"),
  make_option(c("-c", "--c2"), 
              type = "numeric", 
              default = 0.5,  
              help = "c^2", 
              metavar = "numeric"),
  make_option(c("-A", "--alpha"), 
              type = "numeric", 
              default = 0.05,  
              help = "FDR", 
              metavar = "numeric")
  
  
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

### set up parameters ###
seed <- opt$seed
replicates.num <- opt$replicatenum
n <- opt$nsample # sample size
K <- opt$Knum # total phenotypes
scene <- opt$scene # 1, 2, 3, 4, 5, 6
beta <- opt$beta
lambda <- lambdas_setup(beta = beta, scene = scene, K = K)
H1 <- which(sapply(lambda, function(x) any(x != 0))) # the categories whose are associated with x
MAF <- opt$MAF
rho_a <- opt$rhoa
rho_w  <- opt$rhow
c2 <- opt$c2
alpha <- opt$alpha
print(opt)
### run simulation ###
if(!is.na(seed)){
  set.seed(seed+n+K+scene+
             as.numeric(gsub("0\\.", "", as.character(beta)))+
             as.numeric(gsub("0\\.", "", as.character(MAF)))+
             as.numeric(gsub("0\\.", "", as.character(rho_a)))+
             as.numeric(gsub("0\\.", "", as.character(rho_w)))+
           as.numeric(gsub("0\\.", "", as.character(c2))))
}
seeds.replicate <- sample(0:100000, size = replicates.num)

KCLC.P <- matrix(nrow = replicates.num, ncol = length(lambda))
HCLC.P <- matrix(nrow = replicates.num, ncol = length(lambda))
MANOVA.P <- matrix(nrow = replicates.num, ncol = length(lambda))
multiPhen.P <- matrix(nrow = replicates.num, ncol = length(lambda))
TATES.P <- matrix(nrow = replicates.num, ncol = length(lambda))
for(i in 1:replicates.num){
  simulated.data <- simulation(n=n, MAF=MAF, lambda=lambda, rho_a=rho_a, rho_w=rho_w, c = sqrt(c2), seed = seeds.replicate[i])
  
  # KCLC
  Kmeans.clusters <- kmeans_Silhouette(simulated.data$y, simulated.data$cat)
  KCLC.P[i, ] <- CLC(simulated.data$y, simulated.data$x, Kmeans.clusters, simulated.data$cat) 
  
  # HCLC
  HCM.clusters <- HCM(simulated.data$y, simulated.data$cat)
  HCLC.P[i, ] <- CLC(simulated.data$y, simulated.data$x, HCM.clusters, simulated.data$cat)
  
  # MANOVA
  MANOVA.P[i, ] <- MANOVA.test(simulated.data$y, simulated.data$x, simulated.data$cat)
  
  # multiPhen
  multiPhen.P[i, ] <- multiphen.test(simulated.data$y, simulated.data$x, simulated.data$cat)
  
  # TATES
  TATES.P[i, ] <- TATES.test(simulated.data$y, simulated.data$x, simulated.data$cat)
}

outputdir <- "./output"
if (!dir.exists(outputdir)) {
  dir.create(outputdir, recursive = TRUE)
}
# create a subfolder, which is named by parameter used
subfolder <- paste0("scene", scene, 
                    "-RN", replicates.num, 
                    "-N", n, 
                    "-K", K, 
                    "-beta", gsub("0\\.", "", as.character(beta)),
                    "-rhoa", gsub("0\\.", "", as.character(rho_a)),
                    "-rhow", gsub("0\\.", "", as.character(rho_w)),
                    "-c2", gsub("0\\.", "", as.character(c2)),
                    "-timestamp", format(Sys.time(), "%Y_%m_%d_%H_%M_%S"))
subfolderdir <- paste0(outputdir, "/", subfolder)
if (!dir.exists(subfolderdir)) {
  dir.create(subfolderdir, recursive = TRUE)
}

### save P values ###
pdir <- paste0(subfolderdir, "/P")
if (!dir.exists(pdir)) {
  dir.create(pdir, recursive = TRUE)
}
saveRDS(KCLC.P, file = paste0(pdir, "/KCLC.RDS"))
saveRDS(HCLC.P, file = paste0(pdir, "/HCLC.RDS"))
saveRDS(MANOVA.P, file = paste0(pdir, "/MANOVA.RDS"))
saveRDS(multiPhen.P, file = paste0(pdir, "/multiPhen.RDS"))
saveRDS(TATES.P, file = paste0(pdir, "/TATES.RDS"))
writeLines(as.character(H1), paste0(pdir, "/H1.txt"))

### FDR ###
FDR_estimation <- function(P, alpha, H1) {
  fdr <- numeric(nrow(P))
  for(i in 1:nrow(P)){
    t <- T_estimation(P[i, ], alpha)
    stat <- sum(P[i, -H1]<=t)/max(1, sum(P[i, ]<= t))
    fdr[i] <- stat
  }
  return(list(fdr = fdr, mean = mean(fdr), variance = var(fdr)))
}

Power_estimation <- function(P, alpha, H1){
  power <- numeric(nrow(P))
  for(i in 1:dim(P)[1]){
    t <- T_estimation(P[i, ], alpha)
    po <- sum(P[i,Ha]<= t)/length(H1)
    power[i] <- po
  }
  return(list(power = power, mean = mean(power), variance = var(power)))
}


FDR_Power <- list(KCLC=list(fdr=FDR_estimation(KCLC.P, alpha, H1), 
                       power=Power_estimation(KCLC.P, alpha, H1)), 
             HCLC=list(fdr=FDR_estimation(HCLC.P, alpha, H1), 
                       power=Power_estimation(HCLC.P, alpha, H1)),
             MANOVA=list(fdr=FDR_estimation(MANOVA.P, alpha, H1), 
                         power=Power_estimation(MANOVA.P, alpha, H1)), 
             multiPhen=list(fdr=FDR_estimation(multiPhen.P, alpha, H1), 
                            power=Power_estimation(multiPhen.P, alpha, H1)), 
             TATES=list(fdr=FDR_estimation(TATES.P, alpha, H1), 
                        power=Power_estimation(TATES.P, alpha, H1)))
saveRDS(FDR_Power,  paste0(subfolderdir, "/FDR_Power.RDS"))


