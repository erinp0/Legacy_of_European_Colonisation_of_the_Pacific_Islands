library(rstan)
library("gtools")

# Firstly, prepare the data for stan
refpanel_profiles <- read.table("gb_refpanel.txt")
sumtotal=rowSums(refpanel_profiles)[1]

poly_panel <- read.table("gbonly_vectors.txt",header=T,row.names=1)
target=poly_panel
## header.true <- function(df) {
##   names(df) <- as.character(unlist(df[1,]))
##   df[-1,]
## }
## target = header.true(poly_panel)

shipcounts <- read.table("manifestcounts2.txt")

## Name matching
rownames(shipcounts)[rownames(shipcounts)=="Cumbria"]="GB_Cumbria"
rownames(shipcounts)[rownames(shipcounts)=="Devon"]="GB_Devon"
rownames(shipcounts)[rownames(shipcounts)=="Ireland"]="GB_Ireland"
rownames(shipcounts)[rownames(shipcounts)=="Italy"]="BB_italians"
rownames(shipcounts)[rownames(shipcounts)=="Lincoln"]="GB_Lincolnshire"
rownames(shipcounts)[rownames(shipcounts)=="Merseyside"]="GB_Merseyside"
rownames(shipcounts)[rownames(shipcounts)=="NE Scotland"]="GB_NE_Scot"
rownames(shipcounts)[rownames(shipcounts)=="Northumbria"]="GB_Northumberland"
rownames(shipcounts)[rownames(shipcounts)=="NW Scotland"]="GB_NW_Scot"
rownames(shipcounts)[rownames(shipcounts)=="NW Wales"]="GB_NW_Wales"
rownames(shipcounts)[rownames(shipcounts)=="Orkneys"]="GB_Orkney"
rownames(shipcounts)[rownames(shipcounts)=="S England"]="GB_S_England"
rownames(shipcounts)[rownames(shipcounts)=="SC England"]="GB_SC_England"
rownames(shipcounts)[rownames(shipcounts)=="Northumbria"]="GB_Northumberland"
rownames(shipcounts)[rownames(shipcounts)=="SE England"]="GB_SE_England"
rownames(shipcounts)[rownames(shipcounts)=="North America"]="HM3_tsi_usa"
rownames(shipcounts)[rownames(shipcounts)=="Yorkshire"]="GB_N_Yorkshire"
rownames(shipcounts)[rownames(shipcounts)=="Anglia"]="GB_Anglia"
rownames(shipcounts)[rownames(shipcounts)=="SE Wales"]="GB_SE_Wales"
rownames(shipcounts)[rownames(shipcounts)=="France"]="HGDP_french"
rownames(shipcounts)[rownames(shipcounts)=="NI SW Scotland"]="GB_Noi_SWscot"
rownames(shipcounts)[rownames(shipcounts)=="Cornwall"]="GB_Cornwall"
rownames(shipcounts)[rownames(shipcounts)=="Germany"]="BB_german"
rownames(shipcounts)[rownames(shipcounts)=="NC England"]="GB_NC_England"

#Tidying up ship manifest data
shipwide <- t(shipcounts)
colkeep <- colnames(shipwide)[colnames(shipwide) %in% colnames(refpanel_profiles)]
shipwide <- cbind(shipwide, "Pacific"=0)

#Tidying up ref panel data

#Deleting other ancestry groups
groups_of_pacific = c("HGDP_pima","HGDP_japanese","HGDP_cambodians","HGDP_papuan","NAN_Melanesian","HGDP_tujia","HGDP_han","HGDP_Burusho")
## ref_colkeep = c(colkeep,groups_of_pacific)
##refpanel_profiles = refpanel_profiles[ref_colkeep]
## remove=row.names(refpanel_profiles)[!(row.names(refpanel_profiles) %in% ref_colkeep)]


#summing contribution of pacific groups to each of the groups (adding pacific col)
pacificsums=rowSums(refpanel_profiles[,groups_of_pacific])
refpanel_profiles=cbind(refpanel_profiles[,colkeep],"Pacific"=pacificsums)
## refpanel_profiles = refpanel_profiles[(row.names(refpanel_profiles) %in% ref_colkeep),]
refpanel_profiles = refpanel_profiles/rowSums(refpanel_profiles)*sumtotal

#Summing the rows of pacific groups into one row
##pacific_rows <- refpanel_profiles[groups_of_pacific,]
##rowSums(pacific_rows)

source("nnls.R")
targetpacificsums=rowSums(target[,groups_of_pacific])
target=cbind(target[colkeep],"Pacific"=targetpacificsums)
target=target/rowSums(target)*sumtotal
targetA=admix.nnls.all(as.matrix(target),as.matrix(refpanel_profiles))
pacific_makeup=colMeans(targetA[,groups_of_pacific])
pacific_makeup=pacific_makeup / sum(pacific_makeup)

#summing up pacific ancestry to one group
pacificrow = matrix(pacific_makeup,nrow=1) %*% as.matrix(refpanel_profiles[groups_of_pacific,])

refpanel_profiles <-rbind(refpanel_profiles,"Pacific"=pacificrow)
############
#Now delete the individual pacific group rows
refpanel_profiles = refpanel_profiles[!(row.names(refpanel_profiles) %in% groups_of_pacific),]


#Now delete the individual pacific group col in favour of combined one
##refpanel_profiles = refpanel_profiles[,!(names(refpanel_profiles) %in% groups_of_pacific)]
#refpanel_profiles = refpanel_profiles/rowSums(refpanel_profiles)

#add pacific row for shipwide
shipwide = rbind(shipwide, "Pacific" =c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,750))
shiptotalcounts=rowSums(shipwide)
shipwide=shipwide/rowSums(shipwide)


# reorder the rows in refpanel to be same order as col in shipwide
new_ref = refpanel_profiles[colnames(shipwide),]

#Calculate B matrix by multiplying the reference panel for locations
#by the counts of the ship
B = data.matrix(shipwide)%*%data.matrix(new_ref)
B = B/rowSums(B)

#target data must be integer and sum to 128 for multinomial
#so we normalise to 1, times 128 and round. 
#70x128 for morgans in 2^g
target = target/rowSums(target)
target = round(target*70*128)

# ------------ HAVEN@T GOT BEYOND HERE ---------------------------------------


# This is the original frequency (187 crew (44 O, 143 S), 750 pacific)
f_old=shiptotalcounts/sum(shiptotalcounts) # c(44/937,143/937,750/937)

set.seed(2) # create simulated values which are reproducible

# M is number of individuals in our sample X
M = 30



list_of_data = list(K=dim(B)[1], 
      M=M,
      N=dim(B)[2],
      x=target,
      f_old=f_old,
      B=B,
      twototheG=128,
      inverse_sigma=1/1000)

# MODEL FITTED HERE
fit7g <- stan(
  file = "Mod4Stan.stan",  # Stan program
  data = list_of_data,  # named list of data
  chains = 2,              # number of Markov chains
  warmup = 1000,           # number of warmup iterations per chain
  iter = 6000,             # total number of iterations per chain
  cores = 1,               # number of cores (could use one per chain)
)

## We correlate the true B estimates with the observed ones
## and see the model performs well
## Extract all the parameters we will check


# CHECKING ESS RHAT ETC
print(fit7g)
neff = summary(fit7g)$summary[,'n_eff']
neff = neff/4000
print(neff)
pairs(fit7g, pars=c('f_true[1]','f_true[2]','f_true[3]','sigma_fchange'))

## Here is the output from the model

pars7g=as.data.frame(fit7g)

#col means of each param
mpars7g=colMeans(pars7g)

# estimates of indB, indF, sigma_fchange
indBest7g=matrix(mpars7g[grep("indB",names(mpars7g))],nrow=M)

indFest7g=matrix(mpars7g[grep("indf",names(mpars7g))],nrow=M)

sigmaest7g=mpars7g[grep("sigma_fchange",names(mpars7g))]

fest7g=mpars7g[grep("f_true",names(mpars7g))]


## Compare old mixing prop with estimate and observed average
#f_true_est is average of f_true from the model<- what we are interested in
#obs_f_true is average indF from the model
old_true = data.frame(f_old=f_old,
           f_true_estimate=fest7g,
           observed_f_true=colMeans(indFest7g))


#---------------HYPOTHESIS TESTING---------------------


#calculate sample cov matrix for f_true
sample_f_true = pars7g[,c("f_true[1]","f_true[2]","f_true[3]")]
sample_cov = cov(sample_f_true)

#calculate sample mean for f_true
sample_mean = fest7g

#calculate xbar - mu
xbar_mu = sample_mean - f_old

## mahalanobis distance
N=35
G=7
sigma_1.data <- c(f_old[1], 0,0, 0,f_old[2],0,0, 0,f_old[3])
sigma_1 <- matrix(sigma_1.data, nrow=3, ncol =3, byrow=TRUE)
sigma_2.data <- c(f_old[1]^2, f_old[1]*f_old[2],f_old[1]*f_old[3], f_old[1]*f_old[2],f_old[2]^2,f_old[3]*f_old[2],f_old[1]*f_old[3], f_old[3]*f_old[2],f_old[3]^2)
sigma_2 <- matrix(sigma_2.data, nrow=3, ncol =3, byrow=TRUE)


sigma_matrix <- ((1/((N/G)+1))*sigma_1) + (G^2/(G*N + N^2))*sigma_2
inv_sigma = solve(sigma_matrix,tol = 1e-18)
mahala = t(xbar_mu)%*% inv_sigma %*% xbar_mu



#--------------PLOts--------------------
library(Ternary)
TernaryPlot(alab = expression(bold('f_true'[1])), 
            blab = expression(bold('f_true'[2])),
            clab = expression(bold('f_true'[3])))
TernaryPoints(pars7g[c("f_true[1]", "f_true[2]","f_true[3]")] ,
              col='blue', pch=20, cex=0.5)
TernaryPoints(f_old, col='black', pch=4, cex=0.5)







#--------------OTHER ANALYSIS-----------------------

plot(indBest7g,list_of_data$indB)

plot(colMeans(list_of_data$indB),colMeans(indBest7g))

#store each individuals lambda estimate in table
est_of_lambda = data.frame(f_estimate_1 = indFest7g[,1],
                        f_estimate_2 =indFest7g[,2],
                        f_estimate_3 =indFest7g[,3])

#average lambda 
lambda_new_est = data.frame(Officer_lambda = mean(indFest7g[,1]),
                            Sailor_lambda =mean(indFest7g[,2]),
                            Pacific_lambda =mean(indFest7g[,3]))


# Calculate estimated number of gpps from each of the 3 groups (O,S,P)
est_gpps = data.frame(f_estimate_1 = indFest7g[,1]*128,
                      f_estimate_2 =indFest7g[,2]*128,
                      f_estimate_3 =indFest7g[,3]*128)



