# To run this script, you should have the packages "tensor" and
# "pracma" installed. If you save this script e.g. by the name
# "SK_LV_model.R", you can run it from the command line by invoking:
#
# Rscript SK_LV_model.R [outfile] [steps] [dt] [species] [loci]
#
# where the command line arguments are as follows:
# - [outfile]: name of file to save results in
# - [steps]: the number of time steps to run the program for (to reach
#            equilibrium, this may need to be large, e.g., 1e7)
# - [dt]: the size of each time step (e.g., 0.02)
# - [species]: the number of initial species (e.g., 51)
# - [loci]: the number of distinct loci contributing to the
#           quantitative trait (e.g., 25 means there are 25 loci and
#           therefore 2 * 25 + 1 = 51 distinct genotypes)
#
# The output is an (S x G) table, where S is the number of species and
# G is the number of genotypes. The (i,j)th entry is the density of
# species i's genotype G at the end of the simulation.

args <- commandArgs(trailingOnly=TRUE) # read command line arguments
fout <- args[1] # name of output file
numsteps <- as.numeric(args[2]) # total number of time steps
dt <- as.numeric(args[3]) # size of one time step
numsp <- as.numeric(args[4]) # number of species
n <- as.numeric(args[5]) # number of loci

# Load libraries
require(tensor)
require(pracma)

# Parameters
numtr <- 2*n+1 # Number of possible genotypes 
w <- 0.1 # width of the competition kernel
genotype <- seq(from=-0.6, to=0.6, l=numtr) # pre-environment phenotypes
sigma <- runif(numsp, min=0.01, max=0.05) # environmental noise widths
r0 <- outer(sigma, genotype, function(s, f) 1/2*erf((0.5-f)/(sqrt(2)*s)) +
            1/2*erf((0.5+f)/(sqrt(2)*s))) # intrinsic growth; species with
                                          # trait value outside the inverval
                                          # [-0.5, 0.5] have zero growth rate

# Competition coefficients, with environmental noise taken into account
A <- array(dim=c(numsp, numtr, numsp, numtr))
mu <- genotype
for (q in 1:numsp) for (i in 1:numtr) for (s in 1:numsp) {
    A[q,i,s,] <- w/sqrt(2*sigma[q]^2+2*sigma[s]^2+w^2) *
        exp(-(mu[i]-mu)^2/(2*sigma[q]^2+2*sigma[s]^2+w^2))
}

# Shpak-Kondrashov inheritance: approximate diploid model based on
# exact haploid model
haplR <- array(0, dim=rep(n+1,3))
for(i in 0:n) for(j in 0:i) for(k in 0:min(n,(i+j))){ 
    haplR[1+i,1+j,1+k] <- sum(dhyper(max(0, i+j-n):min(i, j), i, n-i, j) *
                              dbinom(k-(max(0, i+j-n):min(i,j)), i+j -
                                     2*(max(0, i+j-n):min(i, j)), prob=0.5))
}
for (k in 0:n) {
    haplR[,,1+k] <- haplR[,,1+k] + t(haplR[,,1+k])
    diag(haplR[,,1+k]) <- diag(haplR[,,1+k])/2
}
indexsum.haplR <- matrix(0, 2*n+1, 2*n+1)
for(k in 0:n){
    for(i in 0:n) indexsum.haplR[1+i,1+k] <- haplR[1+i,1,1+k]
    for(j in 0:n) indexsum.haplR[1+j+n,1+k] <- haplR[1+n,1+j,1+k]
}

R <- array(dim=rep(numtr, 3))
for (i in 0:(2*n)) for (j in 0:(2*n)) for (q in 0:(2*n)) {
    R[1+i,1+j,1+q] <- sum(indexsum.haplR[1+i,1+(0:q)] *
                          indexsum.haplR[1+j,1+q-(0:q)])
}

# Initial Conditions
N <- matrix(0, numsp, numtr)
for(sp in 1:numsp){ 
    tr <- sample(seq(numtr), size=1) 
    N[sp,tr] <- 1
}
N <- N + runif(length(N), 0, 0.1)

# Dynamics
for (simtime in 0:numsteps) {
    if ((simtime%%100)==0) { # Save data and eliminate extinct species
        A <- A[rowSums(N)>0,,rowSums(N)>0,] # every 100 time steps
        r0 <- r0[rowSums(N)>0,]
        N <- N[rowSums(N)>0,]
        N[N<1e-8] <- 0
        txt <- paste0("# time: ", simtime, " * ", dt, " (out of ",
                      numsteps, " * ", dt, ")\n",
                      "# sigmas: ", paste(sigma, collapse=" "))
        if (!any(is.na(N))) {
            write.table(txt, fout, quote=FALSE, row.names=FALSE,
                        col.names=FALSE)
            write.table(N, fout, append=TRUE, row.names=FALSE,
                        col.names=FALSE, quote=FALSE)
        } else stop("Abundances went NA!")
    }
    simtime <- simtime + 1
    
    # Reproduction
    Nprime <- sapply(seq(dim(R)[3]),
                     function(i) rowSums(tensor(N,R,2,1)[,,i]*N))
    Nprime <- N + Nprime/rowSums(Nprime)*rowSums(N)*r0*dt
    
    # Selection
    N <- Nprime*exp(-dt*tensor(A, Nprime, c(3,4), c(1,2)))
}
