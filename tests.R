
# Simulate trait data

n <- 10 # number of indiviuals
m <- 3 # number of traits
trait.names <- c("size","horn","wing")

Z <- matrix(rnorm(n*m), nrow=n, ncol=m) # matrix of trait values
W <- rpois(n,lambda=4) # vector of absolute fitness values
colnames(Z) <- trait.names

# Calculate:
# distributional selection differentials (DSD)
# directional and non-directional components of the DSD (dD and dN)
# distributional selection gradients (delta)
# linear selection differentials (s)
# linear selection gradients (beta)
analyze.selection(Z,W)

# Maximizer functions for each trait
h.matrix <- maximizer(Z,W)

# Same as previous except that trait values are given alongside the corresponding maximizer values
z.h.matrix. <- maximizer(Z,W,include.traits=TRUE)

# Variance-covariance matrix for the maximizers (i.e. for the random variables h_i(Z_i))
H <- cov(h.matrix,h.matrix)

# DSD using the cumulative integral definition
DSD.cum <- DSD(Z,W)$DSD

# DSD using the covariance definition
DSD.cov <- as.vector(cov(W/mean(W),h.matrix))

# Check that DSD values are equal under the cumulative integral and covariance definitions
isTRUE(all.equal(DSD.cum,DSD.cov))

# Matrix product of H and delta
H.times.delta <- as.vector(H %*% delta(Z,W)$delta)

# Check that DSD = H delta
isTRUE(all.equal(DSD.cum,H.times.delta))

# p-values for the DSD, dD, dN, and delta using a permutation test (might take a little while to run)
selection.permutation.test(Z,W,bootlength=10000)

# Subset the trait 'horn'

z=Z[,"horn"]

# Optimal flow for 'horn'
myFlow <- minFlowSelection(z,W)

# Decompose the optimal flow into directional and non-directional components (D and N)
decomp <- flowDecomposition(z,myFlow)

# Compose D and N to recover the original flow
recomp <- flowComposition(z,decomp$D,decomp$N)

# Check that the composition of D and N equals the original flow
isTRUE(all.equal(myFlow,recomp))

# DSD on 'horn' using the Earth mover's definition (with an additional factor of n/(n-1) to reduce finite-sample bias and maintain compatibility with the other definitions)
DSD.EM.horn <- (n/(n-1))*sum(myFlow*abs(outer(z,z,"-")))

# Check that the DSD values are equal under the cumulative integral and Earth mover's definitions
isTRUE(all.equal(DSD(z,W)$DSD,DSD.EM.horn))

# Work done under the component flows D and N
DSD.D.horn <- (n/(n-1))*sum(decomp$D*abs(outer(z,z,"-")))
DSD.N.horn <- (n/(n-1))*sum(decomp$N*abs(outer(z,z,"-")))

# Check that the work done under D and N equals dD and dN respectively
isTRUE(all.equal(DSD.components(z,W)$dD,DSD.D.horn))
isTRUE(all.equal(DSD.components(z,W)$dN,DSD.N.horn))
