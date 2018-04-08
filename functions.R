# This script contains functions for the distributional selection differential (DSD), maximizers (h), and associated functions defined in the paper: # Henshaw JM, Zemel Y (2017). A unified measure of linear and nonlinear selection on quantitative traits. Methods in Ecology and Evolution 8(5): 604-614 (doi:10.1111/2041‑210X.12685)

# Functions for analyzing selection:

# FUNCTION NAME               RETURNS

# maximizer                   maximizers (h)
# DSD                         distributional selection differentials (DSD)
# DSD.components              DSD and their directional and non-directional components (dD and dN)
# delta                       distributional selection gradients (delta)
# s                           linear selection differentials (s)
# beta                        linear selection gradients (beta)
# analyze.selection           values for DSD, dD, dN, delta, s, and beta
# selection.permutation.test  p-values for DSD, dD, dN, and s, based on a simple permutation test
            
# Arguments for the above functions are:
# a vector or matrix of trait values (Z), and
# a vector of fitness values (W).
# If Z is a matrix (i.e. if there are two or more traits), each row should represent an individual and each column should represent a trait.

# Note one difference from Henshaw & Zemel (2017) is that the DSD and its components are calculated with a divisor of (n-1) instead of (n). This reduces finite-sample bias. See: Head ML, Kahn AT, Henshaw JM, Keogh JS, Jennions MD (2017). Sexual selection on male body size, genital length and heterozygosity: consistency across habitats and social settings. Journal of Animal Ecology 86(6): 1458-1468 (doi:10.1111/1365-2656.12742)

# Functions to calculate optimal flows and their components:

# FUNCTION NAME               RETURNS

# minFlowGeneral              optimal flow between two arbitrary trait distributions
# minFlowSelection            optimal flow between the trait distributions before and after selection
# flowDecomposition           breakdown of a flow into directional and non-directional component flows (D and N)
# flowComposition             composition of a directional flow D and non-directional flow N


#####     #####     #####     #####     #####


# 'input.check' returns errors if the trait matrix and fitness vector are not in appropriate formats

input.check<-function(Z,W)

{if(!is.vector(Z)&&!is.matrix(Z)&&!is.data.frame(Z)) stop("Z must be either: (1) a matrix or data frame with rows for each individual and columns for each trait, or (2) a vector of trait values for each individual, in the case where only one trait is being analysed")
  
  if(!is.vector(W)) stop("W must be a vector of fitness values")
  
  # Extract dimensions from Z
  if(is.vector(Z)){
    n <- length(Z)
    ntraits <- 1
  }
  else{
    n <- dim(Z)[1]
    ntraits <- dim(Z)[2]
  }
  
  if(n!=length(W))stop("The number of individuals must be equal for Z and W")

  if(!is.numeric(Z))stop("Trait values Z must be numeric")

  if(!is.numeric(W))stop("Fitness values W must be numeric")

  if(!is.finite(Z)||!is.finite(W))stop("Currently only complete data are supported. This may change in the future.")

  if(any(W<0))stop("Fitness values W must be non-negative")

  if(all(W==0))stop("At least one fitness value must be non-zero")

}

# 'standardize' returns variance-standardised trait values

standardize<-function(Z){(Z-mean(Z))/sd(Z)}


# 'maximizer' returns the maximizer for each trait in Z.
# If standard=TRUE, then traits are variance-standardized before analysis.
# If include.traits=TRUE, then trait values are returned alongside the maximizer values.
# If ordered=FALSE, then maximizers are returned in the same ordering as the inputted traits.
# If ordered=TRUE, then each set of trait values is sorted from smallest to largest, and maximizers are returned in the corresponding ordering (which is different for each trait). Note that rows no longer correspond to individuals in this case.

maximizer <- function(Z,W,standard=FALSE,include.traits=FALSE,ordered=FALSE){
    
  input.check(Z,W)
  
  # Extract dimensions from Z and standardize traits if required
  if(is.vector(Z)){
    n <- length(Z)
    ntraits <- 1
    if(standard){Z <- standardize(Z)}
  }
  else{
    n <- dim(Z)[1]
    ntraits <- dim(Z)[2]
    if(standard){Z<-apply(Z,2,standardize)}
  }
  
  # Set up a data frame h.matrix for the results, with column names given by cols
  if(include.traits){
    cols<-if (is.null(colnames(Z))) {paste(c("z","h.z"), as.vector(mapply(rep,c(1:ntraits),2)), sep = "")} else {paste(c("","h."), as.vector(mapply(rep,colnames(Z),2)), sep = "")}
    h.matrix<-data.frame(matrix(NA,nrow=n,ncol=2*ntraits,dimnames=list(NULL,cols)))
  }
  else{
    cols<-if (is.null(colnames(Z))) {paste("h.z", 1:ntraits, sep = "")} else {paste("h.", colnames(Z), sep = "")}
    h.matrix<-data.frame(matrix(NA,nrow=n,ncol=ntraits,dimnames=list(NULL,cols)))
  }
  
  # Calculate the maximizer for each trait Zi 
  for (i in c(1:ntraits)){
      
    if(is.vector(Z)) {Zi <- Z} else {Zi<-Z[,i]}
    
    # Order trait values from smallest to largest, and order fitness correspondingly
    or<-order(Zi)
    Zi<-Zi[or]
    Wi<-W[or]
    
    # Calculate the maximizer for Zi
    hdiff <- diff(Zi) * sign(mean(Wi)*c(1:(n-1)) - cumsum(Wi)[-n]) # Represents the differences h(Zi[k]) - h(Zi[k-1])
    h <- cumsum(c(0, hdiff)) # Represents h(Zi[k]), with an arbitrary value for h(Zi[1])
    h <- h - mean(h) # Normalize so that the mean of the maximizer is zero
    
    # Put the results in h.matrix
    if(include.traits){
      h.matrix[,2*i-1]<-if (ordered) {Zi} else {Zi[order(or)]}
      h.matrix[,2*i]<-if (ordered) {h} else {h[order(or)]}
    }
    else{
    h.matrix[,i]<-if (ordered) {h} else {h[order(or)]}
    }
  }
    
  h.matrix
    
}

# 'DSD' returns distributional selection differentials for the trait vector/matrix Z and fitness vector W.
# It is based on the cumulative integral definition of the DSD: see Henshaw & Zemel (2017).
# If standard=TRUE, then all traits are variance-standardized before analysis.

DSD<-function(Z,W,standard=FALSE)
  
{input.check(Z,W)
  
  # Extract dimensions from Z and standardize traits if required
  if(is.vector(Z)){
    n <- length(Z)
    ntraits <- 1
    if(standard){Z <- standardize(Z)}
  }
  else{
    n <- dim(Z)[1]
    ntraits <- dim(Z)[2]
    if(standard){Z<-apply(Z,2,standardize)}
  }
  
  # Set up a data frame DSD.vector for the results, with row names given by rows
  rows<-if (is.null(colnames(Z))) {paste("z", 1:ntraits, sep = "")} else {colnames(Z)}
  DSD.vector<-data.frame(DSD=vector(length=ntraits),row.names=rows)
  
  # Calculate the DSD for each trait Zi
  for (i in c(1:ntraits)){
  
    if(is.vector(Z)) {Zi <- Z} else {Zi<-Z[,i]}
    
    # Order trait values from smallest to largest, and order fitness correspondingly
    or<-order(Zi)
    Zi<-Zi[or]
    Wi<-W[or]
    
    # Calculate relative fitness
    w<-Wi/mean(Wi)
    
    # Calculate the DSD for Zi; the formula here is equivalent to that in Henshaw & Zemel (2017), except with a divisor of (n-1) instead of n to reduce finite-sample bias.
    cumulant<-cumsum(1-w)[-n]
    DSD.vector[i,]<-(1/(n-1))*sum(diff(Zi)*abs(cumulant))
    
    }

  DSD.vector

}

# 'DSD.components' returns the DSD along with their directional and non-directional components (dD and dN) for the trait vector/matrix Z and fitness vector W.
# If standard=TRUE, then all traits are variance-standardized before analysis.

DSD.components<-function(Z,W,standard=FALSE){
  
  input.check(Z,W)
  
  # Extract dimensions from Z and standardize traits if required
  if(is.vector(Z)){
    ntraits <- 1
    if(standard){Z <- standardize(Z)}
  }
  else{
    ntraits <- dim(Z)[2]
    if(standard){Z<-apply(Z,2,standardize)}
  }
  
  # Save values of the DSD
  DSDtemp=DSD(Z,W)$DSD
  
  # Row names for the results
  rows<-if (is.null(colnames(Z))) {paste("z", 1:ntraits, sep = "")} else {colnames(Z)}
  
  # Return a data frame with values for d, dD, and dN
  data.frame(DSD=DSDtemp,dD=abs(cov(Z,W/mean(W))),dN=DSDtemp-abs(cov(Z,W/mean(W))),row.names=rows)
  
}

# 'delta' returns distributional selection gradients for the trait vector/matrix Z and fitness vector W.
# If standard=TRUE, then all traits are variance-standardized before analysis.

delta<-function(Z,W,standard=FALSE){
  
  input.check(Z,W)
  
  # Extract dimensions from Z and standardize traits if required
  if(is.vector(Z)){
    ntraits <- 1
    if(standard){Z <- standardize(Z)}
  }
  else{
    ntraits <- dim(Z)[2]
    if(standard){Z<-apply(Z,2,standardize)}
  }
  
  # Row names for the results
  rows<-if (is.null(colnames(Z))) {paste("z", 1:ntraits, sep = "")} else {colnames(Z)}
  
  # Return a data frame with distributional selection gradients (delta).
  # delta values are calculated as coefficients of the regression of relative fitness on the maximizers.
  data.frame(delta=lm(W/mean(W) ~ 1 + data.matrix(maximizer(Z,W)))$coefficients[2:(1+ntraits)],row.names=rows)
  
}


# 's' returns linear selection differentials for the trait vector/matrix Z and fitness vector W.
# If standard=TRUE, then all traits are variance-standardized before analysis.

s<-function(Z,W,standard=FALSE)
  
{input.check(Z,W)
  
  # Extract dimensions from Z and standardize traits if required
  if(is.vector(Z)){
    ntraits <- 1
    if(standard){Z <- standardize(Z)}
  }
  else{
    ntraits <- dim(Z)[2]
    if(standard){Z<-apply(Z,2,standardize)}
  }
  
  # Row names for the results
  rows<-if (is.null(colnames(Z))) {paste("z", 1:ntraits, sep = "")} else {colnames(Z)}
  
  # Return a data frame with linear selection differentials (s)
  data.frame(s=cov(Z,W/mean(W)),row.names=rows)
  
}

# 'beta' returns linear selection gradients for the trait vector/matrix Z and fitness vector W.
# If standard=TRUE, then all traits are variance-standardized before analysis.

beta<-function(Z,W,standard=FALSE)
  
{input.check(Z,W)
  
  # Extract dimensions from Z and standardize traits if required
  if(is.vector(Z)){
    ntraits <- 1
    if(standard){Z <- standardize(Z)}
  }
  else{
    ntraits <- dim(Z)[2]
    if(standard){Z<-apply(Z,2,standardize)}
  }
  
  # Row names for the results
  rows<-if (is.null(colnames(Z))) {paste("z", 1:ntraits, sep = "")} else {colnames(Z)}
  
  # Center traits for regression
  Z.centered <- Z - colMeans(Z)
  
  # Return a data frame with linear selection gradients (beta).
  # beta values are calculated as coefficients of the regression of relative fitness on trait values.
  data.frame(beta=lm(W/mean(W) ~ 1 + Z.centered)$coefficients[2:(1+ntraits)],row.names=rows)
  
}

# 'analyze.selection' merges results from the functions DSD.components, delta, s and beta.
# It returns the DSD, dD, dN, delta, s, and beta.
# If standard=TRUE, then all traits are variance-standardized before analysis.

analyze.selection <- function(Z,W,standard=FALSE){
  
  cbind(DSD.components(Z,W,standard), delta(Z,W,standard), s(Z,W,standard), beta(Z,W,standard))

}

# 'selection.permutation.test' returns p-values for d, dD, dN and s, based on a simple permutation test.
# If standard=TRUE, then all traits are variance-standardized before analysis.
# bootlength determines the number of permutations used; more permutations means more accurate p-values.

selection.permutation.test<-function(Z,W,standard=FALSE,bootlength=10000)
  
{input.check(Z,W)
  
  # Extract dimensions from Z and standardize traits if required
  if(is.vector(Z)){
    ntraits <- 1
    if(standard){Z <- standardize(Z)}
  }
  else{
    ntraits <- dim(Z)[2]
    if(standard){Z<-apply(Z,2,standardize)}
  }
  
  # Row names for the results
  rows<-if (is.null(colnames(Z))) {paste("z", 1:ntraits, sep = "")} else {colnames(Z)}
  
  # Set up a data frame for the results
  results<-data.frame(DSD=vector(length=ntraits),
                      pvalue.DSD=vector(length=ntraits),
                      dD=vector(length=ntraits),
                      pvalue.dD=vector(length=ntraits),
                      dN=vector(length=ntraits),
                      pvalue.dN=vector(length=ntraits),
                      s=vector(length=ntraits),
                      pvalue.s=vector(length=ntraits),
                      row.names=rows)
  
  # Save the observed values of d, dD, dN and s to the results data frame
  DSD.components.temp=DSD.components(Z,W)
  results$DSD<-DSD.components.temp$DSD
  results$dD<-DSD.components.temp$dD
  results$dN<-DSD.components.temp$dN
  results$s<-s(Z,W)$s
  
  #Set up tables that will later contain values of d, dD, dN and s corresponding to each permutation
  boottable.DSD <- matrix(nrow = bootlength, ncol = ntraits)
  boottable.dD <- matrix(nrow = bootlength, ncol = ntraits)
  boottable.dN <- matrix(nrow = bootlength, ncol = ntraits)
  boottable.s <- matrix(nrow = bootlength, ncol = ntraits)
  
  # For each permutation, run this loop
  for(i in 1:bootlength){
    
    # Permute fitness values
    Wtemp<-sample(W)
    
    # Save values of d, dD, dN and s for this permutation to the boottables
    DSD.components.temp=DSD.components(Z,Wtemp)
    boottable.DSD[i,]<-DSD.components.temp$DSD
    boottable.dD[i,]<-DSD.components.temp$dD
    boottable.dN[i,]<-DSD.components.temp$dN
    boottable.s[i,]<-s(Z,Wtemp)$s
  
  }
  
  # Calculate p-values for d, dD, dN and s for each trait as the proportion of boottable values that exceed the observed value (for s, absolute values are used)
  for (i in 1:ntraits){
    
    results$pvalue.DSD[i] <- length(boottable.DSD[,i][boottable.DSD[,i]>results$DSD[i]])/bootlength
    results$pvalue.dD[i]  <- length(boottable.dD[,i][boottable.dD[,i]>results$dD[i]])/bootlength
    results$pvalue.dN[i]  <- length(boottable.dN[,i][boottable.dN[,i]>results$dN[i]])/bootlength
    results$pvalue.s[i]  <- length(boottable.s[,i][abs(boottable.s[,i])>abs(results$s[i])])/bootlength
    
  }
  
  results
  
}

# 'minFlowGeneral' returns an optimal flow between two arbitrary distributions (p and q) over the trait values z.
# The output is a square matrix, where the (i,j)th entry represents the mass moved from ith to the jth trait value.
# Note that z is a vector, not a matrix.
# The distributions (p and q) are vectors of probabilities adding up to one.
# For the particular case of the trait distributions before and after fitness, see the function minFlowSelection below.
# If standard=TRUE, then trait values are variance-standardized before analysis.
# If ordered=FALSE, the flow matrix is returned with rows and columns in the same ordering as the inputted trait values.
# If ordered=TRUE, then trait values are sorted from smallest to largest, and flow matrix is ordered correspondingly.

minFlowGeneral <- function(z, p, q,standard=FALSE,ordered=FALSE)
{
  n<-length(z)
  
  # Check input and return errors if there are any problems
  if(!is.vector(z)||!is.vector(p)||!is.vector(q)||!isTRUE(all.equal(n,length(p),length(q)))) stop("z, p and q must be vectors of the same length")
  if(any(p < 0)||!isTRUE(all.equal(sum(p), 1)))	stop("p should be a vector of non-negative numbers that sum to one")
  if(any(q < 0)||!isTRUE(all.equal(sum(q), 1)))	stop("q should be a vector of non-negative numbers that sum to one")
  
  # Standardize trait values if required
  z <- standardize(z)
  
  # Reorder trait values
  or <- order(z)
  z <- z[or]
  p <- p[or]
  q <- q[or]
  
  # Recast the optimal transportation problem as a linear program.
  # This can be done using a more dedicated algorithm like the Hungarian method.
  constraints <- matrix(0, nrow = 2*n, ncol = n^2)
  for(i in 1:n)
  {
    constraints[i, n*0:(n-1) + i] <- 1
    constraints[i + n, 1:n + (i-1)*n] <- 1
  }
  
  # Remove last redundant constraint
  constraints <- constraints[-2*n, ]
  nconst <- 2*n - 1 # =nrow(constraints)
  
  cost.vec <- abs(outer(z, z, "-"))  # Distance between any pair of trait values (in L1 norm)
  
  res <- lpSolve::lp(direction = "min", objective.in = cost.vec, const.mat = constraints, const.dir = rep("==", nconst), const.rhs = c(p, q[-n]))
  
  flow <- matrix(res$solution, ncol = n, nrow = n, byrow = FALSE)
  
  if(ordered){flow}else{flow[order(or),order(or)]}
}

# 'minFlowSelectionl' returns an optimal flow between the trait distributions before and after selection.
# z is a vector of trait values and W is a vector of fitness values.
# If standard=TRUE, then trait values are variance-standardized before analysis.
# If ordered=FALSE, the flow matrix is returned with rows and columns in the same ordering as the inputted trait values.
# If ordered=TRUE, then trait values are sorted from smallest to largest, and flow matrix is ordered correspondingly.

minFlowSelection <- function(z, W,standard=FALSE,ordered=FALSE){
  
  n = length(z)

  # Check input and return errors if there are any problems
  if(!is.vector(z)||!is.vector(W)||!isTRUE(all.equal(n,length(W)))) stop("Trait values z and fitness values W must be vectors of the same length")
  if(any(W < 0))	stop("Fitness values W must be non-negative")
  
  minFlowGeneral(z,rep(1/n,n),W/sum(W),standard,ordered)
  
}

# 'flowDecomposition' decomposes any flow over trait values z into directional and non-directional component flows (D and N).
# z is a vector of trait values (not a matrix).
# If ordered=FALSE, then D and N are returned with rows and columns in the same ordering as the inputted trait values.
# If ordered=TRUE, then trait values are sorted from smallest to largest, and D and N are ordered correspondingly.

flowDecomposition <- function(z,flow,ordered=FALSE){
  
  n <- length(z)
  
  # Check input and return errors if there are any problems
  if(!is.vector(z)) stop("Trait values z must be a vector")
  if(!is.matrix(flow)||!isTRUE(all.equal(c(n, n), dim(flow)))) stop("flow must be a square matrix of the same dimension as the trait vector z")
  if(any(flow<0)||!isTRUE(all.equal(sum(flow),1)))	stop("flow should be a matrix of non-negative numbers that sum to one")
  
  # Order the trait vector and flow matrix
  or<-order(z)
  z<-z[or]
  flow<-flow[or,or]
  
  # Calculate D and N as per the Supporting Information of Henshaw & Zemel (2017)
  p <- rowSums(flow)
  q <- colSums(flow)
  L <- ! (H <- row(flow) < col(flow))
  l <- sum(   (flow * outer(z, z, "-"))[L]    )
  s <- sum(z * (q - p))
  if(s <= 0)
  {
    if(identical(l, 0))	D <- 0*flow else
    {
      D <- -flow*s/l
      D[H] <- 0
    }
  } else
  {
    h <- l + s
    D <- flow*s/h
    D[L] <- 0
  }
  
  D <- D + diag(p - rowSums(D))
  N <- flow - D
  N <- N + diag(q - colSums(N))
  
  if(!ordered){D<-D[order(or),order(or)]
  N<-N[order(or),order(or)]}
  
  list(D = D, N = N)
}

# 'flowComposition' returns the composition of a directional flow D and a non-directional flow N over trait values z.
# z is a vector, not a matrix.

flowComposition <- function(z,D,N)
{
  # Check input and return errors if there are any problems
  if(!is.matrix(N) || !is.matrix(D) || !is.numeric(D) || !is.numeric(N))  stop("D and N must be numeric matrices")
  n <- dim(D)[1]
  if(!isTRUE(all.equal(c(n, n), dim(D),dim(N))))  stop("D and N must be square matrices of the same dimension")
  if(!is.vector(z)||!is.numeric(z)) stop("z must be a numeric vector of trait values")
  if(any(N < 0) || any(D < 0))	stop("All enties of D and N must be non-negative")
  
  # Reorder trait values and component flows
  or <- order(z)
  D <- D[or,or]
  N <- N[or,or]
  
  if(!isTRUE(all.equal(colSums(D), rowSums(N))))  stop("The column sums of D must equal the row sums of N")
  if(!isTRUE(all.equal(D[lower.tri(D)], rep(0, n*(n-1)/2)))
    && !isTRUE(all.equal(D[upper.tri(D)], rep(0, n*(n-1)/2)))  )  warning("D must be a triangular matrix")
  
  flow <- D + N - diag(colSums(D))
  
  flow[order(or),order(or)]
}