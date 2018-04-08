# DSD
R code for calculating the distributional selection differential (DSD), maximizers (h), and associated functions defined in the paper:

Henshaw JM, Zemel Y (2017). A unified measure of linear and nonlinear selection on quantitative traits. Methods in Ecology and Evolution 8(5): 604-614 (doi:10.1111/2041‑210X.12685)

Requires the R package lpSolve for the calculation of optimal flows.

## Contents
### functions.R
This file contains the following functions:

| Function | Description |
| --- | --- |
| maximizer | calculates the maximizers (h) for a vector/matrix of trait values Z and a vector of fitness values W |
| DSD | calculates the distributional selection differentials (DSD) for a vector/matrix of trait values Z and a vector of fitness values W |
| DSD.components | calculates the DSD along with their directional and non-directional components (dD and dN) |
| delta | calculates the distributional selection gradients (delta) for a vector/matrix of trait values Z and a vector of fitness values W |
| s | calculates the linear selection differentials (s) for a vector/matrix of trait values Z and a vector of fitness values W |
| beta | calculates the linear selection gradients (beta) for a vector/matrix of trait values Z and a vector of fitness values W |
| analyze.selection | calculates DSD, dD, dN, delta, s, and beta for a vector/matrix of trait values Z and a vector of fitness values W |
| selection.permutation.test | calculates p-values for DSD, dD, dN, and s, based on a simple permutation test |
| minFlowGeneral | calculates an optimal flow between two arbitrary trait distributions (p and q) over trait values z |
| minFlowSelection | calculates an optimal flow between the trait distributions before and after selection for a vector of trait values z and a vector of fitness values W|
| flowDecomposition | decomposes any flow over a vector of trait values z into directional and non-directional component flows (D and N)|
| flowComposition | composes a directional flow D and a non-directional flow N|

### tests.R
This file contains some simulations for trying out the functions in the file functions.R

### Copyright (c) 2018 Yoav Zemel and Jonathan M. Henshaw
