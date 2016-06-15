
# 0. Setup ----------------------------------------------------------------

# generates simulated spatial (raster) datasets if not installed, the
# RandomFields library can be installed using this command:
# install.packages("RandomFields")
library(RandomFields)
library(raster)
set.seed(42)

# 1. Generate spp data ----------------------------------------------------

# user defined parameters:
# set the number of rows and columns of the raster:
nr <- 100 # rows
nc <- 100 # cols

# set the number of species:
ns <- 10 # species

# run the following code without modification:
N <- nr * nc
x <- seq(0, 100, length.out = nc)
y <- seq(0, 100, length.out = nr)

# the random field parameters are set here and can be adjusted
model <- "stable"
mean <- 0
variance <- 4
nugget <- 1
scale <- 20
alpha <- 1

# Dynamically create species cache file
spp_cache_file <- file.path("cache", paste0("spp_", nr, "x", nc, ".RData"))

# If species data for this dimension (nr * nc) exists, load that instead of
# recreating it.
if (file.exists(spp_cache_file)) {
  load(spp_cache_file)
} else {
  # species array
  spp <- array(0, dim = c(ns, nr, nc))
  for (i in 1:ns) {
    message("Creating species ", i)
    f <- GaussRF(x = x, y = y, model = model, grid = TRUE,
                 param = c(mean, variance, nugget, scale, alpha))
    f[which(f < 0)] <- 0
    spp[i,,] <- f
  }
  # save the species data as an R data object
  save(spp, file = spp_cache_file)
}

# Make a RasterStack for plotting purposes
spps <- stack()

for (i in 1:ns) {
  spps <- addLayer(spps, raster(spp[i,,]))
}
names(spps) <- paste0("species", 1:ns)

plot(spps, nc = 3, nr = 4)

# 2. Generate cost data ---------------------------------------------------

# reset the random field parameters (can be adjusted):
mean <- 1000
variance <- 100
nugget <- 1
scale <- 20
alpha <- 1.5
x <- seq(0, 100, length.out = nc)
y <- seq(0, 100, length.out = nr)

# Dynamically create cost cache file
cost_cache_file <- file.path("cache", paste0("cost_", nr, "x", nc, ".RData"))

# If cost data for this dimension (nr * nc) exists, load that instead of
# recreating it.
if (file.exists(cost_cache_file)) {
  load(cost_cache_file)
} else {
  cost <- GaussRF(x = x, y = y, model = model, grid = TRUE,
                  param = c(mean, variance, nugget, scale, alpha))
  # Save the species data as an R data object
  save(cost, file = cost_cache_file)
}

plot(raster(cost), main = "Cost")

# Create equal cost surface. This is a cheap operation, no need for caching.
cost_equal <- raster(matrix(rep(1, nr * nc), nrow = nr, ncol = nc))

plot(cost_equal, main = "Cost equal")

# 3. Optimization ---------------------------------------------------------

# Load the gurobi R package (requires that both Gurobi and the Gurobi R package
# have both been installed):
require(gurobi)

solve_for_target <- function(targetlevel) {

  # solve the basic reserve selection problem
  # create the constraint matrix with appropriate dimensions
  constr <- matrix(0, nrow = ns, ncol = N)

  # in our example we have one contraint for each of 10 species so we populate the
  # matrix from our species rasters:
  for (i in 1:ns) {
    constr[i,] <- as.vector(t(spp[i,,]))
  }

  # create a vector of default equalities (these can be adjusted later, though we
  # do not need to adjust them in this example):
  sense <- rep(">=", ns)

  # set the targets (the right hand side of the constraint equations); here, we
  # assume that 25% of the total value of each species raster must be met or
  # exceeded:
  rhs <- rep(0, ns)
  for (i in 1:ns) {
    rhs[i] <- targetlevel * sum(spp[i,,])
  }

  # set up Gurobi model object in R
  model <- list()

  # set this as a minimisation problem:
  model$modelsense <- "min"
  # set all decision variables as binary:
  model$vtype <- "B"
  # vector of state values that are being minimised (costs in our example):
  model$obj <- as.vector(t(cost_equal))
  # assign the constraints matrix object, and the right hand side and sense vectors:
  model$A <- constr
  model$rhs <- rhs
  model$sense <- sense
  # set the parameters that control the algorithm (the algorithm stops when a gap
  # of 0.5% is achieved in this example):
  params <- list(Presolve = 2, MIPGap = 0.005)
  # solve the problem
  result <- gurobi(model, params)
  return(result)
}

solutions <- list()
target_levels <- c(0.02, 0.05, 0.1, 0.25, 0.5, 0.8)

for (targetlevel in target_levels) {
  solutions[[as.character(targetlevel)]] <- solve_for_target(targetlevel)
}

reserve_networks <- stack()

for (solution in solutions) {
  reserve_network <- matrix(solution$x, ncol = nc, nrow = nr)
  reserve_networks <- addLayer(reserve_networks, raster(reserve_network))
}
names(reserve_networks) <- paste0("targetlevel_", target_levels)
sum_raster <- sum(reserve_networks)
z_colors <- c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF", "#E0F3F8", "#91BFDB")
plot(sum_raster, col = rev(z_colors))

# 4. Connectivity ---------------------------------------------------------

# create the planning unit (PU) neighbour matrix
neighbours <- matrix(0, nrow = N, ncol = 4)
# all middle PUs have 4 neighbours:
for (i in 2:(nc - 1)) {
  for (j in 2:(nr - 1)) {
    id <- (j - 1) * nc + i
    neighbours[id, 1] <- id - nc
    neighbours[id, 2] <- id + nc
    neighbours[id, 3] <- id - 1
    neighbours[id, 4] <- id + 1
  }
}

# first and last row, excluding corners
for (i in 2:(nc - 1)) {
  id <- i
  neighbours[id, 1] <- id + nc
  neighbours[id, 2] <- id - 1
  neighbours[id, 3] <- id + 1
  id <- (nr - 1) * nc + i
  neighbours[id, 1] <- id - nc
  neighbours[id, 2] <- id - 1
  neighbours[id, 3] <- id + 1
}

# first and last columns, excluding corners
for (i in 2:(nr - 1)) {
  id <- (i - 1) * nc + 1
  neighbours[id, 1] <- id + nc
  neighbours[id, 2] <- id - nc
  neighbours[id, 3] <- id + 1
  id <- i * nc
  neighbours[id, 1] <- id + nc
  neighbours[id, 2] <- id - nc
  neighbours[id, 3] <- id - 1
}

# corners
id <- 1
neighbours[id, 1] <- id + nc
neighbours[id, 2] <- id + 1
id <- nc
neighbours[id, 1] <- id + nc
neighbours[id, 2] <- id - 1
id <- nc * (nr - 1) + 1
neighbours[id, 1] <- id - nc
neighbours[id, 2] <- id + 1
id <- nc * nr
neighbours[id, 1] <- id - nc
neighbours[id, 2] <- id - 1

# 5. Optimize with connectivity -------------------------------------------

# load the gurobi R package (requires that both Gurobi and the Gurobi R package
# have both been installed):
require(gurobi)
# load the Matrix package for sparse matrices
require(Matrix)

# re-specify the number of rows and columns of the raster if they are not
# already set in R:
nr <- 100 # rows
nc <- 100 # cols
# re-specify the number of species:
ns <- 10 # species
# run the following code without modification:
N <- nr*nc

# set a default penalty factor here:
b <- 50

# determine the number of additional decision variables (M) required:
M <- length(which(neighbours > 0))
# construct the cost vector from the cost raster and with an additional M
# replicates of b for the additional decision variables:
objvals <- c(as.vector(t(cost)), rep(b, M))

# to construct the sparse constraints matrix we need three empty vectors that
# are subsequently populated; it is often difficult to calculate the required
# size a priori, so we create large vectors and adjust their size later;
# # mc, mr and mz refer to the column position, row position and value of each
# cell in the constraint matrix:
mc <- rep(0, 1000000)
mr <- rep(0, 1000000)
mz <- rep(0, 1000000)

# similarly, we create large vectors for the equality relationship and right
# hand side values, and adjust their length later; we set the default sense
# to <= so that we do not have to specify it in the iterative loop later:
sense <- rep("<=", 100000)
rhs <- rep(0, 100000)

# the first 10 rows of the constraint matrix correspond to the species targets,
# which we populate as follows:
mc[1:(N * ns)] <- rep(c(1:N), ns)
mr[1:(N*ns)] <- rep(c(1:ns), each = N)
for (i in 1:ns) {
  pos <- (i - 1) * N + 1
  mz[pos:(pos + N - 1)] <- as.vector(t(spp[i,,]))
}
# we need to adjust the sense matrix for these species targets:
sense[1:ns] <- ">="
# we populate the right hand side values as we did in the previous section:
# set the targets (the right hand side of the constraint equations); here, we
# assume that 25% of the total value of each species raster must be met or
# exceeded:
targetlevel <- 0.25
target <- rep(0, ns)
for (i in 1:ns) {
  target[i] <- targetlevel * sum(spp[i,,])
}
rhs[1:ns] <- target

# now work through all the neighbours, adding two constraints for each one we
# use cr to indicate the current row of the constraint matrix, idx refers to the
# current index of the mc, mr and mz vectors, and z refers to the index of the
# new decision variable (column).
z <- N + 1
idx <- N * ns + 1
cr <- ns + 1

for (i in 1:N) {
  ids <- neighbours[i, which(neighbours[i,] > 0)]
  for (j in 1:length(ids)) {
    # z_ij - x_i <= 0
    mc[idx] <- i
    mr[idx] <- cr
    mz[idx] <- -1
    idx <- idx + 1
    mc[idx] <- z
    mr[idx] <- cr
    mz[idx] <- 1
    idx <- idx + 1
    cr <- cr + 1

    # z_ij - x_j <= 0
    mc[idx] <- ids[j]
    mr[idx] <- cr
    mz[idx] <- -1
    idx <- idx + 1
    mc[idx] <- z
    mr[idx] <- cr
    mz[idx] <- 1
    idx <- idx + 1
    cr <- cr + 1

    z <- z + 1
  }
}

# finally, resize the vectors to the right size:
idx <- idx - 1
mc <- mc[1:idx]
mr <- mr[1:idx]
mz <- mz[1:idx]
cr <- cr - 1
sense <- sense[1:cr]
rhs <- rhs[1:cr]

# construct the sparse matrix:
constr <- sparseMatrix(i = mr, j = mc, x = mz)

# set up Gurobi model object in R
model <- list()
# set this as a minimisation problem:
model$modelsense <- "min"
# set all decision variables as binary:
model$vtype <- "B"
# vector of state values that are being minimised (costs in our example):
model$obj <- objvals
# assign the constraints matrix object, and the right hand side and sense
# vectors:
model$A <- constr
model$rhs <- rhs
model$sense <- sense
# set the parameters that control the algorithm (the algorithm stops when a gap
# of 0.5% is achieved in this example):
params <- list(Presolve = 2,MIPGap = 0.005)
# solve the problem
result <- gurobi(model, params)

# 6. Pareto frontier ------------------------------------------------------

# create a vector of penalty values to use:
pen <- seq(0, 175, 25)
#pen <- seq(0, 50, 25)
# create some empty data objects to store the objective values and solutions for
# visualisation:
sols <- list()
solcost <- rep(0, length(pen))

# loop through the penalties, optimising for each one:
for (p in 1:length(pen)) {
  objvals <- c(as.vector(t(cost)), -rep(pen[p], M))
  model$obj <- objvals
  result <- gurobi(model, params)
  sols[[p]] <- result$x
  solcost[p] <- result$objval
  save(sols, file = "sols.RData")
  save(solcost, file = "solcost.RData")
}

# Make a RasterStack for plotting purposes
solutions <- stack()

sols_mtrxs <- list()

for (i in 1:length(sols)) {
  sols_mtrxs[[i]] <- matrix(sols[[i]], ncol = nc, nrow = nr)
}

for (i in 1:length(sols)) {
  solutions <- addLayer(solutions, raster(sols_mtrxs[[i]]))
}
solutions <- addLayer(solutions, raster(cost))
names(solutions) <- c(paste0("solution", 1:length(sols)), "cost")

plot(solutions, nc = 3, nr = 3)
plot(rev(pen), solcost, type = "l", xlab = "Penalty for aggregation",
     ylab = "Solution cost")

