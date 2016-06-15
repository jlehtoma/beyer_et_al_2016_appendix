

# Solving conservation planning problems with integer linear programming (appendices)

HAWTHORNE L. BEYER (1,3), YANN DUJARDIN (2), MATT WATTS (1) & HUGH P. 
POSSINGHAM(1)  

1 ARC Centre of Excellence for Environmental Decisions, Centre for Biodiversity 
& Conservation Science, University of Queensland, Brisbane QLD 4072, Australia  
2 CSIRO Ecosystem Sciences, Ecosciences Precinct, Dutton Park, QLD 4102, 
Australia  
3 Corresponding author email: hawthorne@spatialecology.com  

Original publication:

Beyer, H. L., Dujardin, Y., Watts, M. E., & Possingham, H. P. (2016). 
Solving conservation planning problems with integer linear programming. 
Ecological Modelling, 328, 14–22. DOI: [10.1016/j.ecolmodel.2016.02.005](http://dx.doi.org/10.1016/j.ecolmodel.2016.02.005)

## A Summary

Section B provides a detailed case study of how to implement the Marxan with 
Zones (Watts et al., 2009) objective functions in an integer linear programming 
(ILP) framework using the linearisation techniques described in the main text.

Section C is a detailed description of implementing ILP optimisation problems 
using R and Gurobi, including documented code for the simulations reported in 
the main text.

## B Linearisation of the Marxan With Zones reserve selection optimisation problem

A variation on the Marxan problem in which each planning unit can take on more 
than two states is called ‘Marxan with Zones’ (Watts et al., 2009). This 
functionality allows several additional types of problems to be addressed 
depending on what the stratification of states represent (management actions, 
zoning types, etc). For example, Klein et al. (2010) use this approach to assign 
marine planning units to one of five zones that regulated the conservation 
status and level of fishing allowed, ranging from complete conservation 
protection to unrestricted commercial fishery.

MarxanWith Zones solves the following optimisation problem:

$$min \sum_{i=1}^{N} \sum_{k=1}^{Z} c_{ik}x_{ik} +  \sum_{k=1}^{Z} \sum_{m=1}^{Z} b_{km} \sum_{i=1}^{N} \sum_{j=1}^{N} x_{ik}x_{jm}v_{ijkm}$$  
  
$$s.t \sum_{k=1}^{Z} x_{ik} = 1, i \in N$$ 
  
$$\sum_{i=1}^{N} \sum_{k=1}^{Z} a_{ij}p_{jk}x_{ij} \geq T_{j},j \in N$$
  
$$\sum_{i=1}^{N} a_{ij}x_{ik} \geq U_{jk},j \in N,k \in K$$

$$x_{i} \in \{0,1\},i \in N$$

where $N$ and $Z$ refer to the number of planning units and zones (actions) 
respectively, $x_{ik}$ is a binary decision variable that is 1 when unit $i$ is 
assigned to zone $k$ and 0 otherwise. In the first part of the objective 
function $c_{ik}$ is the cost of planning unit $i$ when assigned to zone $k$. 
Thus, the costs associated with selecting a planning unit can vary among the 
different zone types and planning units. The second part of the objective 
function provides a way of enforcing or facilitating the assignment of zones 
between pairs of units (often neighbouring units). For example, conflicting 
zones in neighbouring units could be discouraged by assigning a higher penalty 
in the array $v$. It is also possible to specify $v$ to facilitate greater 
clumping of units of a given zone type by assigning a cost to same-type pairs 
that is less than the cost of any different-type pairs. Although $v$ is large (a 
four dimensional array), it may contain many zeros if costs are only assigned to 
neighbouring units, making it more straightforward to linearise this problem 
(because we can ignore any instances of $x_{ik}x_{jm}$ for which 
$v_{ijkm} = 0$). As with the Marxan objective function (Eqn 4), parameter matrix 
$b_{km}$ can be adjusted by the decision maker to control the overall strength of 
the penalties.

The first constraint ensures that each planning unit is assigned to exactly one 
zone. In the second constraint $a_{ij}$ represents the amount or value of 
feature $j$ in unit $i$ and $p_{jk}$ is the proportion of that feature that is 
preserved when the unit is assigned to zone $k$. Thus, summed over all planning 
units and zones, this constraint ensures that minimum representation targets 
($T_{j}$) are met. The final constraint ensures that targets ($U$) for the 
amounts of feature $j$ are represented by specific zone types. For example, if 
the targets refer to species abundances, and the zones refer to multiple land 
use types, then this last constraint could be used to ensure that at least 50% 
of the abundance was represented by fully protected conservation areas.

As with Marxan, the constraints in the Marxan with Zones problem are 
incorporated into the objective function in the form of a shortfall penalty 
function (Watts et al., 2009) that allows solutions to be found even when all 
targets are not met. When all targets are met, however, the shortfall penalty is 
0. For the purposes of implementing Marxan with Zones as an ILP problem we 
choose not to merge the constraints with the objective function as this allows 
for greater flexibility in customising the problem.

The Marxan with Zones problem can be linearised in a similar way as the 
linearisation of the Marxan problem. The first component of the objective 
function is linear with respect to the decision variables $x$ and requires no 
modification. The second component of the objective function contains the 
quadratic expression $x_{ik}x_{jm}$, which can be linearised by replacing it 
with a new decision variable $z_{ijkm}$ and implementing the following 
additional constraints:  
 
$$z_{ijkm} - x_{ik} - x_{jm} \leq -1$$

This constraint ensures that $z_{ijkm} = 1$ when $x_{ik} = x_{jm} = 1$, and 
that $z_{ij} = 0$ when either $x_{ik} = 0$ or $x_{jm} = 0$. Note that this 
constraint differs from the two constraints needed to linearise the Marxan 
objective function because the sign of the quadratic term differs. We can also 
express the objective function using set notation, which makes it clearer that 
the penalty term is only evaluated for neighbours:

$$min \sum_{i=1}^{N} \sum_{k=1}^{Z} c_{ik}x_{ik} +  \sum_{k=1}^{Z} \sum_{m=1}^{Z} b_{km} \sum_{(i,j) \in E} x_{ik}x_{jm}v_{ijkm}$$

where $E$ defines the set of all neighbouring planning units.

In the worst case scenario the linearisation of the Marxan with Zones problem 
could result in approximately $NZ + N^{2}Z^{2}$ decision variables and 
$N^{2}Z^{2}$ constraints in addition to the structural constraints defining the 
targets. However, the linearisation terms can be omitted whenever
$v_{ijkm} = 0$, which typically applies to all non-connected planning units. In 
fact, in most applications the array $v$ is sparse. For example, in the case of 
a hexagon grid a planning unit will have at most 6 neighbours, thereby adding 
less than $6 N Z$ decision variables (note that planning units on the edges of 
the grid have fewer than 6 neighbours) and $6 N Z$ constraints, which is a 
fraction of the size of the worst case scenario. Nevertheless, the dimension of
the problem can still increase rapidly with $N$ and $Z$.

A number of customisations of this problem are possible. It may be appropriate 
to allow no zone assignment to some units in cases where only a subset of unit 
and zone combinations are required to meet targets. This can be achieved by 
changing the first constraint to:

$$\sum_{k=1}^{Z} x_{ik} \leq 1,$$  

thereby relaxing the constraint that forces all planning units to be assigned to 
exactly one zone.

Another problem is forcing some units to be assigned particular zones. Even 
though there is no ‘decision’ to be made for such units it is often necessary to 
include them in optimisation problems because they can affect zone assignments 
to other units via the second component of the optimisation function. For 
example, in a planning problem it may be unrealistic to transform urban areas 
into any other land use type, so they would be fixed as urban. But if a penalty 
has been implemented for placing protected areas adjacent to urban area these 
planning units would need to be included in the problem. Although it is possible 
to force zone assignments by assigning large costs to alternative zones in the 
cost matrix $c_{ik}$, a more direct and definitive approach is to implement a 
constraint $x_{ik} = 1$ for each $i_{k}$ corresponding to existing urban 
planning units. We note that such constraints can also be implemented in Marxan 
with Zones (see Watts et al., 2008).

Similarly, to prevent the assignment of incompatible zones in adjacent planning 
units (e.g. a dredge spoil dump next to a marine protected area) two constraints 
can be implemented:

$$x_{ik} + x_{jm} \leq 1, i \in N, j \in N$$ 
$$x_{im} + x_{jk} \leq 1, i \in N, j \in N$$

where $k$ and $m$ are the incompatible zones.

In Marxan With Zones these relationships will be enforced less explicitly by 
means of the values that are assigned to the cost matrix cik and the penalty 
array $v_{ijkm}$. Indeed, adding additional constraints to enforce these effects 
is detrimental to processing times as it increases the dimension size of the 
problem. The benefit of using the constraint approach, however, is that it is 
definitive. Finding suitable values for $c$ and $v$ to prevent undesirable 
solutions from arising often involves a process of trial and error, particularly 
when parameter $b$ is also being adjusted as this affects the penalties. The 
decision maker must inspect the solutions to determine whether undesirable 
assignments have occurred and this can be difficult to do effectively for very 
large problems. Implemented as a constraint, we are certain that these 
undesirable assignments are precluded, thereby removing the reliance on users to 
appropriately tune the problem formulation.

## C Implementation of linear programming methods to conservation planning

Here, we explain how to solve ILP problems using R and the commercial 
optimisation software Gurobi. This includes code for simulating 10,000 planning 
units with cost and species data for 10 species (Section 2), using this data to 
solve a simple, linear ILP reserve selection problem (Section 3), using this 
data to solve the non-linear Marxan objective function (Section 4), and some 
general advice for solving large ILP problems (Section 5). We recommend that 
readers understand how the simpler reserve selection problem is solved (Section 
3) before moving on to the Marxan problem (Section 4), which is considerably 
more complex.

We note that copying and pasting code from a PDF format can sometimes cause 
issues with how quotation marks are translated into plain text. If you 
experience problems you may wish to first look at whether replacing the single 
and double quotation marks resolves the problem.

### C.1 Alternative linear programming software

Alternatives to Gurobi include the commercial software CPLEX (IBM, Inc., 2009), 
MatLab (MATLAB, 2014) and Mathematica (Mathematica, 2014), and the open source 
software lpSolve (Berkelaar et al., 2005). We have found that processing times 
with lpSolve can be considerably longer than Gurobi. Furthermore, the R 
interface to lpSolve (the `lpSolve` R package) lacks some important features, 
including the ability to accept constraint matrices in the sparse format and 
setting to determine when the solver exits (a time limit or gap statistic). 
Thus, only the smaller ILP problems we evaluate could be implemented though the 
lpSolve R interface and the more complex problems would have to be implemented 
through the lpSolve programming interface (the API) or possibly the command line 
interface. We do not explore this functionality here.

### C.2 Generating spatial datasets

The following R code was used to simulate input data for the case study, 
including 10 raster (matrix) datasets representing the value of each planning 
unit (one raster cell) to each of 10 species of conservation concern, and 
another raster representing the cost of selecting a planning unit for inclusion 
in the reserve system. To run this code it is necessary to first install the 
`RandomFields` library in R.

#### 1. Species data representing the value of each planning unit (raster cell) 
to each species


```r
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

# species array
spp <- array(0, dim = c(ns, nr, nc))
for (i in 1:ns) {
  f <- GaussRF(x = x, y = y, model = model, grid = TRUE,
               param = c(mean, variance, nugget, scale, alpha))
  f[which(f < 0)] <- 0
  spp[i,,] <- f
}

# save the species data as an R data object
# save(spp, file="spp.RData")
# the data can subsequently be reloaded using this command:
# load(file="spp.RData")
```


```r
# Make a RasterStack for plotting purposes
spps <- stack()

for (i in 1:ns) {
  spps <- addLayer(spps, raster(spp[i,,]))
}
names(spps) <- paste0("species", 1:ns)

plot(spps, nc = 3, nr = 4)
```

![](analysis_files/figure-html/plot-species-1.png)<!-- -->

### 2. Simulate the cost associated with selecting each planning unit


```r
# reset the random field parameters (can be adjusted):
mean <- 1000
variance <- 100
nugget <- 1
scale <- 20
alpha <- 1.5
x <- seq(0, 100, length.out = nc)
y <- seq(0, 100, length.out = nr)

cost <- GaussRF(x = x, y = y, model = model, grid = TRUE,
                param = c(mean, variance, nugget, scale, alpha))
# save the species data as an R data object
# save(cost, file="cost.RData")
# the data can subsequently be reloaded using this command:
# load(file="cost.RData")
```



```r
plot(raster(cost), main = "Cost")
```

![](analysis_files/figure-html/plot-cost-1.png)<!-- -->


```r
# Create equal cost surface
cost_equal <- raster(matrix(rep(1, nr * nc), nrow = nr, ncol = nc))
```



## C.3 Solving the basic reserve selection problem using ILP

Here, we demonstrate how to solve the basic reserve selection problem (see Eqn 1 
in the main text) using R and the commercial optimisation software Gurobi (Gurobi 
Optimization, Inc., 2014). This is largely simply a problem of data organisation 
so that the call to the Gurobi interface in R can be made (the `gurobi` function 
in the `gurobi` pacakge in R).

The costs (the values that are being minimised) must be represented in a vector 
in which each item corresponds to a decision variable. Decision variables are 
synonymous with planning units in problems that do not require linearisation. 
Thus, if there are 10,000 planning units, the cost vector will be of length 
10,000 for linear problems. We address the more complex issue of non-linear 
problems below.

The constraints must be organised in a matrix in which the columns correspond to 
decision variables and the rows correspond to individual constraints. If there 
are 10,000 planning units and 15 constraints, the matrix will have dimensions 
10,000 columns by 15 rows for linear problems. This matrix represents the left 
hand side of each constraint equation. Associated with the constraints matrix 
are two vectors representing the right hand side of each constraint equation 
(e.g. the target values) and the equality relationship for each constraint 
equation (<, >,≤,≥, or =). Thus the right hand side is a numerical vector and 
the equality relationship is a character vector, each with a length equal to the 
number of rows of the constraint matrix (15 in the example above).

Note that the order of values in the cost vector and the constrain matrix must 
match. I.e. the 10th item in the cost vector corresponds to a planning unit that 
must also refer to the planning unit in the 10th column of the constraints 
matrix. Often the hardest part of implementing ILP problems is ensuring that the 
order of planning units in these data objects is correct. Note the use of the 
transpose function, `t()`, in the code below to ensure that the order is correct 
when converting the spatial data (a 2D matrix) to a linear vector.

There are also a small number of other parameters that must be set. The 
`modelsense` parameter indicates whether the function should be minimised or 
maximised (`min` or `max` respectively). The `vtype` parameter is used to 
indicate that the decision variables are binary (`B`). Finally, a parameter 
vector can be specified with a range of options that control the algorithm, such 
as indicating when the algorithm should force termination, or the gap that 
should be achieved before termination (details of these settings are available 
in the Gurobi documentation).

This example assumes that the relevant spatial datasets have been loaded in R 
(see the previous section).


```r
# load the gurobi R package (requires that both Gurobi and the Gurobi R package 
# have both been installed): 
require(gurobi)
```

```
## Loading required package: gurobi
```

```
## Loading required package: slam
```

```r
# re-specify the number of rows and columns of the raster if they are not 
# already set in R: 
nr <- 100 # rows 
nc <- 100 # cols 
# re-specify the number of species: 
ns <- 10 # species 
# run the following code without modification: 
N <- nr * nc

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
targetlevel <- 0.05 
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
result <- gurobi(model,params)
```

```
## Warning for adding variables: zero or small (< 1e-13) coefficients, ignored
## Optimize a model with 10 rows, 10000 columns and 50850 nonzeros
## Coefficient statistics:
##   Matrix range    [2e-04, 9e+00]
##   Objective range [1e+00, 1e+00]
##   Bounds range    [1e+00, 1e+00]
##   RHS range       [2e+02, 9e+02]
## Found heuristic solution: objective 729
## Presolve time: 0.05s
## Presolved: 10 rows, 10000 columns, 50850 nonzeros
## Variable types: 0 continuous, 10000 integer (10000 binary)
## Presolved: 10 rows, 10000 columns, 50850 nonzeros
## 
## 
## Root relaxation: objective 2.607932e+02, 51 iterations, 0.03 seconds
## 
##     Nodes    |    Current Node    |     Objective Bounds      |     Work
##  Expl Unexpl |  Obj  Depth IntInf | Incumbent    BestBd   Gap | It/Node Time
## 
##      0     0  260.79323    0   10  729.00000  260.79323  64.2%     -    0s
## H    0     0                     263.0000000  260.79323  0.84%     -    0s
## H    0     0                     262.0000000  260.79323  0.46%     -    0s
## 
## Explored 0 nodes (51 simplex iterations) in 0.27 seconds
## Thread count was 4 (of 4 available processors)
## 
## Optimal solution found (tolerance 5.00e-03)
## Best objective 2.620000000000e+02, best bound 2.610000000000e+02, gap 0.3817%
```

The result object contains several data objects including the objective value 
achieved (`result$objval` = 262) and the vector of decision 
variable values (0 or 1 in our example because the variables are binary).


```r
reserve_network <- matrix(result$x, ncol = nc, nrow = nr)
image(reserve_network, col = c("white", "black"))
```

![](analysis_files/figure-html/plot-optimal-reserve-1.png)<!-- -->

## C.4 Solving the Marxan optimisation problem using ILP

Here, we demonstrate how to solve the basic reserve selection problem (see Eqn 4 
in the main text) using R and the commercial optimisation software Gurobi 
(Gurobi Optimization, Inc., 2014). As with the simple reserve selection problem 
this is largely simply a problem of data organisation. What makes it more 
complicated, however, is that the objective function must be linearised thereby 
adding new decision variables. As discussed in the main text, numerous 
additional decision variables are sometimes needed to linearise a non-linear 
objective function.

The data structures needed are similar to those described in the previous 
section. One difference is that the dimensions of the cost vector and the 
constraint matrix will increase as new decision variables are added. As the 
additional decision variables associated with each linearisation term only 
affect two planning units at a time, one consequence of linearisation is that 
the constraint matrix becomes both sparse and large. As such, it is no longer 
desirable to express the constraint matrix in its full format as this could 
consume a great deal of memory (RAM) even for modestly sized problems. Instead, 
we use the sparse matrix format.

Another important difference is that a new data structure is needed to describe 
the neighbour-relationships. In our example we look at the 4 cardinal neighbours 
for each planning unit (though units on the edge of the raster have less than 4 
neighbours, and this too must be accounted for). There are numerous ways such 
data structures could be implemented. In the case of a regular grid of planning 
units with a fairly consistent number of neighbours we implement the data 
structures as a matrix. In the case of polygons with an irregular number of 
neighbours it would be better to implement this data structure as a list. Thus, 
some adaptation of this code may be needed to apply it to other problems.

### C.4.1 Create data structure defining neighbours

Here, we assume that the planning units are numbered sequentially beginning from 
1. A value of 0 indicated NoData for planning units at the edge of the raster 
with fewer than 4 neighbours. The following code assigns the neighbour ID 
numbers, which is straightforward to calculate arithmetically for a grid of 
planning units. If the planning units were irregular polygons then determining 
neighbours would require the use of GIS functions.


```r
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
```

### C.4.2 Solve the linearised Marxan objective function

In this example the number of decision variables will be one for each planning 
unit plus one for each pair of neighbours as we assume a symmetric cost for not 
selecting each pair of neighbours (the cost for not selecting $A$ if $B$ is 
selected is the same as the cost for not selecting $B$ if $A$ is selected). In 
the case of asymmetric costs, twice as many additional decision variables would 
be required to linearise the objective function.

We also assume an arbitrary ‘cost’ of b units for failing to select neighbouring 
planning units. This value has no significance in our example other than to 
control the degree of aggregation of planning units, and is subjectively 
adjusted by the decision maker. Specifically we assume here that 
$\sf{v_{if}} = 1$ for all neighbours, $\sf{v_{ij}} = 0$ for all non-neighbours, 
and $b$ is adjusted by the decision maker (see Eqn 4 in the main text). If real 
costs have been estimated they could certainly be used here, including different 
values for each pair of planning ($\sf{v_{ij}} = 0$) units if that data is 
available.


```r
# load the gurobi R package (requires that both Gurobi and the Gurobi R package 
# have both been installed): 
require(gurobi) 
# load the Matrix package for sparse matrices 
require(Matrix)
```

```
## Loading required package: Matrix
```

```r
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
```

```
## Warning for adding variables: zero or small (< 1e-13) coefficients, ignored
## Optimize a model with 79210 rows, 49600 columns and 209250 nonzeros
## Coefficient statistics:
##   Matrix range    [2e-04, 9e+00]
##   Objective range [5e+01, 1e+03]
##   Bounds range    [1e+00, 1e+00]
##   RHS range       [1e+03, 5e+03]
## Found heuristic solution: objective 3.41861e+06
## Presolve removed 79200 rows and 39600 columns
## Presolve time: 0.09s
## Presolved: 10 rows, 10000 columns, 50850 nonzeros
## Variable types: 0 continuous, 10000 integer (10000 binary)
## Presolved: 10 rows, 10000 columns, 50850 nonzeros
## 
## 
## Root relaxation: objective 1.572765e+06, 58 iterations, 0.03 seconds
## 
##     Nodes    |    Current Node    |     Objective Bounds      |     Work
##  Expl Unexpl |  Obj  Depth IntInf | Incumbent    BestBd   Gap | It/Node Time
## 
##      0     0 1572765.40    0   10 3418605.08 1572765.40  54.0%     -    0s
## H    0     0                    1574210.0984 1572765.40  0.09%     -    0s
## 
## Explored 0 nodes (58 simplex iterations) in 0.31 seconds
## Thread count was 4 (of 4 available processors)
## 
## Optimal solution found (tolerance 5.00e-03)
## Best objective 1.574210098374e+06, best bound 1.572765404712e+06, gap 0.0918%
```

The result object contains several data objects including the objective value 
achieved (`result$objval`) and the vector of decision variable values (0 or 1 in 
our example because the variables are binary).


```r
reserve_network <- matrix(result$x, ncol = nc, nrow = nr)
image(reserve_network, col = c("white", "black"))
```

![](analysis_files/figure-html/plot-optimal-reserve-2-1.png)<!-- -->

### C.4.3 Calculating the Pareto frontier

Once the basic data structures have been established the Pareto frontier can be 
estimated by iteratively solving the optimisation problem for different values 
of the penalty factor $b$. In our example that simply involves adjusting the cost 
vector. Determining the range of penalty values to explore, and the interval 
between them requires trial and error.


```r
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
  save(solcost, file="solcost.RData")
}
```

Note that the above code saves the sols and solcost results objects within the 
loop. This is precautionary. For some problems the optimisation can take a 
considerable amount of time to run and if the objects are not saved regularly 
they could be lost if there is a problem with the computer.


```r
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
plot(pen, solcost, type="l", xlab="Penalty", ylab="Solution cost")
```
