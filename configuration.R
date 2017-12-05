
source('helpers.R')

burnin <- 10^4
iterations <- 5*10^3

LOG <- FALSE

# List of scenarios to be simulated
scenario.definitions <- list(
  list(
    k=1,
    mu=0,
    kappa=0,
    weight=1,
    name="uniform"
  ),
  list(
    k=1,
    mu=0,
    kappa=10,
    weight=1,
    name="single_component"
  ),
  list(
    k=2,
    mu=c(-pi*(0.33/2), pi*(0.33/2)),
    kappa=c(10, 10),
    weight=c(0.5, 0.5),
    name="two_components"
  ),
  list(
    k=2,
    mu=c(0, pi),
    kappa=c(10, 10),
    weight=c(0.5, 0.5),
    name="two_components_opposite"
  ),
  list(
    k=3,
    mu=c(-pi*0.33, 0, pi*0.33),
    kappa=rep(10, 3),
    weight=rep(1/3, 3),
    name="three_components"
  )
)
