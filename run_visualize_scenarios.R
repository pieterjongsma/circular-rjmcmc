#! /usr/local/bin/Rscript --vanilla

source('configuration.R')
source('visualization.R')

for (definition in scenario.definitions) {
  file <- paste0("figures/scenario_", definition$name, ".pdf")
  pdf(file, width=3, height=5)
  par(mar=c(0, 0, 0, 0))
  plot(circular(NA, rotation="clock", zero=pi/2), xlim=c(-1.2, 1.2), ylim=c(-2, 2))
  CurveForParameters(definition, add=TRUE)
  dev.off()
}
