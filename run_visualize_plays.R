#! /usr/local/bin/Rscript --vanilla

#
# Genres:
# 164 | Deep House, Paris
# 332 | Indie Electronic, Amsterdam
# 263 | Relax, Paris
# 2   | Soul, Amsterdam
# 

plotPlays <- function(plays, bandwidth) {
  par(mfrow=c(1, 1), mar=c(2, 2, 2, 2))
  plays.density <- density(plays, bw=bandwidth)
  frame <- 1.4
  plot(plays.density, xlim=c(-frame, frame), ylim=c(-frame, frame), xlab=NA, ylab=NA, main=NA)
}

load('plays.RData')

for (i in 1:length(plays)) {
  filename <- paste0("figures/genre_", i, ".pdf")
  pdf(file=filename, width=5, height=5)
  plotPlays(plays[[i]], 200)
  dev.off()
}

all.plays <- circular(unlist(plays), template="clock24")

filename <- "figures/genres_joined.pdf"
pdf(file=filename, width=5, height=5)
plotPlays(all.plays, 200)
dev.off()
