#!/usr/bin/env Rscript

# contributed by Bo Li

argv = commandArgs(TRUE)
if (length(argv) != 2) {
   cat("Usage: Rscript plotPRcurves.R input.table output.pdf\n")
   q(status = 1)
}

library(PRROC)

plotPR = function(id, progs, data) {
	idx = data[,1] == progs[id]
	if (id == 1) {
		plot(data[idx,2], data[idx,3], type = 'l', lwd = 3, col = id, lty = id, xlim = c(0, 1), ylim = c(0, 1), xlab = "Recall", ylab = "Precision")
	} else {
		par(new = T)
		plot(data[idx,2], data[idx,3], type = 'l', lwd = 3, col = id, lty = id, xlim = c(0, 1), ylim = c(0, 1), xlab = "", ylab = "")
	}
}

data = read.table(argv[1])
progs = levels(data[,1])

pdf(argv[2])
par(mar = c(5, 4, 8, 2) + 0.1, xpd = TRUE)
a = lapply(1:length(progs), plotPR, progs, data)
legend(x = -0.06, y = 1.3, legend = progs, ncol = 3, lwd = 3, col = 1:length(progs), lty = 1:length(progs), cex = 0.54)
dev.off()
