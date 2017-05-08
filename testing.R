# Load data
data <- read.csv("WoodyClimateHemi.12May.csv")
data$herb <- data$woodiness=="H"

plot.herb <- function(var, data, by, ...){
    bins <- with(data, seq(min(data[,var]),max(data[,var]),by=by))
    x <- with(data, table(cut(data[,var],bins), herb))[-1,]
    weight <- as.numeric(scale(rowSums(x)[-1]))
    weight <- (weight - min(weight) + 1) * 2
    x <- apply(x, 1, function(y) y[2]/sum(y))
    plot(x ~ bins[-1:-2], type="n", axes=FALSE, ylab="Proportion herbaceous", ylim=c(0,1), lwd=weight, ...)
    for(i in seq_along(x))
        lines(bins[(i+2):(i+3)], x[i:(i+1)], lwd=weight[i])
    axis(1)
    axis(2)
}

data <- data[,c("species","herb","tmin.lo")]
data <- na.omit(data)

plot.herb("tmin.lo", data, by=5)
