a <- matrix(rnorm(1000000), nrow=10000)

# (1) Totally unvectorized version of taking the rowSums.
# Shown here to mimic the performance of a completely 
# naive implementation of an arbitrary operation.
system.time({
    out <- numeric(nrow(a))
    for (i in seq_len(nrow(a))) {
        total <- 0
        for (j in seq_len(ncol(a))) { 
            total <- total + a[i,j]
        }
        out[i] <- total
    }
})

# (2) Semi-vectorized version of taking the rowSums.
# Here, the operation inside the loop is vectorized.
system.time({
    out2 <- numeric(nrow(a))
    for (i in seq_len(nrow(a))) {
        out[i] <- sum(a[i,])
    }
})

# (3) Completely vectorized.
system.time(out3 <- rowSums(a))

