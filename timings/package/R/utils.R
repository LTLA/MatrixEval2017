se <- function(x) { sd(x)/sqrt(length(x)) }

timeExprs <- function(..., times=10) {
    output <- numeric(times)
    expr <- substitute(system.time(...))
    for (it in seq_len(times)) {
        output[it] <- eval(expr)["elapsed"]
    }
    return(mean(output)*1e3) # Get to milliseconds
}
