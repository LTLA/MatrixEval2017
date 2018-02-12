se <- function(x) { sd(x)/sqrt(length(x)) }

timeExprs <- function(..., times=10) {
    output <- numeric(times)
    expr <- substitute(system.time(...))
    for (it in seq_len(times)) {
        gc() # Trigger garbage collection to avoid it being in the timing.
        output[it] <- eval(expr)["elapsed"]
    }
    return(mean(output)*1e3) # Get to milliseconds
}

writeToFile <- function(..., timings, file, overwrite) {
    write.table(data.frame(..., Time=mean(timings), SE=se(timings)), file=file, 
                append=!overwrite, col.names=overwrite, sep="\t", quote=FALSE, row.names=FALSE)
    return(invisible(NULL))
}
