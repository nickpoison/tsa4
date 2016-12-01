grid <-
function (nx = NULL, ny = nx, col = gray(.9), lty = 1, 
    lwd = par("lwd"), equilogs = TRUE) 
{
    if (is.null(nx) || (!is.na(nx) && nx >= 1)) {
        log <- par("xlog")
        if (is.null(nx)) {
            ax <- par("xaxp")
            if (log && equilogs && ax[3L] > 0) 
                ax[3L] <- 1
            at <- axTicks(1, axp = ax, log = log)
        }
        else {
            U <- par("usr")
            at <- seq.int(U[1L], U[2L], length.out = nx + 1)
            at <- (if (log) 
                10^at
            else at)[-c(1, nx + 1)]
        }
        abline(v = at, col = col, lty = lty, lwd = lwd)
    }
    if (is.null(ny) || (!is.na(ny) && ny >= 1)) {
        log <- par("ylog")
        if (is.null(ny)) {
            ax <- par("yaxp")
            if (log && equilogs && ax[3L] > 0) 
                ax[3L] <- 1
            at <- axTicks(2, axp = ax, log = log)
        }
        else {
            U <- par("usr")
            at <- seq.int(U[3L], U[4L], length.out = ny + 1)
            at <- (if (log) 
                10^at
            else at)[-c(1, ny + 1)]
        }
        abline(h = at, col = col, lty = lty, lwd = lwd)
    }
}
