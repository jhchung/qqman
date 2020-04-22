#' Creates a Q-Q plot
#' 
#' Creates a quantile-quantile plot from p-values from a GWAS study. 
#' Sub threshold points are sampled to reduce plotting time. Currently uses reservior sampling
#' 
#' @param pvector A numeric vector of p-values.
#' @param p_thresh P-value threshold for sampling. P-values > \code{p_thresh} will be sampled.
#' @param n_sample Integer indicating number of sub-threshold points to plot
#' @param ... Other arguments passed to \code{plot()}
#' 
#' @return A Q-Q plot.
#' 
#' @keywords visualization qq qqplot
#' 
#' @importFrom stats ppoints
#' @import utils
#' @import graphics
#' @import wrswoR
#' 
#' @examples
#' qq_sample(gwasResults$P, p_thresh = 1e-5, n_sample = 10000)
#' 
#' @export
qq_sample <- function (pvector, p_thresh = 1e-5, n_sample = Inf, ...) {
  
  require(wrswoR)
  if (!is.numeric(pvector)) 
    stop("Input must be numeric.")
  pvector <- pvector[!is.na(pvector) & !is.nan(pvector) & !is.null(pvector) & 
                       is.finite(pvector) & pvector < 1 & pvector > 0]
  
  df <- data.frame(o = -log10(sort(pvector, decreasing = FALSE)),
                   e = -log10(ppoints(length(pvector))))
  
  if (n_sample != Inf){
    top_obs <- df[df$o >= -log10(p_thresh), ]
    other_obs <- df[df$o < -log10(p_thresh), ]
    
    # Weighted resampling
    sample_idx <- wrswoR::sample_int_expj(n = nrow(other_obs), 
                                          size = n_sample, 
                                          prob = other_obs$o)
    other_obs <- other_obs[sample_idx, ]
    df <- rbind(top_obs, other_obs)
    
    message("Keep ", nrow(df), " variants")
  }
  
  
  def_args <- list(pch = 20, 
                   xlim = c(0, max(df$e)), 
                   ylim = c(0, max(df$o)), 
                   xlab = expression(Expected ~ ~-log[10](italic(p))), 
                   ylab = expression(Observed ~ ~-log[10](italic(p))))
  dotargs <- list(...)
  tryCatch(do.call("plot", c(list(x = df$e, y = df$o), 
                             def_args[!names(def_args) %in% 
                                        names(dotargs)], dotargs)), warn = stop)
  abline(0, 1, col = "red")
}
