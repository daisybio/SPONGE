
sponge_batch_correction <- function (expression, batches= sample(1:3, nrow(expression), replace=T), mean.only = FALSE)
{
  library(COMBAT)
  library(sva)

  expression_t <- t(expression)

  if(mean.only)
  {
      # non-parametric adjustment, mean-only version
      combat_edata = ComBat(dat=expression_t, batch=batches, mod=NULL, par.prior=FALSE, mean.only=TRUE)
  }
  else{
      # parametric adjustment
      combat_edata = ComBat(dat=expression_t, batch=batches, mod=NULL, par.prior=TRUE, prior.plots=FALSE)
  }


  return(t(combat_edata))
}

