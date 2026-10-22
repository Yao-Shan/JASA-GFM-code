# Check fitting and inference on four n=24 networks.
arg <- grep("^--file=",commandArgs(FALSE),value=TRUE)
root <- dirname(dirname(normalizePath(gsub("~+~", " ", sub("^--file=","",arg), fixed=TRUE),mustWork=TRUE)))
api <- new.env(parent=globalenv()); source(file.path(root,"R","gfm_baseline.R"),local=api)
truth <- api$gfb_generate_truth(24L, alpha_seed=817285L)
rows <- list()
for (model in c("gaussian","logit","probit","poisson")) {
 A <- api$gfb_generate_network(truth,model=model,seed=2001L)
 stopifnot(identical(A,api$gfb_generate_network(truth,model=model,seed=2001L)))
 fit <- api$gfb_fit(A,K=2L,model=model,control=list(tolerance=.05))
 stopifnot(inherits(fit,"gfb_fit"), all(is.finite(fit$U_hat)), all(is.finite(fit$alpha_hat)),
           max(abs(colMeans(fit$U_hat)))<1e-10,
           max(abs(crossprod(fit$U_hat)/24-diag(2)))<1e-10)
 frozen <- serialize(fit,NULL)
 inference <- tryCatch(api$gfb_infer(fit,U_star=truth$U_star,alpha_star=truth$alpha_star,
                      Lambda_star=truth$Lambda_star,nodes=1L,debias=TRUE),error=function(e)e)
 stopifnot(identical(frozen,serialize(fit,NULL))) # inference leaves the fit unchanged
 if(inherits(inference,"error")) {
   status <- conditionMessage(inference)
   # Keep failed inference visible without retries.
 } else {
   stopifnot(nrow(inference$nodes)==1L,
             all(c("z_u1","z_u2","z_alpha_joint") %in% names(inference$nodes)))
   status <- inference$lambda$status
 }
 rows[[model]] <- data.frame(model=model,converged=fit$diagnostics$converged,
                            n_iter=fit$diagnostics$n_iter,inference_status=status)
}
print(do.call(rbind,rows),row.names=FALSE)
cat("PASS: four deterministic-seed small fits, network repeatability and fit reuse checked; no calibration claim.\n")
