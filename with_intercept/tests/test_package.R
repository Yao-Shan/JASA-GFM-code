# Check package settings without fitting; run from any working directory.
arg <- grep("^--file=", commandArgs(FALSE), value=TRUE)
stopifnot(length(arg)==1L)
root <- dirname(dirname(normalizePath(gsub("~+~", " ", sub("^--file=", "", arg), fixed=TRUE), mustWork=TRUE)))
api <- new.env(parent=globalenv())
set.seed(37); before <- .Random.seed
source(file.path(root,"R","gfm_baseline.R"), local=api)
stopifnot(identical(before,.Random.seed))
source(file.path(root,"simulations","paper_settings.R"), local=api)
count <- 0L
check <- function(x) {stopifnot(isTRUE(x)); count <<- count+1L}
check(all(c("gfb_fit","gfb_infer","gfb_run_simulation","gfb_generate_truth") %in% ls(api)))
check(all(vapply(c("R","simulations","plots","docs","tests"),
                 function(x) dir.exists(file.path(root,x)), logical(1))))
# Check that ordinary-R defaults do not require local accelerators.
check(!dir.exists(file.path(root,"runtime")))
check(identical(formals(api$gfb_run_simulation)$use_fast_fit,FALSE))
check(identical(formals(api$gfb_run_simulation)$use_fast_inference,FALSE))
for (group in paste0("G",1:5)) for (scenario in c("equal","mixed")) {
  settings <- api$gfb_paper_settings(group, scenario)
  check(settings$alpha_seed==817285L && settings$design_seed==238687L)
  check(settings$n_sim==1000L && settings$K==2L && identical(settings$nodes,1L))
  check(settings$control$tolerance==.05 && !settings$control$backtracking)
  check(settings$control$max_debias_memory_mb==if(settings$n==2000) 6144 else 2048)
  check(settings$n_cores==1L && settings$inference_workers==1L)
  check(identical(settings$lambda_values, if(scenario=="equal") c(.5,.5) else c(.7,-.5)))
  for (model in settings$models) {
    config <- api$gfb_config(settings$n,settings$K,model,control=settings$control)
    check(config$learning_rate==switch(model,gaussian=.1,poisson=.05,1))
    check(config$max_iter==switch(model,gaussian=10000L,poisson=80000L,40000L))
    check(config$tolerance==.05 && config$M_bound_U==20 && config$M_bound_alpha==20)
  }
}
truth <- api$gfb_generate_truth(24L, alpha_seed=817285L)
check(max(abs(colMeans(truth$U_star)))<1e-12)
check(max(abs(crossprod(truth$U_star)/24-diag(2)))<1e-12)
check(abs(sum(truth$alpha_star^2)-1)<1e-12)
check(identical(before,.Random.seed))
check(identical(truth,api$gfb_generate_truth(24L,alpha_seed=817285L)))
check(api$gfb_config(24L)$tolerance==.001) # library default
# Verify the alpha initializer against a direct matrix solve.
B <- outer(1:12,1:12,function(i,j) sin(i+j)+cos(i-j))
expected <- as.vector(solve(12*diag(12)+matrix(1,12,12),rowSums(B)))
check(max(abs(api$.gfb_additive_initial(B)-expected))<1e-12)
if (requireNamespace("knitr",quietly=TRUE)) {
  for (name in c("node_inference.Rmd","all_nodes.Rmd")) {
    pars <- lapply(knitr::knit_params(readLines(file.path(root,"simulations",name))),function(x)x$value)
    check(identical(pars$run_simulation,FALSE))
    check(!any(c("use_fast_fit","use_fast_inference","runtime_dir") %in% names(pars)))
    text <- paste(readLines(file.path(root,"simulations",name)),collapse="\n")
    check(!grepl("use_fast_fit|use_fast_inference|runtime_dir|runtime/",text))
    check(pars$alpha_seed==817285 && pars$n_cores==1 && pars$n_sim==3)
    # Check paths and defaults using setup and truth generation only.
    script <- tempfile(fileext=".R")
    knitr::purl(file.path(root,"simulations",name), output=script, quiet=TRUE)
    env <- new.env(parent=globalenv())
    old <- getwd(); setwd(root)
    tryCatch(sys.source(script,envir=env),finally=setwd(old))
    unlink(script)
    check(is.null(env$simulation))
  }
}
cat(sprintf("PASS: %d package, settings and deterministic formula checks. No network fitting.\n",count))
