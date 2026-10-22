# Test plotting with small deterministic fixtures.
arg <- grep("^--file=", commandArgs(FALSE), value=TRUE)
root <- dirname(dirname(normalizePath(gsub("~+~", " ", sub("^--file=", "", arg), fixed=TRUE), mustWork=TRUE)))
api <- new.env(parent=globalenv())
set.seed(37); before <- .Random.seed
source(file.path(root, "plots", "plot_individual.R"), local=api)
stopifnot(identical(before, .Random.seed))
work <- tempfile("gfb-plot-test-"); dir.create(work)
fixture <- function(K) {
  B <- 10L; m <- K*(K+1L)/2L
  z <- matrix(sin(seq_len(B*m)), B, m)
  y <- matrix(cos(seq_len(B*(K+1L))), B, K+1L)
  values <- seq(.3, .3+.1*(K-1L), length.out=K)
  true <- c(values/sqrt(2), rep(0,m-K))
  t1 <- data.frame(model="gaussian",replicate=seq_len(B),fit_converged=TRUE,
                   status="ok",node=1L,chi_sq_Y=rowSums(y^2),z_alpha_joint=y[,K+1L],alpha_marginal_z=y[,1L])
  for (j in seq_len(K)) t1[[paste0("z_u",j)]] <- y[,j]
  t2 <- data.frame(model="gaussian",replicate=seq_len(B),fit_converged=TRUE,status="ok",chi_sq_debiased=rowSums(z^2))
  for (j in seq_len(m)) {
    t2[[paste0("z_debiased_",j)]] <- z[,j]
    t2[[paste0("raw_error_",j)]] <- seq_len(B)/1000 + .01*j
    t2[[paste0("debiased_error_",j)]] <- seq_len(B)/2000 - .002*j
    t2[[paste0("lambda_v_hat_",j)]] <- true[j]+t2[[paste0("raw_error_",j)]]
    t2[[paste0("lambda_v_debiased_",j)]] <- true[j]+t2[[paste0("debiased_error_",j)]]
  }
  # Exclude failed or incomplete pairs for both estimators.
  t2$status[9] <- "failed"; t2$lambda_v_hat_1[10] <- NA_real_
  list(schema_version="gfm-baseline-v1",
       settings=list(n=24L,K=K,n_sim=B,lambda_values=values),truth=list(Lambda_star=diag(values,nrow=K)),
       node_summary=t1,lambda_summary=t2)
}
for (K in c(1L,3L)) {
  sim <- fixture(K); path <- file.path(work,paste0("native",K,".rds")); saveRDS(sim,path)
  hash <- unname(tools::md5sum(path))
  out <- api$gfb_plot_saved(path,file.path(work,paste0("figures",K)),font="serif",histogram_breaks=8)
  m <- K*(K+1L)/2L
  stopifnot(nrow(out$figures)==K+3L+m+1L+2L*m, identical(hash,unname(tools::md5sum(path))))
  stopifnot(!any(grepl("z_raw|chi_sq_raw",out$figures$statistic)))
  h <- read.csv(file.path(out$output_dir,"histogram_values.csv"))
  c <- api$.gfb_plot_coordinates(K)
  for (j in seq_len(nrow(c))) {
    d <- h[h$row==c$row[j] & h$col==c$col[j] & h$scale=="matrix_entry",]
    expected <- sim$lambda_summary[[paste0("lambda_v_hat_",j)]][d$replicate] * c$multiplier[j]
    stopifnot(max(abs(d$raw-expected))<1e-12,!any(d$replicate==9L))
    e <- h[h$row==c$row[j] & h$col==c$col[j] & h$scale=="matrix_entry_error",]
    stopifnot(identical(d$replicate,e$replicate),max(abs(d$raw-e$raw-d$truth))<1e-12,
              max(abs(d$debiased-e$debiased-d$truth))<1e-12)
  }
  stopifnot(!any(h$replicate==10L & h$row==1L & h$col==1L))
  bad <- sim; bad$lambda_summary$raw_error_1[1] <- 10
  saveRDS(bad,file.path(work,"bad.rds"))
  stopifnot(inherits(try(api$gfb_read_plot_data(file.path(work,"bad.rds")),silent=TRUE),"try-error"))
  bad <- sim; bad$node_summary$chi_sq_Y[1] <- 100
  saveRDS(bad,file.path(work,"bad.rds"))
  stopifnot(inherits(try(api$gfb_read_plot_data(file.path(work,"bad.rds")),silent=TRUE),"try-error"))
  stopifnot(inherits(try(api$gfb_plot_saved(path,out$output_dir,font="serif"),silent=TRUE),"try-error"))
}
# Check the historical CSV schema using the same fixture values.
sim <- fixture(3L); t1 <- sim$node_summary; t2 <- sim$lambda_summary
for (stage in c("t1","t2")) {
  d <- get(stage); d$n <- 24L; d$K <- 3L; d$scenario <- "lambda_1"
  d$lambda_values <- paste(sim$settings$lambda_values,collapse=",")
  d$converged <- d$fit_converged; d$inference_status <- d$status
  d$status <- ifelse(d$fit_converged,"success","failed")
  if (stage=="t2") for (j in 1:6) {
    for (prefix in c("lambda_v_hat","lambda_v_debiased","raw_error","debiased_error"))
      names(d)[names(d)==paste0(prefix,"_",j)] <- paste0(prefix,"_v",j)
    names(d)[names(d)==paste0("z_debiased_",j)] <- paste0("z_lambda_v",j,"_debiased")
  }
  assign(stage,d)
}
csvdir <- file.path(work,"csv"); dir.create(csvdir)
write.csv(t1,file.path(csvdir,"theorem1_records.csv"),row.names=FALSE)
write.csv(t2,file.path(csvdir,"theorem2_records.csv"),row.names=FALSE)
loaded <- api$gfb_read_plot_data(csvdir)
stopifnot(nrow(loaded$t1)==10L,nrow(loaded$t2)==10L, sum(loaded$t2$valid)==9L)
t1 <- rbind(t1,t1[1,]);write.csv(t1,file.path(csvdir,"theorem1_records.csv"),row.names=FALSE)
stopifnot(inherits(try(api$gfb_read_plot_data(csvdir),silent=TRUE),"try-error"))
stopifnot(identical(before,.Random.seed))
cat("PASS: native/historical plotting schemas, K=1/K=3, paired histogram scales, truth and chi-square validation, failure filters and overwrite protection.\n")
cat("Temporary fixture PDFs:",work,"\n")
