# Define paper settings without fitting models or writing output.
gfb_paper_settings <- function(group = c("G1", "G2", "G3", "G4", "G5"),
                               scenario = c("equal", "mixed"), n_cores = 1L) {
  group <- match.arg(group); scenario <- match.arg(scenario)
  if (length(n_cores) != 1L || !is.numeric(n_cores) || !is.finite(n_cores) ||
      n_cores != floor(n_cores) || n_cores < 1 || n_cores > 7) stop("n_cores must be 1..7.")
  n <- c(G1=1000L, G2=1000L, G3=1000L, G4=2000L, G5=2000L)[[group]]
  base <- c(G1=2000L, G2=12000L, G3=32000L, G4=32000L, G5=2000L)[[group]]
  list(n=n, K=2L, models=c("gaussian","logit","probit","poisson"), n_sim=1000L,
       nodes=1L, lambda_values=if(scenario=="equal") c(.5,.5) else c(.7,-.5),
       sigma=1, design_seed=238687L, alpha_seed=817285L, replication_seed_base=base,
       alpha_min=-2, alpha_max=0, alpha_normalization="l2",
       control=list(tolerance=.05, backtracking=FALSE, M_bound_U=20, M_bound_alpha=20,
                    M_bound_Lambda=20, probability_eps=.001,
                    max_debias_memory_mb=if(n==2000L) 6144 else 2048),
       debias=TRUE, debias_method="pdf_bp_bm", save_fits=FALSE,
       n_cores=as.integer(n_cores), inference_workers=1L,
       worker_backend=if(n_cores==1L) "in_process" else "psock")
}
