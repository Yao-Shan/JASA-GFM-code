# Computational convention and proof limitations

The fitting, utility and inference R modules preserve the calculations in the
local study sources, with shorter source comments. The simulation driver also
defaults its two acceleration switches to `FALSE`. The notebooks remove
acceleration options and runtime paths; this upload package uses ordinary R.
Local runtime sources are retained outside the package.
`core_source_manifest.csv` records source provenance and the packaging changes.
Portable notebook paths and trial defaults differ from the production studies;
report-study settings are supplied explicitly.

## Current default: pdf_bp_bm

With y_i=(u_i,alpha_i), w_i=(Lambda u_i,1) and
V_i=[sum_{j!=i} E(-ell''_ij) w_j w_j']^{-1}, the retained definition is

    BP_yi = -1/2 sum_{j!=i} V_j (u_i,0)
    DeltaY = D D' BP_Y / n^2

It uses the complete (K+1)-dimensional inverse; the alpha component of BP is not
forced to zero. The implementation is not replaced by full constraint curvature.
For ell=log p, define kappa_ij=E(ell'_ij ell''_ij) and s_ij=w_i' V_j w_i:

    BM_yi = sum_{j!=i} kappa_ij (s_ij+s_ji) w_j
    BM_Lambda = sum_i sum_{j!=i} kappa_ij s_ij d_ij

The ordered sum for BM_Lambda has no extra half. Kappa is zero for Gaussian,
Logit and Poisson, generally nonzero for Probit. `correction_score_ell2` and
`correction_BM` are aliases and must not be summed twice.

Let R=P(P I_YY P)^+P, Ieff=I_LL-I_LY R I_YL and
T(B)=Ieff^{-1}(B_Lambda-I_LY R B_Y). The code uses B1 with negative one-half
expected-third-derivative contractions, and adds

    correction_ell3 = T(B1)
    correction_BM = -T(BM)
    correction_constraint = + I_LL^{-1} I_LY D D' BP_Y / n^2

These signs are explicit computational conventions; they do not literally copy
all signs in the draft proof. The exact finite-sample projection is
P=I-D(D'D)^{-1}D', not the unqualified plug-in I-DD'/n. Optional
`constrained_second_order` remains a separate method, not an additional correction.

## Scope of validation

Computational regression checks do not establish the increasing-dimension theorem.
The draft's treatment of BP's column space and the identification of a deterministic
penalty term with a random Hessian term require further theoretical review. This
publication package preserves the computations and reports their convention;
it does not claim to have repaired those proof steps or established that they caused
any observed finite-sample deviation. The draft PDF itself is not required at runtime.

Current truth alpha is the whole normalized Uniform(-2,0) draw; the fitted alpha
is not constrained to unit norm. Full-matrix initialization retains reconstructed
diagonal entries before computing (n I + 11')^{-1} Theta_initial 1. Likelihood and
gradient fitting use observed off-diagonal edges. Probit correction contractions
use expected derivatives/information; observed information standardizes the reported
statistics. See the source for complete definitions.
