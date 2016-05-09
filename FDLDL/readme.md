# Fuzzy Dictionary Learning idea 
* Folder: FDLDL
* Combination of FDDL and Latent Dictionary Learning
* Time created: 4/7/2016 9:13:16 AM

## Cost function contains:
 * __No Reconstruction term ||Y - DX||__
 * The Fidelity term ||Yc - DWcXc||_F^2 + \sum ||DWiXc||_F^2
 * The Fisher term on sparse coefficients (lambda2)
 * Sparsity term (lambda1)
 * Latent DL constraint on `W` and `D`.
 * Nonnegative constrain only. 
 * __No equality constrain__

## Solution:

## Todo:
  - [x] cost function unchanged
  - [ ] updateX
      + Similar to Fuzzy Dictionary
  - [ ] updateD 
  - [ ] updateW


 