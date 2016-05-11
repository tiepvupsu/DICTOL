
# Discriminative dictionary Learning Toolbox for Classification.
_This repository is under construction_

* [__SRC__](#src): Sparse Representation-based Classification [[1]](#fn_src). 

* [__ODL__](#odl): Online dictionary Learning [[2]](#fn_odl).

* [__LCKSVD__](#lcksvd): Label consistent K-SVD[^fn_lcksvd] [[3]](#fn_ldl).

* [__FDDL__](#FDDL): Fisher discrimination dictionary learning [[4]](#fn_fdd).

* [__DLSI__](#DLSI): dictionary Learning with Structured Incoherence [[5]](#fn_dls)

* [__DFDL__](#dfdl): Discriminative Feature-Oriented dictionary Learning [[6]](#fn_dfd)

* [__DLCOPAR__](#dlcopar) [[7]](#fn_cor) 

* [__LDL__](#LDL): Latent dictionary learning [[8]](#fn_ldl).

* [__LRSDL__](#LRSDL): Low-rank Shared dictionary Learning [[9]](#fn_lrs). 

* [__FuzzyDL__](#FuzzyDL):

# Notation
* `Y`: signals.
* `D`: dictionary.
* `X`: sparse coefficient.
* `d`: signal dimension. `d = size(Y, 1)`.
* `C`: number of classes.
* `c`: class index.
* `n_c`: number of training samples in class `c`. Typically, all `n_c` are the same and equal to `n`.
* `N`: total number of training samples.
* `Y_range`: an array storing range of each class, suppose that labels are sorted in a ascending order. 
  Example: If `Y_range = [0, 10, 25]`, then:
     - There are two classes, samples from class 1 range from 1 to 10, from class 2 range from 11 to 25. 
     - In general, samples from class `c` range from `Y_range(c) + 1` to `Y_range(c+1)`
     - We can observe that number of classes `C = numel(Y_range) - 1`.
* `k_c`: number of atoms in class-specific dictionary `c`. Typically, all `n_c` are the same and equal to `k`.
* `k_0`: number of atoms in the shared-dictionary 
* `K`: total number of dictionary atoms. 
* `D_range`: similar to `Y_range` but used for dictionary without the shared dictionary. 

# Support functions

All of the following functions are located in subfolder `utils`.

#### `get_block_col`
* Extract a block of columns from a matrix.
* Syntax: `Mc = get_block_col(M, c, col_range)`
    - `M`: the big matrix `M = [M_1, M_2, ...., M_C]`.
    - `c`: block index.
    - `col_range`: range of samples, see `Y_range` and `D_range` above.
* Example: `M` has 25 columns and `col_range = [0, 10, 25]`, then `get_block_col(M, 1, col_range)` will output the first block of `M`, i.e. `M(:, 1:10)`.

#### `get_block_row`
* Extract a block of rows from a matrix.
* Syntax: `Mc = get_block_row(M, c, row_range)`
    - `M`: the big matrix `M = [M_1; M_2; ....; M_C]`.
    - `c`: block index.
    - `row_range`: range of samples, see `Y_range` and `D_range` above.
* Example: `M` has 40 rows and `row_range = [0, 10, 25, 40]`, then `get_block_row(M, 2, row_range)` will output the second block of `M`, i.e. `M(11:25, :)`.

#### `get_block`
* Extract a submatrix of a matrix 
* Syntax: `Mij = get_block(M, i, j, row_range, col_range)`
    - `M` the big matrix: `M = [ M11, M12, ..., M1m; M21, M22, ..., M2m; ... ; Mn1, Mn2,..., Mnm]`
    - `i`: row block index 
    - `j`: column block index 
    - `row_range`: row range
    - `col_range`: columns range 
* Note: `get_block(M, i, j, row_range, col_range) = get_block_col(get_block_row(M, i, row_range), j, col_range).`


#### `label_to_range`
* Convert from Labels to Ranges
* Example: if `label = [1 1 1 2 2 2 2 3 3]`, then `range = [0, 3, 7, 9]`.
* Syntax: `range = label_to_range(label)`

#### `range_to_label`
* Convert from Ranges to Labels
* Example: if `range = [0, 3, 5]`` then `label = [1 1 1 2 2]``
* Syntax: `label = range_to_label(range)`

#### `norm1`
* Return norm 1 of a matrix, which is sum of absolute value of all element of that matrix.
* Syntax: `res = norm1(X)`

#### `normF2`
* Return square of the Frobenius norm, which is sum of square of all elements in a matrix
* Syntax: `res = normF2(X)`

#### `normc`
* Normalize columns of a matrix: norm 2 of each columns equals to 1. This function is a built-in function in some recent MATLAB version.
* Syntax: `M1 = normc(M1)`

#### `vec`
* Vectorization of a matrix. This function is a built-in function in some recent MATLAB version.
* Syntax: `a = vec(A)`

#### `nuclearnorm`
* Return nuclear norm of a matrix.
* Syntax `res = nuclearnorm(X)`
* `

#### `shrinkage`
* Soft thresholding function.
* Syntax: ` X = shrinkage(U, lambda)`
* Solve the following optimization problem:
  `X = arg min_X 0.5*||X - U||_F^2 + lambda||X||_1`
  where `U` and `X` are matrices with same sizes. `lambda` can be either positive a scalar or a positive matrix (all elements are positive) with same size as `X`. In the latter case, it is a weighted problem.

#### `shrinkage_rank`
* Singular value thresholding algorithm for matrix completion [[10]](#fn_shr).
* Syntax: `Y = shrinkage_rank(D, lambda)` 
* Solve the following optimization problem:
  `X = arg min_X 0.5*||X - U||_F^2 + lambda*||X||_*`
  where `||X||_*` is the nuclear norm.

#### `fista`
* A Fast Iterative Shrinkage-Thresholding Algorithm for Linear Inverse Problems[[11]](#fn_fista).
* Solve the problem: `X = arg min_X F(X) = f(X) + lambda||X||_1` where:
  - `X`: variable, can be a matrix.
  - `f(X)` is a smooth convex function with continuously differentiable with Lipschitz continuous gradient `L(f)` (Lipschitz constant of the gradient of `f`).
* Syntax: `[X, iter] = fista(grad, Xinit, L, lambda, opts, calc_F)` where:
  - INPUT:
    + `grad`: a _function_ calculating gradient of `f(X)` given `X`.
    + `Xinit`: initial guess.
    + `L`: the Lipschitz constant of the gradient of `f(X)`.
    + `lambda`: a regularization parameter, can be either positive a scalar or a weighted matrix.
    + `opts`: a _structure_ variable describing the algorithm.
      * `opts.max_iter`: maximum iterations of the algorithm. Default `300`.
      * `opts.tol`: a tolerance, the algorithm will stop if difference between two successive `X` is smaller than this value. Default `1e-8`.
      * `opts.show_progress`: showing `F(X)` after each iteration or not. Default `false`. 
    + `calc_F`: optional, a _function_ calculating value of `F` at `X` via `feval(calc_F, X)`. 
  - OUTPUT:
    + `X`: solution.
    + `iter`: number of iterations.

#### `lasso_fista`
* Syntax: `[X, iter] = lasso_fista(Y, D, Xinit, lambda, opts)`
* Solving a Lasso problem using FISTA [[11]](#fn_fista): `X = arg min_X 0.5*||Y - DX||_F^2 + lambda||X||_1`. Note that `lambda` can be either a positive scalar or a matrix with positive elements.
  - INPUT:
    + `Y, D, lambda`: as in the problem.
    + `Xinit`: Initial guess 
    + `opts`: options. See also [`fista`](#fista)
  - OUTPUT:
    + `X`: solution.
    + `iter`: number of fistat iterations.
* **Note**:
  - _To see a toy example, un this function without inputs_
  - _Can be used for solving a Weighted Lasso problem_.






# SRC
* Sparse Representation-based classification implementation [[1]](#fn_src).
* Folder `SRC`.

#### `SRC_pred`
* Classification based on SRC.
* Syntax: `[pred, X] = SRC_pred(Y, D, D_range, opts)`
  - INPUT:
    + `Y`: test samples.
    + `D`: the total dictionary. `D = [D_1, D_2, ..., D_C]` with `D_c` being the _c-th_ class-specific dictionary.
    + `D_range`: range of class-specific dictionaries in `D`. See also [Notation](#notation), [`label_to_range`](#label_to_range), [`range_to_label`](#range_to_label), [`get_block_col`](#get_block_col).
    + `opts`: options.
      * `opts.lambda`: `lambda` for the Lasso problem.
      * `opts.max_iter`: maximum iterations of fista algorithm. See also [fista](#fista), [lasso_fista](#lasso_fista).
      * others.
  - OUTPUT:
    + `pred`: predicted labels of test samples.
    + `X`: solution of the lasso problem.

# ODL
* An implementation of the well-known Online Dictionary Learning method [[2]](#fn_old).
* Solving the dictionary learning problem:

   `[D, X] = arg min_{D, X} 0.5||Y - DX||_F^2 + lambda||X||_1` subject to `||d_i||_2 <= 1`.
* Folder `./ODL`.

## Supporting funtions:

### `ODL`
* The algorithm to solve the main problem stated above.
* Syntax: `[D, X] = ODL(Y, k, lambda, opts, sc_method)`
  - INPUT: 
    + `Y`: collection of samples.
    + `k`: number of atoms in the desired dictionary.
    + `lambda`: norm 1 regularization parameter.
    + `opts`: option.
    + `sc_method`: sparse coding method used in the sparse coefficient update. Possible values:
      * `'fista'`: using FISTA algorithm. See also [`fista`](#fista).
      * `'spams'`: using SPAMS toolbox [[12]](#fn_spams). 
  - OUTPUT:
    + `D, X`: as in the problem.

 
### `ODL_cost`
* Calculating cost 
* Syntax `cost = ODL_cost(Y, D, X, lambda)`

### `ODL_updateD`
* The dictionary update algorithm in ODL. 
* Solving the optimization problem:

  `D = arg min_D -2trace(E'*D) + trace(D*F*D')` subject to: `||d_i||_2 <= 1`,  where `F` is a positive semidefinite matrix. 
* Syntax `[D, iter] = ODL_updateD(D, E, F, opts)`
  - INPUT: 
    + `D, E, F` as in the above problem.
    + `opts`. options:
      * `opts.max_iter`: maximum number of iterations.
      * `opts.tol`: when the difference between `D` in two successive iterations less than this value, the algorithm will stop.
  - OUTPUT:
    + `D`: solution.
    + `iter`: number of run iterations.


# LCKSVD


# FDDL
## `FDDL_fidelity`
* Syntax: cost = FDDL_fidelity(Y, Y_range, D, D_range, X)
* Calculating the fidelity term in FDDL[[4]](#fn_fdd):
* $\sum_{c=1}^C \Big(\|Y_c - D_cX^c_c\|_F^2 + \sum_{i \neq c} \|D_c X^c_i\|_F^2\Big)$

## `FDDL_discriminative`
* Syntax: cost = FDDL_discriminative(X, Y_range)
* calculating the discriminative term in FDDL[[4]](#fn_fdd):
* $\|X\|_F^2 + \sum_{c=1}^C (\|Xc - Mc\|_F^2 - \|Mc - M\|_F^2) $

## `FDDL_cost`
* `FDDL_cost(Y, Y_range, D, D_range, X, lambda1, lambda2)`


## `FDDL_updateX`

## `FDDL_updateD`

## `FDDL_pred`

# DLSI
## `DLSI_term`
* Syntax: cost = DLSI_term(D, D_range)
* Calculating the structured incoherence term in DLSI [[5]](#fn_dls).
* $\sum_{c=1}^C \sum_{i \neq c} \|D_i^TD_c\|_F^2$

## `DLSI_cost`
## `DLSI_updateX`
#### `DLSI_updateD`
#### `DLSI_pred`

# DFDL

# DLCOPAR

# LDL

# LRSDL
### `LRSDL_top`
* Syntax `LRSDL_top(dataset, N_train, k, lambda1, lambda2, lambda3)`


# FuzzyDL 

# References

<a name="fn_src">[1]</a>. Wright, John, et al. "Robust face recognition via sparse representation." _Pattern Analysis and Machine Intelligence, IEEE Transactions on_ 31.2 (2009): 210-227. [paper](http://www.columbia.edu/~jw2966/papers/WYGSM09-PAMI.pdf )

<a name="fn_old">[2]</a>. Mairal, Julien, et al. "Online learning for matrix factorization and sparse coding." _The Journal of Machine Learning Research 11_ (2010): 19-60. [[paper]](http://www.di.ens.fr/~fbach/mairal10a.pdf)

<a name="fn_lck">[3]</a>. Jiang, Zhuolin, Zhe Lin, and Larry S. Davis. "Label consistent K-SVD: Learning a discriminative dictionary for recognition." _Pattern Analysis and Machine Intelligence, IEEE Transactions on_ 35.11 (2013): 2651-2664. [[Project page]](http://www.umiacs.umd.edu/~zhuolin/projectlcksvd.html)

<a name="fn_fdd">[4]</a>. Yang, Meng, et al. "Fisher discrimination dictionary learning for sparse representation." _Computer Vision (ICCV), 2011 IEEE International Conference on. IEEE_, 2011. [[paper]](http://www4.comp.polyu.edu.hk/~cslzhang/paper/conf/iccv11/FDDL_ICCV_final.pdf), [[code]](http://www4.comp.polyu.edu.hk/~cslzhang/code/FDDL.zip)

<a name="fn_dls">[5]</a>. Ramirez, Ignacio, Pablo Sprechmann, and Guillermo Sapiro. "Classification and clustering via dictionary learning with structured incoherence and shared features." _Computer Vision and Pattern Recognition (CVPR), 2010 IEEE Conference on. IEEE_, 2010. [[paper]](http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=5539964&tag=1)

<a name="fn_dfd">[6]</a>. Discriminative Feature-Oriented dictionary Learning[^fn_dfdl].
 Tiep H. Vu, H. S. Mousavi, V. Monga, A. U. Rao and G. Rao, "Histopathological Image Classification using Discriminative Feature-Oriented dictionary Learning", _IEEE Transactions on Medical Imaging_ , volume 35, issue 3, pages 738-751, March 2016. [[paper]](http://arxiv.org/pdf/1506.05032v5.pdf) [[Project page]](http://signal.ee.psu.edu/dfdl.html)

<a name="fn_cor">[7]</a>. Kong, Shu, and Donghui Wang. "A dictionary learning approach for classification: separating the particularity and the commonality." _Computer Vision ECCV_ 2012. Springer Berlin Heidelberg, 2012. 186-199. [[paper]] (http://www.cs.zju.edu.cn/people/wangdh/papers/draft_ECCV12_particularity.pdf)

<a name="fn_ldl">[8]</a>. Yang, Meng, et al. "Latent dictionary learning for sparse representation based classification." _Proceedings of the IEEE Conference on Computer Vision and Pattern Recognition_. 2014.[[paper]](https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&cad=rja&uact=8&ved=0ahUKEwjRsIrixPjLAhWE9h4KHc-JDiUQFggiMAA&url=http%3A%2F%2Fwww.cv-foundation.org%2Fopenaccess%2Fcontent_cvpr_2014%2Fpapers%2FYang_Latent_dictionary_Learning_2014_CVPR_paper.pdf&usg=AFQjCNGSzR64f4QkoHn6D-668-wpB2xIcQ&sig2=oKGPUp1o-L9O1Q1XHwhZcg)

<a name="fn_lrs">[9]</a>. Tiep H. Vu, Vishal Monga. "Learning a low-rank shared dictionary for object classification." Submitted to International Conference on Image Processing (ICIP) 2016. [[paper]](http://arxiv.org/abs/1602.00310)

<a name="fn_shr">[10]</a>. A singular value thresholding algorithm for matrix completion." _SIAM Journal on Optimization_ 20.4 (2010): 1956-1982. [[paper]](http://arxiv.org/pdf/0810.3286v1.pdf)

<a name="fn_fista">[11]</a>. Beck, Amir, and Marc Teboulle. "A fast iterative shrinkage-thresholding algorithm for linear inverse problems." _SIAM journal on imaging sciences_ 2.1 (2009): 183-202. [[paper]](http://people.rennes.inria.fr/Cedric.Herzet/Cedric.Herzet/Sparse_Seminar/Entrees/2012/11/12_A_Fast_Iterative_Shrinkage-Thresholding_Algorithmfor_Linear_Inverse_Problems_(A._Beck,_M._Teboulle)_files/Breck_2009.pdf)

<a name="fn_spams"> [12]</a>. [The Sparse Modeling Software](http://spams-devel.gforge.inria.fr/)
