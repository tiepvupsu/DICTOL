
# Discriminative dictionary Learning Toolbox for Classification.
_This repository is under construction_

<!-- MarkdownTOC -->

- [Notation](#notation)
- [Sparse Representation-based classification](#sparse-representation-based-classification)
- [Online Dictionary Learning](#online-dictionary-learning)
- [LCKSVD](#lcksvd)
- [FDDL](#fddl)
- [DLSI](#dlsi)
  - [Cost function](#cost-function)
  - [Training DLSI](#training-dlsi)
  - [DLSI predict new sample](#dlsi-predict-new-sample)
  - [Demo](#demo)
- [DLCOPAR](#dlcopar)
  - [`DLCOPAR`](#dlcopar-1)
  - [`DLCOPAR_cost`](#dlcoparcost)
  - [`DLCOPAR_updateX`](#dlcoparupdatex)
  - [`DLCOPAR_updateD`](#dlcoparupdated)
  - [`DLCOPAR_pred`](#dlcoparpred)
  - [`DLCOPAR_top`](#dlcopartop)
- [LRSDL](#lrsdl)
    - [`LRSDL_top`](#lrsdltop)
- [References](#references)

<!-- /MarkdownTOC -->




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
* `k_c`: number of bases in class-specific dictionary `c`. Typically, all `n_c` are the same and equal to `k`.
* `k_0`: number of bases in the shared-dictionary 
* `K`: total number of dictionary bases. 
* `D_range`: similar to `Y_range` but used for dictionary without the shared dictionary. 




# Sparse Representation-based classification
* Sparse Representation-based classification implementation [[1]](#fn_src).
* Classification based on SRC.
* Syntax: `[pred, X] = SRC_pred(Y, D, D_range, opts)`
  - INPUT:
    + `Y`: test samples.
    + `D`: the total dictionary. `D = [D_1, D_2, ..., D_C]` with `D_c` being the _c-th_ class-specific dictionary.
    + `D_range`: range of class-specific dictionaries in `D`. See also [Notation](#notation).
    + `opts`: options.
      * `opts.lambda`: `lambda` for the Lasso problem. Default: `0.01`.
      * `opts.max_iter`: maximum iterations of fista algorithm. Default:  `100`. [Check this simple implementation of FISTA](https://github.com/tiepvupsu/FISTA) 
  - OUTPUT:
    + `pred`: predicted labels of test samples.
    + `X`: solution of the lasso problem.

# Online Dictionary Learning
* An implementation of the well-known Online Dictionary Learning method [[2]](#fn_old).
* Solving the dictionary learning problem:
<img src = "http://latex2png.com/output//latex_c1709f7ea3f7f523694cf0d6b9a61aa9.png" height = "40"/>    

* Syntax: `[D, X] = ODL(Y, k, lambda, opts, sc_method)`
  - INPUT: 
    + `Y`: collection of samples.
    + `k`: number of bases in the desired dictionary.
    + `lambda`: norm 1 regularization parameter.
    + `opts`: option.
    + `sc_method`: sparse coding method used in the sparse coefficient update. Possible values:
      - `'fista'`: using FISTA algorithm. See also [`fista`](#fista).
      - `'spams'`: using SPAMS toolbox [[12]](#fn_spams). 
  - OUTPUT:
    + `D, X`: as in the problem.

# LCKSVD


# FDDL


# DLSI

## Cost function

<img src = "http://latex2png.com/output//latex_7fdbb1d22a9436aa96c371b0edee496f.png" height = "50"/>

## Training DLSI 
* function `[D, X, rt] = DLSI(Y, Y_range, opts)`
* The main DLSI algorithm 
* INPUT: 
  - `Y, Y_range`: training samples and their labels 
  - `opts`: 
    + `opts.lambda, opts.eta`: `lambda` and `eta` in the cost function 
    + `opts.max_iter`: maximum iterations. 
* OUTPUT:
  - `rt`: total running time of the training process.   
* `[D, X, rt] = argmin_{D, X}(sum 0.5*||Y_c - D_c X_c||_F^2) + lambda*norm1(X) + 0.5*eta * sum_{i \neq c} ||D_i^T D_c||_F^2`
* updating `X`:
  `Xi = arg\min 0.5*||Y_c - D_c X_c||_F^2 + lambda ||X_c||`
* updating `D`:
  `Di = argmin ||Y_c - D_c X_c||_F^2 + \eta ||D_{\c}^T D_c||_F^2`

## DLSI predict new sample 
* function `pred = DLSI_pred(Y, D, opts)`
* predict the label of new input `Y` given the trained dictionary `D` and 
parameters stored in `opts` 

## Demo 
Run `DLSI_top` in Matlab command window.

# DLCOPAR
## `DLCOPAR`
* function `[D, X, rt] = DLCOPAR(Y, Y_range, opts)`
* The main DLCOPAR alg

## `DLCOPAR_cost`
* function cost = DLSI_cost(Y, Y_range, D, D_range, X, opts)        
* Calculating cost function of DLCOPAR with parameters lambda and eta are stored in  `opts.lambda` and `opts.rho`.
* `f(D, X) = 0.5*sum_{c=1}^C 05*||Y - DX||_F^2 + sum_{c=1}^C ( ||Y_c - D_Cp1 X^Cp1_c - D_c X_c^c||F^2 + sum_{i != c}||X^i_c||_F^2) + lambda*||X||_1 + 0.5*eta*sum_{i \neq c}||Di^T*Dc||_F^2`

## `DLCOPAR_updateX`
* function X = DLCOPAR_updateX(Y, Y_range, D, X, opts)
* updating X in DLCOPAR. 

## `DLCOPAR_updateD` 
* function D = DLCOPAR_updateD(Y, Y_range, D, X, opts) 
* update D in DLCOPAR, including both PARTICULAR dictionaries and the COMMON dictionary.
* The algorithm used here is the efficient algorithm presented in LRSDL paper 

## `DLCOPAR_pred`
* function pred = DLCOPAR_pred(Y, D, D_range_ext, opts)
* predict label of the input Y
* INPUT:
   + `opts.classify_mode` = either 'GC' (global coding) or 'LC' (local coding)

## `DLCOPAR_top`
  * function DLCOPAR_top(dataset, N_train, k, k0, lambda, eta)
  * The top function of DLCOPAR 
  * INPUT:
    + `dataset`: name of the dataset stored in `.mat` file in `data` folder. 
      Note that `dataset` is the file name of the `.mat`, excluding `.mat`.
    + `N_train`: number of training samples in each class 
    + `k`: number of bases in EACH PARTICULAR dictionary 
    + `k0`: number of bases in the COMMON dictionary
    + `lambda, eta`: regularization parameters.
  * To run an small example, type `DLCOPAR_top` without input in MATLAB command window.



# LRSDL
### `LRSDL_top`
* Syntax `LRSDL_top(dataset, N_train, k, lambda1, lambda2, lambda3)`



# References

<a name="fn_src">[1]</a>. Wright, John, et al. "Robust face recognition via sparse representation." _Pattern Analysis and Machine Intelligence, IEEE Transactions on_ 31.2 (2009): 210-227. [paper](http://www.columbia.edu/~jw2966/papers/WYGSM09-PAMI.pdf )

<a name="fn_old">[2]</a>. Mairal, Julien, et al. "Online learning for matrix factorization and sparse coding." _The Journal of Machine Learning Research 11_ (2010): 19-60. [[paper]](http://www.di.ens.fr/~fbach/mairal10a.pdf)

<a name="fn_lck">[3]</a>. Jiang, Zhuolin, Zhe Lin, and Larry S. Davis. "Label consistent K-SVD: Learning a discriminative dictionary for recognition." _Pattern Analysis and Machine Intelligence, IEEE Transactions on_ 35.11 (2013): 2651-2664. [[Project page]](http://www.umiacs.umd.edu/~zhuolin/projectlcksvd.html)

<a name="fn_fdd">[4]</a>. Yang, Meng, et al. "Fisher discrimination dictionary learning for sparse representation." _Computer Vision (ICCV), 2011 IEEE International Conference on. IEEE_, 2011. [[paper]](http://www4.comp.polyu.edu.hk/~cslzhang/paper/conf/iccv11/FDDL_ICCV_final.pdf), [[code]](http://www4.comp.polyu.edu.hk/~cslzhang/code/FDDL.zip)

<a name="fn_dls">[5]</a>. Ramirez, Ignacio, Pablo Sprechmann, and Guillermo Sapiro. "Classification and clustering via dictionary learning with structured incoherence and shared features." _Computer Vision and Pattern Recognition (CVPR), 2010 IEEE Conference on. IEEE_, 2010. [[paper]](http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=5539964&tag=1)

<a name="fn_dfd">[6]</a>. Discriminative Feature-Oriented dictionary Learning.
 Tiep H. Vu, H. S. Mousavi, V. Monga, A. U. Rao and G. Rao, "Histopathological Image Classification using Discriminative Feature-Oriented dictionary Learning", _IEEE Transactions on Medical Imaging_ , volume 35, issue 3, pages 738-751, March 2016. [[paper]](http://arxiv.org/pdf/1506.05032v5.pdf) [[Project page]](http://signal.ee.psu.edu/dfdl.html)

<a name="fn_cor">[7]</a>. Kong, Shu, and Donghui Wang. "A dictionary learning approach for classification: separating the particularity and the commonality." _Computer Vision ECCV_ 2012. Springer Berlin Heidelberg, 2012. 186-199. [[paper]] (http://www.cs.zju.edu.cn/people/wangdh/papers/draft_ECCV12_particularity.pdf)

<a name="fn_lrs">[8]</a>. Tiep H. Vu, Vishal Monga. "Learning a low-rank shared dictionary for object classification." Submitted to International Conference on Image Processing (ICIP) 2016. [[paper]](http://arxiv.org/abs/1602.00310)

<a name="fn_shr">[9]</a>. A singular value thresholding algorithm for matrix completion." _SIAM Journal on Optimization_ 20.4 (2010): 1956-1982. [[paper]](http://arxiv.org/pdf/0810.3286v1.pdf)

<a name="fn_fista">[10]</a>. Beck, Amir, and Marc Teboulle. "A fast iterative shrinkage-thresholding algorithm for linear inverse problems." _SIAM journal on imaging sciences_ 2.1 (2009): 183-202. [[paper]](http://people.rennes.inria.fr/Cedric.Herzet/Cedric.Herzet/Sparse_Seminar/Entrees/2012/11/12_A_Fast_Iterative_Shrinkage-Thresholding_Algorithmfor_Linear_Inverse_Problems_(A._Beck,_M._Teboulle)_files/Breck_2009.pdf)

<a name="fn_spams"> [11]</a>. [The Sparse Modeling Software](http://spams-devel.gforge.inria.fr/)
