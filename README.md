# GT-COPR

## Summary of the package
We provide MATLAB codes for GT-CORP algorithm in the following files:

- **GTCOPR.m:** main function of the GTCOPR algorithm
- **tensor_R.m**: convert the pairwise relations to tensors
- **buil_AUX.m:** initialize the auxiliary variables
- **update_AUX.m**: update the auxiliary variables
- **update_collap.m**: compute the derivatives of the first objective which involves tensor collapsing
- **update_graphReg.m**: compute the derivatives of the second objective which involves graph regularization
- **toy.m**: a drive file to test GTCOPR on a random dataset

Note: the tensor framework of GT-COPR implementation is based on [MATLAB Tensor Toolbox](http://www.sandia.gov/~tgkolda/TensorToolbox/), please install/download it first.


## References
1: Zhuliu Li, Wei Zhang, R Stephanie Huang, and Rui Kuang. "Learning a Low-rank Tensor of Pharmacogenomic Multi-relations from Biomedical Networks." IEEE International Conference on Data Mining (ICDM), 2019. 

2: Brett W. Bader, Tamara G. Kolda and others. MATLAB Tensor Toolbox Version 2.6, Available online, February 2015. (URL: http://www.sandia.gov/~tgkolda/TensorToolbox/).
