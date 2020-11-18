**Dirichlet_Process_Mixture_Model_For_RNA_data**
 

This repository contains the code to estimate the models in the paper **On the inference of RNA life cycle kinetic rates from sequencing data by latent Dirichlet Process Mixture Models**.

The code is written in *R* and *C++* and to run it the package *OdePackRT1* is needed, that can be installed by the source archive in the directory **R PACKAGE**. Notice that you need to have *OPENMP* installed.

The posterior sample from the model can be obtained launching the code *PosteriorSamples.R*. The following variables must be specified
```R
  NAME      = ""
  DIR       = ".../Dirichlet_Process_Mixture_Model_For_RNA_data/"
  Nthreads  =
  DoClust   =
```

where *NAME* is the name of the output files, *DIR* is the directory where this repository is downloaded, *Nthreads* is the number of cores that should be used for the computations, and *DoClust* is *TRUE* if the mixture model is used, *FALSE* if not. The results are stored in the file *NAME.Rdata* in the directory **OUT**.

The model output is a list of 18 elements
* **MinK** parameter \beta in the paper. A matrix with number of rows equal to the number of posterior samples, and number of columns equal to the number of parameters. The first 3 columns are the parameters of the k_1, k_2 and k_3 of the first gene. The next three the one of the second...
* **VarK** parameter \sigma^2 in the paper. A matrix with number of rows equal to the number of posterior samples, and number of columns equal to the number of parameters. The first 3 columns are the parameters of the k_1, k_2 and k_3 of the first gene. The next three the one of the second...
* **MeanK** parameter \mu in the paper. A matrix with number of rows equal to the number of posterior samples, and number of columns equal to the number of parameters. The first 3 columns are the parameters of the k_1, k_2 and k_3 of the first gene. The next three the one of the second...
* **EtaK1** parameter \eta_1 in the paper. A matrix with number of rows equal to the number of posterior samples, and number of columns equal to the number of parameters. The first 3 columns are the parameters of the k_1, k_2 and k_3 of the first gene. The next three the one of the second...
* **EtaK2** parameter \eta_2 in the paper. A matrix with number of rows equal to the number of posterior samples, and number of columns equal to the number of parameters. The first 3 columns are the parameters of the k_1, k_2 and k_3 of the first gene. The next three the one of the second...
* **Y** the estimated ODEs. A matrix with number of rows equal to the number of posterior samples, and number of columns equal to the NumberTime*NumberGene*3. The first 3 columns are the ODEs of the first time point for the first gene, the next three are the second time point of the first gene ... row 34, 35, 36 contains the first time point of the secong gene, ...
* **Y0** the estimated first time points of the ODEs. A matrix with number of rows equal to the number of posterior samples, and number of columns equal to NumberGene*3. The first 3 columns are the first time point of the ODEs of the first gene, the next three the ones of the second, ...
* **Beta0** regressive parameter \beta_{j,0} in the paper. A matrix with number of rows equal to the number of posterior samples, and number of columns equal to 3
* **Beta1** regressive parameter \beta_{j,1} in the paper. A matrix with number of rows equal to the number of posterior samples, and number of columns equal to 3
* **SF_1** parameter \rho_{1}(t) in the paper. A matrix with number of rows equal to the number of posterior samples, and number of columns equal to the number of time points
* **SF_2** parameter \rho_{2}(t) in the paper. A matrix with number of rows equal to the number of posterior samples, and number of columns equal to the number of time points
* **SF_3** parameter \rho_{3}(t) in the paper. A matrix with number of rows equal to the number of posterior samples, and number of columns equal to the number of time points
* **ModelOutput_11** Mixture model output for the shape parameters of the first kinetic rates
* **ModelOutput_12** Mixture model output for the starting level of the first kinetic rates
* **ModelOutput_21** Mixture model output for the shape parameters of the second kinetic rates
* **ModelOutput_22** Mixture model output for the starting level of the second kinetic rates
* **ModelOutput_31** Mixture model output for the shape parameters of the third kinetic rates
* **ModelOutput_32** Mixture model output for the starting level of the third kinetic rates

Each *Mixture model output* is composed of a list of 2 elements, where *PostK* is the MAP estimator of the number of clusters for each gene, and
*ModelOutput* is a list of 150 elements, where the k-th element contains the posterior samples of the iterations where the number of clusters is k. *ModelOutput[[k]]* is a list of two elements
* **LikelihoodParameters**
  * *st* parameter z in the paper. A matrix with number of rows equal to the number of posterior samples, and number of columns equal to the number of genes
  * *MAP* MAP estimator of the number of clusters for each gene
  * *pi* parameter pi in the paper. A matrix with number of rows equal to the number of posterior samples, and number of columns equal to k
  * *DPpar* parameter \zeta in the paper.
* **Clustering**
  * *Beta* parameter \Zeta (mean of the normal distribution) in the paper. A matrix with number of rows equal to the number of posterior samples, and number of columns equal to k (if it is the cluster object of the starting level) or k*4 (if it is the cluster object of the shape parameters). The first 1 or 4 elements, are related to the first cluster, the second  1 or 4 to the second, ...
  * *Sigma* parameter \Omega (covariance matrix of the normal distribution) in the paper. A matrix with number of rows equal to the number of posterior samples, and number of columns equal to k (if it is the cluster object of the starting level) or k*4^2 (if it is the cluster object of the shape parameters). The first 1 or 4^2 elements, are related to the first cluster, the second  1 or 4^2 to the second, ...
