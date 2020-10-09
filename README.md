Dirichlet_Process_Mixture_Model_For_RNA_data


This repository contains the code to estimate the models in the paper **On the inference of RNA life cycle kinetic rates from sequencing data by latent Dirichlet Process Mixture Models**.

The code in written in *R* and  to run it the package *OdePackRT1* is needed, that can be installed by the source archive in the directory **R PACKAGE**. Notice that you need to have *OPENMP* installed. 

The posterior sample from the model can be obtained lunching the code *PosteriorSamples.R*. In the code the value of

```R
  NAME      = ""
  DIR       = ".../Dirichlet_Process_Mixture_Model_For_RNA_data/"
  Nthreads  =
  DoClust   =
```

where *NAME* is the name of the output files, *DIR* is the directory where this repository is downloaded, *Nthreads* is the number of cores that should be used for the computations and *DoClust* is *TRUE* if the mixture model is used, *FALSE* if not.
