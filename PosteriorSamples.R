rm(list=ls())
set.seed(100)

#### #### #### #### #### #### #### ####
#### PARAMETERS TO SET
#### #### #### #### #### #### #### ####

NAME  = "" # name of the output files
DIR   = ".../Dirichlet_Process_Mixture_Model_For_RNA_data/" # directory where the github repository is downloaded

Nthreads = 2 # number of cores to be used for the computation
DoClust = T # T if a mixture model is used

#### #### #### #### #### #### #### ####
#### Libraries
#### #### #### #### #### #### #### ####

library(OdePackRT1)

#### #### #### #### #### #### #### ####
#### Data
#### #### #### #### #### #### #### ####

Kmax = 150
DirOUT = paste(DIR,"/out/",sep="")

# The data we want to model
allGenesData  = readRDS(paste(DIR,"data/filteredRpkms.rds",sep=""))

# The kinetic rates estimated by Inspect
InspectRates  = readRDS(paste(DIR,"data/filteredRates.rds",sep=""))

nG    = nrow(allGenesData$premature) # number of genes
nrep  = 3 # number of replicates
nt    = 11 # number of time

## Inspect data
fk = list()
fk[[1]] = InspectRates$synthesis_k1
fk[[2]] = InspectRates$processing_k2
fk[[3]] = InspectRates$degradation_k3

deltat  = 0.5/3
XX = seq(0,16, deltat)

Pmat = matrix(allGenesData$premature[,], ncol=nG, byrow=T)
Mmat = matrix(allGenesData$mature[,], ncol=nG, byrow=T)
Cmat = matrix(allGenesData$nascent[,], ncol=nG, byrow=T)

SS = matrix(unlist(strsplit(colnames(allGenesData$premature),"_")),ncol=2, byrow=T)
MM = data.frame(as.numeric(as.character(SS[,1])),SS[,2])
DATA = matrix(NA, ncol=4, nrow= length(XX)*nG*nrep)
DATA[,4] = XX
MM[MM[,1]==0.2,1] = 1.5
W = round(DATA[,4],2) %in%round(as.numeric(MM[,1]),2)
DATA[W,1:3] =  cbind(c(Pmat), c(Mmat),c(Cmat))

MinKstart = matrix(NA,nrow=nG,ncol=3)
EtaK1start = matrix(NA,nrow=nG,ncol=3)
EtaK2start = matrix(NA,nrow=nG,ncol=3)
MeanKstart = matrix(NA,nrow=nG,ncol=3)
VarKstart = matrix(NA,nrow=nG,ncol=3)
WW = c(1 , 2  ,3  ,4  ,7 ,10, 13, 25, 49,73,97)

# starting values of the algorithm based on the INSPECT data
for(g in 1:nrow(fk[[1]]))
{
  for(kk in 1:3)
  {

    MinKstart[g,kk]  = fk[[kk]][g,1]
    Wm1              = which.max(fk[[kk]][g,])
    Wm2              = which.min(fk[[kk]][g,])

    T1               = (Wm1!=1)&&(Wm1!=11)
    T2               = (Wm2!=1)&&(Wm2!=11)

    if((T1==F)&&(T2==F))
    {
      MM    = abs(fk[[kk]][g,1]+fk[[kk]][g,11])/2
      Wmin  = which.min(abs(fk[[kk]][g,]-MM))
      MeanKstart[g,kk] = XX[WW][Wmin]+ifelse(XX[WW][Wmin]<7,0.2,-0.2)
      EtaK1start[g,kk] = abs(0.5*log(fk[[kk]][g,11]/fk[[kk]][g,1]))
      EtaK2start[g,kk] = abs(0.5*log(fk[[kk]][g,11]/fk[[kk]][g,1]))
      if(fk[[kk]][g,1]>fk[[kk]][g,11])
      {
        EtaK1start[g,kk] = -EtaK1start[g,kk]
      }else{
        EtaK2start[g,kk] = -EtaK2start[g,kk]
      }
      VarKstart[g,kk]  = runif(1,0.001,0.01)

    }else{

      WWWmax  = fk[[kk]][g,which.max(fk[[kk]][g,])]
      WWWmin  = fk[[kk]][g,which.min(fk[[kk]][g,])]

      TT1 = (WWWmax>fk[[kk]][g,1])&&(WWWmax>fk[[kk]][g,11])
      TT2 = (WWWmin<fk[[kk]][g,1])&&(WWWmin<fk[[kk]][g,11])

      if(TT1 && TT2)
      {
        Max1T = max(c(WWWmax-fk[[kk]][g,1],WWWmax-fk[[kk]][g,11]))
        Max2T = max(c(fk[[kk]][g,1]-WWWmin,fk[[kk]][g,11]-WWWmin))

        if(Max1T>Max2T)
        {
          TT3 = T
        }else{
          TT3 = F
        }
      }else{
        TT3 = TT1
      }

      if( TT3 )
      {
        ## max
        WWW  = which.max(fk[[kk]][g,])
        MM   = max(fk[[kk]][g,])
      }else{
        ## min
        WWW  = which.min(fk[[kk]][g,])
        MM   = min(fk[[kk]][g,])
      }

      MeanKstart[g,kk] = XX[WW][WWW]+ifelse(XX[WW][WWW]<7,0.2,-0.2)
      EtaK1start[g,kk] = log(MM/fk[[kk]][g,1])
      EtaK2start[g,kk] = log(MM/fk[[kk]][g,11])
      VarKstart[g,kk]  = runif(1,0.001,0.01)
    }
  }
}

ValInitclust = sample(1:min(Kmax,100),nG*5*3,replace=T)
beta0start  = rep(-4,3)
beta1start  = rep(1.8,3)
MinKstart   = c(t(MinKstart))
EtaK1start  = c(t(EtaK1start))
EtaK2start  = c(t(EtaK2start))
MeanKstart  = c(t(MeanKstart))
VarKstart   = c(t(VarKstart))

SHstart     = c(rep(0.066,33))

DATAApp = matrix(NA, nrow=nG*3*11, ncol=4)
for(g in 1:nG)
{
  ZZ = (1:(97*3)+(g-1)*(97)*3)[c(WW,WW+97,WW+97*2)]
  DATAApp[1:33 +(g-1)*33, ] = DATA[ZZ,]
}

ReCoefStart = list()
SigmaStart = list()
for(i in 1:(Kmax*2))
{
  ReCoefStart[[i]] = rnorm(5,2,1)
  SigmaStart[[i]] = diag(abs(rnorm(5,10,0.2) ))
}

Mod = ODE_MODEL_T2Random(
  NAME        = NAME,
  DIR         = DirOUT,
  Data        = DATA[,1:3],
  nt          = length(XX),
  nrep          = nrep,
  nlatentk      = 3,
  tmin          = 0,
  tmax          = 16,
  nstep         = 2000,
  Kmax         = Kmax,
  priors  = list(

      Beta0 = list(
          type        = "Normal",
          Par1        = rep(0,3),
          Par2        = rep(100,3)
        ),
      Beta1 = list(
        type        = "Normal",
        Par1        = rep(0,3),
        Par2        = rep(100,3)
      ),
      Y0 = list(
        type        = "Uniform",
        Par1        = rep(-5,3),
        Par2        = rep(5,3),
        Par3        = rep(-5.,3),
        Par4        = rep(5.,3),
        Par5        = rep(2,3)
      ),
      MinK   = list(
        type        = "Normal",
        Par1        = rep(0,3),
        Par2        = rep(100,3),
        Par3        = rep(0,3),
        Par4        = rep(200,3),
        Par5        = rep(0,3)
      ),
      EtaK1   = list(
        type        = "Normal",
        Par1        = rep(0,3),
        Par2        = rep(100,3),
        Par3        = rep(-300,3),
        Par4        = rep(300,3),
        Par5        = rep(2,3)
     ),
      EtaK2   = list(
        type        = "Normal",
        Par1        = rep(0,3),
        Par2        = rep(100,3),
        Par3        = rep(-300,3),
        Par4        = rep(300,3),
        Par5        = rep(2,3)
      ),
      MeanK   = list(
        type        = "Normal",
        Par1        = rep(0,3),
        Par2        = rep(100,3),
        Par2        = rep(0,3),
        Par3        = rep(0,3),
        Par4        = rep(16,3),
        Par5        = rep(2,3)
      ),
      VarK   = list(
        type        = "Normal",
        Par1        = rep(0,3),
        Par2        = rep(100,3),
        Par3        = rep(0,3),
        Par4        = rep(20,3),
        Par5        = rep(2,3)
      ),
      RegCoef = list(
          type        = "Normal",
          Par1        = c(0,0,3,8,5),
          Par2        = rep(100,5)
      ),
      sigmaMat = list(
          type        = "InverseWishart",
          Par1        = rep(2,3),
          Par2        = c(diag(1,5))
      ),
      DPpar = list(
          type        = "Gamma",
          Par1        = 1,
          Par2        = 1 ## diagonale Psi
      )
  ),
  start  = list(
    RegCoef   = unlist(ReCoefStart),
    sigmaMat  = unlist(SigmaStart),
    stMat     = ValInitclust,
    piMat     = rep(0.1,Kmax*10),
    DPpar     = rep(1,100),
    Beta0     = beta0start,
    Beta1     = beta1start,
    MinK      = MinKstart,
    VarK      = VarKstart,
    EtaK1     = EtaK1start,
    EtaK2     = EtaK2start,
    MeanK     = MeanKstart,
    SH        = SHstart
  ),
  MCMC = list(
    iter        = 15000000,
    burnin      = 1425000,
    thin        = 30
  ),
  Adapt= list(
    doUpdate    = T,
    start       = 5000,
    end         = 0.9,
    exp         = 0.6,
    acc         = 0.534,
    batch       = 20,
    epsilon     = 0.000000001,
    lambda      = 1,
    sd          = list(
      General     = 0.001,
      SF          = 0.0001
    )
  ),
  IterShow          = 1000,
  DoParallel        = T,
  Nthreads          = Nthreads,
  Verbose           = 1,
  saveODE           = T,
  Rparam            = F,
  Spike             = T,
  LimSpike          = 10,
  LimC              = 0.01,
  sample_Var        = T,
  sample_Ode        = T,
  SampleOdeType     = 1,
  abs_err           = 10^(-8),
  rel_err           = 10^(-4),
  sample_Pi         = T,
  sample_DpPar      = T,
  sample_st         = T,
  sample_RandomPar  = T,
  IterStartRandom   = ifelse(DoClust==T,50000,50000*1000000),
  addStartRandom    = 0,
  Use_pi            = T,
  ExpMolt           = 0.00,
  ExpLambda         = 10,
  Do_MergeSplit     = T,
  TypeK             = 2,
  MoltSF            = 100,
  A                 = NULL
)

# OUTPUT
save.image(paste(DirOUT,NAME,".Rdata",sep=""))
