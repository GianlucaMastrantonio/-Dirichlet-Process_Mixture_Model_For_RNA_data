#### #### #### #### #### #### #### ####
#### Directory
#### #### #### #### #### #### #### ####
rm(list=ls())
set.seed(100)
Nthreads = 2
DoClust = T


NAME  = "ModelOut"
DIR   = "/Users/gianlucamastrantonio/Dropbox (Politecnico di Torino Staff)/lavori/ode/mattia-amir-gianluca-enrico/Submit/Dirichlet_Process_Mixture_Model_For_RNA_data/"


DirOUT = paste(DIR,"/out/",sep="")
#### #### #### #### #### #### #### ####
#### Libraries
#### #### #### #### #### #### #### ####

library(OdePackRT1)

# library(BH)
# library(coda)
# library(MCMCpack)
# library(odeintr)
# library(Rcpp)
# library(deSolve)
# library(sp)
# library(akima)
# library(evd)
# library(truncdist)
#
# library(label.switching)

#### #### #### #### #### #### #### ####
#### Data
#### #### #### #### #### #### #### ####
Kmax = 150
# allGenesData = readRDS(paste(DIRdata,"filteredRpkms.rds",sep=""))
# funckvalue   = readRDS(paste(DIRdata,"filteredRates.rds",sep=""))
# The data we want to model
#load(paste(DIR,"data/Data.Rdata",sep=""))
allGenesData  = readRDS(paste(DIR,"data/filteredRpkms.rds",sep=""))

# The kinetic rates estimated by Inspect
InspectRates  = readRDS(paste(DIR,"data/filteredRates.rds",sep=""))

nG    = nrow(allGenesData$premature) # number of genes
nrep  = 3 # number of replicates
nt    = 11 # number of time


deltat  = 0.5/3
XX = seq(0,16, deltat)





## Inspect data
fk = list()
fk[[1]] = InspectRates$synthesis_k1
fk[[2]] = InspectRates$processing_k2
fk[[3]] = InspectRates$degradation_k3





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
      # EtaK1start[g,kk] = ifelse(fk[[kk]][g,11]>fk[[kk]][g,1],1,-1)*0.5*log(fk[[kk]][g,11]/fk[[kk]][g,1])
      # EtaK2start[g,kk] = ifelse(fk[[kk]][g,11]>fk[[kk]][g,1],-1,1)*0.5*log(fk[[kk]][g,11]/fk[[kk]][g,1])
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
      #runif(1,0.0001,0.0002+1*(8-abs(8-MeanKstart[g,kk]))/16    )

    }
    # plot(XX[WW],fk[[kk]][g,])
    # sss = seq(0,16,by=0.05)
    # lines(sss,BiNormalLogT(sss,start = MinKstart[g,kk], ratio1= EtaK1start[g,kk],ratio2= EtaK2start[g,kk], mu= MeanKstart[g,kk], sigma= VarKstart[g,kk]^0.5 ))
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
        #Par5        = rep(0,3)
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
        #Par5        = rep(2,3)
      ),
      VarK   = list(
        type        = "Normal",
        Par1        = rep(0,3),
        Par2        = rep(100,3),
        Par3        = rep(0,3),
        Par4        = rep(20,3),
        Par5        = rep(2,3)
        #Par5        = rep(2,3)
      ),


      RegCoef = list(
          type        = "Normal",
          Par1        = c(0,0,3,8,5), ## #eta1,2,var media
          Par2        = rep(100,5)
      ),
      sigmaMat = list(
          type        = "InverseWishart",
          Par1        = rep(2,3),
          Par2        = c(diag(1,5)) ## diagonale Psi
      ),

      DPpar = list(
          type        = "Gamma",
          Par1        = 1,
          Par2        = 1 ## diagonale Psi
      )
    ) ,
      start   = list(
      RegCoef = unlist(ReCoefStart),
      sigmaMat = unlist(SigmaStart),
      stMat = ValInitclust,
      piMat = rep(0.1,Kmax*10),
      DPpar = rep(1,100),
      Beta0   = beta0start,
      Beta1   = beta1start,
      MinK = MinKstart,
      VarK =VarKstart,
      EtaK1 = EtaK1start,
      EtaK2 = EtaK2start,
      MeanK = MeanKstart,
      SH    = SHstart

    ),
    # MCMC = list(
    #   iter        = 5000*20,
    #   burnin      = 3500*20,
    #   thin        = 1*20
    # ),
    MCMC = list(
      iter        = 15000000,
      burnin      = 1425000,
      thin        = 30
    ),
    #    MCMC = list(
    #   iter        = 1,
    #   burnin      = 0,
    #   thin        = 1
    # ),
    Adapt= list(
      doUpdate    = T,
      start       = 5000,
      end         = 0.9,
      exp         = 0.6,
      acc        = 0.534,
      batch       = 20,
      epsilon     = 0.000000001,
      lambda      = 1,
      sd          = list(
        General     = 0.001,
        SF          = 0.0001

      )
    ),
    IterShow       = 1000 ,
    DoParallel = T,
    Nthreads = Nthreads,
    Verbose = 1,
    saveODE = T,
    Rparam = F,
     Spike  = T,
    LimSpike = 10,
    LimC   = 0.01,
    sample_Var = T,
    sample_Ode = T,
    SampleOdeType = 1,
    abs_err = 10^(-8),
    rel_err = 10^(-4),
    sample_Pi = T,
    sample_DpPar = T,
    sample_st = T,
    sample_RandomPar = T,
    IterStartRandom = ifelse(DoClust==T,50000*1000000,50000),
    addStartRandom  = 0,
    Use_pi = T,
    ExpMolt = 0.00,
    ExpLambda = 10,
    Do_MergeSplit = T,
    TypeK = 2,
    MoltSF = 100,
    A = NULL
)

str(Mod)

save.image(paste(DirOUT,NAME,".Rdata",sep=""))






  WWW2 = nrow(Mod$MinK)
  MinKstart   = matrix((Mod$MinK[WWW2,]),ncol=3,byrow=T)
  EtaK1start  = matrix((Mod$EtaK1[WWW2,]),ncol=3,byrow=T)
  EtaK2start  = matrix((Mod$EtaK2[WWW2,]),ncol=3,byrow=T)
  MeanKstart  = matrix((Mod$MeanK[WWW2,]),ncol=3,byrow=T)
  VarKstart   = matrix((Mod$VarK[WWW2,]),ncol=3,byrow=T)
  beta0start  = Mod$Beta0[WWW2,]
  beta1start  = Mod$Beta1[WWW2,]

  SHstart     = c(Mod$SF_1[WWW2,],Mod$SF_2[WWW2,],Mod$SF_3[WWW2,])
  #beta0start  = rep(-4,3)
  #beta1start  = rep(1.8,3)


    # ValInitclust = rep(1,nG*5*3)

    # igg = 0

    # ModelOut=Mod$ModelOutput_11
    # k = max(ModelOut$PostK)
    # WNWNrow = nrow(ModelOut$ModelOutput[[k]]$Clustering$st)
    # WNWNcol = ncol(ModelOut$ModelOutput[[k]]$Clustering$st)
    # ValInitclust[igg+1:nG] = ModelOut$ModelOutput[[k]]$Clustering$st[WNWNrow,]
    # igg = igg+nG

    # ModelOut=Mod$ModelOutput_12
    # k = max(ModelOut$PostK)
    # WNWNrow = nrow(ModelOut$ModelOutput[[k]]$Clustering$st)
    # WNWNcol = ncol(ModelOut$ModelOutput[[k]]$Clustering$st)
    # ValInitclust[igg+1:nG] = ModelOut$ModelOutput[[k]]$Clustering$st[WNWNrow,]
    # igg = igg+nG

    # ModelOut=Mod$ModelOutput_21
    # k = max(ModelOut$PostK)
    # WNWNrow = nrow(ModelOut$ModelOutput[[k]]$Clustering$st)
    # WNWNcol = ncol(ModelOut$ModelOutput[[k]]$Clustering$st)
    # ValInitclust[igg+1:nG] = ModelOut$ModelOutput[[k]]$Clustering$st[WNWNrow,]
    # igg = igg+nG

    # ModelOut=Mod$ModelOutput_22
    # k = max(ModelOut$PostK)
    # WNWNrow = nrow(ModelOut$ModelOutput[[k]]$Clustering$st)
    # WNWNcol = ncol(ModelOut$ModelOutput[[k]]$Clustering$st)
    # ValInitclust[igg+1:nG] = ModelOut$ModelOutput[[k]]$Clustering$st[WNWNrow,]
    # igg = igg+nG

    # ModelOut=Mod$ModelOutput_31
    # k = max(ModelOut$PostK)
    # WNWNrow = nrow(ModelOut$ModelOutput[[k]]$Clustering$st)
    # WNWNcol = ncol(ModelOut$ModelOutput[[k]]$Clustering$st)
    # ValInitclust[igg+1:nG] = ModelOut$ModelOutput[[k]]$Clustering$st[WNWNrow,]
    # igg = igg+nG

    # ModelOut=Mod$ModelOutput_32
    # k = max(ModelOut$PostK)
    # WNWNrow = nrow(ModelOut$ModelOutput[[k]]$Clustering$st)
    # WNWNcol = ncol(ModelOut$ModelOutput[[k]]$Clustering$st)
    # ValInitclust[igg+1:nG] = ModelOut$ModelOutput[[k]]$Clustering$st[WNWNrow,]
    # igg = igg+nG

  ret = list(
  MinKstart   = MinKstart,
  EtaK1start  = EtaK1start,
  EtaK2start  = EtaK2start,
  MeanKstart  = MeanKstart,
  VarKstart   = VarKstart,
  beta0start  = beta0start,
  beta1start  = beta1start,
  SHstart     = SHstart
  #ValInitclust = ValInitclust
  )

save(ret, file = paste(DirOUT,NAME,"_INIT.Rdata",sep=""))






  WWW2 = nrow(Mod$MinK)
  MinKstart   = matrix((Mod$MinK[WWW2,]),ncol=3,byrow=T)
  EtaK1start  = matrix((Mod$EtaK1[WWW2,]),ncol=3,byrow=T)
  EtaK2start  = matrix((Mod$EtaK2[WWW2,]),ncol=3,byrow=T)
  MeanKstart  = matrix((Mod$MeanK[WWW2,]),ncol=3,byrow=T)
  VarKstart   = matrix((Mod$VarK[WWW2,]),ncol=3,byrow=T)
  beta0start  = Mod$Beta0[WWW2,]
  beta1start  = Mod$Beta1[WWW2,]
  SHstart     = c(Mod$SF_1[WWW2,],Mod$SF_2[WWW2,],Mod$SF_3[WWW2,])
  #beta0start  = rep(-4,3)
  #beta1start  = rep(1.8,3)


    ValInitclust = rep(1,nG*5*3)

    igg = 0

    ModelOut=Mod$ModelOutput_11
    TT = table(ModelOut$PostK)
    k = as.numeric(names(TT)[which.max(TT)])
    #k = max(ModelOut$PostK)
    WNWNrow = nrow(ModelOut$ModelOutput[[k]]$Clustering$st)
    WNWNcol = ncol(ModelOut$ModelOutput[[k]]$Clustering$st)
    ValInitclust[igg+1:nG] = ModelOut$ModelOutput[[k]]$Clustering$st[WNWNrow,]
    igg = igg+nG

    ModelOut=Mod$ModelOutput_12
    TT = table(ModelOut$PostK)
    k = as.numeric(names(TT)[which.max(TT)])
    #k = max(ModelOut$PostK)
    WNWNrow = nrow(ModelOut$ModelOutput[[k]]$Clustering$st)
    WNWNcol = ncol(ModelOut$ModelOutput[[k]]$Clustering$st)
    ValInitclust[igg+1:nG] = ModelOut$ModelOutput[[k]]$Clustering$st[WNWNrow,]
    igg = igg+nG

    ModelOut=Mod$ModelOutput_21
    TT = table(ModelOut$PostK)
    k = as.numeric(names(TT)[which.max(TT)])
    #k = max(ModelOut$PostK)
    WNWNrow = nrow(ModelOut$ModelOutput[[k]]$Clustering$st)
    WNWNcol = ncol(ModelOut$ModelOutput[[k]]$Clustering$st)
    ValInitclust[igg+1:nG] = ModelOut$ModelOutput[[k]]$Clustering$st[WNWNrow,]
    igg = igg+nG

    #
    #
    #
    #
    #
    # ModelOut=Mod$ModelOutput_21
    # TT = table(ModelOut$PostK)
    # k = as.numeric(names(TT)[which.max(TT)])
    # #k = max(ModelOut$PostK)
    # WNWNrow = nrow(ModelOut$ModelOutput[[k]]$Clustering$st)
    # WNWNcol = ncol(ModelOut$ModelOutput[[k]]$Clustering$st)
    # ValInitclust[igg+1:nG] = ModelOut$ModelOutput[[k]]$Clustering$st[WNWNrow,]
    # igg = igg+nG

    ModelOut=Mod$ModelOutput_22
    TT = table(ModelOut$PostK)
    k = as.numeric(names(TT)[which.max(TT)])
    #k = max(ModelOut$PostK)
    WNWNrow = nrow(ModelOut$ModelOutput[[k]]$Clustering$st)
    WNWNcol = ncol(ModelOut$ModelOutput[[k]]$Clustering$st)
    ValInitclust[igg+1:nG] = ModelOut$ModelOutput[[k]]$Clustering$st[WNWNrow,]
    igg = igg+nG

    ModelOut=Mod$ModelOutput_31
    TT = table(ModelOut$PostK)
    k = as.numeric(names(TT)[which.max(TT)])
    #k = max(ModelOut$PostK)
    WNWNrow = nrow(ModelOut$ModelOutput[[k]]$Clustering$st)
    WNWNcol = ncol(ModelOut$ModelOutput[[k]]$Clustering$st)
    ValInitclust[igg+1:nG] = ModelOut$ModelOutput[[k]]$Clustering$st[WNWNrow,]
    igg = igg+nG

    ModelOut=Mod$ModelOutput_32
    TT = table(ModelOut$PostK)
    k = as.numeric(names(TT)[which.max(TT)])
    #k = max(ModelOut$PostK)
    WNWNrow = nrow(ModelOut$ModelOutput[[k]]$Clustering$st)
    WNWNcol = ncol(ModelOut$ModelOutput[[k]]$Clustering$st)
    ValInitclust[igg+1:nG] = ModelOut$ModelOutput[[k]]$Clustering$st[WNWNrow,]
    igg = igg+nG

  ret = list(
  MinKstart   = MinKstart,
  EtaK1start  = EtaK1start,
  EtaK2start  = EtaK2start,
  MeanKstart  = MeanKstart,
  VarKstart   = VarKstart,
  beta0start  = beta0start,
  beta1start  = beta1start,
  SHstart     = SHstart,
  ValInitclust = ValInitclust
  )

save(ret, file = paste(DirOUT,NAME,"_INIT.Rdata",sep=""))


#load("/Users/gianlucamastrantonio/Dropbox (Politecnico di Torino Staff)/SFTPsave/SFTP_ODE/disma/output/RealAll_ng500.Rdata")
pdf(paste(DirOUT,NAME,".pdf",sep=""))

par(mfrow=c(3,2))

plot(Mod$Beta0[,1], type="l")
#abline(h=BetaP[1],coRl=2)
plot(Mod$Beta1[,1], type="l")
#abline(h=BetaP[2],col=2)

plot(Mod$Beta0[,2], type="l")
#abline(h=BetaM[1],col=2)
plot(Mod$Beta1[,2], type="l")
##abline(h=BetaM[2],col=2)

plot(Mod$Beta0[,3], type="l")
#abline(h=BetaC[1],col=2)
plot(Mod$Beta1[,3], type="l")
#abline(h=BetaC[2],col=2)

par(mfrow=c(4,3))
for(ifd in 1:ncol(Mod$SF_3))
{
  plot(Mod$SF_1[,ifd], type="l")
}
par(mfrow=c(4,3))
for(ifd in 1:ncol(Mod$SF_3))
{
  plot(Mod$SF_2[,ifd], type="l")
}
par(mfrow=c(4,3))
for(ifd in 1:ncol(Mod$SF_3))
{
  plot(Mod$SF_3[,ifd], type="l")
}
par(mfrow=c(3,3))
nNonNA = 11
WW = c(1 , 2  ,3  ,4  ,7 ,10, 13, 25, 49,73,97)
indext = seq(0,16,0.5/3)[WW]
indextOBS = indext
WW = 1:11

Ymean = matrix(colMeans(Mod$Y), ncol=3, byrow=T)
Yq1 = matrix(apply(Mod$Y,2,quantile,probs=0.025), ncol=3, byrow=T)
Yq2 = matrix(apply(Mod$Y,2,quantile,probs=0.975), ncol=3, byrow=T)
WWWW =c(WW,WW+max(WW),WW+max(WW)*2)

Y3app = Mod$Y[,seq(3,ncol(Mod$Y),by=3)]
Y3app_1 = Y3app* matrix(Mod$SF_1,nrow=nrow(Mod$SF_1),ncol=11*(dim(Y3app)[2]/11))
Y3app_2 = Y3app* matrix(Mod$SF_2,nrow=nrow(Mod$SF_1),ncol=11*(dim(Y3app)[2]/11))
Y3app_3 = Y3app* matrix(Mod$SF_3,nrow=nrow(Mod$SF_1),ncol=11*(dim(Y3app)[2]/11))

Y3app_1_mean = colMeans(Y3app_1)
Y3app_2_mean = colMeans(Y3app_2)
Y3app_3_mean = colMeans(Y3app_3)
for(g in 1:nG)
{


  par(mfrow=c(3,3))
  plot(Mod$MinK[, (g-1)*3+1],type="l")
  #abline(h = MinK[g,1], col=2)
  plot(Mod$MinK[,(g-1)*3+2], type="l")
  #abline(h = MinK[g,2], col=2)
  plot(Mod$MinK[,(g-1)*3+3], type="l")
  #abline(h = MinK[g,3], col=2)


  plot(Mod$VarK[,(g-1)*3+1],type="l")
  #abline(h = VarK[g,1], col=2)
  plot(Mod$VarK[,(g-1)*3+2], type="l")
  #abline(h = VarK[g,2], col=2)
  plot(Mod$VarK[,(g-1)*3+3], type="l")
  #abline(h = VarK[g,3], col=2)


  plot(Mod$MeanK[,(g-1)*3+1],type="l")
  #abline(h = MeanK[g,1], col=2)
  plot(Mod$MeanK[,(g-1)*3+2], type="l")
  #abline(h = MeanK[g,2], col=2)
  plot(Mod$MeanK[,(g-1)*3+3], type="l")
  # abline(h = MeanK[g,3], col=2)


  plot(Mod$EtaK1[,(g-1)*3+1],type="l")
  #abline(h = EtaK1[g,1], col=2)
  plot(Mod$EtaK1[,(g-1)*3+2], type="l")
  #abline(h = EtaK1[g,2], col=2)
  plot(Mod$EtaK1[,(g-1)*3+3], type="l")
  #abline(h = EtaK1[g,3], col=2)


  plot(Mod$EtaK2[,(g-1)*3+1],type="l")
  #abline(h = EtaK2[g,1], col=2)
  plot(Mod$EtaK2[,(g-1)*3+2], type="l")
  #abline(h = EtaK2[g,2], col=2)
  plot(Mod$EtaK2[,(g-1)*3+3], type="l")
  #abline(h = EtaK2[g,3], col=2)


  plot(Mod$Y0[,(g-1)*3+1],type="l")
  #abline(h = PinitSIM[g], col=2)
  plot(Mod$Y0[,(g-1)*3+2], type="l")
  #abline(h = MinitSIM[g], col=2)
  plot(Mod$Y0[,(g-1)*3+3], type="l")
  #abline(h = CinitSIM[g], col=2)


  jjj = g
  par(mfrow=c(3,2))
  W = 1:nNonNA+(jjj-1)*nNonNA

  j = 1
  plot(indextOBS,Ymean[W,j], type="l", ylim= c(min(DATA[,j][WWWW+33*(jjj-1)],Yq1[W,j]), max(DATA[,j][WWWW+33*(jjj-1)],Yq2[W,j])  ), main=paste("Gene_",jjj, sep="") )
  lines(indextOBS,Yq1[W,j], lty=2, col=3)
  lines(indextOBS,Yq2[W,j], lty=2, col=3)
  points(indext[c(WW,WW,WW)],DATA[,j][WWWW+33*(jjj-1)],col=2, pch=20, cex=0.7)



  j = 2
  plot(indextOBS,Ymean[W,j], type="l", ylim= c(min(DATA[,j][WWWW+33*(jjj-1)],Yq1[W,j]), max(DATA[,j][WWWW+33*(jjj-1)],Yq2[W,j])  ))
  lines(indextOBS,Yq1[W,j], lty=2, col=3)
  lines(indextOBS,Yq2[W,j], lty=2, col=3)
  points(indext[c(WW,WW,WW)],DATA[,j][WWWW+33*(jjj-1)],col=2, pch=20, cex=0.7)

  j = 3
  TTTT = c( c(Y3app_1_mean[W]),c(Y3app_2_mean[W]),c(Y3app_3_mean[W]) )
  #plot(indextOBS,Ymean[W,j], type="l", ylim= c(min(DATA[,j][WWWW+33*(jjj-1)],Yq1[W,j]), max(DATA[,j][WWWW+33*(jjj-1)],Yq2[W,j])  ))
  plot(indextOBS,Ymean[W,j], type="l", ylim= c(min(DATA[,j][WWWW+33*(jjj-1)]), max(DATA[,j][WWWW+33*(jjj-1)])  ))
  #plot(indextOBS,Ymean[W,j], type="l", ylim= c(min(TTTT), max(TTTT)  ))
  lines(indextOBS,Yq1[W,j], lty=2, col=3)
  lines(indextOBS,Yq2[W,j], lty=2, col=3)
  points(indext[c(WW,WW,WW)],DATA[,j][WWWW+33*(jjj-1)],col=2, pch=20, cex=0.7)
  lines(indextOBS,Y3app_1_mean[W], col=2, lty=2)
  lines(indextOBS,Y3app_2_mean[W], col=3, lty=2)
  lines(indextOBS,Y3app_3_mean[W], col=4, lty=2)


  k1mat  = matrix(nrow=nrow(Mod$MinK), ncol=11)
  k2mat  = matrix(nrow=nrow(Mod$MinK), ncol=11)
  k3mat  = matrix(nrow=nrow(Mod$MinK), ncol=11)
  W = (jjj-1)*3
  for(i in 1:nrow(Mod$MinK))
  {
    for(k in 1:3)
    {
      XX = BiNormalLogT(indextOBS,start = Mod$MinK[i,W+k], ratio1= Mod$EtaK1[i,W+k],ratio2= Mod$EtaK2[i,W+k], mu= Mod$MeanK[i,W+k], sigma= Mod$VarK[i,W+k]^0.5,typek=typek )



      if(k==1)
      {
        k1mat[i,] = XX
      }
      if(k==2)
      {
        k2mat[i,] = XX
      }
      if(k==3)
      {
        k3mat[i,] = XX
      }
    }


  }


  for(k in 1:3)
  {
    if(k==1)
    {
      XX = k1mat
      #ttc = Tk1
    }
    if(k==2)
    {
      XX = k2mat
      #ttc = Tk2
    }
    if(k==3)
    {
      XX = k3mat
      #ttc = Tk3

    }
    XXM = matrix(colMeans(XX), ncol=11, byrow=T)
    XX1 = matrix(apply(XX,2,quantile,probs=0.025), ncol=11, byrow=T)
    XX2 = matrix(apply(XX,2,quantile,probs=0.975), ncol=11, byrow=T)

    plot(indextOBS,XXM[1,], type="l", ylim=c(min(XX1), max(XX2)),main=paste("Gene_",jjj, sep=""))
    lines(indextOBS,XX1[1,], lty=2, col=3)
    lines(indextOBS,XX2[1,], lty=2, col=3)
    points(indextOBS,fk[[k]][jjj,], col=2)
    #lines(indextOBS,ttc[!is.na(ttc[,jjj]),jjj], col=2)
  }



}

dev.off()



#save.image(paste(DirOUT,NAME,".Rdata",sep=""))




for(ic in 1:2)
{
  for(ik in 1:3)
  {
    if(ic==1)
    {
        if(ik==1)
        {
            ModelOut=Mod$ModelOutput_11
            NAME_c = paste(NAME,"_11")
        }
        if(ik==2)
        {
            ModelOut=Mod$ModelOutput_21
            NAME_c = paste(NAME,"_21")
        }
        if(ik==3)
        {
            ModelOut=Mod$ModelOutput_31
            NAME_c = paste(NAME,"_31")
        }
    }
    if(ic==2)
    {
        if(ik==1)
        {
            ModelOut=Mod$ModelOutput_12
            NAME_c = paste(NAME,"_12")
        }
        if(ik==2)
        {
            ModelOut=Mod$ModelOutput_22
            NAME_c = paste(NAME,"_22")
        }
        if(ik==3)
        {
           ModelOut=Mod$ModelOutput_32
           NAME_c = paste(NAME,"_32")
        }


    }


    W = as.numeric(rownames(table(ModelOut$PostK))[(table(ModelOut$PostK)>20)])
    for(k in W)
    {
        pdf(paste(DirOUT,NAME_c,"_K=",k,"NoStand.pdf",sep=""))







        par(mfrow=c(1,2))
        barplot(table(ModelOut$PostK), main="Posterior of K")
        plot(ModelOut$PostK, type="l")
        par(mfrow=c(1,1))
        barplot(table(ModelOut$ModelOutput[[k]]$Clustering$MAP), main="MAP")

        if(ic==1)
        {

          rmnorm=function(n = 1, mean = rep(0, d), varcov)
          {
             d <- if (is.matrix(varcov))
                 ncol(varcov)
             else 1
             z <- matrix(rnorm(n * d), n, d) %*% chol(varcov)
             y <- t(mean + t(z))
             return(y)
          }

          par(mfrow=c(1,2))
          for(kk in 1:k)
          {
            WW = (kk-1)*(4)+1:4
              WW2 = (kk-1)*(4^2)+1:16
            BETA = ModelOut$ModelOutput[[k]]$ LikelihoodParameters$Beta[,WW]
            SIGMA = ModelOut$ModelOutput[[k]]$ LikelihoodParameters$Sigma[,WW2]

            kmat  = matrix(nrow=nrow(ModelOut$ModelOutput[[k]]$ LikelihoodParameters$Beta), ncol=11)
            for(i in 1:nrow(ModelOut$ModelOutput[[k]]$ LikelihoodParameters$Beta))
            {


                SS = rmnorm(1,BETA[i,], matrix(SIGMA[i,], ncol=4))

                XX = BiNormalLogT(indextOBS,start = 0.1, ratio1= SS[1],ratio2= SS[2], mu= exp(SS[4]), sigma= exp(SS[3])^0.5,typek=typek )

                #sdd = ifelse(sd(XX)==0,1,sd(XX))
                #kmat[i,] = (XX-mean(XX))/sdd
                kmat[i,] = XX

            }
            XX = kmat
            #XXM = matrix(colMeans(XX), ncol=11, byrow=T)
            #  XX1 = matrix(apply(XX,2,quantile,probs=0.025), ncol=11, byrow=T)
            #  XX2 = matrix(apply(XX,2,quantile,probs=0.975), ncol=11, byrow=T)

              #plot(indextOBS,XXM[1,], type="l", ylim=c(min(XX1), max(XX2)))
              #lines(indextOBS,XX1[1,], lty=2, col=3)
              #lines(indextOBS,XX2[1,], lty=2, col=3)
          }






          ###






          BB = which(ModelOut$PostK==k)
          MinKm = colMeans(Mod$MinK[BB,])
          VarKm = colMeans(Mod$VarK[BB,])
          MeanKm = colMeans(Mod$MeanK[BB,])
          EtaK1m = colMeans(Mod$EtaK1[BB,])
          EtaK2m = colMeans(Mod$EtaK2[BB,])

          COL = ModelOut$ModelOutput[[k]]$Clustering$MAP
          par(mfrow=c(2,2))
          plot(MinKm[seq(1,nG*3, by=3)+ik-1],VarKm[seq(1,nG*3, by=3)+ik-1], col=COL)
          plot(MinKm[seq(1,nG*3, by=3)+ik-1],MeanKm[seq(1,nG*3, by=3)+ik-1], col=COL)
          plot(MinKm[seq(1,nG*3, by=3)+ik-1],EtaK1m[seq(1,nG*3, by=3)+ik-1], col=COL)
          plot(MinKm[seq(1,nG*3, by=3)+ik-1],EtaK2m[seq(1,nG*3, by=3)+ik-1], col=COL)

          plot(VarKm[seq(1,nG*3, by=3)+ik-1],MeanKm[seq(1,nG*3, by=3)+ik-1], col=COL)
          plot(VarKm[seq(1,nG*3, by=3)+ik-1],EtaK1m[seq(1,nG*3, by=3)+ik-1], col=COL)
          plot(VarKm[seq(1,nG*3, by=3)+ik-1],EtaK2m[seq(1,nG*3, by=3)+ik-1], col=COL)

          plot(MeanKm[seq(1,nG*3, by=3)+ik-1],EtaK1m[seq(1,nG*3, by=3)+ik-1], col=COL)
          plot(MeanKm[seq(1,nG*3, by=3)+ik-1],EtaK2m[seq(1,nG*3, by=3)+ik-1], col=COL)

          plot(EtaK1m[seq(1,nG*3, by=3)+ik-1],EtaK2m[seq(1,nG*3, by=3)+ik-1], col=COL)

          par(mfrow=c(1,1))
          plot(VarKm[seq(1,nG*3, by=3)+ik-1],MeanKm[seq(1,nG*3, by=3)+ik-1], col=COL)
          plot(EtaK1m[seq(1,nG*3, by=3)+ik-1],EtaK2m[seq(1,nG*3, by=3)+ik-1], col=COL)

      KMAT = matrix(0,nrow=ncol(Mod$MinK), ncol=11)
      for(ig in 1:nG)
      {
        for(iapp in 1:length(BB))
        {
          Wapp = (ig-1)*3

          i = BB[iapp]
          KMAT[ig, ] = KMAT[ig, ]+BiNormalLogT(indextOBS,start = 0.1, ratio1= Mod$EtaK1[i,Wapp+ik],ratio2= Mod$EtaK2[i,Wapp+ik], mu= Mod$MeanK[i,Wapp+ik], sigma= Mod$VarK[i,Wapp+ik]^0.5,typek=typek )
        }
      }
      KMAT = KMAT/length(BB)
      # for(ig in 1:nG)
      # {
      #   sdd = ifelse(sd(c(KMAT[ig,]))==0,1,sd(c(KMAT[ig,])))

      #   KMAT[ig,] = (c(KMAT[ig,])-mean(c(KMAT[ig,])))/sdd
      # }

          for(kss in 1:k)
          {
            par(mfrow=c(1,1))

            Ak = which(ModelOut$ModelOutput[[k]]$Clustering$MAP==kss)
            if(length(Ak)>0)
            {



              par(mfrow=c(1,1))
              plot(0,0, xlim=c(0,16), ylim=c(range( c( KMAT[Ak,] ) )), type="n", main=paste("K=",kss))
              for(iii in 1:length(Ak))
              {
                WAk = Ak[iii]
                lines(indextOBS,KMAT[WAk,], col=iii, type="l")
              }
              par(mfrow=c(3,3))
              for(iii in 1:length(Ak))
              {
                WAk = Ak[iii]
                plot(indextOBS,KMAT[WAk,], col=1+kss, type="l", main=paste("Gene ",WAk))
              }
            }




          }

        }








        par(mfrow=c(1,1))
        plot(ModelOut$ModelOutput[[k]]$Clustering$MAP, type="l")

        par(mfrow=c(4,2))
        plot(ModelOut$ModelOutput[[k]]$LikelihoodParameters$Beta)
         par(mfrow=c(4,2))
        plot(ModelOut$ModelOutput[[k]]$LikelihoodParameters$Sigma)
        plot(ModelOut$ModelOutput[[k]]$Clustering$"pi")
        plot(ModelOut$ModelOutput[[k]]$Clustering$DPpar)
        #plot(ModelOut$Clustering$st)



        dev.off()
    }



  }
}

#save.image(paste(DirOUT,NAME,".Rdata",sep=""))





# for(ic in 1:2)
# {
#   for(ik in 1:3)
#   {
#     if(ic==1)
#     {
#         if(ik==1)
#         {
#             ModelOut=Mod$ModelOutput_11
#             NAME_c = paste(NAME,"_11")
#         }
#         if(ik==2)
#         {
#             ModelOut=Mod$ModelOutput_21
#             NAME_c = paste(NAME,"_21")
#         }
#         if(ik==3)
#         {
#             ModelOut=Mod$ModelOutput_31
#             NAME_c = paste(NAME,"_31")
#         }
#     }
#     if(ic==2)
#     {
#         if(ik==1)
#         {
#             ModelOut=Mod$ModelOutput_12
#             NAME_c = paste(NAME,"_12")
#         }
#         if(ik==2)
#         {
#             ModelOut=Mod$ModelOutput_22
#             NAME_c = paste(NAME,"_22")
#         }
#         if(ik==3)
#         {
#            ModelOut=Mod$ModelOutput_32
#            NAME_c = paste(NAME,"_32")
#         }


#     }


#     W = as.numeric(rownames(table(ModelOut$PostK))[(table(ModelOut$PostK)>20)])
#     for(k in W)
#     {
#         pdf(paste(DirOUT,NAME_c,"_K=",k,"Post_LS_NoStand.pdf",sep=""))







#         par(mfrow=c(1,2))
#         barplot(table(ModelOut$PostK), main="Posterior of K")
#         plot(ModelOut$PostK, type="l")
#         par(mfrow=c(1,1))
#         barplot(table(ModelOut$ModelOutput[[k]]$Clustering$MAP), main="MAP")

#         if(ic==1)
#         {

#           rmnorm=function(n = 1, mean = rep(0, d), varcov)
#           {
#              d <- if (is.matrix(varcov))
#                  ncol(varcov)
#              else 1
#              z <- matrix(rnorm(n * d), n, d) %*% chol(varcov)
#              y <- t(mean + t(z))
#              return(y)
#           }

#           par(mfrow=c(1,2))
#           for(kk in 1:k)
#           {
#             WW = (kk-1)*(4)+1:4
#               WW2 = (kk-1)*(4^2)+1:16
#             BETA = ModelOut$ModelOutput[[k]]$ LikelihoodParameters$Beta[,WW]
#             SIGMA = ModelOut$ModelOutput[[k]]$ LikelihoodParameters$Sigma[,WW2]

#             kmat  = matrix(nrow=nrow(ModelOut$ModelOutput[[k]]$ LikelihoodParameters$Beta), ncol=11)
#             for(i in 1:nrow(ModelOut$ModelOutput[[k]]$ LikelihoodParameters$Beta))
#             {


#                 SS = rmnorm(1,BETA[i,], matrix(SIGMA[i,], ncol=4))

#                 XX = BiNormalLogT(indextOBS,start = 0.1, ratio1= SS[1],ratio2= SS[2], mu= exp(SS[4]), sigma= exp(SS[3])^0.5 )

#                 #sdd = ifelse(sd(XX)==0,1,sd(XX))
#                 #kmat[i,] = (XX-mean(XX))/sdd
#                 kmat[i,] = XX

#             }
#             XX = kmat
#             XXM = matrix(colMeans(XX), ncol=11, byrow=T)
#               XX1 = matrix(apply(XX,2,quantile,probs=0.025), ncol=11, byrow=T)
#               XX2 = matrix(apply(XX,2,quantile,probs=0.975), ncol=11, byrow=T)

#               plot(indextOBS,XXM[1,], type="l", ylim=c(min(XX1), max(XX2)))
#               lines(indextOBS,XX1[1,], lty=2, col=3)
#               lines(indextOBS,XX2[1,], lty=2, col=3)
#           }






#           ###






#           BB = which(ModelOut$PostK==k)
#           MinKm = colMeans(Mod$MinK[BB,])
#           VarKm = colMeans(Mod$VarK[BB,])
#           MeanKm = colMeans(Mod$MeanK[BB,])
#           EtaK1m = colMeans(Mod$EtaK1[BB,])
#           EtaK2m = colMeans(Mod$EtaK2[BB,])

#           COL = ModelOut$ModelOutput[[k]]$Clustering$MAP
#           par(mfrow=c(2,2))
#           plot(MinKm[seq(1,nG*3, by=3)+ik-1],VarKm[seq(1,nG*3, by=3)+ik-1], col=COL)
#           plot(MinKm[seq(1,nG*3, by=3)+ik-1],MeanKm[seq(1,nG*3, by=3)+ik-1], col=COL)
#           plot(MinKm[seq(1,nG*3, by=3)+ik-1],EtaK1m[seq(1,nG*3, by=3)+ik-1], col=COL)
#           plot(MinKm[seq(1,nG*3, by=3)+ik-1],EtaK2m[seq(1,nG*3, by=3)+ik-1], col=COL)

#           plot(VarKm[seq(1,nG*3, by=3)+ik-1],MeanKm[seq(1,nG*3, by=3)+ik-1], col=COL)
#           plot(VarKm[seq(1,nG*3, by=3)+ik-1],EtaK1m[seq(1,nG*3, by=3)+ik-1], col=COL)
#           plot(VarKm[seq(1,nG*3, by=3)+ik-1],EtaK2m[seq(1,nG*3, by=3)+ik-1], col=COL)

#           plot(MeanKm[seq(1,nG*3, by=3)+ik-1],EtaK1m[seq(1,nG*3, by=3)+ik-1], col=COL)
#           plot(MeanKm[seq(1,nG*3, by=3)+ik-1],EtaK2m[seq(1,nG*3, by=3)+ik-1], col=COL)

#           plot(EtaK1m[seq(1,nG*3, by=3)+ik-1],EtaK2m[seq(1,nG*3, by=3)+ik-1], col=COL)

#           par(mfrow=c(1,1))
#           plot(VarKm[seq(1,nG*3, by=3)+ik-1],MeanKm[seq(1,nG*3, by=3)+ik-1], col=COL)
#           plot(EtaK1m[seq(1,nG*3, by=3)+ik-1],EtaK2m[seq(1,nG*3, by=3)+ik-1], col=COL)

#       KMAT = matrix(0,nrow=ncol(Mod$MinK), ncol=11)
#       for(ig in 1:nG)
#       {
#         for(iapp in 1:length(BB))
#         {
#           Wapp = (ig-1)*3

#           i = BB[iapp]
#           KMAT[ig, ] = KMAT[ig, ]+BiNormalLogT(indextOBS,start = 0.1, ratio1= Mod$EtaK1[i,Wapp+ik],ratio2= Mod$EtaK2[i,Wapp+ik], mu= Mod$MeanK[i,Wapp+ik], sigma= Mod$VarK[i,Wapp+ik]^0.5 )
#         }
#       }
#       KMAT = KMAT/length(BB)
#       # for(ig in 1:nG)
#       # {
#       #   sdd = ifelse(sd(c(KMAT[ig,]))==0,1,sd(c(KMAT[ig,])))

#       #   KMAT[ig,] = (c(KMAT[ig,])-mean(c(KMAT[ig,])))/sdd
#       # }

#           for(kss in 1:k)
#           {
#             par(mfrow=c(1,1))

#             Ak = which(ModelOut$ModelOutput[[k]]$Clustering$MAP==kss)
#             if(length(Ak)>0)
#             {



#               par(mfrow=c(1,1))
#               plot(0,0, xlim=c(0,16), ylim=c(range( c( KMAT[Ak,] ) )), type="n", main=paste("K=",kss))
#               for(iii in 1:length(Ak))
#               {
#                 WAk = Ak[iii]
#                 lines(indextOBS,KMAT[WAk,], col=iii, type="l")
#               }
#               par(mfrow=c(3,3))
#               for(iii in 1:length(Ak))
#               {
#                 WAk = Ak[iii]
#                 plot(indextOBS,KMAT[WAk,], col=1+kss, type="l", main=paste("Gene ",WAk))
#               }
#             }




#           }

#         }








#         par(mfrow=c(1,1))
#         plot(ModelOut$ModelOutput[[k]]$Clustering$MAP, type="l")

#         par(mfrow=c(4,2))
#         plot(ModelOut$ModelOutput[[k]]$LikelihoodParameters$Beta)
#          par(mfrow=c(4,2))
#         plot(ModelOut$ModelOutput[[k]]$LikelihoodParameters$Sigma)
#         plot(ModelOut$ModelOutput[[k]]$Clustering$"pi")
#         plot(ModelOut$ModelOutput[[k]]$Clustering$DPpar)
#         #plot(ModelOut$Clustering$st)



#         dev.off()
#     }



#   }
# }





# save.image(paste(DirOUT,NAME,".Rdata",sep=""))
