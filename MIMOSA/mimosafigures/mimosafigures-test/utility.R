#'setup the environment
#'
#'load MIMOSA, parallel, filehash,gdata,pracma,reshape and ggplot2
setup<-function(){
  suppressPackageStartupMessages({require(MIMOSA)
                                  require(multicore)
                                  require(parallel)
                                  require(filehash)
                                  require(gdata)
                                  require(pracma)
                                  require(ggplot2)
                                  require(reshape)
                                  require(data.table)
                                  require(SingleCellAssay)
                                  require(gplots)
                                  
  })
}

#'construct an expected vs observed fdr curve
#'
#'Expected vs observed FDr taking the expected fdr and the truth.
fdrCurves<-function(fdr,truth){
  truth<-truth[order(fdr)]
  true.fdr<-cumsum(1-truth)/1:length(truth)
  fdr<-sort(fdr)
  cbind(fdr,true.fdr)
}

#'Construct and ROC curve
#'
#'Given the p-value and the truth, construct an ROC curve.
rocCurves<-function(p,truth){
  #s<-seq(0,1,l=1000)
  #table<-t(sapply(s,function(th){
  #  prop.table(table(test=factor(p<=th,levels=c("FALSE","TRUE")),truth=truth),margin=2)["TRUE",]
  #}))
  #colnames(table)<-c("FPR","TPR")
  #table
  o<-order(p,decreasing=FALSE)
  intervals<-lapply(sapply(sort(unique(p)),function(x)which(as.numeric(p[o])%in%x)),max)
  rates<-lapply(intervals,function(x)data.frame(TPR=sum(truth[o][1:x])/sum(truth),FPR=sum(!truth[o][1:x])/sum(!truth)))
  int.len<-diff(c(0,do.call(c,intervals)))
  do.call(rbind,lapply(1:length(int.len),function(i)data.frame(TPR=rep(rates[[i]][,1],int.len[[i]]),FPR=rep(rates[[i]][,2],int.len[[i]]))))
  #data.frame(TPR=cumsum(truth[o])/sum(truth),FPR=cumsum(!truth[o])/sum(!truth))
}

#'read the HVTN065 data set
readHVTNData<-function(){
  suppressWarnings(data<-fread("data/hvtn065/e065icf_fh_p.csv",stringsAsFactors=FALSE))
  suppressWarnings(rx<-fread("data/hvtn065/rx_v2.csv"))
  setnames(rx,tolower(colnames(rx[,Ptid:=factor(gsub("\\-","",Ptid))])))
  data[,ptid:=factor(ptid)]
  return(merge(rx,data,by=c("protocol","ptid")))
}

#'Construct an expression set
processHVTNData<-function(data){
  setnames(data,toupper(colnames(data)))
  data<-data[VISITNO%in%c("2","12")]
  data[,CYTNUM:=as.integer(mean(CYTNUM)),by=c("ANTIGEN,CYTOKINE,TCELLSUB,ASSAYID,PTID,VISITNO")]
  data[,NSUB:=as.integer(mean(NSUB)),by=c("ANTIGEN,CYTOKINE,TCELLSUB,ASSAYID,PTID,VISITNO")]
  setkeyv(data,c("ANTIGEN","CYTOKINE","TCELLSUB","ASSAYID","PTID","VISITNO"))
  data<-unique(data)
  data<-melt(data,measure=c("CYTNUM","NSUB"),id=c("PTID","PUB_ID","ANTIGEN","CYTOKINE","ASSAYID","VISITNO","RELIABLE","TCELLSUB","RX_CODE"))
  data<-(cast(data,PTID+PUB_ID+ANTIGEN+ASSAYID+RELIABLE+VISITNO+TCELLSUB+RX_CODE~CYTOKINE+variable))
  
  #drop CMV
  data<-subset(data,!ANTIGEN%in%c("CMV","sebctrl"))

  
  #marginalize over IFNg, IL2, and TNFa
  data<-data.table(within(data, cbind(IFNg_CYTNUM<-`IFN+IL2-TNF-_CYTNUM`+`IFN+IL2-TNF+_CYTNUM`+`IFN+IL2+TNF-_CYTNUM`+`IFN+IL2+TNF+_CYTNUM`,
                                    IFNg_NSUB<-`IFN-IL2-TNF+_NSUB`+`IFN-IL2-TNF+_CYTNUM`+`IFN-IL2+TNF-_CYTNUM`+`IFN-IL2+TNF+_CYTNUM`,
                                    IL2_CYTNUM<-`IFN-IL2+TNF-_CYTNUM`+`IFN-IL2+TNF+_CYTNUM`+`IFN+IL2+TNF-_CYTNUM`+`IFN+IL2+TNF+_CYTNUM`,
                                    IL2_NSUB<-`IFN-IL2-TNF+_NSUB`+`IFN-IL2-TNF+_CYTNUM`+`IFN+IL2-TNF-_CYTNUM`+`IFN+IL2-TNF+_CYTNUM`,
                                    TNF_CYTNUM<-`IFN-IL2-TNF+_CYTNUM`+`IFN-IL2+TNF+_CYTNUM`+`IFN+IL2-TNF+_CYTNUM`+`IFN+IL2+TNF+_CYTNUM`,
                                    TNF_NSUB<-`IFN-IL2-TNF+_NSUB`+`IFN-IL2+TNF-_CYTNUM`+`IFN+IL2-TNF-_CYTNUM`+`IFN+IL2+TNF-_CYTNUM`,
                                    IFNgIL2_CYTNUM<-`IFN+IL2+TNF+_CYTNUM`+`IFN+IL2+TNF-_CYTNUM`,
                                    IFNgIL2_NSUB<-`IFN-IL2-TNF+_CYTNUM`+`IFN-IL2-TNF+_NSUB`+`IFN+IL2-TNF+_CYTNUM`+`IFN+IL2-TNF-_CYTNUM`+`IFN-IL2+TNF+_CYTNUM`+`IFN-IL2+TNF-_CYTNUM`,
                                    IL2TNF_CYTNUM<-`IFN-IL2+TNF+_CYTNUM`+`IFN+IL2+TNF+_CYTNUM`,
                                    IL2TNF_NSUB<-`IFN-IL2-TNF+_CYTNUM`+`IFN-IL2-TNF+_NSUB`+`IFN+IL2-TNF+_CYTNUM`+`IFN-IL2+TNF-_CYTNUM`+`IFN+IL2+TNF-_CYTNUM`+`IFN+IL2-TNF-_CYTNUM`,
                                    IFNgTNF_CYTNUM<-`IFN+IL2-TNF+_CYTNUM`+`IFN+IL2+TNF+_CYTNUM`,
                                    IFNgTNF_NSUB<-`IFN-IL2+TNF-_CYTNUM`+`IFN-IL2-TNF+_NSUB`+`IFN-IL2-TNF+_CYTNUM`+`IFN-IL2+TNF+_CYTNUM`+`IFN+IL2+TNF-_CYTNUM`+`IFN+IL2-TNF-_CYTNUM`)))
  
  data<-(melt(data,id=c("PTID","PUB_ID","ANTIGEN","ASSAYID","TCELLSUB","RX_CODE","RELIABLE","VISITNO")))
  
  data<-cbind(data,colsplit(data$variable,"_",names=c("CYTOKINE","component")))
  data$variable<-NULL
  setnames(data,"component","variable")
  data$VISITNO<-factor(data$VISITNO)
  data$PTID<-factor(data$PTID)
  data<-recast(data,id.var=c("PTID","PUB_ID","ANTIGEN","ASSAYID","TCELLSUB","RX_CODE","RELIABLE","VISITNO","CYTOKINE","variable"),PTID+PUB_ID+ANTIGEN+ASSAYID+TCELLSUB+RX_CODE+VISITNO+CYTOKINE+RELIABLE~variable,value="value")
  data<-subset(data,RELIABLE%in%"Y")
  h<-huber(with(data,CYTNUM/NSUB))
  data<-data[((with(data,CYTNUM/NSUB)-h$mu)/h$s)<30,]
  data<-within(data,NSUB<-NSUB-CYTNUM)
  return(data)
}

constructEset<-function(d){
  E<-ConstructMIMOSAExpressionSet(d,reference=ANTIGEN%in%"negctrl",measure.columns=c("CYTNUM","NSUB"),other.annotations=c("ANTIGEN","CYTOKINE","TCELLSUB","VISITNO","ASSAYID","PTID","PUB_ID","RX_CODE"),default.cast.formula=component~ASSAYID+VISITNO+PTID+RX_CODE+PUB_ID+CYTOKINE+TCELLSUB+ANTIGEN+RefTreat,.variables=.(VISITNO,PTID,ASSAYID,TCELLSUB,CYTOKINE),featureCols=1)
  return(E)
}

#extract the results into a table
summarizeMIMOSAFit <- function (result) {
  hvtn.summary<-(do.call(rbind,lapply(result,function(x){if(!inherits(x,"logical")){
    if(class(x@result)=="MDMixResult"){
      colnames(x@z)<-c("z1","z2")
      cbind(Method="EM",NSUB=x@result@data$n.stim[,1],CYTNUM=x@result@data$n.stim[,2],NSUB_REF=x@result@data$n.unstim[,1],CYTNUM_REF=x@result@data$n.unstim[,2],effect_size=prop.table(as.matrix(x@result@data$n.stim),1)[,2]-prop.table(as.matrix(x@result@data$n.unstim),1)[,2],pData(x@result),x@z)
    }else if(class(x@result)=="MCMCResult"){
      colnames(x@z)<-c("z1","z2")
      cbind(Method="MCMC",NSUB=x@result@n.stim[,1],CYTNUM=x@result@n.stim[,2],NSUB_REF=x@result@n.unstim[,1],CYTNUM_REF=x@result@n.unstim[,2],effect_size=prop.table(as.matrix(x@result@n.stim),1)[,2]-prop.table(as.matrix(x@result@n.unstim),1)[,2],pData(x@result),x@z)
    }
  }})))
  return(hvtn.summary)
}

#Computes  Fisher's exact test on the summarized data
FisherFromSummary <- function (mysummary) {
  suppressPackageStartupMessages(require(doParallel))
  registerDoParallel()
  fisher.p<-foreach(i=1:nrow(mysummary),.combine=c)%dopar%{
    fisher.test(rbind(as.matrix(mysummary[i,c("NSUB_REF","CYTNUM_REF")]),as.matrix(mysummary[i,c("NSUB","CYTNUM")])),alternative="greater")$p.value
  }
  mysummary$fisher.p<-fisher.p
  return(mysummary)
}


#quick ROC curve
rocPlot <- function (mysummary,usecol="fisher.p",method="Fisher") {
  mysummary$truth<-(grepl("T",mysummary$RX_CODE)&mysummary$VISITNO%in%"12")&mysummary$effect_size>0
    
  o<-order(mysummary[,usecol],decreasing=FALSE)
  this<-data.frame(mysummary$CYTOKINE,mysummary$TCELLSUB,Method=method,TPR=cumsum(mysummary[o,"truth"])/sum(mysummary$truth),FPR=cumsum(!mysummary[o,"truth"])/sum(!mysummary$truth))
  #if(include.fisher){
  #  o<-order(mysummary$fisher.p,decreasing=FALSE)
  #  fisher<-data.frame(mysummary$CYTOKINE,mysummary$TCELLSUB,Method="Fisher",TPR=cumsum(mysummary[o,"truth"])/sum(mysummary$truth),FPR=cumsum(!mysummary[o,"truth"])/sum(!mysummary$truth))
  #  this<-rbind(fisher,this)
  #}
  return(this)
}


quickROCs <- function (mysummary,usecol="fisher.p",method="Fisher") {
  do.call(rbind,by(mysummary,factor(mysummary$CYTOKINE:mysummary$TCELLSUB),function(x)rocPlot(x,usecol=usecol,method=method)))
}


#compute the FDR curves
fdrCurve<-function(mysummary,method="Fisher",usecol="fisher.p",adjust=p.adjust(usecol,"fdr")){
  T<-mysummary[,usecol]
  usecol<-T
  T<-eval(parse(text=deparse(substitute(adjust))))
  truth<-(grepl("T",mysummary$RX_CODE)&mysummary$VISITNO%in%"12")&mysummary$effect_size>0
  this<-data.frame(Method=method,fdrCurves(T,  truth))
#   if(include.fisher){
#     fisher<-data.frame(Method="Fisher",fdrCurves(p.adjust(mysummary$fisher.p,"fdr"),  (grepl("T",mysummary$RX_CODE)&mysummary$VISITNO%in%"12")&mysummary$effect_size>0))
#     this<-rbind(fisher,this)
#   }
  return(this)
}

#LRT for MIMOSA
LRT <- function (mysummary,alternative="one.sided") {
  mysummary$lrt.p<-apply(mysummary[,c("NSUB","CYTNUM","NSUB_REF","CYTNUM_REF")],1,function(x){
    ps<-prop.table(x[1:2])
    pu<-prop.table(x[3:4])
    la<-dmultinom(x[1:2],prob=ps,log=TRUE)+dmultinom(x[3:4],prob=pu,log=TRUE)
    ln<-dmultinom(x[1:2],prob=pu,log=TRUE)+dmultinom(x[3:4],prob=pu,log=TRUE)
#     browser()
    p<-pchisq(-2*(ln-la),df=1,lower.tail=FALSE)
    if(alternative%in%"one.sided"&&ps[2]<pu[2]){
      p<-1
    }
    p
  })
  mysummary
}



ComputeROCs <- function (mcmc,em) {
  #Compute ROC curves
  ROC<-rbind(quickROCs(subset(mcmc,TCELLSUB%in%"cd3+/cd4+"),usecol="fisher.p",method="Fisher"),
             quickROCs(subset(mcmc,TCELLSUB%in%"cd3+/cd4+"),usecol="z1",method="MIMOSA (MCMC)"),
             quickROCs(subset(em,TCELLSUB%in%"cd3+/cd4+"),usecol="z1",method="MIMOSA (EM)"),
             quickROCs(subset(mcmc,TCELLSUB%in%"cd3+/cd4+"),usecol="lrt.p",method="LRT"),
             quickROCs(subset(mcmc,TCELLSUB%in%"cd3+/cd4+"),usecol="fold",method="RankFoldChange"))
  setnames(ROC,c("mysummary.CYTOKINE","mysummary.TCELLSUB"),c("CYTOKINE","TCELLSUB"))
  ROC
}

ComputeFDRs <- function (mcmc, em) {
  
  #Compute FDR curves
  FDR<-rbind(
    do.call(rbind,
            by(mcmc,factor(mcmc$CYTOKINE:mcmc$TCELLSUB),function(x)data.frame(x$CYTOKINE,x$TCELLSUB,fdrCurve(x,method="Fisher",usecol="fisher.p",adjust=p.adjust(usecol,"fdr"))))),
    do.call(rbind,
            by(mcmc,factor(mcmc$CYTOKINE:mcmc$TCELLSUB),function(x)data.frame(x$CYTOKINE,x$TCELLSUB,fdrCurve(x,method="MIMOSA (MCMC)",usecol=c("z1","z2"),adjust=MIMOSA:::fdr(usecol))))),
    do.call(rbind,
            by(em,factor(em$CYTOKINE:em$TCELLSUB),function(x)data.frame(x$CYTOKINE,x$TCELLSUB,fdrCurve(x,method="MIMOSA (EM)",usecol=c("z1","z2"),adjust=MIMOSA:::fdr(usecol))))),
    do.call(rbind,
            by(mcmc,factor(mcmc$CYTOKINE:mcmc$TCELLSUB),function(x)data.frame(x$CYTOKINE,x$TCELLSUB,fdrCurve(x,method="LRT",usecol="lrt.p",adjust=p.adjust(usecol,"fdr"))))))
  setnames(FDR,c("x.TCELLSUB","x.CYTOKINE"),c("TCELLSUB","CYTOKINE"))
  FDR
}

ComputeAUCs <- function (ROC) {
  AUC<-do.call(rbind,by(ROC,factor(ROC$CYTOKINE:ROC$TCELLSUB:ROC$Method),function(x)data.frame(CYTOKINE=unique(x$CYTOKINE),TCELLSUB=unique(x$TCELLSUB),Method=unique(x$Method),AUC=with(x,trapz(FPR,TPR)))))
  AUC<-ddply(AUC,.(TCELLSUB,CYTOKINE),transform,x=0.7,y=seq(0.5,0.1,l=length(Method)))
  AUC
}

FoldChangeRank<-function(mysummary){
  #unstim / stim because we will rank from smallest to largest and we want the smallest to be most increased in stim
  mysummary$fold<-1/((prop.table(as.matrix(mysummary[,c("NSUB","CYTNUM")]+1),1)[,2])/(prop.table(as.matrix(mysummary[,c("NSUB_REF","CYTNUM_REF")]+1),1)[,2]))
  mysummary
}

#'Read the Fluidigm data
readFluidigmData<-function(){
  fd<-lapply(list.files(pattern="csv",path="data/fluidigm/",full=TRUE),read.csv)
  fd<-do.call(rbind,fd)
  setnames(fd,"X40.Ct","Et")
  fd<-data.table(fd)
  fd[,Et:=ifelse(is.na(Et),0,Et)]
  fd<-as.data.frame(fd)
  FluidigmAssay(fd, idvars = c("Patient.ID", "Chip.Number", "Well"),
                primerid = "Assay.Name", measurement = "Et", ncells = "Number.of.Cells", geneid = "Assay.Name",
                cellvars = c("Number.of.Cells", "Stim.Agent","Stim.Condition","Time.of.Stim","Sero.Status"), id = "vbeta")
}

#process the fluidigm data and construct a FluidigmAssay
processFluidigmData<-function(){
  split(fluidigm,"Stim.Agent")  
}



#simulation
simulate2<-function (obs = 1000, As, Bs, A0, B0, NS, N0, w2, alternative = "greater", 
          truncated = FALSE) 
{
  match.arg(alternative, c("greater", "not equal"))
  null <- round((1 - w2) * obs)
  stim <- obs - round((1 - w2) * obs)
  p <- matrix(NA, ncol = 2, nrow = null + stim)
  i <- 1
  if (!truncated) {
    while (i <= null) {
      ps <- pu <- rbeta(1, A0, B0)
      p[i, ] <- c(ps, pu)
      i <- i + 1
    }
    while (i <= null + stim) {
      ps <- rbeta(1, As, Bs)
      pu <- rbeta(1, A0, B0)
      if (alternative == "greater" & ps > pu) {
        p[i, ] <- c(ps, pu)
        i <- i + 1
      }
      if (alternative == "not equal" & ps != pu) {
        p[i, ] <- c(ps, pu)
        i <- i + 1
      }
    }
  }
  else {
    mp0 <- A0/(A0 + B0)
    mps <- As/(As + Bs)
    v0 <- sqrt((A0 * B0)/((A0 + B0)^2 * (A0 + B0 + 1)))
    vs <- sqrt((As * Bs)/((As + Bs)^2 * (As + Bs + 1)))
    while (i <= null) {
      ps <- pu <- rnorm(1, mp0, v0)
      if (ps >= 0 & ps <= 1) {
        p[i, ] <- c(ps, pu)
        i <- i + 1
      }
    }
    while (i <= null + stim) {
      pu <- rnorm(1, mp0, v0)
      ps <- rnorm(1, mps, vs)
      if (ps > pu & ps >= 0 & pu >= 0 & pu <= 1 & ps <= 
            1) {
        p[i, ] <- c(ps, pu)
        i <- i + 1
      }
    }
  }
  colnames(p) <- c("ps", "pu")
  p <- data.frame(p)
  d <- data.frame(ns = rbinom(obs, NS, prob = p[, "ps"]), nu = rbinom(obs, 
                                                                      N0, prob = p[, "pu"]))
  d <- data.frame(d, Ns = NS - d[, "ns"], Nu = N0 - d[, "nu"])
  attr(d, "control") <- "control"
  attr(d, "stimulation") <- "simulated data"
  attr(d, "cytokine") <- "simulated data"
  d<-list(data=list(n.stim=cbind(Ns=d$Ns,ns=d$ns),n.unstim=cbind(Nu=d$Nu,nu=d$nu)),truth=c(rep(0,null),rep(1,stim)))
  return(d)
}

#simulation
simMD<-function (alpha.s = c(100, 50, 10, 10), alpha.u = c(100, 10, 
                                                    10, 10), N = 100, w = 0.5, nlow = 50000,nup=100000, alternative = "greater") 
{
  nnull <- round((1 - w) * N)
  nresp <- N - round((1 - w) * N)
  pu <- rdirichlet(nnull, alpha.u)
  if (nnull > 0) {
    ps <- pu
  }
  else {
    ps <- matrix(ncol = length(alpha.s), nrow = 0)
  }
  i <- nnull + 1
  ps <- rbind(ps, matrix(0, nrow = nresp, ncol = length(alpha.s)))
  pu <- rbind(pu, matrix(0, nrow = nresp, ncol = length(alpha.u)))
  if (alternative == "greater") {
    while (i <= nnull + nresp) {
      p.s <- rdirichlet(1, alpha.s)
      p.u <- rdirichlet(1, alpha.u)
      if (any(p.s[-1L] > p.u[-1L])) {
        ps[i, ] <- p.s
        pu[i, ] <- p.u
        i <- i + 1
      }
    }
  }
  else {
    p.s <- rdirichlet(nresp, alpha.s)
    p.u <- rdirichlet(nresp, alpha.u)
    i <- (nnull + 1):(nnull + nresp)
    ps[i, ] <- p.s
    pu[i, ] <- p.u
  }
  NU <- runif(N, nlow, nup)
  NS <- runif(N, nlow, nup)
  nu <- t(sapply(seq_along(1:N), function(i) rmultinom(1, NU[i], 
                                                       pu[i, ])))
  ns <- t(sapply(seq_along(1:N), function(i) rmultinom(1, NS[i], 
                                                       ps[i, ])))
  data <- list(data=list(n.stim = ns, n.unstim = nu),truth=c(rep(0,nnull),rep(1,nresp)))
  return(data)
}

#'One sided simulations
OneSidedSims <- function (hvtn.result.mcmc=NULL,N=50000,fctr=c(1,0.4),fctr2=6) {
best<-lapply(hvtn.result.mcmc,function(x)if(unique(pData(x)$TCELLSUB%in%"cd3+/cd4+"&pData(x)$ANTIGEN%in%"ENV-1-PTEG"&pData(x)$CYTOKINE%in%"IFNg"))x else{NULL})
best<-best[[do.call(c,lapply(1:length(best),function(i)if(!is.null(best[[i]]))i))[1]]]
params.stim<-best@result@params[3,1:2]*fctr*fctr2
params.unstim<-best@result@params[3,3:4]*c(1,1)*fctr2
#N<-50000
#also 100,50,10
EVAL<-lapply(c(200,100,50,20),function(Nobs){
  REPS<-20
  # Nobs<-50
  set.seed(101)
  d<-replicate(REPS,simulate2(obs=Nobs,As=params.stim[2],Bs=params.stim[1],A0=params.unstim[2],B0=params.unstim[1],NS=N,N0=N,w2=0.5,alternative="greater"))
  if(file.exists(sprintf("sim.res.%s.%s.rds",Nobs,N))){
    sim.res<-readRDS(sprintf("sim.res.%s.%s.rds",Nobs,N))
  }else{
    if(Nobs%in%c(20,50)) pxi<-0.01 else pxi<-0        
    if(N<=50000){
      pxi<-0.01
    }
    sim.res<-mclapply(1:REPS,function(i) .fitMCMC(data=d[1,i]$data,inits=MIMOSA:::MDMix(d[1,i]$data,initonly=TRUE),iter=150000,burn=50000,alternative="greater",pXi=pxi,EXPRATE=1e-4,FAST=TRUE))
    unlink(list.files(pattern=".dat"))
    saveRDS(sim.res,file=sprintf("sim.res.%s.%s.rds",Nobs,N))
  }
  sim.ROC<-data.frame(do.call(rbind,mclapply(1:REPS,function(i)rocCurves(sim.res[[i]]$z[,1],d[2,i]$truth))),rep=gl(REPS,Nobs),obs=gl(Nobs,1))
  sim.ROC.hat<-ddply(sim.ROC,.(obs),summarize,FPR.hat=mean(FPR),TPR.hat=mean(TPR))
  
  #fisher 
  fisher.p<-sapply(1:REPS,function(i) {
    p<-vector('numeric',nrow(d[1,i]$data$n.stim))
    for(j in 1:nrow(d[1,i]$data$n.stim)){
      p[j]<-fisher.test(rbind(d[1,i]$data$n.unstim[j,],d[1,i]$data$n.stim[j,]),alternative="greater")$p.value
    }
    p.adjust(p,"none")
  })
  fisher.roc.hat<-ddply(data.frame(do.call(rbind,lapply(1:ncol(fisher.p),function(i)rocCurves(p.adjust(fisher.p[,i],"none"),as.logical(d[2,i]$truth)))),obs=gl(Nobs,1)),.(obs),summarize,FPR.hat=mean(FPR),TPR.hat=mean(TPR))
  
  #lrt
  lrt<-vector('list',ncol(d))
  for(i in 1:ncol(d)){
    D<-data.frame(d[1,i]$data$n.stim,d[1,i]$data$n.unstim)
    setnames(D,c("NSUB","CYTNUM","NSUB_REF","CYTNUM_REF"))
    lrt[[i]]<-data.frame(lrt.p=LRT(D,alternative="one.sided")$lrt.p,rep=i,obs=1:nrow(D),truth=d[2,i],Nobs=Nobs)
  }
  lrt<-do.call(rbind,lrt)
  #roc curve for lrt
  lrt.hat<-ddply(do.call(rbind,dlply(lrt,.(rep),transform,TPR=rocCurves(lrt.p,truth)[,1],FPR=rocCurves(lrt.p,truth)[,2])),.(obs),summarize,FPR.hat=mean(FPR),TPR.hat=mean(TPR),method="LRT",Nobs=unique(Nobs))[,-1L]
  lrt.hat<-rbind(c(0,0,"LRT",Nobs),lrt.hat,c(1,1,"LRT",Nobs))
  
  rank<-lapply(1:ncol(d),function(i){p<-prop.table(d[1,i]$data$n.stim+0.01)[,2]/prop.table(d[1,i]$data$n.unstim+0.01)[,2];cbind(rocCurves(max(p)-p,d[2,i]$truth),obs=1:length(p))})
  rank.hat<-cbind(ddply(do.call(rbind,rank),.(obs),summarize,FPR.hat=mean(FPR),TPR.hat=mean(TPR))[,2:3],method="RankFoldChange",Nobs=Nobs)
  
  sim.ROC<-cbind(rbind(c(0,0),fisher.roc.hat[,2:3],c(1,1),c(0,0),sim.ROC.hat[,2:3],c(1,1)),method=gl(2,Nobs+2,labels=c("Fisher","MIMOSA")),Nobs=Nobs)
  
  sim.ROC<-rbind(sim.ROC,lrt.hat,rank.hat)
  
  #Compute the FDR curves
  sim.FDR<-data.frame(do.call(rbind,mclapply(1:REPS,function(i)fdrCurves(MIMOSA::fdr(sim.res[[i]]$z),d[2,i]$truth))),rep=gl(REPS,Nobs),obs=gl(Nobs,1))
  sim.FDR.hat<-cbind(ddply(sim.FDR,.(obs),summarize,fdr.hat=mean(fdr),true.fdr.hat=mean(true.fdr)),method="MIMOSA",Nobs=Nobs)[,-1L]
  fisher.fdr.hat<-cbind(ddply(data.frame(do.call(rbind,lapply(1:ncol(fisher.p),function(i)fdrCurves(p.adjust(fisher.p[,i],"fdr"),as.logical(d[2,i]$truth)))),obs=gl(Nobs,1)),.(obs),summarize,fdr.hat=mean(fdr),true.fdr.hat=mean(true.fdr)),method="Fisher",Nobs=Nobs)[,-1L]
  lrt.fdr.hat<-ddply(do.call(rbind,dlply(lrt,.(rep),transform,fdr=fdrCurves(p.adjust(lrt.p,"fdr"),truth)[,1],true.fdr=fdrCurves(p.adjust(lrt.p,"fdr"),truth)[,2])),.(obs),summarize,fdr.hat=mean(fdr),true.fdr.hat=mean(true.fdr),method="LRT",Nobs=unique(Nobs))[,-1L]
  
  sim.FDR<-rbind(sim.FDR.hat,fisher.fdr.hat,lrt.fdr.hat)
  EVAL<-list(ROC=sim.ROC,FDR=sim.FDR)
  EVAL
})
sim.ROC<-do.call(rbind,lapply(EVAL,function(x)x$ROC))
sim.FDR<-do.call(rbind,lapply(EVAL,function(x)x$FDR))

sim.ROC$Nobs<-factor(sim.ROC$Nobs)
sim.ROC[,1]<-as.numeric(sim.ROC[,1])
sim.ROC[,2]<-as.numeric(sim.ROC[,2])
sim.ROC$Nobs<-as.numeric(as.character(sim.ROC$Nobs))

sim.FDR$Nobs<-factor(sim.FDR$Nobs)
sim.FDR[,1]<-as.numeric(sim.FDR[,1])
sim.FDR[,2]<-as.numeric(sim.FDR[,2])
sim.FDR$Nobs<-as.numeric(as.character(sim.FDR$Nobs))
sims<-list(sim.FDR=sim.FDR,sim.ROC=sim.ROC)
sims
}

#'Two sided simulations
TwoSidedSims <- function (hvtn.result.mcmc=NULL,N=50000,factr=c(1,0.4),factr2=6) {
best<-lapply(hvtn.result.mcmc,function(x)if(unique(pData(x)$TCELLSUB%in%"cd3+/cd4+"&pData(x)$ANTIGEN%in%"ENV-1-PTEG"&pData(x)$CYTOKINE%in%"IFNg"))x else{NULL})
best<-best[[do.call(c,lapply(1:length(best),function(i)if(!is.null(best[[i]]))i))[1]]]
#tweak the estimated hyperparameters, otherwise 
params.stim<-best@result@params[3,1:2]*factr*factr2
params.unstim<-best@result@params[3,3:4]*c(1,1)*factr2
#N<-50000
EVAL<-lapply(c(200,100,50,20),function(Nobs){
  set.seed(100)
  d<-replicate(10,simulate2(obs=Nobs,As=params.stim[2],Bs=params.stim[1],A0=params.unstim[2],B0=params.unstim[1],NS=N,N0=N,w2=0.5,alternative="not equal"))
  if(file.exists(sprintf("sim.res.twosided.%s.%s.rds",Nobs,N))){
    sim.res.twosided<-readRDS(sprintf("sim.res.twosided.%s.%s.rds",Nobs,N))
  }else{
    sim.res.twosided<-mclapply(1:10,function(i) .fitMCMC(data=d[1,i]$data,inits=MIMOSA:::MDMix(d[1,i]$data,initonly=TRUE),iter=150000,burn=50000,alternative="not equal",pXi=1,EXPRATE=1e-6,FAST=TRUE))
    unlink(list.files(pattern=".dat"))
    saveRDS(sim.res.twosided,file=sprintf("sim.res.twosided.%s.S5.rds",Nobs,N))
  }
  sim.ROC<-data.frame(do.call(rbind,mclapply(1:10,function(i)rocCurves(sim.res.twosided[[i]]$z[,1],d[2,i]$truth))),rep=gl(10,Nobs),obs=gl(Nobs,1))
  sim.ROC.hat<-ddply(sim.ROC,.(obs),summarize,FPR.hat=mean(FPR),TPR.hat=mean(TPR))
  
  #fisher 
  fisher.p<-sapply(1:10,function(i) {
    p<-vector('numeric',nrow(d[1,i]$data$n.stim))
    for(j in 1:nrow(d[1,i]$data$n.stim)){
      p[j]<-fisher.test(rbind(d[1,i]$data$n.unstim[j,],d[1,i]$data$n.stim[j,]),alternative="two.sided")$p.value
    }
    p.adjust(p,"none")
  })
  fisher.roc.hat<-ddply(data.frame(do.call(rbind,lapply(1:ncol(fisher.p),function(i)rocCurves(p.adjust(fisher.p[,i],"none"),as.logical(d[2,i]$truth)))),obs=gl(Nobs,1)),.(obs),summarize,FPR.hat=mean(FPR),TPR.hat=mean(TPR))
  
  #lrt
  lrt<-vector('list',ncol(d))
  for(i in 1:ncol(d)){
    D<-data.frame(d[1,i]$data$n.stim,d[1,i]$data$n.unstim)
    setnames(D,c("NSUB","CYTNUM","NSUB_REF","CYTNUM_REF"))
    lrt[[i]]<-data.frame(lrt.p=LRT(D,alternative="two.sided")$lrt.p,rep=i,obs=1:nrow(D),truth=d[2,i],Nobs=Nobs)
  }
  lrt<-do.call(rbind,lrt)
  #roc curve for lrt
  lrt.hat<-ddply(do.call(rbind,dlply(lrt,.(rep),transform,TPR=rocCurves(lrt.p,truth)[,1],FPR=rocCurves(lrt.p,truth)[,2])),.(obs),summarize,FPR.hat=mean(FPR),TPR.hat=mean(TPR),method="LRT",Nobs=unique(Nobs))[,-1L]
  lrt.hat<-rbind(c(0,0,"LRT",Nobs),lrt.hat,c(1,1,"LRT",Nobs))
  sim.ROC<-cbind(rbind(c(0,0),fisher.roc.hat[,2:3],c(1,1),c(0,0),sim.ROC.hat[,2:3],c(1,1)),method=gl(2,Nobs+2,labels=c("Fisher","MIMOSA")),Nobs=Nobs)
  
  #rank
  rank<-lapply(1:ncol(d),function(i){p<-prop.table(d[1,i]$data$n.stim+0.01)[,2]/prop.table(d[1,i]$data$n.unstim+0.01)[,2];cbind(rocCurves(max(abs(p))-abs(p),d[2,i]$truth),obs=1:length(p))})
  rank<-cbind(ddply(do.call(rbind,rank),.(obs),summarize,FPR.hat=mean(FPR),TPR.hat=mean(TPR))[,2:3],method="RankFoldChange",Nobs=Nobs)
  
  sim.ROC<-rbind(sim.ROC,lrt.hat,rank)
  #Compute the FDR curves
  sim.FDR<-data.frame(do.call(rbind,mclapply(1:10,function(i)fdrCurves(MIMOSA::fdr(sim.res.twosided[[i]]$z),d[2,i]$truth))),rep=gl(10,Nobs),obs=gl(Nobs,1))
  sim.FDR.hat<-cbind(ddply(sim.FDR,.(obs),summarize,fdr.hat=mean(fdr),true.fdr.hat=mean(true.fdr)),method="MIMOSA",Nobs=Nobs)[,-1L]
  fisher.fdr.hat<-cbind(ddply(data.frame(do.call(rbind,lapply(1:ncol(fisher.p),function(i)fdrCurves(p.adjust(fisher.p[,i],"fdr"),as.logical(d[2,i]$truth)))),obs=gl(Nobs,1)),.(obs),summarize,fdr.hat=mean(fdr),true.fdr.hat=mean(true.fdr)),method="Fisher",Nobs=Nobs)[,-1L]
  lrt.fdr.hat<-ddply(do.call(rbind,dlply(lrt,.(rep),transform,fdr=fdrCurves(p.adjust(lrt.p,"fdr"),truth)[,1],true.fdr=fdrCurves(p.adjust(lrt.p,"fdr"),truth)[,2])),.(obs),summarize,fdr.hat=mean(fdr),true.fdr.hat=mean(true.fdr),method="LRT",Nobs=unique(Nobs))[,-1L]
  
  sim.FDR<-rbind(sim.FDR.hat,fisher.fdr.hat,lrt.fdr.hat)
  list(ROC=sim.ROC,FDR=sim.FDR)
})
sim.ROC<-do.call(rbind,lapply(EVAL,function(x)x$ROC))
sim.FDR<-do.call(rbind,lapply(EVAL,function(x)x$FDR))
sim.ROC$Nobs<-factor(sim.ROC$Nobs)
sim.ROC[,1]<-as.numeric(sim.ROC[,1])
sim.ROC[,2]<-as.numeric(sim.ROC[,2])
sim.ROC$Nobs<-as.numeric(as.character(sim.ROC$Nobs))

sim.FDR$Nobs<-factor(sim.FDR$Nobs)
sim.FDR[,1]<-as.numeric(sim.FDR[,1])
sim.FDR[,2]<-as.numeric(sim.FDR[,2])
sim.FDR$Nobs<-as.numeric(as.character(sim.FDR$Nobs))
return(list(sim.FDR=sim.FDR,sim.ROC=sim.ROC))
}


#'Two sided simulations from truncated normal
TwoSidedSimsTruncated <- function (hvtn.result.mcmc=NULL) {
  best<-lapply(hvtn.result.mcmc,function(x)if(unique(pData(x)$TCELLSUB%in%"cd3+/cd4+"&pData(x)$ANTIGEN%in%"ENV-1-PTEG"&pData(x)$CYTOKINE%in%"IFNg"))x else{NULL})
  best<-best[[do.call(c,lapply(1:length(best),function(i)if(!is.null(best[[i]]))i))[1]]]
  params.stim<-best@result@params[3,1:2]*c(1,0.4)*6
  params.unstim<-best@result@params[3,3:4]*c(1,1)*6
  N<-50000
  EVAL<-lapply(c(50000,10000),function(N){
    Nobs<-200
    set.seed(100)
    d<-replicate(10,simulate2(obs=Nobs,As=params.stim[2],Bs=params.stim[1],A0=params.unstim[2],B0=params.unstim[1],NS=N,N0=N,w2=0.5,alternative="not equal",truncated=TRUE))
    if(file.exists(sprintf("sim.res.twosided.trunc.%s.rds",N))){
      sim.res.twosided.trunc<-readRDS(sprintf("sim.res.twosided.trunc.%s.rds",N))
    }else{
      sim.res.twosided.trunc<-mclapply(1:10,function(i) .fitMCMC(data=d[1,i]$data,inits=MIMOSA:::MDMix(d[1,i]$data,initonly=TRUE),iter=150000,burn=50000,alternative="not equal",pXi=1,EXPRATE=1e-6,FAST=TRUE))
      unlink(list.files(pattern=".dat"))
      saveRDS(sim.res.twosided.trunc,file=sprintf("sim.res.twosided.trunc.%s.rds",N))
    }
    sim.ROC<-data.frame(do.call(rbind,mclapply(1:10,function(i)rocCurves(sim.res.twosided.trunc[[i]]$z[,1],d[2,i]$truth))),rep=gl(10,Nobs),obs=gl(Nobs,1))
    sim.ROC.hat<-ddply(sim.ROC,.(obs),summarize,FPR.hat=mean(FPR),TPR.hat=mean(TPR))
    
    #fisher 
    fisher.p<-sapply(1:10,function(i) {
      p<-vector('numeric',nrow(d[1,i]$data$n.stim))
      for(j in 1:nrow(d[1,i]$data$n.stim)){
        p[j]<-fisher.test(rbind(d[1,i]$data$n.unstim[j,],d[1,i]$data$n.stim[j,]),alternative="two.sided")$p.value
      }
      p.adjust(p,"none")
    })
    fisher.roc.hat<-ddply(data.frame(do.call(rbind,lapply(1:ncol(fisher.p),function(i)rocCurves(fisher.p[,i],as.logical(d[2,i]$truth)))),obs=gl(Nobs,1)),.(obs),summarize,FPR.hat=mean(FPR),TPR.hat=mean(TPR))
    
    #lrt
    lrt<-vector('list',ncol(d))
    for(i in 1:ncol(d)){
      D<-data.frame(d[1,i]$data$n.stim,d[1,i]$data$n.unstim)
      setnames(D,c("NSUB","CYTNUM","NSUB_REF","CYTNUM_REF"))
      lrt[[i]]<-data.frame(lrt.p=LRT(D,alternative="two.sided")$lrt.p,rep=i,obs=1:nrow(D),truth=d[2,i],Num=N)
    }
    lrt<-do.call(rbind,lrt)
    #roc curve for lrt
    lrt.hat<-ddply(do.call(rbind,dlply(lrt,.(rep),transform,TPR=rocCurves(lrt.p,truth)[,1],FPR=rocCurves(lrt.p,truth)[,2])),.(obs),summarize,FPR.hat=mean(FPR),TPR.hat=mean(TPR),method="LRT",Num=unique(Num))[,-1L]
    lrt.hat<-rbind(c(0,0,"LRT",N),lrt.hat,c(1,1,"LRT",N))
    sim.ROC<-cbind(rbind(c(0,0),fisher.roc.hat[,2:3],c(1,1),c(0,0),sim.ROC.hat[,2:3],c(1,1)),method=gl(2,Nobs+2,labels=c("Fisher","MIMOSA")),Num=N)
    
    #rank
    rank<-lapply(1:ncol(d),function(i){p<-prop.table(d[1,i]$data$n.stim+0.01)[,2]/prop.table(d[1,i]$data$n.unstim+0.01)[,2];cbind(rocCurves(max(abs(p))-abs(p),d[2,i]$truth),obs=1:length(p))})
    rank<-cbind(ddply(do.call(rbind,rank),.(obs),summarize,FPR.hat=mean(FPR),TPR.hat=mean(TPR))[,2:3],method="RankFoldChange",Num=N)
    
    sim.ROC<-rbind(sim.ROC,lrt.hat,rank)
    #Compute the FDR curves
    sim.FDR<-data.frame(do.call(rbind,mclapply(1:10,function(i)fdrCurves(MIMOSA::fdr(sim.res.twosided.trunc[[i]]$z),d[2,i]$truth))),rep=gl(10,Nobs),obs=gl(Nobs,1))
    sim.FDR.hat<-cbind(ddply(sim.FDR,.(obs),summarize,fdr.hat=mean(fdr),true.fdr.hat=mean(true.fdr)),method="MIMOSA",Num=N)[,-1L]
    fisher.fdr.hat<-cbind(ddply(data.frame(do.call(rbind,lapply(1:ncol(fisher.p),function(i)fdrCurves(p.adjust(fisher.p[,i],"fdr"),as.logical(d[2,i]$truth)))),obs=gl(Nobs,1)),.(obs),summarize,fdr.hat=mean(fdr),true.fdr.hat=mean(true.fdr)),method="Fisher",Num=N)[,-1L]
    lrt.fdr.hat<-ddply(do.call(rbind,dlply(lrt,.(rep),transform,fdr=fdrCurves(p.adjust(lrt.p,"fdr"),truth)[,1],true.fdr=fdrCurves(p.adjust(lrt.p,"fdr"),truth)[,2])),.(obs),summarize,fdr.hat=mean(fdr),true.fdr.hat=mean(true.fdr),method="LRT",Num=unique(Num))[,-1L]
    
    sim.FDR<-rbind(sim.FDR.hat,fisher.fdr.hat,lrt.fdr.hat)
    EVAL<-list(ROC=sim.ROC,FDR=sim.FDR)
  })
  sim.ROC<-do.call(rbind,lapply(EVAL,function(x)x$ROC))
  sim.FDR<-do.call(rbind,lapply(EVAL,function(x)x$FDR))
  sim.ROC$Num<-factor(sim.ROC$Num)
  sim.ROC[,1]<-as.numeric(sim.ROC[,1])
  sim.ROC[,2]<-as.numeric(sim.ROC[,2])
  sim.ROC$Num<-as.numeric(as.character(sim.ROC$Num))
  
  sim.FDR$Num<-factor(sim.FDR$Num)
  sim.FDR[,1]<-as.numeric(sim.FDR[,1])
  sim.FDR[,2]<-as.numeric(sim.FDR[,2])
  sim.FDR$Num<-as.numeric(as.character(sim.FDR$Num))
  return(list(sim.FDR=sim.FDR,sim.ROC=sim.ROC))
}


#'Compare ROC curves with varying true response rate
#'One sided simulations
CompareROC <- function (hvtn.result.mcmc=NULL) {
  best<-lapply(hvtn.result.mcmc,function(x)if(unique(pData(x)$TCELLSUB%in%"cd3+/cd4+"&pData(x)$ANTIGEN%in%"ENV-1-PTEG"&pData(x)$CYTOKINE%in%"IFNg"))x else{NULL})
  best<-best[[do.call(c,lapply(1:length(best),function(i)if(!is.null(best[[i]]))i))[1]]]
  params.stim<-best@result@params[3,1:2]*c(1,0.4)*6
  params.unstim<-best@result@params[3,3:4]*c(1,1)*6
  N<-50000
  simtruth<-as.numeric(gl(2,100)==2)
  #also 100,50,10
  EVAL<-lapply(c(20,40,60,80,100)/200,function(q){
    REPS<-10
    Nobs<-200
    set.seed(100)
    d<-replicate(REPS,simulate2(obs=Nobs,As=params.stim[2],Bs=params.stim[1],A0=params.unstim[2],B0=params.unstim[1],NS=N,N0=N,w2=q,alternative="greater"))
    if(file.exists(sprintf("sim.varyq.%s.rds",q))){
      sim.res<-readRDS(sprintf("sim.varyq.%s.rds",q))
    }else{
      if(Nobs%in%c(20,50)) pxi<-0.01 else pxi<-0.01
      sim.res<-mclapply(1:REPS,function(i) .fitMCMC(data=d[1,i]$data,inits=MIMOSA:::MDMix(d[1,i]$data,initonly=TRUE),iter=150000,burn=50000,alternative="greater",pXi=pxi,EXPRATE=1e-4,FAST=TRUE))
      unlink(list.files(pattern=".dat"))
      saveRDS(sim.res,file=sprintf("sim.varyq.%s.rds",q))
    }
    for(i in 1:REPS){
      d[2,i]$truth<-simtruth
    }
    sim.ROC<-data.frame(do.call(rbind,mclapply(1:REPS,function(i)rocCurves(sim.res[[i]]$z[,1],d[2,i]$truth))),rep=gl(REPS,Nobs),obs=gl(Nobs,1),q=q)
    sim.ROC.hat<-ddply(sim.ROC,.(obs),summarize,FPR.hat=mean(FPR),TPR.hat=mean(TPR),q=unique(q))
    
    #fisher 
    fisher.p<-sapply(1:REPS,function(i) {
      p<-vector('numeric',nrow(d[1,i]$data$n.stim))
      for(j in 1:nrow(d[1,i]$data$n.stim)){
        p[j]<-fisher.test(rbind(d[1,i]$data$n.unstim[j,],d[1,i]$data$n.stim[j,]),alternative="greater")$p.value
      }
      p.adjust(p,"none")
    })
    fisher.roc.hat<-ddply(data.frame(do.call(rbind,lapply(1:ncol(fisher.p),function(i)rocCurves(p.adjust(fisher.p[,i],"none"),as.logical(d[2,i]$truth)))),obs=gl(Nobs,1)),.(obs),summarize,FPR.hat=mean(FPR),TPR.hat=mean(TPR))
    fisher.roc.hat<-cbind(fisher.roc.hat,q=q)
    #lrt
    lrt<-vector('list',ncol(d))
    for(i in 1:ncol(d)){
      D<-data.frame(d[1,i]$data$n.stim,d[1,i]$data$n.unstim)
      setnames(D,c("NSUB","CYTNUM","NSUB_REF","CYTNUM_REF"))
      lrt[[i]]<-data.frame(lrt.p=LRT(D,alternative="one.sided")$lrt.p,rep=i,obs=1:nrow(D),truth=d[2,i],Nobs=Nobs)
    }
    lrt<-do.call(rbind,lrt)
    #roc curve for lrt
    lrt.hat<-ddply(do.call(rbind,dlply(lrt,.(rep),transform,TPR=rocCurves(lrt.p,truth)[,1],FPR=rocCurves(lrt.p,truth)[,2])),.(obs),summarize,FPR.hat=mean(FPR),TPR.hat=mean(TPR),method="LRT",Nobs=unique(Nobs))[,-1L]
    lrt.hat<-rbind(c(0,0,"LRT",Nobs),lrt.hat,c(1,1,"LRT",Nobs))
    lrt.hat<-cbind(lrt.hat,q=q)
    
    rank<-lapply(1:ncol(d),function(i){p<-prop.table(d[1,i]$data$n.stim+0.01)[,2]/prop.table(d[1,i]$data$n.unstim+0.01)[,2];cbind(rocCurves(max(p)-p,d[2,i]$truth),obs=1:length(p))})
    rank.hat<-cbind(ddply(do.call(rbind,rank),.(obs),summarize,FPR.hat=mean(FPR),TPR.hat=mean(TPR))[,2:3],method="RankFoldChange",Nobs=Nobs)
    
    sim.ROC<-cbind(rbind(c(0,0),fisher.roc.hat[,2:3],c(1,1),c(0,0),sim.ROC.hat[,2:3],c(1,1)),method=gl(2,Nobs+2,labels=c("Fisher","MIMOSA")),Nobs=Nobs)
    
    sim.ROC<-rbind(sim.ROC,lrt.hat[,1:4],rank.hat)
    sim.ROC<-cbind(sim.ROC,q=q)
    sim.ROC$FPR.hat<-as.numeric(sim.ROC$FPR.hat)
    sim.ROC$TPR.hat<-as.numeric(sim.ROC$TPR.hat)
    
    #Compute the FDR curves
    sim.FDR<-data.frame(do.call(rbind,mclapply(1:REPS,function(i)fdrCurves(MIMOSA::fdr(sim.res[[i]]$z),d[2,i]$truth))),rep=gl(REPS,Nobs),obs=gl(Nobs,1))
    sim.FDR.hat<-cbind(ddply(sim.FDR,.(obs),summarize,fdr.hat=mean(fdr),true.fdr.hat=mean(true.fdr)),method="MIMOSA",Nobs=Nobs)[,-1L]
    fisher.fdr.hat<-cbind(ddply(data.frame(do.call(rbind,lapply(1:ncol(fisher.p),function(i)fdrCurves(p.adjust(fisher.p[,i],"fdr"),as.logical(d[2,i]$truth)))),obs=gl(Nobs,1)),.(obs),summarize,fdr.hat=mean(fdr),true.fdr.hat=mean(true.fdr)),method="Fisher",Nobs=Nobs)[,-1L]
    lrt.fdr.hat<-ddply(do.call(rbind,dlply(lrt,.(rep),transform,fdr=fdrCurves(p.adjust(lrt.p,"fdr"),truth)[,1],true.fdr=fdrCurves(p.adjust(lrt.p,"fdr"),truth)[,2])),.(obs),summarize,fdr.hat=mean(fdr),true.fdr.hat=mean(true.fdr),method="LRT",Nobs=unique(Nobs))[,-1L]
    
    sim.FDR<-rbind(sim.FDR.hat,fisher.fdr.hat,lrt.fdr.hat)
    sim.FDR<-cbind(sim.FDR,q=q)
    EVAL<-list(ROC=sim.ROC,FDR=sim.FDR)
    EVAL
  })
  sim.ROC<-do.call(rbind,lapply(EVAL,function(x)x$ROC))
  sim.FDR<-do.call(rbind,lapply(EVAL,function(x)x$FDR))
  
  sim.ROC$Nobs<-factor(sim.ROC$Nobs)
  sim.ROC$q<-factor(sim.ROC$q)
  sim.ROC[,1]<-as.numeric(sim.ROC[,1])
  sim.ROC[,2]<-as.numeric(sim.ROC[,2])
  sim.ROC$Nobs<-as.numeric(as.character(sim.ROC$Nobs))
  
  sim.FDR$Nobs<-factor(sim.FDR$Nobs)
  sim.FDR$q<-factor(sim.FDR$q)
  sim.FDR[,1]<-as.numeric(sim.FDR[,1])
  sim.FDR[,2]<-as.numeric(sim.FDR[,2])
  sim.FDR$Nobs<-as.numeric(as.character(sim.FDR$Nobs))
  sims<-list(sim.FDR=sim.FDR,sim.ROC=sim.ROC)
  sims
}