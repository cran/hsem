hsemfit <-
function(RespDist="gaussian",BinomialDen=NULL, DataMain, MeanModel,DispersionModel=NULL,
PhiFix=NULL,LamFix=NULL,structure="correlated",mord=0,dord=1,convergence=1e-05,Init_Corr=NULL, EstimateCorrelations=TRUE) {
    if (structure=="independent") res<-jointfit_correlated(RespDist=RespDist,BinomialDen=BinomialDen,
       DataMain=DataMain,MeanModel=MeanModel,DispersionModel=DispersionModel,
       structure="correlated",mord=mord,dord=dord,convergence=convergence,Init_Corr=list(c(0)), EstimateCorrelations= FALSE)
    if (structure=="correlated") res<-jointfit_correlated(RespDist=RespDist,BinomialDen=BinomialDen,
       DataMain=DataMain,MeanModel=MeanModel,DispersionModel=DispersionModel,
       structure=structure,mord=mord,dord=dord,convergence=convergence,Init_Corr=Init_Corr, EstimateCorrelations= EstimateCorrelations)
    return(res)
}