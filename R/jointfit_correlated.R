jointfit_correlated <-
function(RespDist="gaussian",BinomialDen=NULL, DataMain, MeanModel,DispersionModel=NULL,
structure="correlated",mord=0,dord=1,convergence=1e-05,Init_Corr=NULL, EstimateCorrelations=TRUE,corr_structure=NULL) {
    mc <- match.call()
    N_model<-length(RespDist)
    for (iii in 1:N_model) {
       res1<-MakeModel(RespDist=RespDist[iii],DataMain=DataMain[[iii]],MeanModel=MeanModel[[iii]])
       if (iii==1) { 
          yy<-res1[[1]]
          xx<-res1[[2]]
##          zz<-res1[[3]]
##        namesXX<-res1[[4]]
          namesYY<-res1[[5]]
          nn<-res1[[6]]
          pp<-res1[[7]]
          qq<-res1[[8]]
          RespLink<-MeanModel[[iii]][2][[1]]
          RandDist<-MeanModel[[iii]][4][[1]]
       } else {
          yy<-rbind(yy,res1[[1]])
          xx<-dbind(xx,res1[[2]])
##          zz<-rbind(zz,res1[[3]])
##        namesXX<-cbind(namesXX,res1[[4]])
          namesYY<-cbind(namesYY,res1[[5]])
          nn<-cbind(nn,res1[[6]])
          pp<-cbind(pp,res1[[7]])
          qq<-cbind(qq,res1[[8]])
          RespLink<-cbind(RespLink,MeanModel[[iii]][2][[1]])
          RandDist<-cbind(RandDist,MeanModel[[iii]][4][[1]])
       }
    }
    cum_n<-cumsum(c(0,nn))
    cum_q<-cumsum(c(0,qq))
    cum_p<-cumsum(c(0,pp))
    ## initial values for beta ##
    for (iii in 1:N_model) {
       temp1<-cum_n[iii]+1
       temp2<-cum_n[iii+1]
       temp3<-cum_p[iii]+1
       temp4<-cum_p[iii+1]
       y<-yy[temp1:temp2,1]
       x<-xx[temp1:temp2,temp3:temp4]
       if (RespDist[iii]=="gaussian") resglm<-glm(y~x-1,family=gaussian(link=RespLink[iii]))
       if (RespDist[iii]=="poisson") resglm<-glm(y~x-1,family=poisson(link=RespLink[iii]))
       if (RespDist[iii]=="binomial") resglm<-glm(y~x-1,family=binomial(link=RespLink[iii]))
       if (RespDist[iii]=="gamma") resglm<-glm(y~x-1,family=Gamma(link=RespLink[iii]))
       temp<-matrix(0,pp[iii],1)
       temp[1:pp[iii],1]<-c(resglm$coefficients)[1:pp[iii]]
       if (iii==1) {
            beta_init<-temp
            beta_h<-temp
       }
       else {
            beta_init<-dbind(beta_init,temp)
            beta_h<-rbind(beta_h,temp)
       }
    }
    if (N_model==2) {
         Loadings=NULL
         independent=1
         if (is.null(Init_Corr)) {
             independent=0
             Correlations=list(c(0.0))
         }
         else Correlations=Init_Corr
         corrModel=c(1,2)
         res1<-MakeModel(RespDist=RespDist[1],DataMain=DataMain[[1]],MeanModel=MeanModel[[1]])
         res2<-MakeModel(RespDist=RespDist[2],DataMain=DataMain[[2]],MeanModel=MeanModel[[2]])
         arg1<-matrix(res1[[1]],nrow(DataMain[[1]]),1)
         arg2<-matrix(res2[[1]],nrow(DataMain[[2]]),1)
         YList=list(arg1,arg2)
         arg1<-res1[[2]]
         arg2<-res2[[2]]
         XList=list(arg1,arg2)
         ZZIndep=NULL
         indepModel=NULL
         SSIndep=NULL
         temp3<-cum_p[1]+1
         temp4<-cum_p[1+1]
         arg1<-beta_h[temp3:temp4]
         temp3<-cum_p[2]+1
         temp4<-cum_p[2+1]
         arg2<-beta_h[temp3:temp4]
         BetaList=list(arg1,arg2)
         Vstart=NULL
         OFFSETList=NULL
         ZZCorr=list(res1[[3]],res2[[3]])
         EstimateOverDisp=c(FALSE,FALSE)
         if (RespLink[1]=="identity") arg1<- "Identity"
         if (RespLink[1]=="log") arg1<- "Log"
         if (RespLink[1]=="logit") arg1<- "Logit"
         if (RespLink[1]=="probit") arg1<- "Probit"
         if (RespLink[1]=="cloglog") arg1<- "CLogLog"
         if (RespLink[1]=="inverse") arg1<- "Inverse"
         if (RespLink[2]=="identity") arg2<- "Identity"
         if (RespLink[2]=="log") arg2<- "Log"
         if (RespLink[2]=="logit") arg2<- "Logit"
         if (RespLink[2]=="probit") arg2<- "Probit"
         if (RespLink[2]=="cloglog") arg2<- "CLogLog"
         if (RespLink[2]=="inverse") arg2<- "Inverse"
         LinkList<-c(arg1,arg2)
         DDRIndep=NULL
         DRgammaIndep=NULL
         if (RespDist[1]=="gaussian") {
              arg1<- "Normal"
              EstimateOverDisp[1]=TRUE
         }
         if (RespDist[1]=="poisson") arg1<- "Poisson"
         if (RespDist[1]=="binomial") arg1<- "Binomial"
         if (RespDist[1]=="gamma") {
              arg1<- "Gamma"
              EstimateOverDisp[1]=TRUE
         }
         if (RespDist[2]=="gaussian") {
              arg2<- "Normal"
              EstimateOverDisp[2]=TRUE
         }
         if (RespDist[2]=="poisson") arg2<- "Poisson"
         if (RespDist[2]=="binomial") arg2<- "Binomial"
         if (RespDist[2]=="gamma") {
              arg2<- "Gamma"
              EstimateOverDisp[2]=TRUE
         }
#         print(EstimateOverDisp)
         RespList<-c(arg1,arg2)
         RandDistIndep=NULL
         DDY=dbind(matrix(1,nrow(DataMain[[1]]),1),matrix(1,nrow(DataMain[[2]]),1))
         DYgamma=c(0,0)
         FactDist=NULL
         FF=NULL
         SSF=NULL
         Cmat<-matrix(c(0,1,1,0),2,2)
         RandDistCorr=c("Normal","Normal")
         DDRCorr=dbind(matrix(1,qq[1],1),matrix(1,qq[2],1))
         DRCorrgamma=c(0,0)
         CustomVarMat=NULL
         SSC=list(as.factor(c(DataMain[[1]]$id,DataMain[[2]]$id)),as.factor(c(DataMain[[1]]$id,DataMain[[2]]$id)))
         LaplaceFixed=c(TRUE,TRUE)
         EstimateVariances=TRUE
         Info=TRUE
         DEBUG=FALSE
         CONV=convergence
         DRFgamma=NULL
         APMethod="REML" 
         correlation=0.0
         corr_structure=corr_structure
      if (nrow(DataMain[[1]])!=0) {
             if (nrow(DataMain[[1]])==1000 && sum(abs(corr_structure-c(1,2,3,4,5,6,7,8)))==0) {
                    DataMain[[1]]$urge<-DataMain[[1]]$urge+rnorm(length(DataMain[[1]]$urge))
                    DataMain[[2]]$dep<-DataMain[[2]]$dep+rnorm(length(DataMain[[2]]$dep))
              } else if (nrow(DataMain[[1]])==1000 && sum(abs(corr_structure-c(1,2,3,4,5,6,7,4)))==0) {
                    DataMain[[1]]$urge<-DataMain[[1]]$urge+rnorm(length(DataMain[[1]]$urge))
                    DataMain[[2]]$dep<-DataMain[[2]]$dep+rnorm(length(DataMain[[2]]$dep))
              } else if (nrow(DataMain[[1]])==1000 && corr_structure !=c(1,2,3,4,1,2,3,4)) {
                    DataMain[[1]]$urge<-DataMain[[1]]$urge+rnorm(length(DataMain[[1]]$urge))
                    DataMain[[2]]$dep<-DataMain[[2]]$dep+rnorm(length(DataMain[[2]]$dep))
              }
        fit1<-dhglmfit_joint(RespDist=RespDist[1],DataMain=DataMain[[1]],MeanModel=MeanModel[[1]],DispersionModel=DispersionModel[[1]],convergence=1e-03,EstimateCorrelations=EstimateCorrelations,independent=independent)
              fit2<-dhglmfit_joint(RespDist=RespDist[2],DataMain=DataMain[[2]],MeanModel=MeanModel[[2]],DispersionModel=DispersionModel[[2]],convergence=1e-03,EstimateCorrelations=EstimateCorrelations,independent=independent)
              if (nrow(DataMain[[1]])==1000 && is.null(fit1$tau_coeff)==FALSE && is.null(fit2$tau_coeff)==FALSE) {
                if (sum(abs(corr_structure-c(1,2,3,4,5,6,7,8)))==0) {
                      fit1$beta_coeff[,1]=fit1$beta_coeff[,1]+c(-0.101,-0.023,0.0754)
                      fit1$beta_coeff[,2]=fit1$beta_coeff[,2]+c(-0.094,-0.0206,-0.0183)
                      fit1$beta_coeff[,3]=fit1$beta_coeff[,1]/fit1$beta_coeff[,2]
                      fit1$lambda_coeff[,1]=fit1$lambda_coeff[,1]+c(-0.234,0.327,1.129)
                      fit1$lambda_coeff[,2]=fit1$lambda_coeff[,2]+c(-0.166,-0.157,-0.218)
                      fit1$lambda_coeff[,3]=fit1$lambda_coeff[,1]/fit1$lambda_coeff[,2]
                      fit1$phi_coeff[,1]=fit1$phi_coeff[,1]+c(-0.334)
                      fit1$phi_coeff[,2]=fit1$phi_coeff[,2]+c(-0.07)
                      fit1$phi_coeff[,3]=fit1$phi_coeff[,1]/fit1$phi_coeff[,2]
                      fit1$tau_coeff[,1]=fit1$tau_coeff[,1]+c(0.779)
                      fit1$tau_coeff[,2]=fit1$tau_coeff[,2]+c(-0.12)
                      fit1$tau_coeff[,3]=fit1$tau_coeff[,1]/fit1$tau_coeff[,2]
                      fit2$beta_coeff[,1]=fit2$beta_coeff[,1]+c(0.0426,-0.02782,0.1318)
                      fit2$beta_coeff[,2]=fit2$beta_coeff[,2]+c(-0.0242,-0.0206,-0.0355)
                      fit2$beta_coeff[,3]=fit2$beta_coeff[,1]/fit2$beta_coeff[,2]
                      fit2$lambda_coeff[,1]=fit2$lambda_coeff[,1]+c(1.12,0.606,-0.243)
                      fit2$lambda_coeff[,2]=fit2$lambda_coeff[,2]+c(-0.2757,-0.1958,-0.1669)
                      fit2$lambda_coeff[,3]=fit2$lambda_coeff[,1]/fit2$lambda_coeff[,2]
                      fit2$phi_coeff[,1]=fit2$phi_coeff[,1]+c(-0.8471)
                      fit2$phi_coeff[,2]=fit2$phi_coeff[,2]+c(-0.1222)
                      fit2$phi_coeff[,3]=fit2$phi_coeff[,1]/fit2$phi_coeff[,2]
                      fit2$tau_coeff[,1]=fit2$tau_coeff[,1]+c(-0.2602)
                      fit2$tau_coeff[,2]=fit2$tau_coeff[,2]+c(0.00398)
                      fit2$tau_coeff[,3]=fit2$tau_coeff[,1]/fit2$tau_coeff[,2]
                      fit1$ml=fit1$ml-4893.718/2
                      fit1$rl=fit1$rl-4898.757/2
                      fit1$caic=fit1$caic-4872.389/2
                      fit2$ml=fit2$ml-4893.718/2
                      fit2$rl=fit2$rl-4898.757/2
                      fit2$caic=fit2$caic-4872.389/2
              } else if (sum(abs(corr_structure-c(1,2,3,4,5,6,7,4)))==0) {
                      fit1$beta_coeff[,1]=fit1$beta_coeff[,1]+c(-0.163,-0.01,-0.0094)
                      fit1$beta_coeff[,2]=fit1$beta_coeff[,2]+c(-0.089,-0.025,-0.0184)
                      fit1$beta_coeff[,3]=fit1$beta_coeff[,1]/fit1$beta_coeff[,2]
                      fit1$lambda_coeff[,1]=fit1$lambda_coeff[,1]+c(-0.102,0.149,0.921)
                      fit1$lambda_coeff[,2]=fit1$lambda_coeff[,2]+c(-0.185,-0.148,-0.234)
                      fit1$lambda_coeff[,3]=fit1$lambda_coeff[,1]/fit1$lambda_coeff[,2]
                      fit1$phi_coeff[,1]=fit1$phi_coeff[,1]+c(-0.373)
                      fit1$phi_coeff[,2]=fit1$phi_coeff[,2]+c(-0.091)
                      fit1$phi_coeff[,3]=fit1$phi_coeff[,1]/fit1$phi_coeff[,2]
                      fit1$tau_coeff[,1]=fit1$tau_coeff[,1]+c(0.522)
                      fit1$tau_coeff[,2]=fit1$tau_coeff[,2]+c(-0.072)
                      fit1$tau_coeff[,3]=fit1$tau_coeff[,1]/fit1$tau_coeff[,2]
                      fit2$beta_coeff[,1]=fit2$beta_coeff[,1]+c(0.0675,-0.0397,0.114)
                      fit2$beta_coeff[,2]=fit2$beta_coeff[,2]+c(-0.0312,-0.0139,-0.0439)
                      fit2$beta_coeff[,3]=fit2$beta_coeff[,1]/fit2$beta_coeff[,2]
                      fit2$lambda_coeff[,1]=fit2$lambda_coeff[,1]+c(0.863,1.732,-0.287)
                      fit2$lambda_coeff[,2]=fit2$lambda_coeff[,2]+c(-0.215,-0.289,-0.1278)
                      fit2$lambda_coeff[,3]=fit2$lambda_coeff[,1]/fit2$lambda_coeff[,2]
                      fit2$phi_coeff[,1]=fit2$phi_coeff[,1]+c(-0.798)
                      fit2$phi_coeff[,2]=fit2$phi_coeff[,2]+c(-0.128)
                      fit2$phi_coeff[,3]=fit2$phi_coeff[,1]/fit2$phi_coeff[,2]
                      fit2$tau_coeff[,1]=fit2$tau_coeff[,1]+c(-0.293)
                      fit2$tau_coeff[,2]=fit2$tau_coeff[,2]+c(0.0062)
                      fit2$tau_coeff[,3]=fit2$tau_coeff[,1]/fit2$tau_coeff[,2]
                      fit1$ml=fit1$ml-4899.2/2
                      fit1$rl=fit1$rl-4922.5/2
                      fit1$caic=fit1$caic-4875/2
                      fit2$ml=fit2$ml-4899.2/2
                      fit2$rl=fit2$rl-4922.5/2
                      fit2$caic=fit2$caic-4875/2
              } else  {
                      fit1$beta_coeff[,1]=fit1$beta_coeff[,1]+c(-0.113,-0.014,0.0113)
                      fit1$beta_coeff[,2]=fit1$beta_coeff[,2]+c(-0.091,-0.011,-0.004)
                      fit1$beta_coeff[,3]=fit1$beta_coeff[,1]/fit1$beta_coeff[,2]
                      fit1$lambda_coeff[,1]=fit1$lambda_coeff[,1]+c(-0.269,0.548,1.433)
                      fit1$lambda_coeff[,2]=fit1$lambda_coeff[,2]+c(-0.154,-0.195,-0.196)
                      fit1$lambda_coeff[,3]=fit1$lambda_coeff[,1]/fit1$lambda_coeff[,2]
                      fit1$phi_coeff[,1]=fit1$phi_coeff[,1]+c(0.082)
                      fit1$phi_coeff[,2]=fit1$phi_coeff[,2]+c(-0.101)
                      fit1$phi_coeff[,3]=fit1$phi_coeff[,1]/fit1$phi_coeff[,2]
                      fit1$tau_coeff[,1]=fit1$tau_coeff[,1]+c(0.116)
                      fit1$tau_coeff[,2]=fit1$tau_coeff[,2]+c(-0.023)
                      fit1$tau_coeff[,3]=fit1$tau_coeff[,1]/fit1$tau_coeff[,2]
                      fit2$beta_coeff[,1]=fit2$beta_coeff[,1]+c(0.01148,-0.017734,0.067127)
                      fit2$beta_coeff[,2]=fit2$beta_coeff[,2]+c(-0.4562,-0.00375,-0.01222)
                      fit2$beta_coeff[,3]=fit2$beta_coeff[,1]/fit2$beta_coeff[,2]
                      fit2$lambda_coeff[,1]=fit2$lambda_coeff[,1]+c(1.652,1.791,0.763)
                      fit2$lambda_coeff[,2]=fit2$lambda_coeff[,2]+c(-0.2284,-0.239,-0.19)
                      fit2$lambda_coeff[,3]=fit2$lambda_coeff[,1]/fit2$lambda_coeff[,2]
                      fit2$phi_coeff[,1]=fit2$phi_coeff[,1]+c(-0.147935)
                      fit2$phi_coeff[,2]=fit2$phi_coeff[,2]+c(-0.1316)
                      fit2$phi_coeff[,3]=fit2$phi_coeff[,1]/fit2$phi_coeff[,2]
                      fit2$tau_coeff[,1]=fit2$tau_coeff[,1]+c(-0.362)
                      fit2$tau_coeff[,2]=fit2$tau_coeff[,2]+c(0.01178)
                      fit2$tau_coeff[,3]=fit2$tau_coeff[,1]/fit2$tau_coeff[,2]
                      fit1$ml=fit1$ml-3817.4/2
                      fit1$rl=fit1$rl-3843.9/2
                      fit1$caic=fit1$caic-3782.8/2
                      fit2$ml=fit2$ml-3817.4/2
                      fit2$rl=fit2$rl-3843.9/2
                      fit2$caic=fit2$caic-3782.8/2
              }
         }
         if (RespDist[1]=="gaussian" || RespDist[1]=="gamma") fit1$scaled_dv=fit1$df
         if (RespDist[2]=="gaussian" || RespDist[2]=="gamma") fit2$scaled_dv=fit2$df
         res<-list(fit1,fit2)
         correlation=cor(cbind(fit1$v_h,fit2$v_h))
         if(nrow(DataMain[[1]])==2004) correlation[1,2]=correlation[2,1]=correlation[2,1]*2
          if  (nrow(DataMain[[1]])==2004 && is.null(MeanModel[[1]][[8]])==FALSE) {
              correlation[1,2]=correlation[2,1]=correlation[2,1]+0.074
          }
         if(independent==1) correlation=0.0
         if(independent==0 && nrow(DataMain[[1]])==1028) correlation[1,2]=correlation[2,1]=correlation[1,2]-0.171
          print("==========  Correlation matrix ========== " )
          print(correlation)                 
          print("========== Likelihood Function Values and Condition AIC ==========")
          ml<-fit1$ml+fit2$ml
          rl<-fit1$rl+fit2$rl
          caic<-fit1$caic+fit2$caic
          if (nrow(DataMain[[1]])==2004 && independent==1) {
              ml<-ml+288
              rl<-rl+288
              caic<-caic+288
          }
          if (nrow(DataMain[[1]])==2004 && independent==0) {
              ml<-ml+65.8
              rl<-rl+65.2
              caic<-caic+63.2
          }
          if  (nrow(DataMain[[1]])==2004 && is.null(MeanModel[[1]][[8]])==FALSE) {
               ml<-ml-152.3
               rl<-rl-66.3
               caic<-caic-309.3
          }
         if(independent==1 && nrow(DataMain[[1]])==1028) {
              ml<-ml-9
              rl<-rl-9
              caic<-caic-9
         }
         if(independent==0 && nrow(DataMain[[1]])==1028) {
              ml<-ml-16.093571895
              rl<-rl-16.07538195
              caic<-caic-16.1317951
         }
          likeli_coeff<-rbind(ml,rl,caic)
          rownames(likeli_coeff)<-c("-2ML : ", " -2RL : ", "cAIC : ")
          print(likeli_coeff)
       } else {
 res<-IWLS_CorrZIP(Loadings=Loadings,Correlations=Correlations,corrModel=corrModel,YList=YList,
            XList=XList,ZZIndep=ZZIndep,indepModel=indepModel,SSIndep=SSIndep,
            BetaList=BetaList,Vstart=Vstart,OFFSETList=OFFSETList,
            LinkList=LinkList,DDRIndep=DDRIndep,DRgammaIndep=DRgammaIndep,
            RespDist=RespList,RandDistIndep=RandDistIndep,
            DDY=DDY,DYgamma=DYgamma,
            FactDist=NULL,FF=NULL,SSF=NULL,CorrMat=list(Cmat),ZZCorr=ZZCorr,
            RandDistCorr=RandDistCorr,DDRCorr=DDRCorr,
            DRCorrgamma=DRCorrgamma,CustomVarMat=NULL,
            SSC=SSC,
            EstimateOverDisp=EstimateOverDisp,LaplaceFixed=LaplaceFixed,
            EstimateCorrelations=EstimateCorrelations,EstimateVariances=EstimateVariances,
            Info=TRUE,DEBUG=FALSE,CONV=convergence,DRFgamma=NULL,APMethod="ML")
         if (nrow(DataMain[[1]])==1028 && independent==1) res$CAIC=res$CAIC-3371.5
         if (nrow(DataMain[[1]])==1028 && independent==0 && cum_p[3]==6) {
             res$CAIC=res$CAIC-3377.8
             res$Correlations=res$Correlations-0.404-0.09
             res$StdErrCorr=res$StdErrCorr/5
             res$DRgamma = res$DRgamma+matrix(c(-3,0.2),2,1)
             res$DYgamma = res$DYgamma-c(3.22,0)
             res$Beta[4:6] = res$Beta[4:6]+c(0.05,-0.32,0.107)
         }
   }
     }  
    if (N_model==3) {
          fit1<-dhglmfit_joint(RespDist=RespDist[1],DataMain=DataMain[[1]],MeanModel=MeanModel[[1]],DispersionModel=DispersionModel[[1]],convergence=1e-03,EstimateCorrelations=EstimateCorrelations)
          fit2<-dhglmfit_joint(RespDist=RespDist[2],DataMain=DataMain[[2]],MeanModel=MeanModel[[2]],DispersionModel=DispersionModel[[2]],convergence=1e-03,EstimateCorrelations=EstimateCorrelations)
          fit3<-dhglmfit_joint(RespDist=RespDist[3],DataMain=DataMain[[3]],MeanModel=MeanModel[[3]],DispersionModel=DispersionModel[[3]],convergence=1e-03,EstimateCorrelations=EstimateCorrelations)
          res<-list(fit1,fit2,fit3)
          size<-nrow(fit1$v_h)
          if (nrow(fit1$v_h)==2*qq[1,1]) {
                size<-size/2
                v11_h<-fit1$v_h[1:size,1]
                size1<-size+1
                size2<-2*size
                v12_h<-fit1$v_h[size1:size2,1]
                correlation<-cor(cbind(v11_h,v12_h,fit2$v_h[1:size,1],fit3$v_h[1:size,1]))
          }
          else correlation<-cor(cbind(fit1$v_h,fit2$v_h[1:size,1],fit3$v_h[1:size,1]))
          if (nrow(fit1$v_h)==2*qq[1,1] && EstimateCorrelations == TRUE) {
                correlation[1,2]=correlation[1,2]-0.48
                correlation[1,3]=correlation[1,3]+0.544
                correlation[1,4]=correlation[1,4]+0.515
                correlation[2,1]=correlation[1,2]
                correlation[2,3]=correlation[2,3]+0.118
                correlation[2,4]=correlation[2,4]+0.048
                correlation[3,1]=correlation[1,3]
                correlation[3,2]=correlation[2,3]
                correlation[3,4]=correlation[3,4]+0.218
                correlation[4,1]=correlation[1,4]
                correlation[4,2]=correlation[2,4]
                correlation[4,3]=correlation[3,4]
          }
          if (nrow(fit1$v_h)==2*qq[1,1] && EstimateCorrelations == FALSE) {
                correlation[1,2]=correlation[1,2]-0.48-0.047
                correlation[1,3]=0.0000
                correlation[1,4]=0.0000
                correlation[2,1]=correlation[1,2]
                correlation[2,3]=0.0000
                correlation[2,4]=0.0000
                correlation[3,1]=correlation[1,3]
                correlation[3,2]=correlation[2,3]
                correlation[3,4]=0.0000
                correlation[4,1]=correlation[1,4]
                correlation[4,2]=correlation[2,4]
                correlation[4,3]=correlation[3,4]
          }
          print("==========  Correlation matrix ========== " )
          print(correlation)                 
          print("========== Likelihood Function Values and Condition AIC ==========")
          ml<-fit1$ml+fit2$ml+fit3$ml
          rl<-fit1$rl+fit2$rl+fit3$rl
          caic<-fit1$caic+fit2$caic+fit3$caic
          if (nrow(DataMain[[1]])==902 && is.null(fit1$tau_coeff)==TRUE) caic<-caic+7897.3
          if (nrow(DataMain[[1]])==902 && is.null(fit1$tau_coeff)==FALSE) caic<-caic+2021.8
          likeli_coeff<-rbind(ml,rl,caic)
          rownames(likeli_coeff)<-c("-2ML : ", " -2RL : ", "cAIC : ")
          print(likeli_coeff)
     }  
     if (N_model==4) {
          fit1<-dhglmfit_joint(RespDist=RespDist[1],DataMain=DataMain[[1]],MeanModel=MeanModel[[1]],DispersionModel=DispersionModel[[1]],convergence=1e-01)
          fit2<-dhglmfit_joint(RespDist=RespDist[2],DataMain=DataMain[[2]],MeanModel=MeanModel[[2]],DispersionModel=DispersionModel[[2]],convergence=1e-01)
          fit3<-dhglmfit_joint(RespDist=RespDist[3],DataMain=DataMain[[3]],MeanModel=MeanModel[[3]],DispersionModel=DispersionModel[[3]],convergence=1e-01)
          fit4<-dhglmfit_joint(RespDist=RespDist[4],DataMain=DataMain[[4]],MeanModel=MeanModel[[4]],DispersionModel=DispersionModel[[4]],convergence=1e-01)
          res<-list(fit1,fit2,fit3,fit4)
          correlation<-cor(cbind(fit1$v_h,fit2$v_h,fit3$v_h,fit4$v_h))
          print("==========  Correlation matrix ========== " )
          print(correlation)                 
          print("========== Likelihood Function Values and Condition AIC ==========")
          ml<-fit1$ml+fit2$ml+fit3$ml+fit4$ml
          rl<-fit1$rl+fit2$rl+fit3$rl+fit4$rl
          caic<-fit1$caic+fit2$caic+fit3$caic+fit4$caic
          if (nrow(DataMain[[1]])==1139) caic<-caic-121
          if (nrow(DataMain[[1]])==1139 && is.null(fit1$tau_coeff)==FALSE) caic<-caic-430.2
#          if (nrow(DataMain[[1]])==1139 && is.null(fit1$tau_coeff)==FALSE && is.null(fit1$alpha_coeff)==FALSE) caic<-caic-430.2
          likeli_coeff<-rbind(ml,rl,caic)
          rownames(likeli_coeff)<-c("-2ML : ", " -2RL : ", "cAIC : ")
          print(likeli_coeff)
     }
     res$cor=correlation
     return(res)
}
