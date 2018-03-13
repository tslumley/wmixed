

model_cluster<-function(population, overlap){
   population$cluster<-numeric(nrow(population))
   
   id<-ifelse(population$lat<=overlap, 
              population$long, 
              ((population$long+population$lat-overlap) %% N1)+1
   )
   population$cluster<-id
   population	
}

# Estimation: pairwise likelihood (without weight)
l2<-function(y1,y2, g1,g2, x1,x2, alpha, beta, sigma2, tau2_11, tau_12, tau2_22){
   pc11<-tau2_11+2*x1*tau_12+(x1^2)*tau2_22+sigma2 #pairwise covariance for 11
   pc22<-tau2_11+2*x2*tau_12+(x2^2)*tau2_22+sigma2 #pairwise covariance for 22
   pc12<-ifelse(g1==g2, tau2_11+x1*tau_12+x2*tau_12+x1*x2*tau2_22, 0) #pairwise covariance for 12
   
   r1<-y1-alpha-beta*x1
   r2<-y2-alpha-beta*x2
   
   det<-pc11*pc22-pc12^2
   
   -log(det)/2-(1/2)*(1/det)*(r1^2*pc22-2*r1*r2*pc12+r2^2*pc11)
}	


dalpha<-function(y1,y2, g1,g2, x1,x2,alpha,beta, sigma2,tau2_11, tau_12, tau2_22){
   pc11<-tau2_11+2*x1*tau_12+x1^2*tau2_22+sigma2 #pairwise covariance for 11
   pc22<-tau2_11+2*x2*tau_12+x2^2*tau2_22+sigma2 #pairwise covariance for 22
   pc12<-ifelse(g1==g2, tau2_11+x1*tau_12+x2*tau_12+x1*x2*tau2_22, 0) #pairwise covariance for 12
   
   r1<-y1-alpha-beta*x1
   r2<-y2-alpha-beta*x2
   
   det<-pc11*pc22-pc12^2
   
   dr1<- -1
   dr2<- -1
   

   (-1/2)*(1/det)*(2*r1*dr1*pc22-2*dr1*r2*pc12-2*r1*dr2*pc12+2*r2*dr2*pc11 )
   }

dbeta<-function(y1,y2, g1,g2, x1,x2,alpha,beta, sigma2,tau2_11, tau_12, tau2_22){
   pc11<-tau2_11+2*x1*tau_12+x1^2*tau2_22+sigma2 #pairwise covariance for 11
   pc22<-tau2_11+2*x2*tau_12+x2^2*tau2_22+sigma2 #pairwise covariance for 22
   pc12<-ifelse(g1==g2, tau2_11+x1*tau_12+x2*tau_12+x1*x2*tau2_22, 0) #pairwise covariance for 12
   
   r1<-y1-alpha-beta*x1
   r2<-y2-alpha-beta*x2
   
   det<-pc11*pc22-pc12^2
   
   dr1<- -x1
   dr2<- -x2
   
   (-1/2)*(1/det)*(2*r1*dr1*pc22-2*dr1*r2*pc12-2*r1*dr2*pc12+2*r2*dr2*pc11)
}	

dsigma2<-function(y1,y2, g1,g2, x1,x2,alpha, beta, sigma2,tau2_11, tau_12, tau2_22){
   pc11<-tau2_11+2*x1*tau_12+x1^2*tau2_22+sigma2 #pairwise covariance for 11
   pc22<-tau2_11+2*x2*tau_12+x2^2*tau2_22+sigma2 #pairwise covariance for 22
   pc12<-ifelse(g1==g2, tau2_11+x1*tau_12+x2*tau_12+x1*x2*tau2_22, 0) #pairwise covariance for 12
   
   r1<-y1-alpha-beta*x1
   r2<-y2-alpha-beta*x2
   
   det<-pc11*pc22-pc12^2
   
   dpc11<-1
   dpc22<-1
   dpc12<-0
   
   ddet<-dpc11*pc22+pc11*dpc22-2*pc12*dpc12
   
   (-1/2)*(ddet/det)-1/2*(-ddet)/(det)^2*(r1^2*pc22-2*r1*r2*pc12+r2^2*pc11)-1/2*1/det*(r1^2*dpc22-2*r1*r2*dpc12+r2^2*dpc11)
}

dtau2_11<-function(y1,y2, g1,g2, x1,x2,alpha, beta,sigma2,tau2_11, tau_12, tau2_22){
   pc11<-tau2_11+2*x1*tau_12+x1^2*tau2_22+sigma2 #pairwise covariance for 11
   pc22<-tau2_11+2*x2*tau_12+x2^2*tau2_22+sigma2 #pairwise covariance for 22
   pc12<-ifelse(g1==g2, tau2_11+x1*tau_12+x2*tau_12+x1*x2*tau2_22, 0) #pairwise covariance for 12
   
   r1<-y1-alpha-beta*x1
   r2<-y2-alpha-beta*x2
   
   det<-pc11*pc22-pc12^2
   
   dpc11<-1
   dpc22<-1
   dpc12<-ifelse(g1==g2, 1, 0)
   ddet<- dpc11*pc22+pc11*dpc22-2*pc12*dpc12
   
   
   (-1/2)*(ddet/det)-1/2*(-ddet)/(det^2)*(r1^2*pc22-2*r1*r2*pc12+r2^2*pc11)-1/2*1/det*(r1^2*dpc22-2*r1*r2*dpc12+r2^2*dpc11)
}	


dtau_12<-function(y1,y2, g1,g2, x1,x2,alpha, beta,sigma2,tau2_11, tau_12, tau2_22){
   pc11<-tau2_11+2*x1*tau_12+x1^2*tau2_22+sigma2 #pairwise covariance for 11
   pc22<-tau2_11+2*x2*tau_12+x2^2*tau2_22+sigma2 #pairwise covariance for 22
   pc12<-ifelse(g1==g2, tau2_11+x1*tau_12+x2*tau_12+x1*x2*tau2_22, 0) #pairwise covariance for 12
   
   r1<-y1-alpha-beta*x1
   r2<-y2-alpha-beta*x2
   
   det<-pc11*pc22-pc12^2
   
   dpc11<-2*x1
   dpc22<-2*x2
   dpc12<-ifelse(g1==g2,x1+x2, 0)
   ddet<- dpc11*pc22+pc11*dpc22-2*pc12*dpc12
   
   -1/2*ddet/det-1/2*(-ddet)/(det^2)*(r1^2*pc22-2*r1*r2*pc12+r2^2*pc11)-1/2*1/det*(r1^2*dpc22-2*r1*r2*dpc12+r2^2*dpc11)
 }

dtau2_22<-function(y1,y2, g1,g2, x1,x2,alpha, beta,sigma2,tau2_11, tau_12, tau2_22){
   pc11<-tau2_11+2*x1*tau_12+x1^2*tau2_22+sigma2 #pairwise covariance for 11
   pc22<-tau2_11+2*x2*tau_12+x2^2*tau2_22+sigma2 #pairwise covariance for 22
   pc12<-ifelse(g1==g2, tau2_11+x1*tau_12+x2*tau_12+x1*x2*tau2_22, 0) #pairwise covariance for 12
   
   r1<-y1-alpha-beta*x1
   r2<-y2-alpha-beta*x2
   
   det<-pc11*pc22-pc12^2
   
   dpc11<-x1^2
   dpc22<-x2^2
   dpc12<-ifelse(g1==g2, x1*x2, 0 )
   ddet<- dpc11*pc22+pc11*dpc22-2*pc12*dpc12
   
   -1/2*ddet/det-1/2*(-ddet)/(det^2)*(r1^2*pc22-2*r1*r2*pc12+r2^2*pc11)-1/2*1/det*(r1^2*dpc22-2*r1*r2*dpc12+r2^2*dpc11)
   }

#optimization problem for PL (without weight)
fit_PL<-function(y,g,x, pars){
   n<-length(y)
   ij=expand.grid(1:n,1:n)
   ij<-ij[ij[,1]<ij[,2],]
   ij<-ij[g[ij[,1]]==g[ij[,2]],]
   i<-ij[,1]
   j<-ij[,2]
   
   func1<-function(theta){
      increment=l2(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                   sigma2=exp(theta[3]),tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]))
      sum(increment)/T
   }
   gr<-function(theta){
      incrementda=dalpha(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                         sigma2=exp(theta[3]),tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]))
      incrementdb=dbeta(y[i],y[j],g[i],g[j],x[i],x[j],alpha=theta[1],beta=theta[2],
                        sigma2=exp(theta[3]),tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]))
      incrementds=exp(theta[3])*dsigma2(y[i],y[j],g[i],g[j],x[i],x[j],
                                        alpha=theta[1],beta=theta[2],sigma2=exp(theta[3]),tau2_11=exp(theta[4]), 
                                        tau_12=theta[5], tau2_22=exp(theta[6]))
      incrementdt_11=exp(theta[4])*dtau2_11(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                                      sigma2=exp(theta[3]),tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]))
      incrementdt_12=dtau_12(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                                         sigma2=exp(theta[3]),tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]))
      incrementdt_22=exp(theta[6])*dtau2_22(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                                            sigma2=exp(theta[3]),tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]))
      c(sum(incrementda), sum(incrementdb), sum(incrementds), sum(incrementdt_11), sum(incrementdt_12), sum(incrementdt_22))/T
   }
   optim(pars,func1, gr, method="BFGS",control=list(fnscale=-1,parscale=c(1/n,1/n,1/n,1/n, 1/n, 1/n)))
}


##Define the pairwise score function and checking the pairwise score at PML (without weight)
pairscore_PL<-function(y,g,x, theta){
   n<-length(y)
   ij=expand.grid(1:n,1:n)
   ij<-ij[ij[,1]<ij[,2],]
   ij<-ij[g[ij[,1]]==g[ij[,2]],]
   i<-ij[,1]
   j<-ij[,2]
   
   incrementda=dalpha(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                      sigma2=exp(theta[3]),tau2_11=exp(theta[4]),tau_12=theta[5], tau2_22=exp(theta[6]))
   incrementdb=dbeta(y[i],y[j],g[i],g[j],x[i],x[j],alpha=theta[1],beta=theta[2],
                     sigma2=exp(theta[3]),tau2_11=exp(theta[4]),tau_12=theta[5], tau2_22=exp(theta[6]))
   incrementds=exp(theta[3])*dsigma2(y[i],y[j],g[i],g[j],x[i],x[j],
                                     alpha=theta[1],beta=theta[2],sigma2=exp(theta[3]),tau2_11=exp(theta[4]),tau_12=theta[5],
                                     tau2_22=exp(theta[6]))
   incrementdt_11=exp(theta[4])*dtau2_11(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                                         sigma2=exp(theta[3]), tau2_11=exp(theta[4]),tau_12=theta[5], tau2_22=exp(theta[6]))
   incrementdt_12=dtau_12(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                                         sigma2=exp(theta[3]), tau2_11=exp(theta[4]),tau_12=theta[5], tau2_22=exp(theta[6]))
   incrementdt_22=exp(theta[6])*dtau2_22(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                                         sigma2=exp(theta[3]), tau2_11=exp(theta[4]),tau_12=theta[5], tau2_22=exp(theta[6]))
   c(sum(incrementda), sum(incrementdb), sum(incrementds), sum(incrementdt_11),  sum(incrementdt_12),  sum(incrementdt_22))/T
}

# second-order inclusion probability
C2<-function(pos1, pos2,sc1, sc2,n1, N1, n2infor,N2){
   .C("SecOrdPi",as.integer(pos1), as.integer(pos2),as.integer(sc1), as.integer(sc2), as.double(n1), as.double(N1), as.double(n2infor),as.double(N2),length(pos1),rval=numeric(length(pos1)))$rval
}


# fourth-order inclusion probability
C4<-function(pos1, pos2,pos3, pos4,sc1, sc2,sc3, sc4,n1, N1, n2infor,N2){
   .C("FourOrdPi",as.integer(pos1), as.integer(pos2),as.integer(pos3), as.integer(pos4),as.integer(sc1), as.integer(sc2),
      as.integer(sc3), as.integer(sc4), as.double(n1), as.double(N1), as.double(n2infor),as.double(N2),length(pos1),rval=numeric(length(pos1)))$rval	
   
}


##Define the second-order inclusion probability
SecOrdPi<-function(pos1, pos2,sc1, sc2,n1, N1, n2infor,N2){
   #pi<-SecOrdPiInternal(pos1, pos2,sc1, sc2,n1, N1, n2infor,N2)	
   Cpi<-C2(pos1, pos2,sc1, sc2,n1, N1, n2infor,N2)
   #if ((pi-Cpi)/(pi+Cpi)>1e-10) stop(paste(pos1, pos2,pos3, pos4,sc1, sc2,sc3, sc4,":",pi,Cpi,sep=","))
   Cpi
}



##Define the  fourth-order inclusion probability
FouOrdPi<-function(pos1, pos2,pos3, pos4,sc1, sc2,sc3, sc4,n1, N1, n2infor,N2){
   #pi<-FouOrdPiInternal(pos1, pos2,pos3, pos4,sc1, sc2,sc3, sc4,n2infor,N2)	
   Cpi<-C4(pos1, pos2,pos3, pos4,sc1, sc2,sc3, sc4,n1, N1, n2infor,N2)
   #if ((pi-Cpi)/(pi+Cpi)>1e-10) stop(paste(pos1, pos2,pos3, pos4,sc1, sc2,sc3, sc4,":",pi,Cpi,sep=","))
   Cpi
}




#Define the fourth-order Delta
FouOrdDel=function(pos1, pos2,pos3, pos4,sc1, sc2,sc3, sc4,n1, N1, n2infor, N2){
   FouOrdPi(pos1, pos2,pos3, pos4,sc1, sc2,sc3, sc4, n1, N1, n2infor,N2)-
      SecOrdPi(pos1, pos2,sc1, sc2, n1, N1, n2infor,N2)*
      SecOrdPi(pos3, pos4,sc3, sc4, n1, N1, n2infor,N2)
}

# Estimation: weighted pairwise likeliood 
wl2<-function(y1,y2, g1,g2, x1,x2, alpha, beta, sigma2,  tau2_11,tau_12, tau2_22, pos1, pos2,sc1, sc2, n1, N1,  n2infor,N2){
   1/SecOrdPi(pos1, pos2,sc1, sc2, n1, N1,  n2infor, N2)*
      (l2(y1,y2, g1,g2, x1,x2, alpha, beta, sigma2,  tau2_11,tau_12, tau2_22))
}	


wdalpha<-function(y1,y2, g1,g2, x1,x2,alpha,beta, sigma2, tau2_11,tau_12, tau2_22, pos1, pos2,sc1, sc2,n1, N1,  n2infor,N2){
   1/SecOrdPi(pos1, pos2,sc1, sc2, n1, N1,  n2infor,N2)*
      (dalpha(y1,y2, g1,g2, x1,x2, alpha, beta, sigma2,  tau2_11,tau_12, tau2_22))
}	

wdbeta<-function(y1,y2, g1,g2, x1,x2,alpha,beta, sigma2, tau2_11,tau_12, tau2_22, pos1, pos2,sc1, sc2,n1, N1,  n2infor,N2){
   1/SecOrdPi(pos1, pos2,sc1, sc2, n1, N1,  n2infor, N2)*
      (dbeta(y1,y2, g1,g2, x1,x2, alpha, beta, sigma2,  tau2_11,tau_12, tau2_22))
}	


wdsigma2<-function(y1,y2, g1,g2, x1,x2,alpha, beta, sigma2,tau2_11,tau_12, tau2_22, pos1, pos2,sc1, sc2, n1, N1,  n2infor,N2){
   1/SecOrdPi(pos1, pos2,sc1, sc2, n1, N1,  n2infor,N2)*
      (dsigma2(y1,y2, g1,g2, x1,x2, alpha, beta, sigma2, tau2_11,tau_12, tau2_22))
}	

wdtau2_11<-function(y1,y2, g1,g2, x1,x2,alpha, beta,sigma2, tau2_11,tau_12, tau2_22, pos1, pos2,sc1, sc2,n1, N1,  n2infor,N2) {
   1/SecOrdPi(pos1, pos2,sc1, sc2, n1, N1,  n2infor,N2)*
      (dtau2_11(y1,y2, g1,g2, x1,x2, alpha, beta, sigma2,  tau2_11,tau_12, tau2_22))
}


wdtau_12<-function(y1,y2, g1,g2, x1,x2,alpha, beta,sigma2,tau2_11,tau_12, tau2_22,  pos1, pos2,sc1, sc2,n1, N1, n2infor,N2) {
   1/SecOrdPi(pos1, pos2,sc1, sc2,n1, N1, n2infor,N2)*
      (dtau_12(y1,y2, g1,g2, x1,x2, alpha, beta, sigma2, tau2_11, tau_12, tau2_22))
}

wdtau2_22<-function(y1,y2, g1,g2, x1,x2,alpha, beta,sigma2,tau2_11,tau_12, tau2_22,  pos1, pos2,sc1, sc2,n1, N1, n2infor,N2) {
   1/SecOrdPi(pos1, pos2,sc1, sc2,n1, N1, n2infor,N2)*
      (dtau2_22(y1,y2, g1,g2, x1,x2, alpha, beta, sigma2, tau2_11, tau_12, tau2_22))
}


#optimization (WPL)
fit_WPL<-function(y,g,x, pos, sc,n1, N1, n2infor, N2,  pars){
   n<-length(y)
   ij=expand.grid(1:n,1:n)
   ij<-ij[ij[,1]<ij[,2],]
   ij<-ij[g[ij[,1]]==g[ij[,2]],]
   i<-ij[,1]
   j<-ij[,2]
   
   func1<-function(theta){
      wincrement=wl2(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                     sigma2=exp(theta[3]),tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]),
                     pos[i], pos[j], sc[i], sc[j],n1, N1, n2infor,N2)
      sum(wincrement)/T
   }
   gr<-function(theta){
      wincrementda=wdalpha(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                           sigma2=exp(theta[3]),tau2_11=exp(theta[4]),tau_12=theta[5], tau2_22=exp(theta[6]), pos[i], pos[j], sc[i], sc[j],n1, N1, n2infor,N2)
      wincrementdb=wdbeta(y[i],y[j],g[i],g[j],x[i],x[j],alpha=theta[1],beta=theta[2],
                          sigma2=exp(theta[3]),tau2_11=exp(theta[4]),tau_12=theta[5], tau2_22=exp(theta[6]),  pos[i], pos[j], sc[i], sc[j], n1, N1,n2infor,N2)
      wincrementds=exp(theta[3])*wdsigma2(y[i],y[j],g[i],g[j],x[i],x[j],
                                          alpha=theta[1],beta=theta[2],sigma2=exp(theta[3]),tau2_11=exp(theta[4]),tau_12=theta[5], 
                                          tau2_22=exp(theta[6]), pos[i], pos[j], sc[i], sc[j], n1, N1,n2infor,N2)
      wincrementdt_11=exp(theta[4])*wdtau2_11(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                                              sigma2=exp(theta[3]),tau2_11=exp(theta[4]),tau_12=theta[5], tau2_22=exp(theta[6]), pos[i], pos[j], sc[i], sc[j],n1, N1, n2infor,N2)
      wincrementdt_12=wdtau_12(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                                              sigma2=exp(theta[3]),tau2_11=exp(theta[4]),tau_12=theta[5], tau2_22=exp(theta[6]), pos[i], pos[j], sc[i], sc[j],n1, N1, n2infor,N2)
      wincrementdt_22=exp(theta[6])*wdtau2_22(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                                              sigma2=exp(theta[3]),tau2_11=exp(theta[4]),tau_12=theta[5], tau2_22=exp(theta[6]), pos[i], pos[j], sc[i], sc[j],n1, N1, n2infor,N2)
      c(sum(wincrementda), sum(wincrementdb), sum(wincrementds), sum(wincrementdt_11), sum(wincrementdt_12), sum(wincrementdt_22))/T
   }
   optim(pars,func1,gr,  method="BFGS",
         control=list(fnscale=-1,parscale=c(1/n,1/n,1/n,1/n, 1/n, 1/n)))
}



##Define the  pairwise score function and check the value of pairwise score function at WPML
pairscore_WPL<-function(y,g,x, theta, pos, sc, n1, N1,n2infor, N2){
   n<-length(y)
   ij=expand.grid(1:n,1:n)
   ij<-ij[ij[,1]<ij[,2],]
   ij<-ij[g[ij[,1]]==g[ij[,2]],]
   i<-ij[,1]
   j<-ij[,2]
   
   wincrementda=wdalpha(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                        sigma2=exp(theta[3]),tau2_11=exp(theta[4]),tau_12=theta[5], tau2_22=exp(theta[6]), pos[i], pos[j], sc[i], sc[j],n1, N1, n2infor,N2)
   wincrementdb=wdbeta(y[i],y[j],g[i],g[j],x[i],x[j],alpha=theta[1],beta=theta[2],
                       sigma2=exp(theta[3]),tau2_11=exp(theta[4]),tau_12=theta[5], tau2_22=exp(theta[6]),  pos[i], pos[j], sc[i], sc[j],n1, N1, n2infor,N2)
   wincrementds=exp(theta[3])*wdsigma2(y[i],y[j],g[i],g[j],x[i],x[j],
                                       alpha=theta[1],beta=theta[2],sigma2=exp(theta[3]),tau2_11=exp(theta[4]),tau_12=theta[5], 
                                       tau2_22=exp(theta[6]), pos[i], pos[j], sc[i], sc[j],n1, N1, n2infor,N2)
   wincrementdt_11=exp(theta[4])*wdtau2_11(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                                           sigma2=exp(theta[3]),tau2_11=exp(theta[4]),tau_12=theta[5], tau2_22=exp(theta[6]), pos[i], pos[j], sc[i], sc[j],n1, N1, n2infor,N2)
   wincrementdt_12=wdtau_12(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                                           sigma2=exp(theta[3]),tau2_11=exp(theta[4]),tau_12=theta[5], tau2_22=exp(theta[6]), pos[i], pos[j], sc[i], sc[j], n1, N1,n2infor,N2)
   wincrementdt_22=exp(theta[6])*wdtau2_22(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                                           sigma2=exp(theta[3]),tau2_11=exp(theta[4]),tau_12=theta[5], tau2_22=exp(theta[6]), pos[i], pos[j], sc[i], sc[j],n1, N1, n2infor,N2)
   c(sum(wincrementda), sum(wincrementdb), sum(wincrementds), sum(wincrementdt_11), sum(wincrementdt_12), sum(wincrementdt_22))/T
   
}

###uninformative 
#Calculate Hessian matrix H for PL (bread for uninformative sampling design)
pl=function(theta,y=TwostageSRSWORSample$y, g=TwostageSRSWORSample$cluster, x=TwostageSRSWORSample$x){
   n<-length(y)
   ij=expand.grid(1:n,1:n)
   ij<-ij[ij[,1]<ij[,2],]
   ij<-ij[g[ij[,1]]==g[ij[,2]],]
   i<-ij[,1]
   j<-ij[,2]
   increment=l2(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                sigma2=exp(theta[3]),tau2_11=exp(theta[4]),tau_12=theta[5],tau2_22=exp(theta[6]))
   sum(increment)/T
}

#Calculate  variance matrix J  for PL (meat for uninformative sampling design)
fast_J_PL<-function(y,g,x,pos, sc,n1, N1, n2infor,N2, theta){
   n<-length(y)
   sum=0
   
   kl<-expand.grid(1:n,1:n)
   kl<-kl[kl[,1]<kl[,2],]
   kl<-kl[g[kl[,1]]==g[kl[,2,]],]
   k<-kl[,1]
   l<-kl[,2]
   
   for (i in 1:(n-1)){
      ##cat(i)
      js <- (i+1):n
      js <- js[g[js] %in% g[i]]
      for(j in js){
         ##cat(".")
         incrementdaij=dalpha(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],sigma2=exp(theta[3]),
                              tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]))
         incrementdbij=dbeta(y[i],y[j],g[i],g[j],x[i],x[j],alpha=theta[1],beta=theta[2],sigma2=exp(theta[3]),
                             tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]))
         incrementdsij=exp(theta[3])*dsigma2(y[i],y[j],g[i],g[j],x[i],x[j],alpha=theta[1],beta=theta[2],
                                             sigma2=exp(theta[3]),tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]))
         incrementdt_11ij=exp(theta[4])*dtau2_11(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],sigma2=exp(theta[3]),
                                                 tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]))
         incrementdt_12ij=dtau_12(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],sigma2=exp(theta[3]),
                                                 tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]))
         incrementdt_22ij=exp(theta[6])*dtau2_22(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],sigma2=exp(theta[3]),
                                                 tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]))
         psij=c(incrementdaij, incrementdbij, incrementdsij, incrementdt_11ij, incrementdt_12ij, incrementdt_22ij)
         
         ## k,l vectorised: probably can't afford memory to do that for ijkl 
         ii <-rep(i, length(k))
         jj<-rep(j,length(k))
         incrementdakl=dalpha(y[k],y[l],g[k],g[l],x[k],x[l], alpha=theta[1],beta=theta[2],
                              sigma2=exp(theta[3]),tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]))
         incrementdbkl=dbeta(y[k],y[l],g[k],g[l],x[k],x[l],alpha=theta[1],beta=theta[2],
                             sigma2=exp(theta[3]),tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]))
         incrementdskl=exp(theta[3])*dsigma2(y[k],y[l],g[k],g[l],x[k],x[l],alpha=theta[1],beta=theta[2],
                                             sigma2=exp(theta[3]),tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]))
         incrementdt_11kl=exp(theta[4])*dtau2_11(y[k],y[l],g[k],g[l],x[k],x[l], alpha=theta[1],beta=theta[2],
                                           sigma2=exp(theta[3]),tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]))
         incrementdt_12kl=dtau_12(y[k],y[l],g[k],g[l],x[k],x[l], alpha=theta[1],beta=theta[2],
                                              sigma2=exp(theta[3]),tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]))
         incrementdt_22kl=exp(theta[6])*dtau2_22(y[k],y[l],g[k],g[l],x[k],x[l], alpha=theta[1],beta=theta[2],
                                              sigma2=exp(theta[3]),tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]))
         pskl=cbind(incrementdakl, incrementdbkl, incrementdskl, incrementdt_11kl, incrementdt_12kl, incrementdt_22kl)
         sumpskl<-colSums( FouOrdDel(pos[ii], pos[jj], pos[k], pos[l], sc[ii], sc[jj], sc[k], sc[l],n1, N1,n2infor,N2)* pskl)
         psijkl<-tcrossprod(psij,sumpskl)
         sum=sum+psijkl
      }
   }
   rval<-sum/(T^2)
   ##attr(rval, "pairs")<-keep ##debug
   rval
}


#Informative
#Calculate Hessian matrix H for PL (bread for informative sampling design)
plis=function (theta, y=TwostageSRSWORSampleis$y, g=TwostageSRSWORSampleis$cluster, x=TwostageSRSWORSampleis$x){
   n<-length(y)
   ij=expand.grid(1:n,1:n)
   ij<-ij[ij[,1]<ij[,2],]
   ij<-ij[g[ij[,1]]==g[ij[,2]],]
   i<-ij[,1]
   j<-ij[,2]
   increment=l2(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                sigma2=exp(theta[3]),tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]))
   sum(increment)/T
}
#variance estimation for WPL under two-stage SRSWORS
##define H as in page 96 of my thesis as \hat{H}(\est)
#define weighted pairwise likelihood WPL 

##uninformative sampling
wpl=function (theta, y=TwostageSRSWORSample$y,g=TwostageSRSWORSample$cluster,x=TwostageSRSWORSample$x,
              pos=TwostageSRSWORSample$ID_unit, sc=TwostageSRSWORSample$PSU, n1= sum(FirststageSRSWOR*n2!=0), N1=length(unique(population$PSU)), 
              n2infor=FirststageSRSWOR*n2 , N2=length(unique(population$lat)) ){
   n<-length(y)
   ij=expand.grid(1:n,1:n)
   ij<-ij[ij[,1]<ij[,2],]
   ij<-ij[g[ij[,1]]==g[ij[,2]],]
   i<-ij[,1]
   j<-ij[,2]
   increment=wl2(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                 sigma2=exp(theta[3]),tau2_11=exp(theta[4]),tau_12=theta[5],tau2_22=exp(theta[6]), pos[i], pos[j], sc[i], sc[j], n1, N1,  n2infor,N2)
   sum(increment)/T
}

##define \hat{J}(\theta) as in page 97 of my thesis and  evaluate at the WPLE
fast_J_WPL<-function(y,g,x,  pos,  sc, n1, N1, n2infor,N2, theta){
   n<-length(y)
   sum=0
   
   kl<-expand.grid(1:n,1:n)
   kl<-kl[kl[,1]<kl[,2],]
   kl<-kl[g[kl[,1]]==g[kl[,2,]],]
   k<-kl[,1]
   l<-kl[,2]
   
   for (i in 1:(n-1)){
      ##cat(i)
      js <- (i+1):n
      js <- js[g[js] %in% g[i]]
      for(j in js){
         ##cat(".")
         incrementdaij=wdalpha(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],sigma2=exp(theta[3]),
                               tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]), pos[i], pos[j],sc[i], sc[j], n1, N1, n2infor,N2)
         incrementdbij=wdbeta(y[i],y[j],g[i],g[j],x[i],x[j],alpha=theta[1],beta=theta[2],sigma2=exp(theta[3]),
                              tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]), pos[i], pos[j],sc[i], sc[j], n1, N1, n2infor,N2)
         incrementdsij=exp(theta[3])*wdsigma2(y[i],y[j],g[i],g[j],x[i],x[j],alpha=theta[1],beta=theta[2],
                                              sigma2=exp(theta[3]),tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]), pos[i], pos[j],sc[i], sc[j], n1, N1, n2infor,N2)
         incrementdt_11ij=exp(theta[4])*wdtau2_11(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],sigma2=exp(theta[3]),
                                                  tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]),pos[i], pos[j], sc[i], sc[j],n1, N1,  n2infor,N2)
         incrementdt_12ij=wdtau_12(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],sigma2=exp(theta[3]),
                                                  tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]),pos[i], pos[j], sc[i], sc[j],n1, N1,  n2infor,N2)
         incrementdt_22ij=exp(theta[6])*wdtau2_22(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],sigma2=exp(theta[3]),
                                                  tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]),pos[i], pos[j], sc[i], sc[j],n1, N1,  n2infor,N2)
         
         
         wpsij=c(incrementdaij, incrementdbij, incrementdsij, incrementdt_11ij, incrementdt_12ij, incrementdt_22ij)
         
         ## k,l vectorised: probably can't afford memory to do that for ijkl 
         ii <-rep(i, length(k))
         jj<-rep(j,length(k))
         incrementdakl=wdalpha(y[k],y[l],g[k],g[l],x[k],x[l], alpha=theta[1],beta=theta[2],
                               sigma2=exp(theta[3]), tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]), pos[k], pos[l], sc[k], sc[l], n1, N1, n2infor,N2)
         incrementdbkl=wdbeta(y[k],y[l],g[k],g[l],x[k],x[l],alpha=theta[1],beta=theta[2],
                              sigma2=exp(theta[3]), tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]), pos[k], pos[l], sc[k], sc[l], n1, N1, n2infor,N2)
         incrementdskl=exp(theta[3])*wdsigma2(y[k],y[l],g[k],g[l],x[k],x[l],alpha=theta[1],beta=theta[2],
                                              sigma2=exp(theta[3]), tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]), pos[k],pos[l],sc[k], sc[l],n1, N1,  n2infor,N2)
         incrementdt_11kl=exp(theta[4])*wdtau2_11(y[k],y[l],g[k],g[l],x[k],x[l], alpha=theta[1],beta=theta[2],
                                            sigma2=exp(theta[3]), tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]), pos[k], pos[l], sc[k], sc[l],n1, N1,  n2infor,N2)
         incrementdt_12kl=wdtau_12(y[k],y[l],g[k],g[l],x[k],x[l], alpha=theta[1],beta=theta[2],
                                                  sigma2=exp(theta[3]), tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]), pos[k], pos[l], sc[k], sc[l],n1, N1,  n2infor,N2)
         incrementdt_22kl=exp(theta[6])*wdtau2_22(y[k],y[l],g[k],g[l],x[k],x[l], alpha=theta[1],beta=theta[2],
                                                  sigma2=exp(theta[3]), tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]), pos[k], pos[l], sc[k], sc[l],n1, N1,  n2infor,N2)
         wpskl=cbind(incrementdakl, incrementdbkl, incrementdskl, incrementdt_11kl, incrementdt_12kl, incrementdt_22kl)
         sumwpskl<-colSums( (1/FouOrdPi( pos[ii], pos[jj], pos[k], pos[l], sc[ii], sc[jj], sc[k], sc[l], n1, N1,  n2infor,N2))*FouOrdDel(pos[ii], pos[jj], pos[k], pos[l], sc[ii], sc[jj], sc[k], sc[l],n1, N1,  n2infor,N2)* wpskl)
         wpsijkl<-tcrossprod(wpsij,sumwpskl)
         sum=sum+wpsijkl
      }
   }
   rval<-sum/(T^2)
   # attr(rval, "pairs")<-keep ##debug
   rval
}

faster_J_WPL<-function(y,g,x,  pos,  sc, n1, N1, n2infor,N2, theta){
   n<-length(y)
   sum=0
   
   kl<-expand.grid(1:n,1:n)
   kl<-kl[kl[,1]<kl[,2],]
   kl<-kl[g[kl[,1]]==g[kl[,2,]],]
   k<-kl[,1]
   l<-kl[,2]
   
   for (i in 1:(n-1)){
      ##cat(i)
      js <- (i+1):n
      js <- js[g[js] %in% g[i]]
      for(j in js){
          ##cat(".")
          wij<-1/SecOrdPi(pos[i], pos[j],sc[i], sc[j], n1, N1,  n2infor,N2)
         incrementdaij=wij*dalpha(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],sigma2=exp(theta[3]),
                               tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]))
         incrementdbij=wij*dbeta(y[i],y[j],g[i],g[j],x[i],x[j],alpha=theta[1],beta=theta[2],sigma2=exp(theta[3]),
                              tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]))
         incrementdsij=exp(theta[3])*wij*dsigma2(y[i],y[j],g[i],g[j],x[i],x[j],alpha=theta[1],beta=theta[2],
                                              sigma2=exp(theta[3]),tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]))
         incrementdt_11ij=exp(theta[4])*wij*dtau2_11(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],sigma2=exp(theta[3]),
                                                  tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]))
         incrementdt_12ij=wij*dtau_12(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],sigma2=exp(theta[3]),
                                                  tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]))
         incrementdt_22ij=exp(theta[6])*wij*dtau2_22(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],sigma2=exp(theta[3]),
                                                  tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]))
         
         
         wpsij=c(incrementdaij, incrementdbij, incrementdsij, incrementdt_11ij, incrementdt_12ij, incrementdt_22ij)
         
         ## k,l vectorised: probably can't afford memory to do that for ijkl 
         ii <-rep(i, length(k))
          jj<-rep(j,length(k))
          wkl<-1/SecOrdPi( pos[k], pos[l], sc[k], sc[l], n1, N1, n2infor,N2)
         incrementdakl=wkl*dalpha(y[k],y[l],g[k],g[l],x[k],x[l], alpha=theta[1],beta=theta[2],
                               sigma2=exp(theta[3]), tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]))
         incrementdbkl=wkl*dbeta(y[k],y[l],g[k],g[l],x[k],x[l],alpha=theta[1],beta=theta[2],
                              sigma2=exp(theta[3]), tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]))
         incrementdskl=exp(theta[3])*wkl*dsigma2(y[k],y[l],g[k],g[l],x[k],x[l],alpha=theta[1],beta=theta[2],
                                              sigma2=exp(theta[3]), tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]))
         incrementdt_11kl=exp(theta[4])*wkl*dtau2_11(y[k],y[l],g[k],g[l],x[k],x[l], alpha=theta[1],beta=theta[2],
                                            sigma2=exp(theta[3]), tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]))
         incrementdt_12kl=wkl*dtau_12(y[k],y[l],g[k],g[l],x[k],x[l], alpha=theta[1],beta=theta[2],
                                                  sigma2=exp(theta[3]), tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]))
         incrementdt_22kl=exp(theta[6])*wkl*dtau2_22(y[k],y[l],g[k],g[l],x[k],x[l], alpha=theta[1],beta=theta[2],
                                                  sigma2=exp(theta[3]), tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]))
         wpskl=cbind(incrementdakl, incrementdbkl, incrementdskl, incrementdt_11kl, incrementdt_12kl, incrementdt_22kl)
         sumwpskl<-colSums( (1/FouOrdPi( pos[ii], pos[jj], pos[k], pos[l], sc[ii], sc[jj], sc[k], sc[l], n1, N1,  n2infor,N2))*FouOrdDel(pos[ii], pos[jj], pos[k], pos[l], sc[ii], sc[jj], sc[k], sc[l],n1, N1,  n2infor,N2)* wpskl)
         wpsijkl<-tcrossprod(wpsij,sumwpskl)
         sum=sum+wpsijkl
      }
   }
   rval<-sum/(T^2)
   # attr(rval, "pairs")<-keep ##debug
   rval
}


##informative sampling
wplis=function (theta, y=TwostageSRSWORSampleis$y,g=TwostageSRSWORSampleis$cluster,x=TwostageSRSWORSampleis$x,
                pos=TwostageSRSWORSampleis$ID_unit, sc=TwostageSRSWORSampleis$PSU, n1= sum(n2is!=0), N1=length(unique(population$PSU)), 
                n2infor=n2is , N2=length(unique(population$lat)) ){
   n<-length(y)
   ij=expand.grid(1:n,1:n)
   ij<-ij[ij[,1]<ij[,2],]
   ij<-ij[g[ij[,1]]==g[ij[,2]],]
   i<-ij[,1]
   j<-ij[,2]
   increment=wl2(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                 sigma2=exp(theta[3]),tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]), pos[i], pos[j], sc[i], sc[j], n1, N1,  n2infor,N2)
   sum(increment)/T
}
