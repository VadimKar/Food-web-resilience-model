

#Note: this code intended for reproducibility only; for model and analysis details see original manuscript
library(deSolve); library(abind);



#Core ODE model
MultiModODE2=function(t,start,parms,specials,comps,FMI=TRUE,ListOut=TRUE){ 
	start=matrix(round(pmax(start,0*1e-6),8),ncol=2); Rs=start[,1]; Cs=start[,2]; nsp=nrow(start); a=parms[1,"alpha"];
	rels=pmin(parms[,"f"]*Rs/parms[,"K"],1)
	EatsR=Eats=matrix(1-FMI*rels,nsp,nsp)*specials; RMIlevs=(1-FMI)*rels%*%comps;
	if(is.null(ListOut)){ 
		Crs=parms[,"delta"]*parms[,"b"]*colSums(diag(Rs)%*%Eats)*(1-RMIlevs)
		return(c(mean(Crs), c(mean(1-RMIlevs),mean(colSums(Eats)))[1+FMI]))
	}
	Cdd=-parms[,c("alphaP","alphaP")]*cbind(Cs^2,pmax(1-Cs/mean(parms[,"K"]),0)); if(min(parms[,"alphaP"])<0) Cdd[,1]=0 else Cdd[,2]=0; #cbind(Ps^2,exp(-Ps/mean(parms[,"K"]))); Cdd=matrix(0,1,2); 
	dRs = Rs*(parms[,"r"]*(1-((diag(1-a,nsp)+a)%*%Rs)/parms[,"K"]) - EatsR%*%(Cs*parms[,"delta"]))
	dCs = Cs*((parms[,"delta"]*parms[,"b"]*(Rs%*%Eats) + Cdd[,2])*(1-RMIlevs) - parms[,"m"]) + Cdd[,1]
	if(ListOut) return(list(c(dRs,dCs))); return(c(dRs,dCs));
}
#Expanded ODE model to account for stage structure
#Now f=g=consumption of adults relative to juves,  and  alphaP = gamma0 = juvenile maturation rate
MultiModODEmech=function(t,start,parms,specials,comps,FMI=TRUE,ListOut=TRUE){ 
	start=matrix(round(pmax(start,1e-6),8),ncol=3); Rs=start[,1]; Cs=start[,2]; Js=start[,3]; nsp=nrow(start);
	rels=pmin(parms[,"f"]*Rs/parms[,"K"],1); RMIlevs=(1-FMI)*rels%*%comps;
	specialsA=matrix(parms[,"f"],nsp,nsp)*specials; #Eats=matrix(1-FMI*rels,nsp,nsp)*specials;

	a=0.025; #Preserve same resource competition level as in main analysis
	Cdd=-parms[,c("alphaP","alphaP")]*cbind(Cs^2,pmax(1-Cs/mean(parms[,"K"]),0)); if(min(parms[,"alphaP"])<0) Cdd[,1]=0 else Cdd[,2]=0; 
	if(!FMI){
	dRs = Rs*(parms[,"r"]*(1-((diag(1-a,nsp)+a)%*%Rs)/parms[,"K"]) - specials%*%(Cs*parms[,"delta"]))
	dJs = Cs*(parms[,"delta"]*parms[,"b"]*(Rs%*%specials) + Cdd[,2])*(1-RMIlevs) - Js*parms[,"alpha"]
	dCs = Js*parms[,"alpha"] - Cs*parms[,"m"] + Cdd[,1]; }
	if(FMI){
	dJs = Rs*parms[,"r"]*(1-((diag(1-a,nsp)+a)%*%Rs)/parms[,"K"]) - Js*(specials%*%(Cs*parms[,"delta"])) - Js*parms[,"alpha"]
	dRs = Js*parms[,"alpha"] - Rs*(specialsA%*%(Cs*parms[,"delta"]))
	dCs = Cs*(parms[,"delta"]*parms[,"b"]*(Js%*%specials + Rs%*%specialsA) + Cdd[,2] - parms[,"m"]) + Cdd[,1]; }	
	if(ListOut) return(list(c(dRs,dCs,dJs))); return(c(dRs,dCs,dJs));
}



#Functions to select and implement ODE solver
Finerun=function(FUN,start,Trun,PARMS,FinAvg=round(Trun/5),NSPi=nrow(PARMS$Parms)){ 
	FUN=MultiModODE2; Juves=FALSE; if(mean(PARMS$Parms[,"alpha"])>0.05){ Juves=TRUE; FUN=MultiModODEmech; };
	colMeans(tail(lsoda(c(start,rep(0.5,NSPi*Juves)),1:Trun,FUN,parms=PARMS$Parms,specials=PARMS$Specials,comps=PARMS$Comps,FMI=PARMS$FMI,ListOut=TRUE),FinAvg))[1+(1:(2*NSPi))]
}
FinerunBEF=function(FUN,start,Trun,PARMS,FinAvg=round(Trun/5)){
	tmp=lsoda(start,1:Trun,FUN,parms=PARMS$Parms,specials=PARMS$Specials,comps=PARMS$Comps,FMI=PARMS$FMI,ListOut=TRUE)[,-1]; Fin=colMeans(tail(tmp,FinAvg));
	Conds=apply(rbind(tmp[TimesBEF,],Fin), 1, function(x) MultiModODE2(t=1,start=x,parms=PARMS$Parms,specials=PARMS$Specials,comps=PARMS$Comps,FMI=PARMS$FMI,ListOut=NULL))
	Conds[1,1:4]=NA; return(c(Fin,Conds[1,],Conds[2,]));
}


#Functions to run bifurcation analyses
gradcollapse=function(FUN,PARMS,parmVar="m",parmseq,transrun=1e3){
	NSP=nrow(PARMS$Specials); temp=matrix(as.numeric(),ncol=1+2*NSP); start=rep(c(0.3,0.8),each=NSP); Min=min(parmseq); 
	for(i in parmseq[parmseq>=Min]){
		if(i>Min) start=temp[nrow(temp),-1]; if(max(tail(start,-NSP))<0.05) break; PARMS$Parms[,parmVar]=i;
		Fins=Finerun(FUN,start,transrun,PARMS); temp=rbind(temp,c(i,pmax(Fins,0.025))); 
	}
	for(j in rev(parmseq[parmseq<i])){
		start=temp[nrow(temp),-1]; if(min(tail(start,-NSP))>2.4) break; PARMS$Parms[,parmVar]=j;
		Fins=Finerun(FUN,start,transrun,PARMS); temp=rbind(temp,c(j,pmax(Fins,0.025)));
	}
	x=parmseq; n=length(x); out=matrix(NA,2*n,2*NSP); ifin=parmseq[which(parmseq==i)-1];
	out[c(which(x==Min):which(x==ifin),n+match(tail(temp[,1],-which(temp[,1]==ifin)[1]),rev(x))),]=temp[,-1]
	return(out);
}
gradcollapseBEF=function(FUN,PARMS,parmVar="m",parmseq,transrun=1e3){
	NSP=nrow(PARMS$Specials); Min=min(parmseq); if(PARMS$FMI==1) Min=min(parmseq);
	Present=PARMS$Parms[,"delta"]!=0.05; start=rep(c(0.3,0.8),each=NSP); temp=matrix(as.numeric(), ncol=(2*NSP)-2 + 2*(length(TimesBEF)+1) + MaxAbs*Nterms + 2);
	for(i in parmseq[parmseq>=Min]){	if(i>Min) start=temp[nrow(temp),1:(2*NSP)]; if(max(tail(start,-NSP))<0.05) break; 
		PARMS$Parms[,parmVar]=i; Fins=pmax(FinerunBEF(FUN,start,transrun,PARMS),0.025); startTmp=head(Fins,2*NSP); 
		CandsI=which(tail(startTmp,NSP)>0.025 & Present); MaxAbsI=length(CandsI)-1; Fins=c(Fins,rep(NA,Nterms*(MaxAbs-MaxAbsI))); #CANT REMOVE DEAD SPPS!!!
		if(MaxAbsI<1){ temp=rbind(temp,Fins); next; };
		for(nCabsTmp in 1:MaxAbsI){ for(j in 1:removalExptReps){ startTmpJ=startTmp; 
				if(befORD==0) startTmpJ[NSP+sample(CandsI,nCabsTmp)]=0
				if(befORD==1) startTmpJ[NSP+CandsI[order(startTmpJ[CandsI])[1:nCabsTmp]]]=0
				Fins=c(Fins,tail(FinerunBEF(FUN,startTmpJ,transrun,PARMS),-2*NSP)); 
		}; };
		temp=rbind(temp,Fins)
	}
	mshort=1+length(parmseq)-which(parmseq==i); if(mshort>0) out=rbind(temp,matrix(NA,nrow=mshort,ncol=ncol(temp))); 
	return(out[,-(NSP+which(!Present))]);
}




#Functions to analyze and plot simulation results
NAzro=function(x){ x[is.na(x)]=0; return(x); }
distinctnss=function(dat,NI=TRUE,Q=0.25,divers=FALSE){
	n=nrow(dat)/2; if(NI) dat=dat[,-(1:4),]; if(divers) dat=dat>0.05; mags=apply(dat,c(1,3),sum,na.rm=TRUE); 
	dmags=head(mags,n)-apply(tail(mags,n),2,rev); dmagMax=matrix(apply(dmags,2,max,na.rm=TRUE),ncol=nreps,byrow=TRUE);
	return(t(apply(dmagMax,1,function(x) c(mean(x),quantile(x,c(Q,1-Q))))))
}
changeanSimp=function(dat,parmVar=mSeq,nreps=nreps,thresh=0.1,endsOnly=TRUE,erlyPersist=NULL,Q=0.2,SppVar=TRUE){
	n=dim(dat)[1]/2; pv2=c(parmVar,rev(parmVar))[-c(n,2*n)]; above=(1:(2*n-2))<(n-1); tr=function(x) if(length(x)==0) NA else x;
	bump=apply(dat, c(2,3), function(x) c(diff(head(x,n))<(-thresh),diff(tail(x,-n))>thresh));
	if(endsOnly) for(i in 1:dim(bump)[3]) bump[,,i]=apply(bump[,,i], 2, function(x){ x[x>0 & !((1:length(x))%in%c(tail(which(x==1 & above),1),which(x==1 & !above)[1]))]=0; return(x); }) 
	bumpSum=apply(bump,c(1,3),sum); if(!is.null(erlyPersist)) return(supsmu(rep(gamSeq,each=30),12-colSums(bumpSum[1:erlyPersist,]))); #bumpSum[1:8,]=0
	nGrads=apply(bump[above,,], 3, function(x) sum(colSums(x,na.rm=TRUE)==0))
	LastFirst=apply(bumpSum,2,function(x) c(sum((pv2*x)[above],na.rm=T)/sum(x[above],na.rm=T),sum((pv2*x)[!above],na.rm=T)/sum(x[!above],na.rm=T)))
	if(SppVar) LastFirstQs=apply(bumpSum,2,function(x) c(quantile(rep(pv2[above],NAzro(x[above])),c(Q,1-Q),na.rm=TRUE),quantile(rep(pv2[!above],NAzro(x[!above])),c(Q,1-Q),na.rm=TRUE)))
	
	nColRec=apply(bumpSum,2,function(x) c(length(na.omit(x[above & x>0])),length(na.omit(x[!above & x>0]))))
	ColSyncs=apply(bumpSum[above,], 2, function(x){ m=max(x,na.rm=T); if(m<2) return(c(0,m)); return(c(m,sum(x,na.rm=T)-m)); })
	temp=cbind(LastCol=LastFirst[1,],FirstRec=LastFirst[2,],nCol=nColRec[1,],nRec=nColRec[2,],nGrads,nColBig=ColSyncs[1,],nColNonbig=ColSyncs[2,])
	if(SppVar){ temp=cbind(temp,t(LastFirstQs)); ret=apply(temp,2,function(x) colMeans(matrix(x,nrow=nreps),na.rm=TRUE)); } else {
		Qts=cbind(t(apply(matrix(temp[,1],nrow=nreps),2,quantile,c(Q,1-Q),na.rm=TRUE)),t(apply(matrix(temp[,2],nrow=nreps),2,quantile,c(Q,1-Q),na.rm=TRUE)))
		ret=cbind(apply(temp,2,function(x) colMeans(matrix(x,nrow=nreps),na.rm=TRUE)),Qts); }
	ret[is.na(ret)]=0; return(ret);
}
changeplotSimp=function(datRMI,datFMI,Fvars=NULL,Qwant=c(0.25,0.025),give=FALSE,xlm=NULL,ylm=NULL){
	if(is.null(Fvars)){ xax=gamSeq; Xlab="Connectance"; } else { xax=Fvars; Xlab="Diversity Factor"; }
	FMIan=FMIanU=changeanSimp(datFMI,mSeq,nreps,Q=Qwant[1]); RMIan=RMIanU=changeanSimp(datRMI,mSeq,nreps,Q=Qwant[1]);
	polyPlot=function(ypts,COL) polygon(cbind(c(xax,rev(xax)),c(ypts[,1]-0.003,rev(ypts[,2])+0.003)),col=c(rgb(red=0,green=0,blue=0,alpha=0.1),rgb(red=1,green=0,blue=0,alpha=0.1))[COL],border=0)
	HystPlotD=function(x,COL){ if(is.null(ylm)) ylm=c(quantile(x[,3:6],0.02,na.rm=TRUE),0.01+max(x[,1:2],na.rm=TRUE)); 
		matplot(xax,x[,1:2],type="l",lty=1:2,lwd=3,col=COL,ylab="Consumer mortality",xlab=Xlab,xlim=xlm,ylim=ylm); polyPlot(x[,3:4],COL); polyPlot(x[,5:6],COL); }

	par(mfrow=c(2,2))
	HystPlotD(FMIan[,c(1:2,8:11)],2); if(length(Qwant)==2){ FMIanU=changeanSimp(datFMI,mSeq,nreps,Q=Qwant[2]); RMIanU=changeanSimp(datRMI,mSeq,nreps,Q=Qwant[2]); }; Uplt=(min(FMIan[,8:11]==FMIanU[,8:11],na.rm=TRUE)==0);
	matlines(xax,FMIanU[,8:9],lty=1,col=2*Uplt); HystPlotD(RMIan[,c(1:2,8:11)],1); matlines(xax,RMIanU[,8:9],lty=1,col=1*Uplt); 
	MagDivers=cbind(distinctnss(datFMI,NI=FALSE,Q=Qwant[1],divers=TRUE),distinctnss(datRMI,NI=FALSE,Q=Qwant[1],divers=TRUE));
	MagAbunds=cbind(distinctnss(datFMI,NI=FALSE,Q=Qwant[1]),distinctnss(datRMI,NI=FALSE,Q=Qwant[1]))
	matplot(xax,MagDivers[,c(1,4)],ylim=c(0,max(MagDivers[,c(1,4)])),type="l",lwd=2,col=2:1,lty=1,xlim=xlm,xlab="State distinctiveness in consumer species richness"); polyPlot(MagDivers[,2:3],2); polyPlot(MagDivers[,5:6],1);
	matplot(xax,MagAbunds[,c(1,4)],ylim=c(0,max(MagAbunds[,c(1,4)])),type="l",lwd=2,col=2:1,lty=1,xlab="State distinctiveness in consumer abundance",xlim=xlm); polyPlot(MagAbunds[,2:3],2); polyPlot(MagAbunds[,5:6],1);
	if(give) return(list(FMI=FMIan,RMI=RMIan))
}
BfnVisSimp=function(FMIdat,RMIdat,trl=3,gamLvl=10,xlm=c(0.05,0.35),Lwd=c(2,3),lty=1:2,Cols=sample(1:8)){
	FMItrans=FMIdat[,,nreps*gamLvl+trl]; RMItrans=RMIdat[,,nreps*gamLvl+trl];
	ymxs=c(max(rbind(FMItrans,RMItrans),na.rm=TRUE),max(rowSums(rbind(FMItrans,RMItrans),na.rm=TRUE)))
	mplot=function(dat){
		matplot(mSeq,dat[1:ntests[2],],xlim=xlm,ylim=c(0,ymxs[1]),type="l",col=Cols,lty=lty[1],lwd=Lwd[1],xlab="Consumer mortality",ylab="Consumer abundances")
		matlines(mSeq,dat[(2*ntests[2]):(ntests[2]+1),],col=Cols,lty=lty[2],lwd=Lwd[1])
		matplot(mSeq,rowSums(dat[1:ntests[2],]),xlim=xlm,ylim=c(0,ymxs[2]),type="l",col=1,lty=lty[1],lwd=Lwd[2],xlab="",ylab="Total Consumer abundance")
		matlines(mSeq,rowSums(dat[(2*ntests[2]):(ntests[2]+1),]),col=4,lty=lty[2],lwd=Lwd[2])
	}
	par(mfcol=c(2,2)); mplot(FMItrans); mplot(RMItrans);
}
bef2ColsFun=function(dat,TAU){
	x=dat[,1:10,]; L=dim(x)[1]; mDown=t(apply(x[1:L,,], c(2,3), bumpfind, thresh=TH, x=1:L))
	maxCol=apply(mDown,1,function(x) which.max(tail(binsum(x,breaks=L,rng=c(0,L)),L))); maxCol[maxCol==1]=round(rowMeans(mDown,na.rm=TRUE))[maxCol==1]; #this will be NA if no hystereses
	if(TAU>1e2) TAU=maxCol-(TAU-1e2); 
	
	sumfun=function(x) c(c(apply(x,2,median,na.rm=TRUE)[-1],median(x[,2]/x[,1],na.rm=TRUE))[c(1:3,5,4)],quantile(x[,5],c(0.25,0.75),na.rm=TRUE),quantile(x[,2],c(0.25,0.75),na.rm=TRUE),quantile(x[,2]/x[,1],c(0.25,0.75),na.rm=TRUE))
	ats0=t(sapply(11:dim(dat)[2],function(y) diag(dat[pmax(maxCol-TAU,1),y,]))); ats=ats0[-(1:5),];
	summ=sumfun(t(ats0[1:5,])); for(nAbs in 1:9) summ=rbind(sumfun(matrix(as.numeric(ats[(nAbs-1)*5*removalExptReps + 1:(5*removalExptReps),]),ncol=5,byrow=TRUE)),summ);
	colnames(summ)=c("Csns","Cinter","Cintra","FeedRel","Fun","FQL","FQU","SQL","SQU","RQL","RQU"); return(summ);
}
changeplotDbef2Simp=function(datbef2FMI,datbef2RMI,gamLvl=ntests[1]-1,Qwant=0.3,NI=TRUE,ylm=c(0,0.85)){
	#Entries in dim(befXMI)[2] are Finals w/o removals and c(mean(Cgrowth) over timesBEF+Trun,mean(Function) over timesBEF+Trun) iterated over (1) replicates and (2) extinction levels 1:9.
	#bef2 has 3 timesBEF and replaces mean(Cgrowth) over timesBEF+Trun with Rsns,Csns,Cinter,Cintra at Trun
	KEEPS=c(rep(c(FALSE,TRUE),c(12,10)),  rep(rep(c(TRUE,FALSE,TRUE),c(4,3,1)),1+removalExptReps*9)); #Toss Rs and mean(Function) over timesBEF
	TRLs=(gamLvl-1)*nreps+(1:nreps); #select connectance level
	FMIan=bef2ColsFun(datbef2FMI[,KEEPS,TRLs],head(TAU,1)); RMIan=bef2ColsFun(datbef2RMI[,KEEPS,TRLs],tail(TAU,1));
	matplot(cbind(FMIan[,5:7],RMIan[,5:7]),type="l",lty=c(1,3,3,1,3,3),lwd=c(3,2,2,3,2,2),col=c(2,8,8,1,8,8),ylim=ylm);
	if(NI==2){ matplot(FMIan[,c(4,10:11)],col=2,lwd=c(3,2,2),lty=c(1,3,3),type="l"); par(new=TRUE); matplot(RMIan[,c(4,10:11)],col=1,lwd=c(3,2,2),lty=c(1,3,3),type="l",xaxt="n",yaxt="n",xlab="",ylab=""); axis(4); }
	if(NI==3) matplot(cbind(FMIan[,c(1,8:9)],RMIan[,c(1,8:9)]),type="l",lty=c(1,3,3,1,3,3),lwd=c(3,2,2,3,2,2),col=c(2,8,8,1,8,8))
	Es=round(seq(1,8,len=20),2); print(Es[c(which.min(sapply(Es,function(y) cv(diff(nrm(RMIan[,5])^y)))),which.min(sapply(Es,function(y) cv(diff(nrm(FMIan[,5])^y)))))]);
	return(list(FMI=FMIan,RMI=RMIan))
}




#Functions setting up food webs
parmsvar=c("r","K","b","delta","alphaP")
rlnormT=function(n,mu,sig) rlnorm(n, meanlog=log(mu/sqrt(1+(sig/mu)^2)), sdlog=sqrt(log(1+(sig/mu)^2)) )
parmgen=function(parms0,FMI,NSP=12,Seed=1,extras=c(nCabs=1,Feedvar=0.3,ISvar=0.25,QvarTrans=0.125)){
	Parms=t(matrix(parms0,length(parms0),NSP)); set.seed(Seed); colnames(Parms)=names(parms0);
	ParmsDevs=parms0["fvar"]*matrix(runif(NSP*length(parms0),-1,1),nrow=NSP)
	Parms[,parmsvar]=(1+ParmsDevs[,1:length(parmsvar)])%*%diag(parms0[parmsvar])
	Parms[,"f"]=parms0["f"]*(1+extras["QvarTrans"]*ParmsDevs[,length(parms0)])
	set.seed(Seed+1); Cabs=sample(1:NSP,extras["nCabs"]); Parms[Cabs,"delta"]=0.05;
	set.seed(Seed+2); Fs=t(VarWeb(1-extras["Feedvar"],NSP,NSP)); set.seed(Seed+3); 
	Ss=VarWeb3(NSP,C=Parms[1,"gamma"],precis=0.05); set.seed(Seed+4); Ss=apply((Ss>0)*matrix(rlnormT(NSP^2,1,extras["ISvar"]),NSP,NSP),2,function(x) x/sum(x));
	return(list(Parms=round(Parms,3),Specials=round(Ss,5),Comps=round(Fs,5),FMI=FMI,Cabs=Cabs))
}
parmtest=function(PARMS,NSP=nrow(PARMS$Parms)){
	PARMS$Parms[,"m"]=0.02; His=(tail(Finerun(MultiModODE2,rep(c(0.5,0.25),each=NSP),4e2,PARMS),NSP)>0.25)[-PARMS$Cabs];
	return(sum(His))
}
parmgenValidSimp=function(parms0,FMI,NSP=12,Seed=1,triesMax=5e2,extras=c(nCabs=1,Feedvar=0.3,ISvar=0.25,QvarTrans=0.125),outseed=FALSE){
	set.seed(Seed); SeedsGen=rnorm(triesMax,1e6,1e6); tryn=1; nCpres=1[-1]; statisfied=FALSE;
	while(tryn<=triesMax & !statisfied){
		NP=tryCatch(parmtest(parmgen(Seed=SeedsGen[tryn],parms0,FMI,NSP,extras)), error=function(e) 0)
		nCpres=c(nCpres, NP); statisfied=NP==(NSP-extras["nCabs"]); tryn=tryn+1;
	}; 
	SeedBest=SeedsGen[tryn-1]; if(tryn>triesMax) SeedBest=SeedsGen[which.max(nCpres)];
	if(outseed) return(SeedBest); return(parmgen(parms0,FMI,NSP,SeedBest,extras))
}
#Generate matrix of recrutiment interactions
VarWeb=function(Y,k,N=k,smudge=0.025){ #Upper bound on Y is 0.915
	Yn=pmin(Y,.915); scl=1/pmax(2.6^((0.743-Yn)/0.149)-0.3,1e-3);
	out=matrix(nrow=0,ncol=k); for(i in 1:1e3){ 
		if(nrow(out)>N) return(out[1:N,])
		temp0=matrix(rgamma(k*1e3,scl),ncol=k,byrow=TRUE); temp=temp0/as.vector(temp0%*%rep(1,k));
		gens=1-apply(temp,1,max); out=rbind(out, temp[which(gens>(Yn-smudge) & gens<(Yn+smudge)),]);
	}
}
#Generate trophic connectance matrix
VarWeb3=function(nsp,Y=NULL,C=pmin(Y+1/nsp,1),precis=0.025){ #Need either C or Y; Y superseeds.
	if(C==1/nsp) return(diag(nsp)[sample(1:nsp),]);
	cmat=matrix(0,nsp,nsp); n=0; while(n<1e5 & (min(c(colSums(cmat),rowSums(cmat)))==0 | abs(mean(cmat)-C)>precis)){ n=n+1; cmat=matrix(rbinom(nsp^2,1,C),nsp,nsp); };
	return(apply(cmat,2,function(x) x/sum(x)));
}


runBif=function(name,Torun,Reps=1:120,triesMax=200){
	model=floor(Torun["mType"]); traitDiversity=model==2; stageStructured=model==3; removalExpt=model==4;
	parms0[c("alphaP","delta","f","K")]=Torun[c("B","delta","f","Ki")];
	if(stageStructured) parms0[c("alpha","alphaP")]=c(Torun["B"], 0.075)
	NSP=12; FMI=Torun["Feedvar"]==0; set.seed(11); Seeds=round(rnorm(max(Reps),1e6,1e6)); 

	Fvars=Torun["Fvar"]; if(traitDiversity){ Fvars=FvarsDiversity; gamSeq=gamDiversity; };
	Outs=Outs0=array(dim=c(2*length(mSeq),NSP-Torun["nCabs"],0)); rownames(Outs)=c(mSeq,rev(mSeq));
	if(removalExpt){ 
		assign("MaxAbs",NSP-Torun["nCabs"]-1,envir=.GlobalEnv); assign("TimesBEF",c(5,20,40),envir=.GlobalEnv); 
		assign("Nterms",removalExptReps*2*(length(TimesBEF)+1),envir=.GlobalEnv);
		Outs=Outs0=array(dim=c(length(mSeq), (2*NSP)-2 + 2*(length(TimesBEF)+1) + MaxAbs*Nterms, 0)); 
	}

	#Implement and save:
	errfun=function(x) tryCatch(x, error=function(e) matrix(NA,nrow(Outs),ncol(Outs)))
	PARTS=length(Reps)>20 #If many replicates, reduce memory use by splitting up results and them reassembling at end
	for(j in Fvars) for(i in gamSeq) for(k in Reps){
		parms0["gamma"]=i; parms0["fvar"]=j; extras=c(nCabs=as.numeric(Torun["nCabs"]),Feedvar=as.numeric(Torun["Feedvar"]),ISvar=as.numeric(parms0["ISvar"]),QvarTrans=as.numeric(Torun["QVT"]));
		PARMS=parmgenValidSimp(parms0,FMI,NSP,Seed=Seeds[k],extras,triesMax=triesMax);
		collapsefun=list(gradcollapse,gradcollapseBEF)[[(1+removalExpt)]]
		res=errfun(collapsefun(MultiModODE2,parmVar="m",parmseq=mSeq,transrun=Trun,PARMS=PARMS))
		if(!removalExpt) res=res[,-c(1:NSP,NSP+PARMS$Cabs)]; Outs=abind(Outs,res,along=3);
		print(c(i,k,j)); if(PARTS & k==max(Reps)){ saveRDS(Outs, paste0("DivAlts0.5",name,"_",i,"_",j,".rds")); Outs=Outs0; };
	}

	if(PARTS){ Outs=Outs0; for(z in Fvars) for(y in gamSeq) Outs=abind(Outs,readRDS(paste0("DivAlts0.5",name,"_",y,"_",z,".rds"))); }
	saveRDS(Outs, paste0("DivAlts0.5",name,".rds"))
	assign(name,Outs,envir=.GlobalEnv)
}




#Parameters used in simulations
parms0=c(K=1.35,r=1,delta=1.1,b=1,m=0.05,gamma=0,alpha=0.025,alphaP=0.15,f=1.2,fvar=0.15,ISvar=0.5)
ntests=c(14,160); 
nreps=120; #Number of replicate food webs to examine
triesMax=500; #In setting up each replicate, max number of candidate food web parameterizations to examine before settling for food web parameterization with the most consumers present at low mortality. Rarely used but controls run time.
Trun=200 #Duration of each ODE implementation

#Note: under default settings each simulation taks 1-3 days on a 2-core CPU.
#For faster times suggest using: nreps=15; ntests[1]=6; triesMax=35; removalExptReps=2;


gamSeq=seq(0.1,0.85,length=ntests[1]); #Connectance levels to examine
mSeq=seq(0.025,1.65,length=ntests[2]); #Range and resolution of mortality values to scan
FvarsDiversity=c(0.01,seq(0.05,0.55,len=ntests[1]-1)); gamDiversity=0.25;  #Diversity levels to examine
removalExptReps=10; befORD=0;  #befORD: 0 identity of removed consumer chosen randomly extinctions, for 1 species removed is the least-abundant consumer






#Run simulations and figures:

#####Figure 3: connectance effect
runBif("RMIconn",c(Feedvar=0.35,Fvar=0.15,delta=1.1,B=0.12,f=1.2,Ki=1.35,nCabs=2,QVT=0.33, mType=1), Reps=1:nreps,triesMax=triesMax)
runBif("FMIconn",c(Feedvar=0,   Fvar=0.15,delta=1.1,B=0.12,f=1.2,Ki=1.35,nCabs=2,QVT=0.33, mType=1), Reps=1:nreps,triesMax=triesMax)
changeplotSimp(Qwant=c(0.25,0.025),xlm=c(0.1,0.83), datRMI=RMIconn, datFMI=FMIconn)


#####Figure 2: example transitions
BfnVisSimp(FMIconn,RMIconn, 37,11, c(0.025,0.28),Lwd=c(1,3)) #Note: due to seed differences, replication of this figure is not exact


#####Figure S4: connectance effect without consumer density dependence
runBif("RMIconnNoDD",c(Feedvar=0.35,Fvar=0.1,delta=1.1,B=0,f=1.25,Ki=1.35,nCabs=2,QVT=0.33, mType=1), Reps=1:nreps,triesMax=triesMax)
runBif("FMIconnNoDD",c(Feedvar=0,   Fvar=0.1,delta=1.1,B=0,f=1.25,Ki=1.35,nCabs=2,QVT=0.33, mType=1), Reps=1:nreps,triesMax=triesMax)
changeplotSimp(Qwant=c(0.25),xlm=c(0.1,0.83), datRMI=RMIconnNoDD, datFMI=FMIconnNoDD)


#####Figure S2: connectance effect in stage-explicit model
runBif("RMIconnSS",c(Feedvar=0.35,Fvar=0.4,delta=1.15,B=0.3,f=1.5,Ki=1.5,nCabs=4,QVT=0.5, mType=3), Reps=1:nreps,triesMax=triesMax)
runBif("FMIconnSS",c(Feedvar=0,   Fvar=0.4,delta=1.2,B=0.5,f=0.02,Ki=1.5,nCabs=4,QVT=0.5, mType=3), Reps=1:nreps,triesMax=triesMax)
changeplotSimp(Qwant=c(0.25),xlm=c(0.1,0.83), datRMI=RMIconnSS, datFMI=FMIconnSS)


#####Figure 4: trait diversity effect
runBif("RMIdivers",c(Feedvar=0.35,Fvar=0.15,delta=1.1,B=0.15,f=1.25,Ki=1.35,nCabs=2,QVT=0.33, mType=2), Reps=1:nreps,triesMax=triesMax)
runBif("FMIdivers",c(Feedvar=0,   Fvar=0.15,delta=1.1,B=0.15,f=1.25,Ki=1.35,nCabs=2,QVT=0.33, mType=2), Reps=1:nreps,triesMax=triesMax)
changeplotSimp(Qwant=0.25,Fvars=FvarsDiversity, datRMI=RMIdivers, datFMI=FMIdivers)


#####Figure 5a: consumer species removal effect on guild-wide recruitment, resource palatability
runBif("RMIbef",c(Feedvar=0.35,Fvar=0.35,delta=1.1,B=0.12,f=1.15,Ki=1.35,nCabs=2,QVT=0.33, mType=4), Reps=1:nreps,triesMax=triesMax)
runBif("FMIbef",c(Feedvar=0,   Fvar=0.35,delta=1.1,B=0.12,f=1.15,Ki=1.35,nCabs=2,QVT=0.33, mType=4), Reps=1:nreps,triesMax=triesMax)
TAU=101; x=changeplotDbef2Simp(Qwant=0.1,gamLvl=12, datbef2RMI=RMIbef, datbef2FMI=FMIbef); 
1-1/c((x[[1]][10,5]/x[[1]][5,5]),  (x[[2]][10,5]/x[[2]][5,5])) #Relative effects of losing 5 consumer species


#####Figure 5b: visualizing collapse cascades
fallFun=function(SEED=11253,FMI=FALSE,simTrun=600,onlyCs=FALSE){
	parms0["gamma"]=0.5; Feedvar=0.35*(!FMI); finess=400; mSeq=seq(0.015,1.2,length=finess); NSP=12;
	PARMS=parmgenValidSimp(parms0,FMI,NSP,Seed=SEED,extras=c(nCabs=2,Feedvar=Feedvar,ISvar=0.5,QvarTrans=0.33),triesMax=100);
	Outs0=gradcollapse(MultiModODE2,parmVar="m",parmseq=mSeq,transrun=Trun,PARMS); dat=Outs0[,-c(1:NSP,NSP+PARMS$Cabs)];
	Th=0.1; ASS=(head(dat,finess)-apply(tail(dat,finess),2,rev))>Th; 
	if(max(ASS,na.rm=TRUE)==0) return(fallFun(33*SEED,FMI,simTrun,onlyCs)); #Re-run analysis in rare FMI cases when a seed produces no ASS
	BC=tail(which(rowSums(ASS)>1),1); DATP=PARMS; DATP$Parms[,"m"]=mSeq[BC]+diff(mSeq)[1]; 
	TS=lsoda(pmax(Outs0[BC,],0.025),1:simTrun,MultiModODE2,parms=DATP$Parms,specials=DATP$Specials,comps=DATP$Comps,FMI=DATP$FMI,ListOut=TRUE)[,-1]
	TS[,1:NSP]=TS[,1:NSP]%*%diag(1/DATP$Parms[,"K"]); return(TS[,-c((1:NSP)*onlyCs,NSP+DATP$Cabs)]);
}
fallAn=function(dc,Th=0.05,tmax=nrow(dc),trng=FALSE){ 
	colT=apply(dc<Th,2,function(x) c(which(x),tmax)[1]); colT[colT==1 | colT==tmax]=NA;    
	if(trng) return(diff(range(colT,na.rm=TRUE)));    
	c(sort(colT-min(colT,na.rm=TRUE)),rep(NA,sum(is.na(colT)))); 
}
set.seed(11); Seeds=runif(2e2,1,1e5); simTrun=600; library(doParallel); library(foreach); registerDoParallel(cores=6);   #storeXMIhigamSave; storeXMISave
storeFMIsave=foreach(i=Seeds,.combine=function(a,b) abind(a,b,along=3),.packages=c("abind","deSolve"),.export=ls()) %dopar% fallFun(i,FMI=TRUE,simTrun)
storeRMIsave=foreach(i=Seeds,.combine=function(a,b) abind(a,b,along=3),.packages=c("abind","deSolve"),.export=ls()) %dopar% fallFun(i,FMI=FALSE,simTrun)
stopImplicitCluster(); TH=0.028; FT=apply(storeFMIsave[,-(1:NSP),],3,fallAn,Th=TH); RT=apply(storeRMIsave[,-(1:NSP),],3,fallAn,Th=TH);
show=2:9; matplot(show, t(rbind(apply(log(FT),1,quantile,c(0.5,0.25,0.75),na.rm=TRUE),apply(log(RT),1,quantile,c(0.5,0.25,0.75),na.rm=TRUE)))[show,], type="l",col=c(2,8,8,1,8,8),lwd=c(3,2,2,3,2,2),lty=c(1,2,2,1,3,3), xlab="Order in extinction",ylab="Extinction time (since 1st extinct)")
matplot((show)-0.1,log(FT[show,]),col=2,pch=1); matpoints((show)+0.1,log(RT[show,]),col=1,pch=1)


#####Figure 5c,d: Plotting example transitions:
trl=1; xmx=200; par(mfrow=1:2); matplot(storeFMIsave[1:xmx,,trl],type="l",col=rep(3:4,c(12,10))); matplot(storeRMIsave[1:xmx,,trl],type="l",col=rep(3:4,c(12,10)));














