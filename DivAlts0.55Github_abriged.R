

#Note: this code intended for reproducibility only; for model and analysis details see original manuscript
# Throughout, the script terms specialized feedbacks as FMI (Feeding-Mediating Interactions) and aggregate feedbacks as RMI (Recruitment-Mediating Interactions)
# Thus, FMI=TRUE denotes specialized feedback model and FMI=FALSE denotes aggregate feedback model
# Also note that models here distinguish consumer diet specialization matrix ("specials") from baseline consumer species-specific consumption rate delta
# Throughout, parameter alphaP regulating consumer density dependence corresponds to parameter beta in the main text
library(deSolve); library(abind);



#Core ODE model.
MultiModODE2=function(t,start,parms,specials,comps,FMI=TRUE,ListOut=TRUE){
	#Ensure non-negative inputs
	start=matrix(round(pmax(start,0),8),ncol=2); Rs=start[,1]; Cs=start[,2]; nsp=nrow(start); a=parms[1,"alpha"];
	#Track each resource's abundance effect on edibility/recruitment
	rels=parms[,"f"]*Rs/parms[,"K"]
	Eats=matrix(1-FMI*rels,nsp,nsp)*specials; RMIlevs=(1-FMI)*rels%*%comps;
	#For Figure 5ab, return community functioning metrics (consumer recruitment & resource uptake)
	if(is.null(ListOut)){ 
		Crs=parms[,"delta"]*parms[,"b"]*colSums(diag(Rs)%*%Eats)*(1-RMIlevs)
		return(c(mean(Crs), c(mean(1-RMIlevs),mean(colSums(Eats)))[1+FMI]))
	}
	dRs = Rs*(parms[,"r"]*(1-((diag(1-a,nsp)+a)%*%Rs)/parms[,"K"]) - Eats%*%(Cs*parms[,"delta"]))
	dCs = Cs*((parms[,"delta"]*parms[,"b"]*(Rs%*%Eats))*(1-RMIlevs) - parms[,"m"]) - parms[,"alphaP"]*Cs^2
	if(ListOut) return(list(c(dRs,dCs))); return(c(dRs,dCs));
}
#Expanded ODE model to account for stage structure
#For this model only, f in code = g in text = consumption of adults relative to juveniles,
#and  alpha in code = gamma0 in text = juvenile maturation rate
MultiModODEmech=function(t,start,parms,specials,comps,FMI=TRUE,ListOut=TRUE){ 
	start=matrix(round(pmax(start,1e-6),8),ncol=3); Rs=start[,1]; Cs=start[,2]; Js=start[,3]; nsp=nrow(start);
	rels=pmin(parms[,"f"]*Rs/parms[,"K"],1); RMIlevs=(1-FMI)*rels%*%comps;
	specialsA=matrix(parms[,"f"],nsp,nsp)*specials;
	a=0.025; #Preserve same resource competition level as in main analysis
	if(!FMI){
	dRs = Rs*(parms[,"r"]*(1-((diag(1-a,nsp)+a)%*%Rs)/parms[,"K"]) - specials%*%(Cs*parms[,"delta"]))
	dJs = Cs*(parms[,"delta"]*parms[,"b"]*(Rs%*%specials))*(1-RMIlevs) - Js*parms[,"alpha"]
	dCs = Js*parms[,"alpha"] - Cs*parms[,"m"] - parms[,"alphaP"]*Cs^2; }
	if(FMI){
	dJs = Rs*parms[,"r"]*(1-((diag(1-a,nsp)+a)%*%Rs)/parms[,"K"]) - Js*(specials%*%(Cs*parms[,"delta"])) - Js*parms[,"alpha"]
	dRs = Js*parms[,"alpha"] - Rs*(specialsA%*%(Cs*parms[,"delta"]))
	dCs = Cs*(parms[,"delta"]*parms[,"b"]*(Js%*%specials + Rs%*%specialsA) - parms[,"m"]) - parms[,"alphaP"]*Cs^2; }	
	if(ListOut) return(list(c(dRs,dCs,dJs))); return(c(dRs,dCs,dJs));
}



#Functions to select and implement ODE solver
Finerun=function(FUN,start,Trun,PARMS,FinAvg=round(Trun/5),NSPi=nrow(PARMS$Parms)){ 
	Juves=FALSE; if(mean(PARMS$Parms[,"alpha"])>0.05){ Juves=TRUE; FUN=MultiModODEmech; };
	colMeans(tail(lsoda(c(start,rep(0.5,NSPi*Juves)),1:Trun,FUN,parms=PARMS$Parms,specials=PARMS$Specials,comps=PARMS$Comps,FMI=PARMS$FMI,ListOut=TRUE),FinAvg))[1+(1:(2*NSPi))]
}
#For Figure 5ab analyses, additionally track community functioning metrics (consumer recruitment & resource uptake)
FinerunBEF=function(FUN,start,Trun,PARMS,FinAvg=round(Trun/5)){
	tmp=lsoda(start,1:Trun,FUN,parms=PARMS$Parms,specials=PARMS$Specials,comps=PARMS$Comps,FMI=PARMS$FMI,ListOut=TRUE)[,-1]; Fin=colMeans(tail(tmp,FinAvg));
	Conds=apply(rbind(tmp[TimesBEF,],Fin), 1, function(x) MultiModODE2(t=1,start=x,parms=PARMS$Parms,specials=PARMS$Specials,comps=PARMS$Comps,FMI=PARMS$FMI,ListOut=NULL))
	Conds[1,1:4]=NA; return(c(Fin,Conds[1,],Conds[2,]));
}


#Functions to run bifurcation analyses
#1 - ramp consumer collapse up & down to detect collapse and recovery points
#parmVar is bifurcation parameter (m throughout manuscript), parmseq are the levels of parmVar to examine, 
#and transrun is the simulation length at each parameter level
gradcollapse=function(FUN,PARMS,parmVar="m",parmseq,transrun=1e3){
	NSP=nrow(PARMS$Specials); temp=matrix(as.numeric(),ncol=1+2*NSP); start=rep(c(0.3,0.8),each=NSP); Min=min(parmseq); 
	#Gradually ramp up consumer mortality until see consumer collapse
	for(i in parmseq[parmseq>=Min]){
		#Stop consumer mortality ramp-up if all consumer densities <0.05
		if(i>Min) start=temp[nrow(temp),-1]; if(max(tail(start,-NSP))<0.05) break;
		PARMS$Parms[,parmVar]=i; Fins=Finerun(FUN,start,transrun,PARMS); temp=rbind(temp,c(i,pmax(Fins,0.025))); 
	}
	#Gradually ramp down consumer mortality
	for(j in rev(parmseq[parmseq<i])){
		start=temp[nrow(temp),-1]; if(min(tail(start,-NSP))>2.4) break; PARMS$Parms[,parmVar]=j;
		Fins=Finerun(FUN,start,transrun,PARMS); temp=rbind(temp,c(j,pmax(Fins,0.025)));
	}
	x=parmseq; n=length(x); out=matrix(NA,2*n,2*NSP); ifin=parmseq[which(parmseq==i)-1];
	out[c(which(x==Min):which(x==ifin),n+match(tail(temp[,1],-which(temp[,1]==ifin)[1]),rev(x))),]=temp[,-1]
	return(out);
}
#2 - repeats ramp-up analysis in gradcollapse for a range of consumer species richness levels
gradcollapseBEF=function(FUN,PARMS,parmVar="m",parmseq,transrun=1e3){
	NSP=nrow(PARMS$Specials); Min=min(parmseq); if(PARMS$FMI==1) Min=min(parmseq);
	start=rep(c(0.3,0.8),each=NSP); temp=matrix(as.numeric(), ncol=(2*NSP)-2 + 2*(length(TimesBEF)+1) + MaxAbs*Nterms + 2);
	#Track any consumers absent throughout analysis (i.e., k=ncol(PARMS$Parms)-sum(Present)):
	Present=PARMS$Parms[,"delta"]!=0.05; 
	for(i in parmseq[parmseq>=Min]){	
		#Stop consumer mortality ramp-up if all consumer densities <0.05
		if(i>Min) start=temp[nrow(temp),1:(2*NSP)]; if(max(tail(start,-NSP))<0.05) break; 
		PARMS$Parms[,parmVar]=i; Fins=pmax(FinerunBEF(FUN,start,transrun,PARMS),0.025); startTmp=head(Fins,2*NSP); 
		#Find which consumers are present and candidantes for removal to test effect of richness. Ensure at least 1 consumer species remains.
		CandsI=which(tail(startTmp,NSP)>0.025 & Present); MaxAbsI=length(CandsI)-1; Fins=c(Fins,rep(NA,Nterms*(MaxAbs-MaxAbsI)));
		if(MaxAbsI<1){ temp=rbind(temp,Fins); next; };
		#Iterate across consumer species richness levels
		for(nCabsTmp in 1:MaxAbsI){
			#Iterate across consumer removal replicates removalExptReps. Remove either randomly (befORD=0) or rarer consumers first (befORD=1)
			for(j in 1:removalExptReps){ startTmpJ=startTmp; 
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
#Measure distinctiveness of states in terms of consumer abundance (divers=FALSE) or species richness (divers=TRUE)
distinctnss=function(dat,NI=TRUE,Q=0.25,divers=FALSE){
	n=nrow(dat)/2; if(NI) dat=dat[,-(1:4),]; if(divers) dat=dat>0.05; mags=apply(dat,c(1,3),sum,na.rm=TRUE); 
	dmags=head(mags,n)-apply(tail(mags,n),2,rev); dmagMax=matrix(apply(dmags,2,max,na.rm=TRUE),ncol=nreps,byrow=TRUE);
	return(t(apply(dmagMax,1,function(x) c(mean(x),quantile(x,c(Q,1-Q))))))
}
#We analyze a series of bifurcation analyses (made via gradcollapse) across levels of trait diversity or connectivity. These are provided in the dat object.
#changeanSimp analyzes these results to quantify state transitions and distinctiveness across levels of diversity or connectivity.
changeanSimp=function(dat,parmVar=mSeq,nreps=nreps,thresh=0.1,Q=0.2){
	#In dat, bifurcation analyses are laid out with (1) a ramp-up in mortality to the point where all consumers are absent (rows 1:n) 
	#followed by (2) a ramp-down in mortality (rows n:2*n). above=TRUE denotes the simulation is from the ramp-up portion preceeding loss of all remaining consumers.
	n=dim(dat)[1]/2; pv2=c(parmVar,rev(parmVar))[-c(n,2*n)]; above=(1:(2*n-2))<(n-1); tr=function(x) if(length(x)==0) NA else x;

	#Find consumer collapses and recoveries.
	#For this, first find "bumps" - any changes in abundance > thresh in one step of parameter change.
	bump=apply(dat, c(2,3), function(x) c(diff(head(x,n))<(-thresh),diff(tail(x,-n))>thresh))
	#Sudden changes in abundance are not exclusive to species' extinction / recovery events. Therefore:
	# - For collapses as mortality ramps up (rows above), remove (make zero) any bumps that happen before species' final bump
	# - For recoveries as mortality ramps down (rows !above), remove (make zero) any bumps that happen after species' first bump
	for(i in 1:dim(bump)[3]) bump[,,i]=apply(bump[,,i], 2, function(x){ x[x>0 & !((1:length(x))%in%c(tail(which(x==1 & above),1),which(x==1 & !above)[1]))]=0; return(x); }) 

	#Count number of species collapsing and recovering at every mortality level change. Species that do neiether went extinct / recovered gradually (=nGrads).
	bumpSum=apply(bump,c(1,3),sum); nGrads=apply(bump[above,,], 3, function(x) sum(colSums(x,na.rm=TRUE)==0));
	#Track mean and quantile mortality levels at which (1) consumers collapsed and (2) consumers recovered.
	LastFirst=apply(bumpSum,2,function(x) c(sum((pv2*x)[above],na.rm=T)/sum(x[above],na.rm=T),sum((pv2*x)[!above],na.rm=T)/sum(x[!above],na.rm=T)))
	LastFirstQs=apply(bumpSum,2,function(x) c(quantile(rep(pv2[above],NAzro(x[above])),c(Q,1-Q),na.rm=TRUE),quantile(rep(pv2[!above],NAzro(x[!above])),c(Q,1-Q),na.rm=TRUE)))
	
	nColRec=apply(bumpSum,2,function(x) c(length(na.omit(x[above & x>0])),length(na.omit(x[!above & x>0]))))
	ColSyncs=apply(bumpSum[above,], 2, function(x){ m=max(x,na.rm=T); if(m<2) return(c(0,m)); return(c(m,sum(x,na.rm=T)-m)); })
	temp=cbind(LastCol=LastFirst[1,],FirstRec=LastFirst[2,],nCol=nColRec[1,],nRec=nColRec[2,],nGrads,nColBig=ColSyncs[1,],nColNonbig=ColSyncs[2,])
	temp=cbind(temp,t(LastFirstQs)); ret=apply(temp,2,function(x) colMeans(matrix(x,nrow=nreps),na.rm=TRUE));
	ret[is.na(ret)]=0; return(ret);
}
#Make Figure 3-4 plots of results from changeanSimp
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
#Visualize bifurcation transitions at a single level of connectance and trait diversity
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

#This part of the code is for Figure 5ab. Some elements ("NI", sns, inter, intra) refer to net species interactions metrics nor used in the final MS.
bumpfind=function(y,x=mSeq,xlag=0,fwd=TRUE,thresh){
	if(is.na(y[2])) return(NA)
	theOne=function(x) c(head(x,1),tail(x,1))[1+fwd]
	d=(1-2*fwd)*diff(y)>thresh;  if(sum(d,na.rm=TRUE)>0) return(x[theOne(which(d)-xlag)]);
	return(x[which.max((1-2*fwd)*diff(y))-xlag]);
}
binsum=function(x,freqs=rep(1,length(x)),breaks=50,rng=rngbin0){
	if(sum(!is.na(x))==0) return(double(2*length(x)));
	nzro=freqs>0; if(sum(nzro,na.rm=TRUE)==0) return(double(2*breaks)); if(length(unique(round(x[nzro],4)))==1) return(c(x[nzro][1],double(breaks-1),1,double(breaks-1)));
	if(!is.null(rng)){ x=c(rng[1],pmin(x,rng[2]),rng[2]); freqs=c(0,freqs,0); }
	tmp=unlist(lapply(split(freqs,cut(x,breaks,labels=FALSE)),sum));
	vals=double(breaks); vals[as.numeric(names(tmp))]=tmp/sum(tmp);
	return(c(min(x,na.rm=TRUE)+diff(range(x,na.rm=TRUE))*(1:breaks)/breaks,vals)); 
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
#generate lognormal RVs
rlnormT=function(n,mu,sig) rlnorm(n, meanlog=log(mu/sqrt(1+(sig/mu)^2)), sdlog=sqrt(log(1+(sig/mu)^2)) )
#generate a food web parameterization - a specific set of species traits and interactions
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
#Not all food web parameterizations are viable: at high connectance or diversity, some consumers can go extinct.
#This almost always happens due to competition that peaks at low mortality levels
#This function tests whether a parameter set allows consumers to coexist at a low mortality level m=0.025
parmtest=function(PARMS,NSP=nrow(PARMS$Parms)){
	PARMS$Parms[,"m"]=0.025; His=(tail(Finerun(MultiModODE2,rep(c(0.5,0.15),each=NSP),4e2,PARMS),NSP)>0.1)[-PARMS$Cabs];
	return(sum(His))
}
#keep generating food web parameterizations until find one where all consumers persist.
#If we don't find such a parameter set, use food web parameterization that allowed the most consumers to persist.
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
#Generate trophic matrix with a given connectance C. Y = elictivity for non-primary resource; not used in analysis.
VarWeb3=function(nsp,Y=NULL,C=pmin(Y+1/nsp,1),precis=0.025){ #Need either C or Y; Y superseeds.
	if(C==1/nsp) return(diag(nsp)[sample(1:nsp),]);
	cmat=matrix(0,nsp,nsp); n=0; while(n<1e5 & (min(c(colSums(cmat),rowSums(cmat)))==0 | abs(mean(cmat)-C)>precis)){ n=n+1; cmat=matrix(rbinom(nsp^2,1,C),nsp,nsp); };
	return(apply(cmat,2,function(x) x/sum(x)));
}



#Implement analysis of food web responses to consumer mortality. 
#This function iterates across levels of connectance or diversity and across replicate food webs.
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
FvarsDiversity=seq(0.025,0.55,len=ntests[1]); gamDiversity=0.25;  #Diversity levels to examine
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





