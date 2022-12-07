
STAARpipeline_pvalue_fig<-function(result,cond_result=NULL){

result2<-unlist(result)

tests<-colnames(result)[5:length(result)]
pvals<-as.numeric(result2[5:length(result2)])



trans_pvals<- -log10(pvals)
range(trans_pvals)
plot(x=c(1:length(trans_pvals)),y=trans_pvals,ylab="-log10(p-value)",xlab="",xlim=c(1,length(trans_pvals)),ylim=c(0,max(trans_pvals+2)),xaxt="n",pch=21,col="black",bg="grey50",cex=2)
axis(side=1,at=c(1:length(trans_pvals)),labels=tests,las=2)

if(!is.null(cond_result)){


result2<-unlist(cond_result)
cond_pvals<-as.numeric(result2[5:length(result2)])
trans_cond_pvals<- -log10(cond_pvals)

points(x=c(1:length(trans_cond_pvals)),y=trans_cond_pvals,pch=21,col="black",bg="indianred",cex=2)

legend("topleft",legend=c("Original","Conditional"),pch=c(21,21),pt.bg=c("grey50","indianred"),col=c("black","black"),bty="n",cex=2)

}else{
legend("topleft",legend=c("Original"),pch=c(21),pt.bg=c("grey50"),col=c("black"),bty="n",cex=1.5)
}

}




STAARpipeline_lovo_pvalue_fig<-function(result){

rmvar<-unlist(result[,"rm_variant"])

staaro<-unlist(result[,ncol(result)])
pvals<-as.numeric(staaro)



trans_pvals<- -log10(pvals)
range(trans_pvals)
plot(x=c(1:length(trans_pvals)),y=trans_pvals,ylab="-log10(p-value)",xlab="",xlim=c(1,length(trans_pvals)),ylim=c(0,max(trans_pvals+2)),xaxt="n",pch=21,col="black",bg="grey50",cex=1.2)
points(x=1,y=trans_pvals[1],pch=23,col="black",bg="indianred",cex=1.5)
axis(side=1,at=c(1:length(trans_pvals)),labels=rmvar,las=2,cex=0.8)
legend("topleft",legend=c("Original","Leave one variant out"),pch=c(23,21),pt.bg=c("indianred","grey50"),col=c("black","black"),bty="n",cex=1.2)

}
