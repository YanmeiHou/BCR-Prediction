# 1.7 TCGA495 fixed RF OOB
# TCGA primary tumors only; fixed original RF features and hyperparameters.
# No added validation records, artificial duplication, seed search or AUC targeting.
# RF-stage OOB conditional on the historical feature selection; not full-pipeline validation.
options(rf.cores=1,stringsAsFactors=FALSE)
Sys.setenv(OMP_NUM_THREADS='1')
required<-c('Biobase','randomForestSRC','pROC','ggplot2','patchwork')
for(pkg in required)stopifnot(requireNamespace(pkg,quietly=TRUE))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(pROC))
out<-project_path('outputs','tcga495')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
for(sub in c('input','tables','figures','models'))dir.create(file.path(out,sub),showWarnings=FALSE)
wr<-function(x,name)write.csv(x,file.path(out,'tables',name),row.names=FALSE,na='')
source_file<-project_path('data','external','bulk','TCGA-PRAD_eSet.RDS')
input_file<-data_path('tcga495','TCGA_495_primary_tumors.csv')
gene_map<-c(ENSG00000124721='DNAH8',ENSG00000131747='TOP2A',ENSG00000101938='CHRDL1',ENSG00000154330='PGM5')
if(!file.exists(input_file)){
 eset<-readRDS(source_file);pd<-Biobase::pData(eset);ex<-Biobase::exprs(eset)
 stopifnot(identical(colnames(ex),rownames(pd)),all(names(gene_map)%in%rownames(ex)))
 x<-as.data.frame(t(ex[names(gene_map),,drop=FALSE]));names(x)<-unname(gene_map)
 all<-data.frame(sample_id=rownames(pd),patient_id=pd$patient_id,sample_type=pd$sample_type,
 BCR=pd$bcr_status,x)
 all$exclusion_reason<-ifelse(all$sample_type!='Primary','Non-primary tissue','')
 bad<-is.na(all$BCR)|!all$BCR%in%0:1|is.na(all$patient_id)|all$patient_id==''|
      !apply(is.finite(as.matrix(all[,unname(gene_map)])),1,all)
 all$exclusion_reason[bad]<-'Missing or invalid endpoint, patient ID or expression'
 wr(all[all$exclusion_reason!='',],'excluded_samples.csv')
 d<-all[all$exclusion_reason=='',];d$exclusion_reason<-NULL
 write.csv(d,input_file,row.names=FALSE)
 wr(data.frame(path=source_file,MD5=unname(tools::md5sum(source_file))),'source_file_hash.csv')
 rm(eset,ex,pd,x,all);gc()
}else d<-read.csv(input_file)
stopifnot(nrow(d)==495,!anyDuplicated(d$sample_id),!anyDuplicated(d$patient_id),
 all(d$sample_type=='Primary'),sum(d$BCR==1)==99,sum(d$BCR==0)==396)
genes<-c('DNAH8','TOP2A','CHRDL1')
# Fixed historical features; PGM5 remains in source data but is not silently added to RF.
center<-colMeans(d[,genes]);scale_value<-sapply(d[,genes],sd)
xx<-as.data.frame(sweep(sweep(as.matrix(d[,genes]),2,center,'-'),2,scale_value,'/'))
wr(data.frame(gene=genes,center=center,scale=scale_value),'training_scaling.csv')
settings<-data.frame(parameter=c('analysis','n_samples','n_patients','BCR','Non_BCR','features',
 'ntree','mtry','nodesize','nsplit','splitrule','bootstrap','samptype','seed','threads','primary_threshold'),
 value=c('Fixed-feature RF-stage OOB on TCGA primary tumors',495,495,99,396,
 paste(genes,collapse=';'),1000,2,5,10,'gini','by.root','swor',777,1,.5))
wr(settings,'model_settings.csv')
cat('Training one fixed RF on 495 unique primary-tumor patients\n');flush.console()
set.seed(777)
fit<-randomForestSRC::rfsrc(BCR~.,data=data.frame(BCR=factor(d$BCR,levels=c(0,1)),xx),
 ntree=1000,mtry=2,nodesize=5,nsplit=10,splitrule='gini',bootstrap='by.root',samptype='swor',
 seed=-777,importance=FALSE,forest=TRUE,block.size=10)
pred<-d[,c('sample_id','patient_id','sample_type','BCR')]
pred$training_probability<-as.numeric(fit$predicted[,'1'])
pred$OOB_probability<-as.numeric(fit$predicted.oob[,'1'])
stopifnot(nrow(fit$xvar)==495,all(is.finite(pred$OOB_probability)),
 all(pred$OOB_probability>=0 & pred$OOB_probability<=1))
saveRDS(list(forest=fit,genes=genes,center=center,scale=scale_value,
 scope='RF-stage OOB conditional on fixed historical features',settings=settings),
 file.path(out,'models','TCGA_495_fixed_feature_RF.rds'))
wr(pred,'all_495_patient_predictions.csv')
calc<-function(y,p,t){
 tp<-sum(p>=t&y==1);tn<-sum(p<t&y==0);fp<-sum(p>=t&y==0);fn<-sum(p<t&y==1)
 data.frame(threshold=t,TP=tp,TN=tn,FP=fp,FN=fn,sensitivity=tp/(tp+fn),specificity=tn/(tn+fp),
 accuracy=(tp+tn)/length(y),balanced_accuracy=(tp/(tp+fn)+tn/(tn+fp))/2,
 PPV=if(tp+fp>0)tp/(tp+fp) else NA_real_,NPV=if(tn+fn>0)tn/(tn+fn) else NA_real_,
 F1=2*tp/(2*tp+fp+fn),Brier=mean((p-y)^2))
}
roc_train<-roc(pred$BCR,pred$training_probability,levels=c(0,1),direction='<',quiet=TRUE)
roc_oob<-roc(pred$BCR,pred$OOB_probability,levels=c(0,1),direction='<',quiet=TRUE)
rocs<-list(Training_fit=roc_train,OOB_RF_stage=roc_oob)
mt<-do.call(rbind,lapply(names(rocs),function(nm){
 rr<-rocs[[nm]];ci<-as.numeric(ci.auc(rr,method='delong'))
 p<-if(nm=='Training_fit')pred$training_probability else pred$OOB_probability
 data.frame(evaluation=nm,n=495,events=99,AUC=as.numeric(auc(rr)),CI95_low=ci[1],CI95_high=ci[3],
 CI_method='DeLong: conditional fixed-score interval; excludes model/feature selection uncertainty',calc(pred$BCR,p,.5))
}))
wr(mt,'performance_metrics.csv')
rocdata<-do.call(rbind,lapply(names(rocs),function(nm){rr<-rocs[[nm]]
 data.frame(evaluation=nm,FPR=1-rr$specificities,TPR=rr$sensitivities,threshold=rr$thresholds)}))
wr(rocdata,'ROC_coordinates.csv')
grid<-do.call(rbind,lapply(sort(unique(c(0,pred$OOB_probability,.5,1))),function(t)calc(pred$BCR,pred$OOB_probability,t)))
# Descriptive only: no OOB-optimized threshold is claimed as independently validated.
wr(grid,'OOB_threshold_sensitivity.csv')
er<-as.data.frame(fit$err.rate,check.names=FALSE)
stopifnot(all(c('all','0','1')%in%names(er)))
er$tree<-seq_len(nrow(er));wr(er,'OOB_error_by_tree.csv')
o<-mt[mt$evaluation=='OOB_RF_stage',]
last<-er[nrow(er),]
stopifnot(abs(last$all-(o$FP+o$FN)/495)<1e-10)
cm<-data.frame(x=c(1,2,1,2),y=c(2,2,1,1),cell=c('TN','FP','FN','TP'),
 count=c(o$TN,o$FP,o$FN,o$TP),row_percent=100*c(o$specificity,1-o$specificity,1-o$sensitivity,o$sensitivity),
 fill=c('#D7EAF4','#FCEEE5','#F2C4A7','#EDF5F9'))
wr(cm,'OOB_confusion_matrix_050.csv')
cuts<-unique(quantile(pred$OOB_probability,seq(0,1,.1)))
bins<-cut(pred$OOB_probability,cuts,include.lowest=TRUE)
cal<-do.call(rbind,lapply(levels(bins),function(b){ii<-bins==b;n<-sum(ii);v<-mean(pred$BCR[ii]);z<-qnorm(.975)
 den<-1+z^2/n;middle<-(v+z^2/(2*n))/den;half<-z*sqrt(v*(1-v)/n+z^2/(4*n^2))/den
 data.frame(bin=b,n=n,mean_prediction=mean(pred$OOB_probability[ii]),observed_BCR=v,
 lower=middle-half,upper=middle+half)}))
stopifnot(sum(cal$n)==495)
wr(cal,'OOB_calibration_bins.csv')
dist<-data.frame(pred,Group=factor(pred$BCR,levels=0:1,labels=c('Non-BCR','BCR')))
wr(do.call(rbind,lapply(0:1,function(y){p<-pred$OOB_probability[pred$BCR==y]
 data.frame(BCR=y,n=length(p),mean=mean(p),SD=sd(p),median=median(p),Q1=quantile(p,.25),Q3=quantile(p,.75),min=min(p),max=max(p))})),
 'OOB_probability_summary.csv')
theme_set(theme_classic(base_size=11,base_family='Arial')+theme(text=element_text(color='black'),
 axis.text=element_text(color='black'),plot.title=element_text(hjust=.5,face='bold'),
 plot.subtitle=element_text(hjust=.5,size=9),legend.title=element_blank()))
palette<-c(Training_fit='#0072B2',OOB_RF_stage='#D55E00')
labs_roc<-setNames(sprintf('%s: %.3f (95%% CI %.3f-%.3f)',c('Training fit','OOB'),mt$AUC,mt$CI95_low,mt$CI95_high),mt$evaluation)
pA<-ggplot(rocdata,aes(FPR,TPR,color=evaluation))+geom_abline(intercept=0,slope=1,linetype=2,color='grey70')+
 geom_path(linewidth=.85)+scale_color_manual(values=palette,breaks=names(labs_roc),labels=labs_roc)+
 coord_equal(xlim=c(0,1),ylim=c(0,1))+labs(title='Training and OOB ROC curves',
 subtitle='Same RF model; all 495 primary-tumor patients',x='False-positive rate',y='True-positive rate')+
 theme(legend.position='bottom',legend.text=element_text(size=8))+guides(color=guide_legend(ncol=1))
elong<-do.call(rbind,lapply(c('all','0','1'),function(k)data.frame(tree=er$tree,error=er[[k]],
 group=switch(k,all='Overall','0'='Non-BCR','1'='BCR'))))
elong<-elong[is.finite(elong$error),]
pB<-ggplot(elong,aes(tree,error,color=group))+geom_line(linewidth=.7)+
 scale_color_manual(values=c(Overall='black','Non-BCR'='#0072B2',BCR='#D55E00'))+
 scale_y_continuous(limits=c(0,1),labels=scales::percent)+
 labs(title='OOB error by tree count',subtitle=sprintf('Final overall error %.1f%%',100*last$all),
 x='Number of trees',y='Misclassification rate')+theme(legend.position='bottom')
pC<-ggplot(cal,aes(mean_prediction,observed_BCR))+geom_abline(slope=1,intercept=0,linetype=2,color='grey60')+
 geom_errorbar(aes(ymin=lower,ymax=upper),width=.012,color='#0072B2',alpha=.5)+
 geom_line(color='#0072B2')+geom_point(color='#0072B2',size=2.5)+coord_equal(xlim=c(0,1),ylim=c(0,1))+
 labs(title='OOB calibration',subtitle=sprintf('Brier %.3f; descriptive Wilson intervals',o$Brier),
 x='Mean predicted probability',y='Observed BCR proportion')
pD<-ggplot(dist,aes(Group,OOB_probability,color=Group,fill=Group))+geom_boxplot(width=.5,alpha=.15,outlier.shape=NA)+
 geom_point(position=position_jitter(width=.13,seed=123),alpha=.4,size=.9)+geom_hline(yintercept=.5,linetype=2,color='grey60')+
 scale_color_manual(values=c('#0072B2','#D55E00'))+scale_fill_manual(values=c('#0072B2','#D55E00'))+
 scale_x_discrete(labels=c('Non-BCR\nn = 396','BCR\nn = 99'))+coord_cartesian(ylim=c(0,1))+
 labs(title='OOB prediction distributions',subtitle='All 495 patients; dashed line = 0.50',x=NULL,y='Predicted BCR probability')+
 theme(legend.position='none')
pE<-ggplot(cm,aes(x,y))+geom_tile(aes(fill=fill),width=.96,height=.96)+scale_fill_identity()+
 geom_text(aes(label=count),size=7,fontface='bold',vjust=-.3)+
 geom_text(aes(label=sprintf('%s | %.1f%%',cell,row_percent)),size=3.5,vjust=2)+
 scale_x_continuous(breaks=1:2,labels=c('Non-BCR','BCR'))+scale_y_continuous(breaks=1:2,labels=c('BCR','Non-BCR'))+
 coord_equal()+labs(title='OOB confusion matrix at threshold 0.50',
 subtitle=sprintf('n = 495 | Sensitivity %.1f%% | Specificity %.1f%%',100*o$sensitivity,100*o$specificity),
 x='Predicted class',y='Observed class')+theme(axis.line=element_blank(),axis.ticks=element_blank())
tlong<-rbind(data.frame(threshold=grid$threshold,value=grid$sensitivity,measure='Sensitivity'),
 data.frame(threshold=grid$threshold,value=grid$specificity,measure='Specificity'))
pF<-ggplot(tlong,aes(threshold,value,color=measure))+geom_line(linewidth=.8)+geom_vline(xintercept=.5,linetype=2,color='grey60')+
 scale_color_manual(values=c(Sensitivity='#D55E00',Specificity='#0072B2'))+
 scale_y_continuous(limits=c(0,1),labels=scales::percent)+labs(title='OOB threshold sensitivity analysis',
 subtitle='Exploratory curve; primary threshold fixed at 0.50',x='Probability threshold',y='Rate')+theme(legend.position='bottom')
plots<-list(A_ROC=pA,B_OOB_error=pB,C_OOB_calibration=pC,D_OOB_distribution=pD,E_OOB_confusion=pE,F_OOB_threshold=pF)
for(nm in names(plots)){
 ggsave(file.path(out,'figures',paste0(nm,'.png')),plots[[nm]],width=6.5,height=5.5,dpi=600)
 ggsave(file.path(out,'figures',paste0(nm,'.pdf')),plots[[nm]],width=6.5,height=5.5,device=cairo_pdf)
}
combined<-wrap_plots(plots,ncol=3,widths=c(1,1,1),heights=c(1,1))+plot_annotation(tag_levels='A',
 title='TCGA: same-model training and OOB evaluation in 495 primary-tumor patients',
 caption=paste('495 unique patients: 99 BCR and 396 Non-BCR. All 52 normal-tissue samples excluded; no artificial duplication.',
 'Fixed historical RF features: DNAH8, TOP2A and CHRDL1. OOB evaluates the RF stage conditional on prior feature selection.',
 'AUC confidence intervals condition on saved scores; they do not include full model-selection uncertainty. No external-validation claim.',sep='\n'),
 theme=theme(plot.title=element_text(hjust=.5,face='bold',size=16),plot.caption=element_text(hjust=0,size=9)))
ggsave(file.path(out,'figures','TCGA_495_training_OOB_six_panels.png'),combined,width=19.5,height=12,dpi=600)
ggsave(file.path(out,'figures','TCGA_495_training_OOB_six_panels.pdf'),combined,width=19.5,height=12,device=cairo_pdf)
ggsave(file.path(out,'figures','TCGA_495_training_OOB_six_panels.tiff'),combined,width=19.5,height=12,dpi=600,compression='lzw')
ggsave(file.path(out,'figures','TCGA_495_training_OOB_preview.png'),combined,width=19.5,height=12,dpi=110)
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
print(mt)
