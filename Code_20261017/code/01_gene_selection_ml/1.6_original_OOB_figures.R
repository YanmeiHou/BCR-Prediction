# ==== 1.6 Original saved-model OOB and AUC figures ====
# Optional; not full-pipeline validation.
# Set the working directory to this extracted package, then source with encoding=UTF-8.
options(digits=17,warn=1)
suppressPackageStartupMessages(library(ggplot2));suppressPackageStartupMessages(library(pROC));suppressPackageStartupMessages(library(patchwork))
src<-project_path('data','derived','oob')
out<-project_path('outputs','oob')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
rd<-function(n)read.csv(file.path(src,n),check.names=FALSE,fileEncoding='UTF-8-BOM')
pred<-rd('Publication_TCGA拟合及OOB逐样本预测.csv');ext<-rd('Publication_外部队列逐样本预测_保存模型推理.csv');err<-rd('Publication_RF保存误差曲线.csv')
oob<-pred[pred$prediction_type=='OOB_RF_stage',];train<-pred[pred$prediction_type=='training_fitted',]
stopifnot(nrow(oob)==547,!anyDuplicated(oob$sample_id),all(is.finite(oob$prediction)))
train<-train[match(oob$sample_id,train$sample_id),];stopifnot(identical(train$BCR,oob$BCR),all(oob$prediction>=0&oob$prediction<=1))
write_data<-function(x,n)write.csv(x,file.path(out,n),row.names=FALSE,fileEncoding='UTF-8',na='')
groups<-list('TCGA fitted'=train,'TCGA OOB'=oob,'DKFZ'=ext[ext$cohort=='DKFZ',],'Taylor'=ext[ext$cohort=='Taylor',],'CancerMap'=ext[ext$cohort=='CancerMap',])
rocs<-list();metrics<-list()
for(n in names(groups)){
 d<-groups[[n]];rr<-roc(d$BCR,d$prediction,levels=c(0,1),direction='<',quiet=TRUE);cc<-ci.auc(rr,method='delong',conf.level=.95)
 rocs[[n]]<-rr
 pp<-as.integer(d$prediction>=.5);tp<-sum(pp==1&d$BCR==1);tn<-sum(pp==0&d$BCR==0);fp<-sum(pp==1&d$BCR==0);fn<-sum(pp==0&d$BCR==1)
 metrics[[n]]<-data.frame(dataset=n,n=nrow(d),events=sum(d$BCR==1),AUC=as.numeric(auc(rr)),AUC_CI95_low=cc[1],AUC_CI95_high=cc[3],Brier=mean((d$prediction-d$BCR)^2),threshold=.5,TP=tp,TN=tn,FP=fp,FN=fn,sensitivity=tp/(tp+fn),specificity=tn/(tn+fp),accuracy=(tp+tn)/nrow(d),CI_method='DeLong; sample-level; conditional on saved predictions; not full-pipeline or patient-cluster-adjusted')
}
mt<-do.call(rbind,metrics);rownames(mt)<-NULL
rocdata<-do.call(rbind,lapply(names(rocs),function(n)data.frame(dataset=n,FPR=1-rocs[[n]]$specificities,TPR=rocs[[n]]$sensitivities,threshold=rocs[[n]]$thresholds)))
write_data(mt,'OOB及对照性能指标.csv');write_data(rocdata,'ROC曲线数据.csv')
options(digits=7)
o<-mt[mt$dataset=='TCGA OOB',];tr<-mt[mt$dataset=='TCGA fitted',]
pal<-c('TCGA fitted'='#314C6F','TCGA OOB'='#008B8B','DKFZ'='#8B9CAF','Taylor'='#8B9CAF','CancerMap'='#8B9CAF')
theme_set(theme_classic(base_family='Arial',base_size=11)+theme(plot.title=element_text(face='bold',size=13,margin=margin(b=6)),plot.subtitle=element_text(size=10,color='#52606D',margin=margin(b=10)),axis.title=element_text(size=10),axis.text=element_text(color='#263442',size=9),panel.grid.major.y=element_line(color='#E8EDF1',linewidth=.3),legend.title=element_blank(),legend.position='bottom',legend.text=element_text(size=9),plot.margin=margin(10,12,10,10)))
p1<-ggplot(rocdata[rocdata$dataset%in%c('TCGA fitted','TCGA OOB'),],aes(FPR,TPR,color=dataset))+geom_abline(slope=1,intercept=0,color='#B4BCC4',linetype=2,linewidth=.5)+geom_path(linewidth=.9)+scale_color_manual(values=pal,breaks=c('TCGA fitted','TCGA OOB'))+coord_equal(xlim=c(0,1),ylim=c(0,1),expand=FALSE)+labs(title='ROC: fitted vs OOB',subtitle=sprintf('Fitted AUC %.3f  |  OOB AUC %.3f',tr$AUC,o$AUC),x='False-positive rate',y='True-positive rate')
order<-c('TCGA fitted','TCGA OOB','DKFZ','Taylor','CancerMap');mt$dataset<-factor(mt$dataset,levels=rev(order))
p2<-ggplot(mt,aes(AUC,dataset,color=dataset))+geom_vline(xintercept=.5,color='#B4BCC4',linetype=2,linewidth=.5)+geom_segment(aes(x=AUC_CI95_low,xend=AUC_CI95_high,yend=dataset),linewidth=.85)+geom_point(size=2.5)+geom_text(aes(x=1.02,label=sprintf('%.3f [%.3f, %.3f]',AUC,AUC_CI95_low,AUC_CI95_high)),hjust=0,size=3,color='#263442')+scale_color_manual(values=pal)+scale_x_continuous(limits=c(.35,1.53),breaks=c(.4,.6,.8,1),expand=expansion(mult=0)) +labs(title='AUC with 95% confidence intervals',subtitle='DeLong intervals; saved model and predictions',x='Area under the ROC curve',y=NULL)+theme(legend.position='none',panel.grid.major.y=element_blank(),axis.text.y=element_text(size=9))
errlong<-do.call(rbind,lapply(c('all','0','1'),function(k)data.frame(trees=err$tree_index,error=err[[k]],type=switch(k,'all'='Overall','0'='Non-BCR','1'='BCR'))));errlong<-errlong[is.finite(errlong$error),]
stopifnot(nrow(errlong)>3)
# Verify the final saved error against classification at threshold 0.5.
saved_final<-err[err$tree_index==max(err$tree_index),]
stopifnot(abs(saved_final$all-(o$FP+o$FN)/o$n)<1e-12,abs(saved_final[['0']]-(1-o$specificity))<1e-12,abs(saved_final[['1']]-(1-o$sensitivity))<1e-12)
p3<-ggplot(errlong,aes(trees,error,color=type))+geom_line(linewidth=.75)+scale_color_manual(values=c('Overall'='#314C6F','Non-BCR'='#008B8B','BCR'='#CC664D'),breaks=c('Overall','Non-BCR','BCR'))+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.25))+labs(title='Saved OOB error by tree count',subtitle='Class-specific error reveals class imbalance',x='Number of trees',y='Misclassification rate')
write_data(errlong,'保存的OOB误差曲线数据.csv')
# Calibration: quantile bins without splitting tied probabilities; Wilson binomial intervals.
cuts<-unique(quantile(oob$prediction,probs=seq(0,1,.1),type=7));oob$bin<-cut(oob$prediction,breaks=cuts,include.lowest=TRUE,labels=FALSE)
cal<-do.call(rbind,lapply(sort(unique(oob$bin)),function(bn){d<-oob[oob$bin==bn,];n<-nrow(d);obs<-mean(d$BCR);z<-qnorm(.975);den<-1+z^2/n;ctr<-(obs+z^2/(2*n))/den;hw<-z*sqrt(obs*(1-obs)/n+z^2/(4*n^2))/den;data.frame(bin=bn,n=n,mean_predicted=mean(d$prediction),observed_BCR=obs,lower=ctr-hw,upper=ctr+hw)}))
stopifnot(sum(cal$n)==547)
p4<-ggplot(cal,aes(mean_predicted,observed_BCR))+geom_abline(slope=1,intercept=0,color='#B4BCC4',linetype=2,linewidth=.5)+geom_errorbar(aes(ymin=lower,ymax=upper),width=.012,color='#008B8B',linewidth=.55)+geom_line(color='#008B8B',linewidth=.75)+geom_point(color='#008B8B',size=2.1)+coord_equal(xlim=c(0,.7),ylim=c(0,.7),expand=FALSE)+labs(title='OOB calibration',subtitle=sprintf('Probability bins; Wilson 95%% intervals | Brier %.3f',o$Brier),x='Mean predicted probability',y='Observed BCR proportion')
write_data(cal,'OOB校准分箱数据.csv')
oob$label<-factor(oob$BCR,levels=c(0,1),labels=c('Non-BCR (n = 448)','BCR (n = 99)'))
p5<-ggplot(oob,aes(label,prediction,fill=label))+geom_violin(width=.75,trim=TRUE,alpha=.22,color=NA)+geom_boxplot(width=.17,outlier.shape=NA,alpha=.7,linewidth=.45)+geom_point(position=position_jitter(width=.12,height=0,seed=777),size=.85,alpha=.45,color='#384B60')+geom_hline(yintercept=.5,linetype=2,color='#CC664D',linewidth=.6)+annotate('text',x=1.45,y=.525,label='Decision threshold = 0.50',hjust=.5,size=3,color='#A84A34')+scale_fill_manual(values=c('#008B8B','#CC664D'))+scale_y_continuous(limits=c(0,1))+labs(title='OOB predictions by true label',subtitle='Each point represents one sample',x=NULL,y='Predicted BCR probability')+theme(legend.position='none')
cm<-data.frame(actual=factor(c('Non-BCR','Non-BCR','BCR','BCR'),levels=c('BCR','Non-BCR')),predicted=factor(c('Non-BCR','BCR','Non-BCR','BCR'),levels=c('Non-BCR','BCR')),count=c(o$TN,o$FP,o$FN,o$TP),cell=c('TN','FP','FN','TP'))
p6<-ggplot(cm,aes(predicted,actual,fill=count))+geom_tile(color='white',linewidth=2)+geom_text(aes(label=paste0(cell,'\n',count)),size=4.2,color='#132C42',lineheight=1.15)+scale_fill_gradient(low='#F3F7F8',high='#8BC9C6')+coord_equal()+labs(title='OOB classification at threshold 0.50',subtitle=sprintf('Sensitivity %.1f%%  |  Specificity %.1f%%',100*o$sensitivity,100*o$specificity),x='Predicted class',y='Observed class')+theme(legend.position='none',panel.grid=element_blank())
write_data(cm,'OOB混淆矩阵.csv');write_data(oob,'OOB逐样本作图数据.csv')
footer<-'RF-stage OOB after feature/model selection; not nested validation of the full pipeline.\nOriginal TCGA sample set includes 52 normal tissues and repeated patient IDs; confidence intervals are sample-level.'
combo<-(p1+p2+p3)/(p4+p5+p6)+plot_layout(heights=c(1,1))+plot_annotation(title='RF-stage out-of-bag evaluation',subtitle='Stepglm[both] + RF  |  TCGA-PRAD: 547 samples, 99 BCR events  |  Saved model; no retraining',caption=footer,tag_levels='A',theme=theme(plot.title=element_text(face='bold',size=21,family='Arial'),plot.subtitle=element_text(size=12,color='#52606D',family='Arial'),plot.caption=element_text(size=10,hjust=0,color='#52606D',family='Arial'),plot.tag=element_text(face='bold',size=16),plot.margin=margin(15,15,12,15)))
saveplot<-function(p,name,w,h){ggsave(file.path(out,paste0(name,'.png')),p,width=w,height=h,dpi=300,device=ragg::agg_png,bg='white');ggsave(file.path(out,paste0(name,'.pdf')),p,width=w,height=h,device=cairo_pdf,bg='white')}
saveplot(combo,'Fig_OOB_综合六联图',17,11.5)
ggsave(file.path(out,'Fig_OOB_预览.png'),combo,width=17,height=11.5,dpi=110,device=ragg::agg_png,bg='white')
ps<-list(p1,p2,p3,p4,p5,p6);nms<-c('Fig_OOB_01_ROC','Fig_OOB_02_AUC对照','Fig_OOB_03_误差曲线','Fig_OOB_04_校准','Fig_OOB_05_预测分布','Fig_OOB_06_混淆矩阵')
for(i in seq_along(ps))saveplot(ps[[i]],nms[i],if(i==2)8.5 else 6.7,5.5)
capture.output(sessionInfo(),file=file.path(out,'sessionInfo_绘图.txt'))
print(mt[,c('dataset','n','AUC','AUC_CI95_low','AUC_CI95_high','Brier','sensitivity','specificity')]);cat('nonmissing error points:',nrow(errlong),'\n')
