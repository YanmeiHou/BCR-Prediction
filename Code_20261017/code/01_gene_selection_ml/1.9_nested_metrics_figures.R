# 1.9 nested metrics figures
options(stringsAsFactors=FALSE)
library(ggplot2);library(patchwork)
root<-project_path('outputs','nested_tuning')
rd<-function(n)read.csv(file.path(root,'tables',n))
wr<-function(d,n)write.csv(d,file.path(root,'tables',n),row.names=FALSE)
d<-rd('outer_predictions.csv');fin<-rd('final_model_predictions.csv')
roc<-function(y,p)pROC::roc(y,p,levels=c(0,1),direction='<',quiet=TRUE)
auc<-function(y,p)as.numeric(pROC::auc(roc(y,p)))
metrics<-function(z,cl){
 rr<-roc(z$BCR,z$prediction);ci<-as.numeric(pROC::ci.auc(rr,method='delong'))
 tp<-sum(cl==1&z$BCR==1);tn<-sum(cl==0&z$BCR==0);fp<-sum(cl==1&z$BCR==0);fn<-sum(cl==0&z$BCR==1)
 data.frame(n=nrow(z),BCR=sum(z$BCR),AUC=as.numeric(pROC::auc(rr)),CI95_low=ci[1],CI95_high=ci[3],
 Brier=mean((z$prediction-z$BCR)^2),TP=tp,TN=tn,FP=fp,FN=fn,sensitivity=tp/(tp+fn),specificity=tn/(tn+fp),
 accuracy=(tp+tn)/nrow(z),balanced_accuracy=(tp/(tp+fn)+tn/(tn+fp))/2)
}
ms<-list();fm<-list()
for(s in unique(d$strategy)){
 z<-d[d$strategy==s,]
 for(t in c('Default_0.50','Inner_selected'))ms[[paste(s,t)]]<-data.frame(strategy=s,threshold_rule=t,metrics(z,if(t=='Default_0.50')z$class_default else z$class_inner_selected))
 for(k in 1:5){zz<-z[z$outer_fold==k,];fm[[paste(s,k)]]<-data.frame(strategy=s,outer_fold=k,metrics(zz,zz$class_inner_selected))}
}
m<-do.call(rbind,ms);f<-do.call(rbind,fm);wr(m,'outer_performance.csv');wr(f,'outer_fold_performance.csv')
summ<-do.call(rbind,lapply(split(f,f$strategy),function(z)data.frame(strategy=z$strategy[1],mean_fold_AUC=mean(z$AUC),SD_fold_AUC=sd(z$AUC),mean_fold_Brier=mean(z$Brier))))
wr(summ,'outer_fold_summary.csv')
fs<-list()
for(s in unique(fin$strategy))for(e in unique(fin$evaluation)){
 z<-fin[fin$strategy==s&fin$evaluation==e,]
 for(t in c('Default_0.50','Inner_selected'))fs[[paste(s,e,t)]]<-data.frame(strategy=s,evaluation=e,threshold_rule=t,threshold=if(t=='Default_0.50').5 else z$threshold[1],metrics(z,as.integer(z$prediction>=if(t=='Default_0.50').5 else z$threshold)))
}
ff<-do.call(rbind,fs);wr(ff,'final_diagnostic_performance.csv')
# Paired outcome-stratified bootstrap of fixed predictions, conditional on this fitted CV run.
base<-d[d$strategy=='Original_RF',];base<-base[order(base$patient_id),];set.seed(20265000)
ixpos<-which(base$BCR==1);ixneg<-which(base$BCR==0)
bootix<-replicate(2000,c(sample(ixpos,replace=TRUE),sample(ixneg,replace=TRUE)),simplify=FALSE)
dd<-list()
for(s in c('Tuned_RF','Regularized_logistic','Selected_pipeline')){
 z<-d[d$strategy==s,];z<-z[order(z$patient_id),];stopifnot(identical(z$patient_id,base$patient_id))
 del<-vapply(bootix,function(i)auc(z$BCR[i],z$prediction[i])-auc(base$BCR[i],base$prediction[i]),numeric(1))
 dd[[s]]<-data.frame(strategy=s,reference='Original_RF',delta_pooled_AUC=auc(z$BCR,z$prediction)-auc(base$BCR,base$prediction),CI95_low=quantile(del,.025),CI95_high=quantile(del,.975),bootstrap_replicates=2000)
}
wr(do.call(rbind,dd),'paired_AUC_differences_conditional.csv')
names_short<-c(Original_RF='Original RF',Tuned_RF='Tuned RF',Regularized_logistic='Regularized logistic',Selected_pipeline='Selected pipeline')
cols<-c('Original RF'='#237DAC','Tuned RF'='#DE7B28','Regularized logistic'='#48866B')
theme_set(theme_classic(base_size=13)+theme(text=element_text(colour='black'),axis.text=element_text(colour='black'),plot.title=element_text(hjust=.5,face='bold'),legend.position='bottom',legend.title=element_blank()))
rocdata<-do.call(rbind,lapply(names(cols),function(s){z<-d[d$strategy==names(names_short)[names_short==s],];r<-roc(z$BCR,z$prediction);data.frame(FPR=1-r$specificities,TPR=r$sensitivities,Model=s)}))
p1<-ggplot(rocdata,aes(FPR,TPR,colour=Model))+geom_abline(slope=1,intercept=0,linetype=2,colour='grey70')+geom_path(linewidth=.8)+scale_colour_manual(values=cols)+coord_equal()+labs(title='Outer-fold ROC curves',x='False-positive rate',y='True-positive rate')
fs2<-f[f$strategy!='Selected_pipeline',];fs2$Model<-unname(names_short[fs2$strategy])
p2<-ggplot(fs2,aes(factor(outer_fold),AUC,group=Model,colour=Model))+geom_line()+geom_point(size=2)+scale_colour_manual(values=cols)+coord_cartesian(ylim=c(.4,1))+labs(title='AUC within each outer fold',x='Outer fold',y='AUC')
z<-d[d$strategy=='Tuned_RF',];z$Group<-ifelse(z$BCR==1,'BCR','Non-BCR')
p3<-ggplot(z,aes(Group,prediction,fill=Group))+geom_boxplot(width=.5,outlier.shape=NA,alpha=.25)+geom_jitter(width=.16,height=0,size=.9,alpha=.35)+scale_fill_manual(values=c('BCR'='#DE7B28','Non-BCR'='#237DAC'))+coord_cartesian(ylim=c(0,1))+labs(title='Tuned RF: outer-fold predictions',x=NULL,y='Predicted BCR probability')+theme(legend.position='none')
cmplot<-function(cl,title){
 tab<-as.data.frame(table(Observed=factor(z$BCR,levels=c(0,1),labels=c('Non-BCR','BCR')),Predicted=factor(cl,levels=c(0,1),labels=c('Non-BCR','BCR'))))
 tab$percent<-tab$Freq/ave(tab$Freq,tab$Observed,FUN=sum)*100
 ggplot(tab,aes(Predicted,Observed,fill=Freq))+geom_tile(colour='white',linewidth=2)+geom_text(aes(label=sprintf('%d\n%.1f%%',Freq,percent)),size=5,colour='black')+scale_fill_gradient(low='#EDF4F7',high='#72AFCC')+scale_y_discrete(limits=c('BCR','Non-BCR'))+coord_equal()+labs(title=title,subtitle='495 patients; percentages within observed class')+theme(legend.position='none',plot.subtitle=element_text(hjust=.5))
}
p4<-cmplot(z$class_default,'Tuned RF: threshold 0.50')
p5<-cmplot(z$class_inner_selected,'Tuned RF: inner-selected thresholds')
rates<-m[m$strategy%in%c('Original_RF','Tuned_RF'),];rr<-rbind(data.frame(rates[,c('strategy','threshold_rule')],Metric='Sensitivity',Value=rates$sensitivity),data.frame(rates[,c('strategy','threshold_rule')],Metric='Specificity',Value=rates$specificity))
rr$Setting<-paste(names_short[rr$strategy],ifelse(rr$threshold_rule=='Default_0.50','0.50','Inner-selected'),sep='\n')
p6<-ggplot(rr,aes(Setting,Value,fill=Metric))+geom_col(position=position_dodge(.8),width=.7)+geom_text(aes(label=sprintf('%.1f%%',100*Value)),position=position_dodge(.8),vjust=-.4,size=3.3)+scale_fill_manual(values=c(Sensitivity='#DE7B28',Specificity='#237DAC'))+scale_y_continuous(limits=c(0,1.08),breaks=seq(0,1,.25),labels=function(x)paste0(100*x,'%'))+labs(title='Sensitivity-specificity trade-off',x=NULL,y='Rate')
plots<-list(A_ROC=p1,B_fold_AUC=p2,C_distribution=p3,D_confusion_050=p4,E_confusion_inner=p5,F_classification=p6)
for(n in names(plots)){
 ggsave(file.path(root,'figures',paste0(n,'.png')),plots[[n]],width=7,height=5.5,dpi=600,bg='white')
 ggsave(file.path(root,'figures',paste0(n,'.pdf')),plots[[n]],width=7,height=5.5,device=cairo_pdf)
}
comb<-wrap_plots(plots,ncol=3)+plot_annotation(title='TCGA: 495 unique tumor patients | nested cross-validation',subtitle='All panels use outer-fold predictions. Selected pipeline equals tuned RF in all five folds.',tag_levels='A')
ggsave(file.path(root,'figures','Nested_CV_six_panels.png'),comb,width=20,height=12,dpi=600,bg='white')
ggsave(file.path(root,'figures','Nested_CV_six_panels.pdf'),comb,width=20,height=12,device=cairo_pdf)
ggsave(file.path(root,'figures','Nested_CV_preview.png'),comb,width=20,height=12,dpi=100,bg='white')
# Final-fit diagnostics are deliberately separate from the outer-fold validation figure.
zf<-fin[fin$strategy=='Final_tuned_RF',];fr<-do.call(rbind,lapply(split(zf,zf$evaluation),function(q){r<-roc(q$BCR,q$prediction);data.frame(FPR=1-r$specificities,TPR=r$sensitivities,Evaluation=if(q$evaluation[1]=='Apparent_training')'Apparent training' else 'Post-tuning OOB diagnostic')}))
pf<-ggplot(fr,aes(FPR,TPR,colour=Evaluation))+geom_abline(slope=1,intercept=0,linetype=2,colour='grey70')+geom_path(linewidth=.8)+coord_equal()+scale_colour_manual(values=c('#237DAC','#DE7B28'))+labs(title='Final RF: training / OOB diagnostic',subtitle='Post-selection diagnostic; primary validation is nested CV',x='False-positive rate',y='True-positive rate')
ggsave(file.path(root,'figures','Final_RF_training_OOB_diagnostic.png'),pf,width=7,height=6,dpi=600,bg='white')
ggsave(file.path(root,'figures','Final_RF_training_OOB_diagnostic.pdf'),pf,width=7,height=6,device=cairo_pdf)
print(m);print(summ);print(do.call(rbind,dd));print(ff)
