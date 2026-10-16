# 1.8 nested model tuning
options(rf.cores=1,stringsAsFactors=FALSE)
Sys.setenv(OMP_NUM_THREADS='1')
for(p in c('randomForestSRC','glmnet','pROC'))stopifnot(requireNamespace(p,quietly=TRUE))
root<-project_path('outputs','nested_tuning')
dir.create(root,recursive=TRUE,showWarnings=FALSE)
for(s in c('input','tables','models','figures'))dir.create(file.path(root,s),showWarnings=FALSE)
src<-data_path('tcga495','TCGA_495_primary_tumors.csv')
dest<-file.path(root,'input','TCGA_495_primary_tumors.csv')
if(!file.exists(dest))file.copy(src,dest)
d<-read.csv(dest);genes<-c('DNAH8','TOP2A','CHRDL1')
stopifnot(nrow(d)==495,!anyDuplicated(d$patient_id),all(d$sample_type=='Primary'),sum(d$BCR)==99)
wr<-function(x,name)write.csv(x,file.path(root,'tables',name),row.names=FALSE,na='')
auc<-function(y,p)as.numeric(pROC::auc(pROC::roc(y,p,levels=c(0,1),direction='<',quiet=TRUE)))
folds<-function(y,k,seed){set.seed(seed);f<-integer(length(y));for(v in 0:1){ix<-which(y==v);f[ix]<-sample(rep(1:k,length.out=length(ix)))};f}
rf<-expand.grid(mtry=c(1,2,3),nodesize=c(5,15,30,60),sampling=c('ordinary','balanced'),stringsAsFactors=FALSE)
rf$family<-'RF';rf$alpha<-NA_real_;rf$lambda<-NA_real_
lg<-expand.grid(alpha=c(0,.5,1),lambda=c(.001,.01,.1,1));lg$family<-'GLMNET';lg$mtry<-NA_integer_;lg$nodesize<-NA_integer_;lg$sampling<-'none'
cols<-c('family','mtry','nodesize','sampling','alpha','lambda')
candidates<-rbind(rf[,cols],lg[,cols]);candidates$id<-sprintf('C%02d',seq_len(nrow(candidates)))
base_id<-candidates$id[candidates$family=='RF' & !is.na(candidates$mtry) & candidates$mtry==2 & candidates$nodesize==5 & candidates$sampling=='ordinary']
base_id<-base_id[!is.na(base_id)];stopifnot(length(base_id)==1)
wr(candidates,'candidate_grid.csv')
fit_model<-function(tr,te,cfg,seed,retain=FALSE){
 mu<-colMeans(tr[,genes]);ss<-sapply(tr[,genes],sd)
 trans<-function(x)sweep(sweep(as.matrix(x[,genes]),2,mu,'-'),2,ss,'/')
 xx<-trans(tr);tx<-trans(te)
 if(cfg$family=='RF'){
  arg<-list(formula=BCR~.,data=data.frame(BCR=factor(tr$BCR,levels=c(0,1)),xx),ntree=1000,
   mtry=cfg$mtry,nodesize=cfg$nodesize,nsplit=10,splitrule='gini',importance=FALSE,forest=TRUE,
   block.size=10,seed=-as.integer(seed))
  set.seed(seed)
  if(cfg$sampling=='balanced'){
   ip<-which(tr$BCR==1);ineg<-which(tr$BCR==0);n_take<-max(1,floor(.632*min(length(ip),length(ineg))))
   samp<-matrix(0L,nrow=nrow(tr),ncol=1000)
   for(b in 1:1000)samp[c(sample(ip,n_take,replace=FALSE),sample(ineg,n_take,replace=FALSE)),b]<-1L
   stopifnot(all(colSums(samp)==2*n_take))
   arg$bootstrap<-'by.user';arg$samp<-samp
  }else{arg$bootstrap<-'by.root';arg$samptype<-'swor'}
  model<-do.call(randomForestSRC::rfsrc,arg)
  pp<-as.numeric(predict(model,newdata=as.data.frame(tx))$predicted[,'1'])
 }else{
  model<-glmnet::glmnet(xx,tr$BCR,family='binomial',alpha=cfg$alpha,lambda=cfg$lambda,
   standardize=FALSE,intercept=TRUE,thresh=1e-8)
  pp<-as.numeric(predict(model,newx=tx,type='response',s=cfg$lambda))
 }
 stopifnot(length(pp)==nrow(te),all(is.finite(pp)),all(pp>=0&pp<=1))
 ans<-list(prediction=pp)
 if(retain)ans<-c(ans,list(model=model,center=mu,scale=ss,config=cfg,genes=genes))
 ans
}
cls_metrics<-function(y,p,t){
 cl<-p>=t;tp<-sum(cl&y==1);tn<-sum(!cl&y==0);fp<-sum(cl&y==0);fn<-sum(!cl&y==1)
 c(TP=tp,TN=tn,FP=fp,FN=fn,sensitivity=tp/(tp+fn),specificity=tn/(tn+fp),
 accuracy=(tp+tn)/length(y),balanced_accuracy=(tp/(tp+fn)+tn/(tn+fp))/2)
}
select_threshold<-function(y,p){
 th<-sort(unique(c(0,p,.5,1)))
 z<-as.data.frame(t(vapply(th,function(t)cls_metrics(y,p,t),numeric(8))))
 z$threshold<-th;z$minimum_rate<-pmin(z$sensitivity,z$specificity)
 i<-order(-z$minimum_rate,-z$balanced_accuracy,-z$threshold)[1]
 z[i,]
}
logs<-list();all_inner_predictions<-list();inner_assignments<-list();outer_selection<-list();outer_predictions<-list();counter<-0L
inner_search<-function(dat,scope,seed){
 f<-folds(dat$BCR,4,seed)
 inner_assignments[[scope]]<<-data.frame(scope=scope,patient_id=dat$patient_id,inner_fold=f)
 score<-list();preds<-list()
 for(j in seq_len(nrow(candidates))){
  cfg<-candidates[j,];pr<-rep(NA_real_,nrow(dat));fold_auc<-fold_brier<-numeric(4)
  for(k in 1:4){
   ix<-f==k
   r<-fit_model(dat[!ix,],dat[ix,],cfg,seed=seed+100+k)
   pr[ix]<-r$prediction
   fold_auc[k]<-auc(dat$BCR[ix],pr[ix]);fold_brier[k]<-mean((pr[ix]-dat$BCR[ix])^2)
   counter<<-counter+1L
   logs[[counter]]<<-data.frame(scope=scope,candidate_id=cfg$id,family=cfg$family,inner_fold=k,
       n_train=sum(!ix),n_valid=sum(ix),AUC=fold_auc[k],Brier=fold_brier[k])
  }
  th<-select_threshold(dat$BCR,pr)
  score[[j]]<-data.frame(scope=scope,candidate_id=cfg$id,family=cfg$family,
   mean_fold_AUC=mean(fold_auc),SD_fold_AUC=sd(fold_auc),mean_Brier=mean(fold_brier),pooled_AUC=auc(dat$BCR,pr),
   selected_threshold=th$threshold,inner_sensitivity=th$sensitivity,inner_specificity=th$specificity)
  preds[[cfg$id]]<-data.frame(scope=scope,candidate_id=cfg$id,patient_id=dat$patient_id,
        inner_fold=f,BCR=dat$BCR,prediction=pr)
  if(j%%6==0){cat(scope,': completed',j,'/',nrow(candidates),'candidates\n');flush.console()}
 }
 scores<-do.call(rbind,score)
 all_inner_predictions[[scope]]<<-do.call(rbind,preds)
 list(scores=scores,preds=preds)
}
pick<-function(scores,family=NULL){
 if(!is.null(family))scores<-scores[scores$family==family,]
 scores[order(-scores$mean_fold_AUC,scores$mean_Brier,scores$candidate_id),][1,]
}
d$outer_fold<-folds(d$BCR,5,20261016)
wr(d[,c('sample_id','patient_id','BCR','outer_fold')],'outer_fold_assignments.csv')
scores_all<-list()
for(k in 1:5){
 cat('\nOUTER FOLD',k,'/ 5\n');flush.console()
 train<-d[d$outer_fold!=k,];test<-d[d$outer_fold==k,]
 search<-inner_search(train,paste0('Outer',k),20262000+k)
 scores_all[[k]]<-search$scores
 choices<-list(Original_RF=search$scores[search$scores$candidate_id==base_id,],
  Tuned_RF=pick(search$scores,'RF'),Regularized_logistic=pick(search$scores,'GLMNET'),
  Selected_pipeline=pick(search$scores))
 fitted<-list()
 for(nm in names(choices)){
  s<-choices[[nm]];cfg<-candidates[candidates$id==s$candidate_id,]
  if(is.null(fitted[[cfg$id]]))fitted[[cfg$id]]<-fit_model(train,test,cfg,20263000+k,retain=TRUE)
  r<-fitted[[cfg$id]]
  outer_selection[[paste(k,nm)]]<-data.frame(outer_fold=k,strategy=nm,s)
  outer_predictions[[paste(k,nm)]]<-data.frame(strategy=nm,outer_fold=k,candidate_id=cfg$id,
   sample_id=test$sample_id,patient_id=test$patient_id,BCR=test$BCR,prediction=r$prediction,
   threshold=s$selected_threshold,class_default=as.integer(r$prediction>=.5),
   class_inner_selected=as.integer(r$prediction>=s$selected_threshold))
 }
 saveRDS(fitted,file.path(root,'models',paste0('outer_fold_',k,'.rds')))
 wr(do.call(rbind,scores_all),'inner_candidate_scores.csv')
 wr(do.call(rbind,outer_predictions),'outer_predictions.csv')
 wr(do.call(rbind,outer_selection),'outer_selected_configurations.csv')
}
cat('\nFULL-DATA TRAINING-SIDE SELECTION (external-fold estimates already fixed)\n');flush.console()
full<-inner_search(d,'Final_full_data',20264000)
scores_all[[6]]<-full$scores
full_choices<-list(Final_tuned_RF=pick(full$scores,'RF'),Final_selected_pipeline=pick(full$scores))
final_predictions<-list();final_cfg<-list();cache<-list()
for(nm in names(full_choices)){
 s<-full_choices[[nm]];cfg<-candidates[candidates$id==s$candidate_id,]
 if(is.null(cache[[cfg$id]]))cache[[cfg$id]]<-fit_model(d,d,cfg,777,retain=TRUE)
 r<-cache[[cfg$id]];r$selected_threshold<-s$selected_threshold
 saveRDS(r,file.path(root,'models',paste0(nm,'.rds')))
 final_cfg[[nm]]<-data.frame(strategy=nm,s,cfg[,setdiff(names(cfg),c('id','family'))])
 final_predictions[[paste0(nm,'_fit')]]<-data.frame(strategy=nm,evaluation='Apparent_training',
   sample_id=d$sample_id,patient_id=d$patient_id,BCR=d$BCR,prediction=r$prediction,threshold=s$selected_threshold)
 if(cfg$family=='RF')final_predictions[[paste0(nm,'_oob')]]<-data.frame(strategy=nm,evaluation='Post_tuning_RF_OOB_diagnostic',
   sample_id=d$sample_id,patient_id=d$patient_id,BCR=d$BCR,prediction=as.numeric(r$model$predicted.oob[,'1']),threshold=s$selected_threshold)
}
wr(do.call(rbind,final_cfg),'final_selected_configurations.csv')
wr(do.call(rbind,final_predictions),'final_model_predictions.csv')
wr(do.call(rbind,scores_all),'inner_candidate_scores.csv')
wr(do.call(rbind,logs),'all_inner_fold_metrics.csv')
wr(do.call(rbind,inner_assignments),'inner_fold_assignments.csv')
wr(do.call(rbind,all_inner_predictions),'all_inner_predictions.csv')
capture.output(sessionInfo(),file=file.path(root,'sessionInfo.txt'))
cat('\nCOMPLETE\n');print(do.call(rbind,final_cfg))
