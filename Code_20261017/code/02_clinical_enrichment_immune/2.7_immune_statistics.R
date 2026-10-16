# ==== 2.7 免疫组间比较、Meta分析与BH校正 ====
# Execution scope: 已复算
# See CHAPTERS.md for order, dependencies and entry points.

# Recalculate statistics from supplied CIBERSORT fractions; no deconvolution rerun.
# Original and primary-tumor populations are retained as separate analyses.
# Both use sample-level tests, without claiming patient-independence correction.
# ---- 2.7.1 immune_statistics ----
immune_statistics <- function(d,samples) {
  require_packages("metafor")
  assert_unique(d,c("cohort","sample_id","cell_type"))
  idx <- match(d$model_sample_id,samples$model_sample_id)
  stopifnot(!anyNA(idx),all(d$BCR==samples$BCR[idx]),all(d$cohort==samples$cohort[idx]),
    all(d$primary_tumor_sensitivity==samples$primary_tumor_sensitivity[idx]),all(d$patient_id==samples$patient_id[idx]),
    all(is.finite(d$fraction)),all(d$fraction>=0&d$fraction<=1))
  cells <- unique(d$cell_type);stopifnot(length(cells)==22)
  stopifnot(all(table(d$model_sample_id)==22),max(abs(tapply(d$fraction,d$model_sample_id,sum)-1))<1e-8)
  populations<-list(original_1025=d,primary_tumor_only=d[d$primary_tumor_sensitivity,])
  tests<-list();effects<-list()
  for(pop in names(populations)) {
    population<-populations[[pop]]
    for(co in c("Pooled",unique(samples$cohort)))for(cell in cells) {
      x<-population[population$cell_type==cell,]
      if(co!="Pooled")x<-x[x$cohort==co,]
      a<-x$fraction[x$BCR==1];b<-x$fraction[x$BCR==0];n1<-length(a);n0<-length(b)
      warnings<-character()
      safe<-function(expr)tryCatch(withCallingHandlers(expr,warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")}),error=function(e){warnings<<-c(warnings,conditionMessage(e));NULL})
      tt<-safe(t.test(a,b,var.equal=FALSE,paired=FALSE,alternative="two.sided",conf.level=.95))
      ww<-safe(wilcox.test(a,b,paired=FALSE,alternative="two.sided",exact=FALSE,correct=TRUE))
      pv<-function(obj)if(is.null(obj)||!is.finite(obj$p.value))NA_real_ else obj$p.value
      tests[[length(tests)+1]]<-data.frame(population=pop,cohort=co,cell_type=cell,n_BCR=n1,n_nonBCR=n0,
        mean_BCR=mean(a),mean_nonBCR=mean(b),sd_BCR=sd(a),sd_nonBCR=sd(b),mean_difference=mean(a)-mean(b),
        welch_p_raw=pv(tt),wilcoxon_p_raw=pv(ww),warnings=paste(unique(warnings),collapse=";"))
      if(co!="Pooled") {
        sp<-sqrt(((n1-1)*var(a)+(n0-1)*var(b))/(n1+n0-2));valid<-n1>1&&n0>1&&is.finite(sp)&&sp>0
        delta<-if(valid)(mean(a)-mean(b))/sp else NA_real_
        variance<-if(valid)1/n1+1/n0+delta^2/(2*(n1+n0)) else NA_real_
        hg<-if(valid)metafor::escalc(measure="SMD",m1i=mean(a),m2i=mean(b),sd1i=sd(a),sd2i=sd(b),n1i=n1,n2i=n0,vtype="LS",correct=TRUE)else NULL
        effects[[length(effects)+1]]<-data.frame(population=pop,cohort=co,cell_type=cell,n_BCR=n1,n_nonBCR=n0,
          cohen_d=delta,variance=variance,hedges_g=if(is.null(hg))NA_real_ else hg$yi,
          hedges_variance=if(is.null(hg))NA_real_ else hg$vi,status=if(valid)"included"else"undefined_zero_pooled_SD")
      }
    }
  }
  tests<-do.call(rbind,tests);effects<-do.call(rbind,effects)
  tests$welch_p_BH22<-tests$wilcoxon_p_BH22<-NA_real_
  for(pop in names(populations))for(co in unique(tests$cohort)) {
    k<-which(tests$population==pop&tests$cohort==co)
    for(method in c("welch","wilcoxon"))tests[k,paste0(method,"_p_BH22")]<-p.adjust(tests[k,paste0(method,"_p_raw")],method="BH",n=22)
  }
  meta<-list()
  for(pop in names(populations))for(cell in cells)for(effect in c("Cohen_d","Hedges_g_sensitivity")) {
    e<-effects[effects$population==pop&effects$cell_type==cell,]
    yi<-if(effect=="Cohen_d")e$cohen_d else e$hedges_g
    vi<-if(effect=="Cohen_d")e$variance else e$hedges_variance
    ok<-is.finite(yi)&is.finite(vi)&vi>0
    if(sum(ok)<2)stop("Insufficient valid cohorts for ",cell)
    fit<-metafor::rma.uni(yi=yi[ok],vi=vi[ok],method="REML",test="z")
    meta[[length(meta)+1]]<-data.frame(population=pop,cell_type=cell,effect_type=effect,k_cohorts=sum(ok),
      cohorts_used=paste(e$cohort[ok],collapse=";"),cohorts_excluded=paste(e$cohort[!ok],collapse=";"),
      pooled_effect=as.numeric(fit$b),SE=fit$se,CI95_low=fit$ci.lb,CI95_high=fit$ci.ub,p_raw=fit$pval,I2_percent=fit$I2,tau2=fit$tau2)
  }
  meta<-do.call(rbind,meta);meta$p_BH22<-NA_real_
  for(pop in names(populations))for(effect in unique(meta$effect_type)) {
    k<-which(meta$population==pop&meta$effect_type==effect);meta$p_BH22[k]<-p.adjust(meta$p_raw[k],method="BH",n=22)
  }
  write_result(tests,"immune_tests.csv");write_result(effects,"immune_cohort_effects.csv");write_result(meta,"immune_meta.csv")
  invisible(list(tests=tests,meta=meta))
}
