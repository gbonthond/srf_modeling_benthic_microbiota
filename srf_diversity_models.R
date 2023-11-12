#### PACKAGES ####
library("mobr")      # mobr_2.0.2
library("mgcv")      # mgcv_1.8-42
library("mgcViz")    # mgcViz_0.1.9
library("MuMIn")     # MuMIn_1.47.5
library("vegan")     # vegan_2.6-4
options(digits = 8)
logit <- function(x) log(x/(1-x))
#################



#### DATA PREPARATION ####

# data files
srf_var <- read.csv("srf_var.csv",row.names = 1,header=T,stringsAsFactors=T);dim(srf_var)
srf_tax <- read.csv("srf_tax.csv", row.names = 1,header=T);dim(srf_tax)
srf_otu <- read.csv("srf_otu.csv",row.names = 1,header=T);dim(srf_otu)
srf_ko  <- read.csv("srf_ko.csv", row.names = 1,header=T);dim(srf_ko)
unique(row.names(srf_otu)==row.names(srf_var) & row.names(srf_ko)==row.names(srf_var))

# tables for higher taxonomic ranks
srf_gen <- aggregate(t(srf_otu),by=srf_tax[,1:6],FUN="sum");dim(srf_gen);rownames(srf_gen) <- paste(srf_gen$phylum,srf_gen$genus,1:nrow(srf_gen),sep="_")
srf_fam <- aggregate(t(srf_otu),by=srf_tax[,1:5],FUN="sum");dim(srf_fam);rownames(srf_fam) <- paste(srf_fam$phylum,srf_fam$family,1:nrow(srf_fam),sep="_")
srf_ord <- aggregate(t(srf_otu),by=srf_tax[,1:4],FUN="sum");dim(srf_ord);rownames(srf_ord) <- paste(srf_ord$phylum,srf_ord$order,1:nrow(srf_ord),sep="_")
srf_cls <- aggregate(t(srf_otu),by=srf_tax[,1:3],FUN="sum");dim(srf_cls);rownames(srf_cls) <- paste(srf_cls$phylum,srf_cls$class,1:nrow(srf_cls),sep="_")
srf_phl <- aggregate(t(srf_otu),by=srf_tax[,1:2],FUN="sum");dim(srf_phl);rownames(srf_phl) <- paste(srf_phl$phylum,srf_phl$phl,1:nrow(srf_phl),sep="_")

# transformations predictors 
srf_var[,"sand_percent"] <- with(srf_var,vf_sand_percent + f_sand_percent + m_sand_percent + c_sand_percent + vc_sand_percent)
srf_var[,"pebble_percent"] <- with(srf_var,vf_pebbles_percent + f_pebbles_percent + pebbles_percent)

# diversity indices
srf_var[,"otu_eff_s"] <- calc_biodiv(srf_otu,index = "S_PIE",groups = rownames(srf_otu))[4]
srf_var[,"gen_eff_s"] <- calc_biodiv(t(srf_gen[,-(1:6)]),index = "S_PIE",groups = rownames(srf_otu))[4]
srf_var[,"ko_eff_s"] <- calc_biodiv(srf_ko,index = "S_PIE",groups = rownames(srf_ko))[4]

#########################



#### EFFECTIVE SPECIES NUMBER ####

## otus ##
{ 
  uGamm_otu_eff_s <- uGamm(otu_eff_s ~ log(seq_depth) + s(grain_phi, k=4, bs="cr") + s(logit(toc_proportion), k=4, bs="cr") + s(shear_stress, k=4, bs="cr") + s(sqrt(mud_percent), k=4,bs="cr") +  s(temperature, k=4, bs="cr") + s(vms_ospar, k=4,bs="cr") + s(station, bs="re"), correlation=corExp(form= ~ x + y_adj),method = "ML", data = srf_var, family = gaussian(link="identity"),na.action=na.fail)
  drd_uGamm_otu_eff_s <- dredge(uGamm_otu_eff_s,trace = 2, fixed = ~ s(station,bs = "re"))
  data.frame("RI"=sw(subset(drd_uGamm_otu_eff_s, delta <= 4)))
  gamm_otu_eff_s <- gamm(get.models(drd_uGamm_otu_eff_s,1)[[1]]$call$formula, 
                     correlation=corExp(form= ~ x + y_adj),
                     method = "REML", data = srf_var)
  summary(gamm_otu_eff_s$gam)
  print(plot.gamViz(getViz(gamm_otu_eff_s$gam), allTerms = F) + 
          l_ciPoly(level=0.95,fill="grey",alpha=0.75) +
          l_ciLine(level=0.95, colour = "black", linetype = 1,size=0.1) +
          l_fitLine(colour = "black",size=.5) + l_rug(mapping = aes(x=x, y=y)) +
          l_points(shape = 19, size = 1, alpha = 0.1) + theme_bw(), pages = 1)
  }

## genera ##
{ 
  uGamm_gen_eff_s <- uGamm(gen_eff_s ~ log(seq_depth) + s(grain_phi, k=4, bs="cr") + s(logit(toc_proportion), k=4, bs="cr") + s(shear_stress, k=4, bs="cr") + s(sqrt(mud_percent), k=4,bs="cr") +  s(temperature, k=4, bs="cr") + s(vms_ospar, k=4,bs="cr") + s(station, bs="re"), correlation=corExp(form= ~ x + y_adj),method = "ML", data = srf_var, family = gaussian(link="identity"),na.action=na.fail)
  drd_uGamm_gen_eff_s <- dredge(uGamm_gen_eff_s,trace = 2, fixed = ~ s(station,bs = "re"))
  data.frame("RI"=sw(subset(drd_uGamm_gen_eff_s, delta <= 4)))
  gamm_gen_eff_s <- gamm(get.models(drd_uGamm_gen_eff_s,1)[[1]]$call$formula, 
                         correlation=corExp(form= ~ x + y_adj),
                         method = "REML", data = srf_var)
  summary(gamm_gen_eff_s$gam)
  print(plot.gamViz(getViz(gamm_gen_eff_s$gam), allTerms = F) + 
          l_ciPoly(level=0.95,fill="grey",alpha=0.75) +
          l_ciLine(level=0.95, colour = "black", linetype = 1,size=0.1) +
          l_fitLine(colour = "black",size=.5) + l_rug(mapping = aes(x=x, y=y)) +
          l_points(shape = 19, size = 1, alpha = 0.1) + theme_bw(), pages = 1)
}

## kos ##
{ 
  uGamm_ko_eff_s <- uGamm(ko_eff_s ~ log(seq_depth) + s(grain_phi, k=4, bs="cr") + s(logit(toc_proportion), k=4, bs="cr") + s(shear_stress, k=4, bs="cr") + s(sqrt(mud_percent), k=4,bs="cr") +  s(temperature, k=4, bs="cr") + s(vms_ospar, k=4,bs="cr") + s(station, bs="re"), correlation=corExp(form= ~ x + y_adj),method = "ML", data = srf_var, family = gaussian(link="identity"),na.action=na.fail)
  drd_uGamm_ko_eff_s <- dredge(uGamm_ko_eff_s,trace = 2, fixed = ~ s(station,bs = "re"))
  data.frame("RI"=sw(subset(drd_uGamm_ko_eff_s, delta <= 4)))
  gamm_ko_eff_s <- gamm(get.models(drd_uGamm_ko_eff_s,1)[[1]]$call$formula, 
                         correlation=corExp(form= ~ x + y_adj),
                         method = "REML", data = srf_var)
  summary(gamm_ko_eff_s$gam)
  print(plot.gamViz(getViz(gamm_ko_eff_s$gam), allTerms = F) + 
          l_ciPoly(level=0.95,fill="grey",alpha=0.75) +
          l_ciLine(level=0.95, colour = "black", linetype = 1,size=0.1) +
          l_fitLine(colour = "black",size=.5) + l_rug(mapping = aes(x=x, y=y)) +
          l_points(shape = 19, size = 1, alpha = 0.1) + theme_bw(), pages = 1)
}

#################################



#### BETA DIVERSITY ####

#function to produce new dataframes with distances
beta_diversity_calculator <- function(X,var) {
  
  dis_mat_ait <- as.matrix(vegdist(X,method = "robust.aitchison"))
  dis_mat_ait[lower.tri(dis_mat_ait)] <- NA
  diag(dis_mat_ait) <- NA
  dis_mat_jac <- as.matrix(vegdist(X,method = "jaccard"))
  dis_mat_jac[lower.tri(dis_mat_jac)] <- NA
  diag(dis_mat_jac) <- NA
  
  df <- data.frame("row"=rownames(dis_mat_ait)[row(dis_mat_ait)], 
                   "col"=colnames(dis_mat_ait)[col(dis_mat_ait)])
  
  for(i in 1:nrow(df)) df[i,"aitchison"]<- dis_mat_ait[rownames(dis_mat_ait) == df[i,"row"],colnames(dis_mat_ait) == df[i,"col"]]
  for(i in 1:nrow(df)) df[i,"jaccard"]  <- dis_mat_jac[rownames(dis_mat_jac) == df[i,"row"],colnames(dis_mat_jac) == df[i,"col"]]
  
  v <- data.frame(var)
  v$row <- rownames(v)
  v$col <- rownames(v)
  
  merged_tables <- merge(merge(na.omit(df),subset(v,select=-col), by="row"),subset(v,select=-row), by="col")
  
  merged_tables$y_new <- with(merged_tables,
                              ifelse((core.x=="a" & core.y=="b"),y_adj.x,
                                     ifelse((core.x=="b" & core.y=="a"),y_adj.y,
                                            ifelse((core.x=="a" & core.y=="c"),y_adj.y,
                                                   ifelse((core.x=="c" & core.y=="a"),y_adj.x,
                                                          ifelse((core.x=="c" & core.y=="b"),y_adj.y,y_adj.x))))))
  
  
  # subsetting to distances within stations only & remove doubled variables
  trimmed_merged <- merged_tables[merged_tables$station.x == merged_tables$station.y,
                                  c("aitchison","jaccard","station.x","x.x","y_new","seq_depth.x","seq_depth.y","grain_phi.x","mud_percent.x","toc_proportion.x","temperature.x","shear_stress.x","vms_ospar.x")]
  
  colnames(trimmed_merged) <- c("aitchison","jaccard","station","x","y_new","seq_depth.x","seq_depth.y","grain_phi","mud_percent","toc_proportion","temperature","shear_stress","vms_ospar")
  invisible(trimmed_merged)
}

srf_otu_beta <- beta_diversity_calculator(srf_otu,srf_var)
srf_gen_beta <- beta_diversity_calculator(srf_gen,srf_var)
srf_ko_beta  <- beta_diversity_calculator(srf_ko,srf_var)

## otu beta diversity ##
{
  uGamm_otu_beta <- uGamm(aitchison ~ log(seq_depth.x) + log(seq_depth.y) + s(grain_phi, k=4, bs="cr") + s(logit(toc_proportion), k=4, bs="cr") + s(shear_stress, k=4, bs="cr") + s(sqrt(mud_percent), k=4,bs="cr") +  s(temperature, k=4, bs="cr") + s(vms_ospar, k=4,bs="cr") + s(station, bs="re"), correlation=corExp(form= ~ x + y_new),method = "ML", data = srf_otu_beta, family = gaussian(link="identity"),na.action=na.fail)
  drd_uGamm_otu_beta <- dredge(uGamm_otu_beta,trace = 2, fixed = ~ s(station,bs = "re"))
  data.frame("RI"=sw(subset(drd_uGamm_otu_beta, delta <= 4)))
  gamm_otu_beta <- gamm(get.models(drd_uGamm_otu_beta,subset = 1)[[1]]$call$formula, 
                       correlation=corExp(form= ~ x + y_new),
                       method = "REML", data = srf_otu_beta)
  summary(gamm_otu_beta$gam)
  print(plot.gamViz(getViz(gamm_otu_beta$gam), allTerms = T) + 
          l_ciPoly(level=0.95,fill="grey",alpha=0.75) +
          l_ciLine(level=0.95, colour = "black", linetype = 1,size=0.1) +
          l_fitLine(colour = "red",size=.5) + l_rug(mapping = aes(x=x, y=y)) +
          l_points(shape = 19, size = 1, alpha = 0.1) + theme_bw(), pages = 1)
  }

## genus beta diversity ##
{
  uGamm_gen_beta <- uGamm(aitchison ~ log(seq_depth.x) + log(seq_depth.y) + s(grain_phi, k=4, bs="cr") + s(logit(toc_proportion), k=4, bs="cr") + s(shear_stress, k=4, bs="cr") + s(sqrt(mud_percent), k=4,bs="cr") +  s(temperature, k=4, bs="cr") + s(vms_ospar, k=4,bs="cr") + s(station, bs="re"), correlation=corExp(form= ~ x + y_new),method = "ML", data = srf_gen_beta, family = gaussian(link="identity"),na.action=na.fail)
  drd_uGamm_gen_beta <- dredge(uGamm_gen_beta,trace = 2, fixed = ~ s(station,bs = "re"))
  data.frame("RI"=sw(subset(drd_uGamm_gen_beta, delta <= 4)))
  gamm_gen_beta <- gamm(get.models(drd_uGamm_gen_beta,subset = 1)[[1]]$call$formula, 
                        correlation=corExp(form= ~ x + y_new),
                        method = "REML", data = srf_gen_beta)
  summary(gamm_gen_beta$gam)
  print(plot.gamViz(getViz(gamm_gen_beta$gam), allTerms = T) + 
          l_ciPoly(level=0.95,fill="grey",alpha=0.75) +
          l_ciLine(level=0.95, colour = "black", linetype = 1,size=0.1) +
          l_fitLine(colour = "red",size=.5) + l_rug(mapping = aes(x=x, y=y)) +
          l_points(shape = 19, size = 1, alpha = 0.1) + theme_bw(), pages = 1)
  }

## functional beta diversity ##
{
  uGamm_ko_beta <- uGamm(aitchison ~ log(seq_depth.x) + log(seq_depth.y) + s(grain_phi, k=4, bs="cr") + s(logit(toc_proportion), k=4, bs="cr") + s(shear_stress, k=4, bs="cr") + s(sqrt(mud_percent), k=4,bs="cr") +  s(temperature, k=4, bs="cr") + s(vms_ospar, k=4,bs="cr") + s(station, bs="re"), correlation=corExp(form= ~ x + y_new),method = "ML", data = srf_ko_beta, family = gaussian(link="identity"),na.action=na.fail)
  drd_uGamm_ko_beta <- dredge(uGamm_ko_beta,trace = 2, fixed = ~ s(station,bs = "re"))
  data.frame("RI"=sw(subset(drd_uGamm_ko_beta, delta <= 4)))
  gamm_ko_beta <- gamm(get.models(drd_uGamm_ko_beta,subset = 1)[[1]]$call$formula, 
                        correlation=corExp(form= ~ x + y_new),
                        method = "REML", data = srf_ko_beta)
  summary(gamm_ko_beta$gam)
  print(plot.gamViz(getViz(gamm_ko_beta$gam), allTerms = T) + 
          l_ciPoly(level=0.95,fill="grey",alpha=0.75) +
          l_ciLine(level=0.95, colour = "black", linetype = 1,size=0.1) +
          l_fitLine(colour = "red",size=.5) + l_rug(mapping = aes(x=x, y=y)) +
          l_points(shape = 19, size = 1, alpha = 0.1) + theme_bw(), pages = 1)
}

#######################