#### PACKAGES ####
library("adespatial") # adespatial_0.3-21
library("psych")      # psych_2.3.3
library("corrplot")   # corrplot_0.92
library("vegan")      # vegan_2.6-4 
library("mvabund")    # mvabund_4.2.1
library("ggplot2")
library("egg")
options(digits = 8)

##################



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

# MEMs for spatial autocorrelation
srf_mems <- dbmem(srf_var[,c("x","y_adj")], silent = F);summary(srf_mems)
srf_var[,c("mem1","mem2")] <- srf_mems
#########################



#### CORRELATIONS ####
col2 <- colorRampPalette(c("#6b0000","#d40000","#ff7909","#fff2b5","#FFFFFF","#FFFFFF","#FFFFFF","#d5f6ff","#b4d4e7","#2166AC","#053062"))

srf_cor_mat <- corr.test(srf_var[,c("grain_mm","grain_slope","pebble_percent","sand_percent","mud_percent","temperature","toc_percent","vms_ospar","shear_stress")],
                         method="spearman",
                         adjust = "none",ci=F)

corrplot(srf_cor_mat[[1]],p.mat=srf_cor_mat[[4]],
         method="color",addCoef.col = "black",
         tl.col="black",tl.cex=1,addgrid.col="white",
         diag=F,type="lower",sig.level = 1,col = col2(200))

#####################



#### PERMANOVAS ####
srf_perms <- with(srf_var, 
                  how(nperm = 9999, 
                      blocks = station))
  
srf_permanova_otu <- adonis2(srf_otu ~ log(seq_depth) + mem1 + mem2 + vms_ospar + temperature + grain_phi + logit(toc_proportion) + sqrt(mud_percent) + shear_stress,
                             data=srf_var, permutations=srf_perms,
                             by = "margin",sqrt.dist = T,
                             method = "robust.aitchison");srf_permanova_otu

srf_permanova_gen <- adonis2(t(srf_gen[,-c(1:6)]) ~ log(seq_depth) + mem1 + mem2 + vms_ospar + temperature + grain_phi + logit(toc_proportion) + sqrt(mud_percent) + shear_stress,
                             data=srf_var, permutations=srf_perms,
                             by = "margin",sqrt.dist = T,
                             method = "robust.aitchison");srf_permanova_gen
srf_permanova_ko <- adonis2(srf_ko ~ log(seq_depth) + mem1 + mem2 + vms_ospar + temperature + grain_phi + logit(toc_proportion) + sqrt(mud_percent) + shear_stress,
                            data=srf_var, permutations=srf_perms,
                            by = "margin",sqrt.dist = T,
                            method = "robust.aitchison");srf_permanova_ko

##################



#### MGLMS ####
# (recommended to run summary.manyglm() to get deviances on an hpc )
mGLM_otu <- manyglm(mvabund(srf_otu) ~ log(seq_depth) + mem1 + mem2 + vms_ospar + temperature + grain_phi + logit(toc_proportion) + sqrt(mud_percent) + shear_stress,
                    show.coef=T,show.fitted=T,show.residuals=T,model=T,x=T,y=T,qr=T,data=srf_var,family = "negative.binomial")
s_mGLM_otu <- summary.manyglm(mGLM_otu,nBoot = 1,block=srf_var$station,show.time = "all",resamp = "case",test = "LR",show.cor = F,show.est = F,show.residuals = F)
s_mGLM_otu$coefficients

mGLM_gen <- manyglm(mvabund(t(srf_gen[,-c(1:6)])) ~ log(seq_depth) + mem1 + mem2 + vms_ospar + temperature + grain_phi + logit(toc_proportion) + sqrt(mud_percent) + shear_stress,
                    show.coef=T,show.fitted=T,show.residuals=T,model=T,x=T,y=T,qr=T,data=srf_var,family = "negative.binomial")
s_mGLM_gen <- summary.manyglm(mGLM_gen,nBoot = 1,block=srf_var$station,show.time = "all",resamp = "case",test = "LR",show.cor = T,show.est = T,show.residuals = T)
s_mGLM_gen$coefficients

mGLM_ko <- manyglm(mvabund(srf_ko) ~ log(seq_depth) + mem1 + mem2 + vms_ospar + temperature + grain_phi + logit(toc_proportion) + sqrt(mud_percent) + shear_stress,
                   show.coef=T,show.fitted=T,show.residuals=T,model=T,x=T,y=T,qr=T,data=srf_var,family = "negative.binomial")
s_mGLM_ko <- summary.manyglm(mGLM_ko,nBoot = 1,block=srf_var$station,show.time = "all",resamp = "case",test = "LR",show.cor = T,show.est = T,show.residuals = T)
s_mGLM_ko$coefficients

###################



#### PARTIAL RESIDUALS FOR NMDS PLOTS ####
# a function to correct for the effect of the sequencing depth
partial_resids_calculator <- function(df,model) {
  rmat <- df
  for(i in 1:ncol(df)){
    rmat[,i] <- df[,i] - 
      exp(model$coefficients["(Intercept)",i] +
            model$coefficients["log(seq_depth)",i]*model$x[,"log(seq_depth)"])
  }
  rmat[rmat < 0] <- 0
  invisible(rmat)
}

rmat_otu <- partial_resids_calculator(df=srf_otu, model=mGLM_otu);min(round(rmat_otu));max(round(rmat_otu))
rmat_gen <- partial_resids_calculator(df=t(srf_gen[,-c(1:6)]), model = mGLM_gen);min(round(rmat_gen));max(round(rmat_gen))
rmat_ko  <- partial_resids_calculator(df=srf_ko, model = mGLM_ko);min(round(rmat_ko));max(round(rmat_ko))

#########################################



#### nMDS ####
srf_mglm_nmds_otu1 <- metaMDS(rmat_otu, distance = "robust.aitchison",trymax = 2000,try = 1000,autotransform = F,k=2)
srf_mglm_nmds_otu2 <- metaMDS(rmat_otu, distance = "robust.aitchison",trymax = 2000,try = 1000,autotransform = F,k=2,previous.best = srf_mglm_nmds1)

srf_mglm_nmds_ko1 <- metaMDS(rmat_ko, distance = "robust.aitchison",trymax = 2000,try = 1000,autotransform = F,k=2)
srf_mglm_nmds_ko2 <- metaMDS(rmat_ko, distance = "robust.aitchison",trymax = 2000,try = 1000,autotransform = F,k=2,previous.best = srf_mglm_nmds_ko1)

srf_mglm_nmds_gen1 <- metaMDS(rmat_gen, distance = "robust.aitchison",trymax = 2000,try = 1000,autotransform = F,k=2)
srf_mglm_nmds_gen2 <- metaMDS(rmat_gen, distance = "robust.aitchison",trymax = 2000,try = 1000,autotransform = F,k=2,previous.best = srf_mglm_nmds_gen1)


# ordination vectors
envs <- with(srf_var,data.frame("grain_phi"=grain_phi,"mud"=sqrt(mud_percent),"toc"=logit(toc_proportion),"temperature"=temperature,"shear"=shear_stress,"SAR"=vms_ospar))
srf_mglm_en_otu <- envfit(srf_mglm_nmds_otu2,envs,permutations = 9999, add=T, na.rm = F)
srf_mglm_en_gen <- envfit(srf_mglm_nmds_gen2,envs,permutations = 9999, add=T, na.rm = F)
srf_mglm_en_ko  <- envfit(srf_mglm_nmds_ko2,envs,permutations = 9999, add=T, na.rm = F)

nmds_ggplot_output <- function(nmds, vars, envfit,arrow){
  df <- vars
  df$NMDS1 <- nmds$points[,1]
  df$NMDS2 <- nmds$points[,2]
  info <- as.data.frame(scores(envfit, display = "vectors",arrow.mul=arrow))
  info$p <- info$vectors$pvals
  info$vars <- rownames(info)
  stress <- nmds$stress
  invisible(list("df"=df,"info"=info,"stress"=stress))
}
ggplot_nmds_otu <- nmds_ggplot_output(srf_mglm_nmds_otu2,srf_var,srf_mglm_en_otu,arrow = 25)
ggplot_nmds_gen <- nmds_ggplot_output(srf_mglm_nmds_gen2,srf_var,srf_mglm_en_gen,arrow = 25)
ggplot_nmds_ko  <- nmds_ggplot_output(srf_mglm_nmds_ko2,srf_var,srf_mglm_en_ko,arrow = 75)

srf_mglm_vms_plot <- ggplot() + 
  geom_point(data=ggplot_nmds_otu$df,aes(x=NMDS1,y=NMDS2,color=vms_ospar),size=3,stroke=0) + 
  scale_color_distiller(palette="Spectral") + theme_bw() + 
  geom_segment(data = ggplot_nmds_otu$info,aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),arrow = arrow(length = unit(5, "mm")), colour = "black") +
  geom_text(data = ggplot_nmds_otu$info, aes(x = NMDS1, y = NMDS2, label = vars),size = 5,nudge_x = 0.0,nudge_y = -0.005) +
  annotate(col="black",geom="text", x=60, y=35, label=paste("stress =",round(ggplot_nmds_otu$stress,3))) +
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12), 
        legend.position = c(.1,.2),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        plot.background = element_blank())
srf_mglm_vms_plot

srf_mglm_vms_plot_gen <- ggplot() + 
  geom_point(data=ggplot_nmds_gen$df,aes(x=NMDS1,y=NMDS2,color=vms_ospar),size=3,stroke=0) + 
  scale_color_distiller(palette="Spectral") + theme_bw() + 
  geom_segment(data = ggplot_nmds_gen$info,aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),arrow = arrow(length = unit(5, "mm")), colour = "black") +
  geom_text(data = ggplot_nmds_gen$info, aes(x = NMDS1, y = NMDS2, label = vars),size = 5,nudge_x = 0.0,nudge_y = -0.005) +
  annotate(col="black",geom="text", x=20, y=20, label=paste("stress =",round(ggplot_nmds_gen$stress,3))) +
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12), 
        legend.position = c(.1,.2),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        plot.background = element_blank())
srf_mglm_vms_plot_gen

srf_mglm_vms_plot_ko <- ggplot() + 
  geom_point(data=ggplot_nmds_ko$df,aes(x=NMDS1,y=NMDS2,color=vms_ospar),size=3,stroke=0) + 
  scale_color_distiller(palette="Spectral") + theme_bw() + 
  geom_segment(data = ggplot_nmds_ko$info,aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),arrow = arrow(length = unit(5, "mm")), colour = "black") +
  geom_text(data = ggplot_nmds_otu$info, aes(x = NMDS1, y = NMDS2, label = vars),size = 5,nudge_x = 0.0,nudge_y = -0.005) +
  annotate(col="black",geom="text", x=70, y=50, label=paste("stress =",round(ggplot_nmds_ko$stress,3))) +
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12), 
        legend.position = c(.1,.2),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        plot.background = element_blank())
srf_mglm_vms_plot_ko

#############



#### DOMINANCE ######
gen_dm <- data.frame(row.names=colnames(srf_gen[,-c(1:6)]));
for(i in 1:nrow(gen_dm)){
  gen_dm[i,"name"] <- rownames(srf_gen[srf_gen[,i+6] == max(srf_gen[,i+6]),])
  gen_dm[i,c("phylum","class","order","family","genus")] <- srf_gen[rownames(srf_gen) == gen_dm[i,"name"] ,c("phylum","class","order","family","genus")]
}
for(i in 1:nrow(gen_dm)){
  gen_dm[i,"dominance"] <- nrow(gen_dm[gen_dm$name == gen_dm$name[i],])
}
gen_dm

##################



#### DEVIANCES ####
devs <- data.frame("names" =c(rownames(s_mGLM_otu$coefficients[-(1:4),]),
                              rownames(s_mGLM_gen$coefficients[-(1:4),]),
                              rownames(s_mGLM_ko$coefficients[-(1:4),])
),
"group"=c(rep("otu",6),
          rep("gen",6),
          rep("ko",6)
),
"deviance"=c(s_mGLM_otu$coefficients[-(1:4),1]/sum(s_mGLM_otu$coefficients[-(1:4),1])*100,
             s_mGLM_gen$coefficients[-(1:4),1]/sum(s_mGLM_gen$coefficients[-(1:4),1])*100,
             s_mGLM_ko$coefficients[-(1:4),1]/sum(s_mGLM_ko$coefficients[-(1:4),1])*100))

ggplot(devs, aes(fill=group, y=deviance, x=reorder(names,deviance))) + 
  geom_bar(position="dodge", stat="identity") +
  coord_flip() +
  geom_hline(yintercept = 0,colour = "black",size=unit(1,"mm")) +
  theme(legend.title = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent",colour = "black",size=1),
        axis.ticks.length = unit(1, "mm"),
        axis.ticks.x = element_line(colour="#333333",size=1),
        axis.ticks.y = element_line(colour="#333333",size=1)) +
  ggtitle("number of kos varying with predictors")

##################



#### DIFFERENTIAL ABUNDANCES ####
gimper <- function(model,vars,ord="coefs"){
  
  coefs <- data.frame(t(model$coefficients)[,vars])
  colnames(coefs) <- vars
  abund <- colSums(model$y)
  stdrs <- data.frame(t(model$stderr.coefficients)[,vars])
  lowers<- coefs-2.576*stdrs
  uppers<- coefs+2.576*stdrs
  
  sumr <- list()
  for(i in c(vars)) {
    sumr[[i]] <- data.frame(row.names = rownames(coefs),
                            "coefs"  = coefs[,i],
                            "abund"  = abund,
                            "lowers" = lowers[,i],
                            "uppers" = uppers[,i])
  }
  
  sumr_negative <- list()
  sumr_positive <- list()
  sumr_total <- list()
  
  for(i in c(vars)) {
    sumr_negative[[i]] <- sumr[[i]][sumr[[i]]$coefs < 0 & sumr[[i]]$lowers < 0 & sumr[[i]]$uppers < 0,]
    sumr_positive[[i]] <- sumr[[i]][sumr[[i]]$coefs > 0 & sumr[[i]]$lowers > 0 & sumr[[i]]$uppers > 0,]
    if(ord=="coefs"){
      sumr_total[[i]] <- rbind(sumr_positive[[i]][order(sumr_positive[[i]]$coefs,decreasing = T),],
                               sumr_negative[[i]][order(sumr_negative[[i]]$coefs,decreasing = T),])
    } else {
      sumr_total[[i]] <- rbind(sumr_positive[[i]][order(sumr_positive[[i]]$abund,decreasing = T),],
                               sumr_negative[[i]][order(sumr_negative[[i]]$abund,decreasing = T),])}
  }
  
  overview <- data.frame(row.names = vars)
  for(i in c(vars)) {
    overview[i,"neg_OTUs"] <- nrow(sumr_negative[[i]])
    overview[i,"pos_OTUs"] <- nrow(sumr_positive[[i]])
    overview[i,"total"] <- nrow(sumr_total[[i]])
  }
  final <- list("tables" = sumr_total,
                "summary" = overview[order(overview$total,decreasing=F),])
  invisible(final)
  
} # function to extract differentially abundant taxa/KOs from mGLMs

## otus ##
{
  dif_otus <- gimper(model = mGLM_otu, vars = c("log(seq_depth)","mem1","mem2","grain_phi","sqrt(mud_percent)","temperature","logit(toc_proportion)","shear_stress","vms_ospar"))
  dif_otus$summary_perc <- dif_otus$summary/ncol(srf_otu)*100
  
  srf_mGLM_otu_difs <- ggplot(dif_otus$summary_perc, aes(x=rownames(dif_otus$summary_perc))) + 
    geom_bar(data = dif_otus$summary_perc,aes(x=reorder(rownames(dif_otus$summary_perc), total), y = total),fill ="#b00000", stat = "identity", position = "dodge") +
    geom_bar(data = dif_otus$summary_perc,aes(x=reorder(rownames(dif_otus$summary_perc), total), y = pos_OTUs),fill ="#0044aa", stat = "identity", position = "dodge") +
    coord_flip() + 
    geom_hline(yintercept = 0,colour = "black",linewidth=unit(1,"mm")) +
    theme(legend.title = element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.position = "none",
          panel.background = element_rect(fill = "transparent"),
          panel.border = element_rect(fill = "transparent",colour = "black",size=1),
          axis.ticks.length = unit(1, "mm"),
          axis.ticks.x = element_line(colour="#333333",size=1),
          axis.ticks.y = element_line(colour="#333333",size=1)) +
    ggtitle("number of OTUs varying with predictors")
  srf_mGLM_otu_difs
  
  data.frame("deviance" =s_mGLM_otu$coefficients[c("log(seq_depth)","mem1","mem2","grain_phi","sqrt(mud_percent)","temperature","logit(toc_proportion)","shear_stress","vms_ospar"),1],
             dif_otus$summary[c("log(seq_depth)","mem1","mem2","grain_phi","sqrt(mud_percent)","temperature","logit(toc_proportion)","shear_stress","vms_ospar"),
                              c("pos_OTUs","neg_OTUs")])
  }

## kos ##
{
  dif_kos <- gimper(model = mGLM_ko, vars = c("log(seq_depth)","mem1","mem2","grain_phi","sqrt(mud_percent)","temperature","logit(toc_proportion)","shear_stress","vms_ospar"))
  dif_kos$summary_perc <- dif_kos$summary/ncol(srf_ko)*100
  
  srf_mGLM_ko_difs <- ggplot(dif_kos$summary_perc, aes(x=rownames(dif_kos$summary_perc))) + 
    geom_bar(data = dif_kos$summary_perc,aes(x=reorder(rownames(dif_kos$summary_perc), total), y = total),fill ="#b00000", stat = "identity", position = "dodge") +
    geom_bar(data = dif_kos$summary_perc,aes(x=reorder(rownames(dif_kos$summary_perc), total), y = pos_OTUs),fill ="#0044aa", stat = "identity", position = "dodge") +
    coord_flip() + 
    geom_hline(yintercept = 0,colour = "black",size=unit(1,"mm")) +
    theme(legend.title = element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.position = "none",
          panel.background = element_rect(fill = "transparent"),
          panel.border = element_rect(fill = "transparent",colour = "black",size=1),
          axis.ticks.length = unit(1, "mm"),
          axis.ticks.x = element_line(colour="#333333",size=1),
          axis.ticks.y = element_line(colour="#333333",size=1)) +
    ggtitle("number of kos varying with predictors")
  srf_mGLM_ko_difs
  }

## genera ##
{
  dif_gens <- gimper(model = mGLM_gen, vars = c("log(seq_depth)","mem1","mem2","grain_phi","sqrt(mud_percent)","temperature","logit(toc_proportion)","shear_stress","vms_ospar"))
  dif_gens$summary_perc <- round(dif_gens$summary/ncol(mGLM_gen$y)*100,2)
  
  srf_mGLM_gen_difs <- ggplot(dif_gens$summary_perc, aes(x=rownames(dif_gens$summary_perc))) + 
    geom_bar(data = dif_gens$summary_perc,aes(x=reorder(rownames(dif_gens$summary_perc), total), y = total),fill ="#b00000", stat = "identity", position = "dodge") +
    geom_bar(data = dif_gens$summary_perc,aes(x=reorder(rownames(dif_gens$summary_perc), total), y = pos_OTUs),fill ="#0044aa", stat = "identity", position = "dodge") +
    coord_flip() + 
    geom_hline(yintercept = 0,colour = "black",size=unit(1,"mm")) +
    theme(legend.title = element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.position = "none",
          panel.background = element_rect(fill = "transparent"),
          panel.border = element_rect(fill = "transparent",colour = "black",size=1),
          axis.ticks.length = unit(1, "mm"),
          axis.ticks.x = element_line(colour="#333333",size=1),
          axis.ticks.y = element_line(colour="#333333",size=1)) +
    ggtitle("number of gens varying with predictors")
  srf_mGLM_gen_difs
  
  # function to generate subsets based on CIs
  generate_subset <- function(table_name, coef_col) {
    subset_negative <- coef_col < 0
    subset_positive <- coef_col > 0
    return(rbind(
      table_name[subset_negative, ][order(table_name[subset_negative, "abund"], decreasing = TRUE), ][1:10, ],
      table_name[subset_positive, ][order(table_name[subset_positive, "abund"], decreasing = TRUE), ][1:10, ]
    ))
  }
  
  # create subsets for each table
  dif_gens_grain_subset <- generate_subset(dif_gens$tables[[4]], dif_gens$tables[[4]]$coefs)
  dif_gens_mud_subset <- generate_subset(dif_gens$tables[[5]], dif_gens$tables[[5]]$coefs)
  dif_gens_tom_subset <- generate_subset(dif_gens$tables[[6]], dif_gens$tables[[6]]$coefs)
  dif_gens_temp_subset <- generate_subset(dif_gens$tables[[7]], dif_gens$tables[[7]]$coefs)
  dif_gens_shear_subset <- generate_subset(dif_gens$tables[[8]], dif_gens$tables[[8]]$coefs)
  dif_gens_vms_subset <- generate_subset(dif_gens$tables[[9]], dif_gens$tables[[9]]$coefs)
  
  generate_ggplot <- function(data, title) {
    ggplot(data, aes(x = reorder(rownames(data), -coefs), y = coefs)) +
      geom_errorbar(aes(ymin = lowers, ymax = uppers), width = 0.15, size = 1) +
      geom_point(size = 3) +
      coord_flip() +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      ggtitle(title)
  }
  
  # Create a list of data frames and titles
  data_list <- list(dif_gens_grain_subset,dif_gens_mud_subset,dif_gens_tom_subset,dif_gens_temp_subset,dif_gens_shear_subset,dif_gens_vms_subset)
  titles <- c("grain", "mud", "tom", "temp", "shear", "vms")
  
  # Generate ggplot objects using lapply
  ggplots_list <- lapply(seq_along(data_list), function(i) {
    generate_ggplot(data_list[[i]], titles[i])
  })
  
  # Arrange the ggplot objects using ggarrange
  srf_dif_abund_composite_plot <- ggarrange(ggplots_list[[1]],
                                            ggplots_list[[2]],
                                            ggplots_list[[3]],
                                            ggplots_list[[4]],
                                            ggplots_list[[5]],
                                            ggplots_list[[6]],
                                            nrow = 3,ncol=2)
  }






################################