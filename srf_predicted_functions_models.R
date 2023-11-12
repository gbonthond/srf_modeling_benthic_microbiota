#### PACKAGES #### 
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
srf_ko  <- read.csv("srf_ko.csv", row.names = 1,header=T);dim(srf_ko)
unique(row.names(srf_ko)==row.names(srf_var))

# transformations predictors 
srf_var[,"sand_percent"] <- with(srf_var,vf_sand_percent + f_sand_percent + m_sand_percent + c_sand_percent + vc_sand_percent)
srf_var[,"pebble_percent"] <- with(srf_var,vf_pebbles_percent + f_pebbles_percent + pebbles_percent)


#########################



#### KEGG MODULES ####

# M00155 respiration
srf_var[, "respiration"] <- rowSums(srf_ko[, colnames(srf_ko) %in% c("K02276","K02275","K02274","K15408","K02277")])

# M00528 nitrification
srf_var[, "nitrification"] <- rowSums(srf_ko[, colnames(srf_ko) %in% c("K10944","K10945","K10946","K10535")])

# M00530.1 dissimilatory nitrate reduction
srf_var[, "nitrate_reduction"] <- rowSums(srf_ko[, colnames(srf_ko) %in% c("K02567","K02568","K00370","K00371","K00374")])

# M00530.2 dissimilatory nitrite reduction
srf_var[, "nitrite_reduction"] <- rowSums(srf_ko[, colnames(srf_ko) %in% c("K00362","K00363","K03385","K15876")])

# M00529 denitrification
srf_var[, "denitrification"] <- rowSums(srf_ko[, colnames(srf_ko) %in% c("K02567","K02568","K00370","K00371","K00374","K00368","K15864","K04561","K02305","K00376")])

# M00596 dissimlatory sulfate reduction
srf_var[, "sulfate_reduction"] <- rowSums(srf_ko[, colnames(srf_ko) %in% c("K00958","K00394","K00395","K11180","K11181")])

# M00595 Thiosulfate oxidation
srf_var[, "sulfur_oxidation"] <- rowSums(srf_ko[, colnames(srf_ko) %in% c("K17222","K17223","K17224","K17225","K17226","K17227")])

# M00174 methane oxidation
srf_var[, "methane_oxidation"] <- rowSums(srf_ko[, colnames(srf_ko) %in% c("K10944","K10945","K10946","K16157","K16158","K16159","K16161","K16160","K16162","K14028","K14029","K23995")])

# M00567 CO2 reduction
srf_var[, "CO2_reduction"] <- rowSums(srf_ko[, colnames(srf_ko) %in% c("K00200","K00201","K00202","K00203","K11261","K00205","K11260","K00204","K00672","K01499","K00319","K13942","K00320","K00577","K00578","K00579","K00580","K00581","K00582","K00583","K00584","K00399","K00401","K00402","K08264","K08265","K22480","K22481","K22482","K03388","K03389","K03390","K14127","K14126","K14128","K22516","K00125")])

#####################



#### metabolic functions ####

## respiration ##
{ 
  uGamm_respiration <- uGamm(log(respiration) ~ log(seq_depth) + s(grain_phi, k=4, bs="cr") + s(logit(toc_proportion), k=4, bs="cr") + s(shear_stress, k=4, bs="cr") + s(sqrt(mud_percent), k=4,bs="cr") +  s(temperature, k=4, bs="cr") + s(vms_ospar, k=4,bs="cr") + s(station, bs="re"), correlation=corExp(form= ~ x + y_adj),method = "ML", data = srf_var, family = gaussian(link="identity"),na.action=na.fail)
  drd_uGamm_respiration <- dredge(uGamm_respiration,trace = 2, fixed = ~ log(seq_depth) + s(station,bs = "re"))
  data.frame("RI"=sw(subset(drd_uGamm_respiration, delta <= 4)))
  gamm_respiration <- gamm(get.models(drd_uGamm_respiration,1)[[1]]$call$formula, 
                         correlation=corExp(form= ~ x + y_adj),
                         method = "REML", data = srf_var)
  summary(gamm_respiration$gam)
  print(plot.gamViz(getViz(gamm_respiration$gam), allTerms = F) + 
          l_ciPoly(level=0.95,fill="grey",alpha=0.75) +
          l_ciLine(level=0.95, colour = "black", linetype = 1,size=0.1) +
          l_fitLine(colour = "black",size=.5) + l_rug(mapping = aes(x=x, y=y)) +
          l_points(shape = 19, size = 1, alpha = 0.1) + theme_bw(), pages = 1)
}

## nitrification ##
{ 
  uGamm_nitrification <- uGamm(log1p(nitrification) ~ log(seq_depth) + s(grain_phi, k=4, bs="cr") + s(logit(toc_proportion), k=4, bs="cr") + s(shear_stress, k=4, bs="cr") + s(sqrt(mud_percent), k=4,bs="cr") +  s(temperature, k=4, bs="cr") + s(vms_ospar, k=4,bs="cr") + s(station, bs="re"), correlation=corExp(form= ~ x + y_adj),method = "ML", data = srf_var, family = gaussian(link="identity"),na.action=na.fail)
  drd_uGamm_nitrification <- dredge(uGamm_nitrification,trace = 2, fixed = ~ log(seq_depth) + s(station,bs = "re"))
  data.frame("RI"=sw(subset(drd_uGamm_nitrification, delta <= 4)))
  gamm_nitrification <- gamm(get.models(drd_uGamm_nitrification,1)[[1]]$call$formula, 
                           correlation=corExp(form= ~ x + y_adj),
                           method = "REML", data = srf_var)
  summary(gamm_nitrification$gam)
  print(plot.gamViz(getViz(gamm_nitrification$gam), allTerms = F) + 
          l_ciPoly(level=0.95,fill="grey",alpha=0.75) +
          l_ciLine(level=0.95, colour = "black", linetype = 1,size=0.1) +
          l_fitLine(colour = "black",size=.5) + l_rug(mapping = aes(x=x, y=y)) +
          l_points(shape = 19, size = 1, alpha = 0.1) + theme_bw(), pages = 1)
}

## nitrate_reduction ##
{ 
  uGamm_nitrate_reduction <- uGamm(log(nitrate_reduction) ~ log(seq_depth) + s(grain_phi, k=4, bs="cr") + s(logit(toc_proportion), k=4, bs="cr") + s(shear_stress, k=4, bs="cr") + s(sqrt(mud_percent), k=4,bs="cr") +  s(temperature, k=4, bs="cr") + s(vms_ospar, k=4,bs="cr") + s(station, bs="re"), correlation=corExp(form= ~ x + y_adj),method = "ML", data = srf_var, family = gaussian(link="identity"),na.action=na.fail)
  drd_uGamm_nitrate_reduction <- dredge(uGamm_nitrate_reduction,trace = 2, fixed = ~ log(seq_depth) + s(station,bs = "re"))
  data.frame("RI"=sw(subset(drd_uGamm_nitrate_reduction, delta <= 4)))
  gamm_nitrate_reduction <- gamm(get.models(drd_uGamm_nitrate_reduction,1)[[1]]$call$formula, 
                           correlation=corExp(form= ~ x + y_adj),
                           method = "REML", data = srf_var)
  summary(gamm_nitrate_reduction$gam)
  print(plot.gamViz(getViz(gamm_nitrate_reduction$gam), allTerms = F) + 
          l_ciPoly(level=0.95,fill="grey",alpha=0.75) +
          l_ciLine(level=0.95, colour = "black", linetype = 1,size=0.1) +
          l_fitLine(colour = "black",size=.5) + l_rug(mapping = aes(x=x, y=y)) +
          l_points(shape = 19, size = 1, alpha = 0.1) + theme_bw(), pages = 1)
}

## nitrite_reduction ##
{ 
  uGamm_nitrite_reduction <- uGamm(log(nitrite_reduction) ~ log(seq_depth) + s(grain_phi, k=4, bs="cr") + s(logit(toc_proportion), k=4, bs="cr") + s(shear_stress, k=4, bs="cr") + s(sqrt(mud_percent), k=4,bs="cr") +  s(temperature, k=4, bs="cr") + s(vms_ospar, k=4,bs="cr") + s(station, bs="re"), correlation=corExp(form= ~ x + y_adj),method = "ML", data = srf_var, family = gaussian(link="identity"),na.action=na.fail)
  drd_uGamm_nitrite_reduction <- dredge(uGamm_nitrite_reduction,trace = 2, fixed = ~ log(seq_depth) + s(station,bs = "re"))
  data.frame("RI"=sw(subset(drd_uGamm_nitrite_reduction, delta <= 4)))
  gamm_nitrite_reduction <- gamm(get.models(drd_uGamm_nitrite_reduction,1)[[1]]$call$formula, 
                           correlation=corExp(form= ~ x + y_adj),
                           method = "REML", data = srf_var)
  summary(gamm_nitrite_reduction$gam)
  print(plot.gamViz(getViz(gamm_nitrite_reduction$gam), allTerms = F) + 
          l_ciPoly(level=0.95,fill="grey",alpha=0.75) +
          l_ciLine(level=0.95, colour = "black", linetype = 1,size=0.1) +
          l_fitLine(colour = "black",size=.5) + l_rug(mapping = aes(x=x, y=y)) +
          l_points(shape = 19, size = 1, alpha = 0.1) + theme_bw(), pages = 1)
}

## denitrification ##
{ 
  uGamm_denitrification <- uGamm(log(denitrification) ~ log(seq_depth) + s(grain_phi, k=4, bs="cr") + s(logit(toc_proportion), k=4, bs="cr") + s(shear_stress, k=4, bs="cr") + s(sqrt(mud_percent), k=4,bs="cr") +  s(temperature, k=4, bs="cr") + s(vms_ospar, k=4,bs="cr") + s(station, bs="re"), correlation=corExp(form= ~ x + y_adj),method = "ML", data = srf_var, family = gaussian(link="identity"),na.action=na.fail)
  drd_uGamm_denitrification <- dredge(uGamm_denitrification,trace = 2, fixed = ~ log(seq_depth) + s(station,bs = "re"))
  data.frame("RI"=sw(subset(drd_uGamm_denitrification, delta <= 4)))
  gamm_denitrification <- gamm(get.models(drd_uGamm_denitrification,1)[[1]]$call$formula, 
                           correlation=corExp(form= ~ x + y_adj),
                           method = "REML", data = srf_var)
  summary(gamm_denitrification$gam)
  print(plot.gamViz(getViz(gamm_denitrification$gam), allTerms = F) + 
          l_ciPoly(level=0.95,fill="grey",alpha=0.75) +
          l_ciLine(level=0.95, colour = "black", linetype = 1,size=0.1) +
          l_fitLine(colour = "black",size=.5) + l_rug(mapping = aes(x=x, y=y)) +
          l_points(shape = 19, size = 1, alpha = 0.1) + theme_bw(), pages = 1)
}

## sulfate_reduction ##
{ 
  uGamm_sulfate_reduction <- uGamm(log(sulfate_reduction) ~ log(seq_depth) + s(grain_phi, k=4, bs="cr") + s(logit(toc_proportion), k=4, bs="cr") + s(shear_stress, k=4, bs="cr") + s(sqrt(mud_percent), k=4,bs="cr") +  s(temperature, k=4, bs="cr") + s(vms_ospar, k=4,bs="cr") + s(station, bs="re"), correlation=corExp(form= ~ x + y_adj),method = "ML", data = srf_var, family = gaussian(link="identity"),na.action=na.fail)
  drd_uGamm_sulfate_reduction <- dredge(uGamm_sulfate_reduction,trace = 2, fixed = ~ log(seq_depth) + s(station,bs = "re"))
  data.frame("RI"=sw(subset(drd_uGamm_sulfate_reduction, delta <= 4)))
  gamm_sulfate_reduction <- gamm(get.models(drd_uGamm_sulfate_reduction,1)[[1]]$call$formula, 
                           correlation=corExp(form= ~ x + y_adj),
                           method = "REML", data = srf_var)
  summary(gamm_sulfate_reduction$gam)
  print(plot.gamViz(getViz(gamm_sulfate_reduction$gam), allTerms = F) + 
          l_ciPoly(level=0.95,fill="grey",alpha=0.75) +
          l_ciLine(level=0.95, colour = "black", linetype = 1,size=0.1) +
          l_fitLine(colour = "black",size=.5) + l_rug(mapping = aes(x=x, y=y)) +
          l_points(shape = 19, size = 1, alpha = 0.1) + theme_bw(), pages = 1)
}

## sulfur_oxidation ##
{ 
  uGamm_sulfur_oxidation <- uGamm(log(sulfur_oxidation) ~ log(seq_depth) + s(grain_phi, k=4, bs="cr") + s(logit(toc_proportion), k=4, bs="cr") + s(shear_stress, k=4, bs="cr") + s(sqrt(mud_percent), k=4,bs="cr") +  s(temperature, k=4, bs="cr") + s(vms_ospar, k=4,bs="cr") + s(station, bs="re"), correlation=corExp(form= ~ x + y_adj),method = "ML", data = srf_var, family = gaussian(link="identity"),na.action=na.fail)
  drd_uGamm_sulfur_oxidation <- dredge(uGamm_sulfur_oxidation,trace = 2, fixed = ~ log(seq_depth) + s(station,bs = "re"))
  data.frame("RI"=sw(subset(drd_uGamm_sulfur_oxidation, delta <= 4)))
  gamm_sulfur_oxidation <- gamm(get.models(drd_uGamm_sulfur_oxidation,1)[[1]]$call$formula, 
                           correlation=corExp(form= ~ x + y_adj),
                           method = "REML", data = srf_var)
  summary(gamm_sulfur_oxidation$gam)
  print(plot.gamViz(getViz(gamm_sulfur_oxidation$gam), allTerms = F) + 
          l_ciPoly(level=0.95,fill="grey",alpha=0.75) +
          l_ciLine(level=0.95, colour = "black", linetype = 1,size=0.1) +
          l_fitLine(colour = "black",size=.5) + l_rug(mapping = aes(x=x, y=y)) +
          l_points(shape = 19, size = 1, alpha = 0.1) + theme_bw(), pages = 1)
}

## methane_oxidation ##
{ 
  uGamm_methane_oxidation <- uGamm(log(methane_oxidation) ~ log(seq_depth) + s(grain_phi, k=4, bs="cr") + s(logit(toc_proportion), k=4, bs="cr") + s(shear_stress, k=4, bs="cr") + s(sqrt(mud_percent), k=4,bs="cr") +  s(temperature, k=4, bs="cr") + s(vms_ospar, k=4,bs="cr") + s(station, bs="re"), correlation=corExp(form= ~ x + y_adj),method = "ML", data = srf_var, family = gaussian(link="identity"),na.action=na.fail)
  drd_uGamm_methane_oxidation <- dredge(uGamm_methane_oxidation,trace = 2, fixed = ~ log(seq_depth) + s(station,bs = "re"))
  data.frame("RI"=sw(subset(drd_uGamm_methane_oxidation, delta <= 4)))
  gamm_methane_oxidation <- gamm(get.models(drd_uGamm_methane_oxidation,1)[[1]]$call$formula, 
                           correlation=corExp(form= ~ x + y_adj),
                           method = "REML", data = srf_var)
  summary(gamm_methane_oxidation$gam)
  print(plot.gamViz(getViz(gamm_methane_oxidation$gam), allTerms = F) + 
          l_ciPoly(level=0.95,fill="grey",alpha=0.75) +
          l_ciLine(level=0.95, colour = "black", linetype = 1,size=0.1) +
          l_fitLine(colour = "black",size=.5) + l_rug(mapping = aes(x=x, y=y)) +
          l_points(shape = 19, size = 1, alpha = 0.1) + theme_bw(), pages = 1)
}

## CO2_reduction ##
{ 
  uGamm_CO2_reduction <- uGamm(log(CO2_reduction) ~ log(seq_depth) + s(grain_phi, k=4, bs="cr") + s(logit(toc_proportion), k=4, bs="cr") + s(shear_stress, k=4, bs="cr") + s(sqrt(mud_percent), k=4,bs="cr") +  s(temperature, k=4, bs="cr") + s(vms_ospar, k=4,bs="cr") + s(station, bs="re"), correlation=corExp(form= ~ x + y_adj),method = "ML", data = srf_var, family = gaussian(link="identity"),na.action=na.fail)
  drd_uGamm_CO2_reduction <- dredge(uGamm_CO2_reduction,trace = 2, fixed = ~ log(seq_depth) + s(station,bs = "re"))
  data.frame("RI"=sw(subset(drd_uGamm_CO2_reduction, delta <= 4)))
  gamm_CO2_reduction <- gamm(get.models(drd_uGamm_CO2_reduction,1)[[1]]$call$formula, 
                           correlation=corExp(form= ~ x + y_adj),
                           method = "REML", data = srf_var)
  summary(gamm_CO2_reduction$gam)
  print(plot.gamViz(getViz(gamm_CO2_reduction$gam), allTerms = F) + 
          l_ciPoly(level=0.95,fill="grey",alpha=0.75) +
          l_ciLine(level=0.95, colour = "black", linetype = 1,size=0.1) +
          l_fitLine(colour = "black",size=.5) + l_rug(mapping = aes(x=x, y=y)) +
          l_points(shape = 19, size = 1, alpha = 0.1) + theme_bw(), pages = 1)
}


############################