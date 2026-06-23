#Plots and stats for paper 2 - SUBMISSION 2

library(dplyr)
library(ggplot2)


#Bar chart aesthetics - width=0.8, linewidth=1, error bar width=0.3, theme_classic, base_size = 35
#dotplot - dotsize = 2% of ymax on axis

# Figure 1B LCV Rab5 area plots DONE-----------------------------------------------
#Re-run CP analysis using smoothed LCV pipeline - helps with edge definition and avoiding splitting of bug
setwd("~/Desktop/Paper 2 figures/Figure 1 morphology + host regulators/20240829 LCV area SMOOTHED untransfected/")

image<- read.csv(file = "20240829_alltransfected_Rab5_WTddotAImage.csv")
legi<- read.csv(file = "20240829_alltransfected_Rab5_WTddotAmergedlegi.csv")
marker<-read.csv(file = "20240829_alltransfected_Rab5_WTddotAAF488.csv")
image2<-read.csv("20241219_MET049_morphImage.csv")
legi2<-read.csv("20241219_MET049_morphmergedlegi.csv")
marker2<-read.csv("20241219_MET049_morphAF488.csv")
codes<-read.csv("codes.csv")

filenames<-data.frame("Filename"=image$FileName_lp, "ImageNumber"=image$ImageNumber)
filenames2<-data.frame("Filename"=image2$FileName_lp,"ImageNumber"=image2$ImageNumber+1000)
filenames<-rbind(filenames,filenames2)
filenames$code<-toupper(substring(filenames$Filename, regexpr("MET|ERB",filenames$Filename) + 7,regexpr("MET|ERB",filenames$Filename)+10))
filenames$code<-gsub(" LP0","LP03",filenames$code)
filenames[which(filenames$code=="RAB5"),"code"]<-substring(filenames[which(filenames$code=="RAB5"),"Filename"], regexpr("MET|ERB",filenames[which(filenames$code=="RAB5"),"Filename"]) + 13,regexpr("MET|ERB",filenames[which(filenames$code=="RAB5"),"Filename"])+16)
filenames$experiment<-substring(filenames$Filename, regexpr("MET|ERB",filenames$Filename),regexpr("MET|ERB",filenames$Filename)+5)
codes$code<-toupper(substring(codes$code,1,4))
filenames<-merge(filenames,codes, all = T)

legi<-data.frame("ImageNumber"=legi$ImageNumber,"LegiNumber"=legi$ObjectNumber,"LegiArea"=legi$AreaShape_Area)
legi2<-data.frame("ImageNumber"=legi2$ImageNumber+1000,"LegiNumber"=legi2$ObjectNumber,"LegiArea"=legi2$AreaShape_Area)
legi<-rbind(legi,legi2)

marker<-data.frame("ImageNumber" = marker$ImageNumber, "MarkerNumber" = marker$ObjectNumber, "MarkerArea" =marker$AreaShape_Area)
marker2<-data.frame("ImageNumber" = marker2$ImageNumber + 1000, "MarkerNumber" = marker2$ObjectNumber, "MarkerArea" =marker2$AreaShape_Area)
marker<-rbind(marker,marker2)

alldata<-merge(legi, marker)
alldata<-merge(alldata, filenames)

#Find area ratio
alldata$ratio<-alldata$MarkerArea/alldata$LegiArea
# 
#trim high outliers - top 0.5%
alldata<-alldata[which(alldata$ratio<quantile(alldata$ratio,probs=0.995)),]

alldata$strain<-factor(alldata$strain,levels = c("WT","ddotA"))
alldata<-alldata[order(alldata$strain),]

#Summarize
alldata_summ<-summarise(group_by(alldata,experiment,strain), meanratio = mean(ratio), n = length(ratio))
alldata_xrep<-summarise(group_by(alldata,strain), meanratio = mean(ratio), n = length(ratio))

pdf(file = "Rab5Area_jitter_wider.pdf",height = 8,width = 5)
ggplot() + geom_jitter(data = alldata, aes(x = strain, y = ratio, fill = experiment), pch = 21, color = "black", width = 0.3, size = 3) + 
  geom_jitter(data = alldata_summ, aes(x = strain, y = meanratio, fill = experiment), pch = 21, color = "black", width = 0.3, size = 8) +
  scale_fill_manual(values = c("slateblue","navy","cornflowerblue","lightblue","grey95"))+
  geom_errorbar(data = alldata_xrep,aes(x = strain, ymax = meanratio, ymin = meanratio), width = 0.7, linewidth = 2)+
  theme_bw(base_size = 30) + theme(legend.position = "none")+ ylim(0,8)+
  labs(x = "Strain", y = "Rab5 area/L.p. area")
dev.off()


sink(file = "Rab5_area_WTvddotA_test.txt")
t.test(meanratio~strain, data = alldata_summ)
sink

write.csv(alldata_summ,"LCVcountsRab5Area.csv", row.names = F)

# Figure 1D UI cell endosome area DONE ----------------------------------------
#Endosome size analysis in UI cells from ddotA samples - masked based on GFP signal using best cyto detection pipeline 11/18
setwd("~/Desktop/Paper 2 figures/Figure 1 morphology + host regulators/20240808 endosome size analysis rerun rabex5 rabaptin5/")

image<-read.csv("20241118_MaskedRab5ASize_FINALImage.csv")
rab5<-read.csv("20241118_MaskedRab5ASize_FINALMaskedRab5.csv")
codes<-read.csv("codes.csv")

filenames<-data.frame("ImageNumber" = image$ImageNumber, "filename" = image$FileName_EGFP)
filenames$code<-toupper(substring(filenames$filename, regexpr("MET|ERB", filenames$filename) + 7, regexpr("MET|ERB", filenames$filename) + 9))
filenames$experiment<-substring(filenames$filename, regexpr("MET|ERB", filenames$filename), regexpr("MET|ERB", filenames$filename) + 5)
codes$code<-toupper(substring(codes$code, 1,3))
filenames<-merge(filenames, codes,all = T)

rab5min<-data.frame("ImageNumber" = rab5$ImageNumber, "Rab5Number" = rab5$ObjectNumber, "Rab5Area" = rab5$AreaShape_Area)

alldata<-merge(filenames, rab5min, all = T)
alldata$Rab5Area_um<-alldata$Rab5Area/5.2902^2
alldata<-alldata[!is.na(alldata$Rab5Area),]

#for the sake of plotting, collapse these wild reps into groups
alldata$experiment<-gsub("ERB128|ERB186", "rep1", alldata$experiment)
alldata$experiment<-gsub("MET029|ERB187", "rep2", alldata$experiment)
alldata$experiment<-gsub("MET039|ERB188", "rep3", alldata$experiment)


#calc per image means and averages across reps
alldata_PIM<-summarise(group_by(alldata, ImageNumber,experiment, marker), avg = mean(Rab5Area_um))
alldata_summ<-summarise(group_by(alldata,experiment, marker), avg = mean(Rab5Area_um), cells = length(unique(ImageNumber)))
alldata_xrep<-summarise(group_by(alldata,marker), avg = mean(Rab5Area_um))


#plot!

pdf("Rab5_Area_UI_masked_wider.pdf",width = 5,height = 8)
ggplot() + geom_jitter(data = alldata_PIM, aes(x = marker, y = avg, fill = experiment), pch = 21, color = "black", width = 0.3, size = 3)+
  geom_jitter(data = alldata_summ, aes(x = marker, y = avg, fill = experiment), pch = 21, color = "black", width = 0.3, size = 8)+
  scale_fill_manual(values = c("white","magenta","grey80")) + 
  geom_errorbar(data = alldata_xrep,aes(x = marker, ymax = avg, ymin = avg),width = 0.7, linewidth = 2)+
  theme_bw(base_size = 30) + theme(legend.position = "none")+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  labs(x = "Construct", y = "Average endosome area per cell (um^2)") + ylim(0,1.2)
dev.off()

#Do stats
bartlett.test(avg~marker,data = alldata_summ) #p>0.05, assume equal variance
model<-aov(avg~marker,data = alldata_summ)

sink(file = "Rab5EndoSize_ANOVA.txt")
bartlett.test(avg~marker,data = alldata_summ)
print("Bartlett p>0.05, assume equal variance")
print("Summary - one way ANOVA")
summary(model)
print("Tukey Kramer post-hoc test results, confidence interval 0.95, alpha=0.05")
TukeyHSD(model, conf.level = 0.95)
sink()

write.csv(alldata_summ, "EndosomeArea.csv",row.names = F)


# Figure 1E Rabex5 Rabaptin5 rec data DONE---------------------------------------------
setwd("/Volumes/AMS_ex/August 2024 trip analysis/Paper 2 figures/Figure 1 morphology + host regulators/")
bar<-read.csv("Endosome_marker_for_plot.csv")
bar$Condition<-paste(bar$Strain,bar$Condition)
bar$pct<-bar$N_pos/bar$N_LCVs*100
bar_r5<-bar[which(bar$Marker=="Rab5A"),]
bar_rest<-bar[which(bar$Marker!="Rab5A"),]

r5_summ<-summarise(group_by(bar_r5, Condition), mpos = mean(pct), sd = sd(pct), strain = unique(Strain))
rest_summ<-summarise(group_by(bar_rest, Condition), mpos = mean(pct), sd = sd(pct),strain = unique(Strain))

pdf(file = "Rab5Percentages_Endo.pdf", height = 8, width = 7)
ggplot() + geom_errorbar(data = r5_summ, aes(y = mpos, x = Condition, ymin = mpos-1, ymax = mpos + sd), width = 0.3, linewidth = 1)+ 
  geom_bar(data = r5_summ, aes(x =Condition, y = mpos, fill = strain), stat = "identity",color="black",width = 0.8, linewidth=1) + 
  geom_dotplot(data = bar_r5, aes(y = pct, x = Condition), stackdir = "center", binaxis = "y", binwidth = 1) + 
  scale_fill_manual(values = c("#BCD2EE","#6A5ACD"))+
  labs(y = "Percent Rab5+ LCVs", x = "Condition") + 
  theme_classic(base_size = 35)+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))  + scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0,50))
dev.off()



# Figure 1F LCV morphology Rabex5 Rabaptin5 OE DONE----------------------------------------------
#re-ran area analysis using smoothed LCV pipeline - cleans up bug object ID significantly
setwd("~/Desktop/Paper 2 figures/Figure 1 morphology + host regulators/20240827 LCV area SMOOTHED rabex5 rabaptin5/")
rep1<-read.csv("Rep1_RabexRabaptin5_Rab5Area.csv")
rep2<-read.csv("Rep2_RabexRabaptin5_Rab5Area.csv")
rep3<-read.csv("Rep3_RabexRabaptin5_Rab5Area.csv")

alldata<-rbind(rep1,rep2,rep3)
alldata<-alldata[which(alldata$ratio<quantile(alldata$ratio,probs = 0.995)),]

alldata$condition<-factor(alldata$condition, levels = c("egfp WT", "rabaptin5 WT", "rabex5 WT", "egfp ddotA","rabaptin5 ddotA","rabex5 ddotA"))
alldata<-alldata[order(alldata$condition),]

alldata_summ<-summarise(group_by(alldata,replicate,condition), meanratio = mean(ratio), n = length(ratio), strain = unique(strain))
alldata_xrep<-summarise(group_by(alldata,condition), meanratio = mean(ratio), n = length(ratio))


#plot

pdf("Rab5LCVareaRbxRbpt_wider.pdf", height = 8, width = 6)
ggplot() + geom_jitter(data = alldata, aes(x = condition, y = ratio, fill = replicate), pch = 21, color = "black", width = 0.3, size = 3)+
  geom_jitter(data = alldata_summ, aes(x = condition, y = meanratio, fill = replicate), pch = 21, color = "black", width = 0.3, size = 8)+
  scale_fill_manual(values = c("navy","cornflowerblue","lightblue"))+
  geom_errorbar(data = alldata_xrep,aes(x = condition, ymax = meanratio, ymin = meanratio),width = 0.7, linewidth = 2)+
  theme_bw(base_size = 30) + theme(legend.position = "none")+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  labs(x = "Construct", y = "Rab5+ area/L.p. area") + ylim(0,6.5)
dev.off()

#stats
bartlett.test(meanratio~condition,data = alldata_summ) #p<0.05, need Games-Howell, bleh
games.howell <- function(grp, obs) {
  
  #Create combinations
  combs <- combn(unique(grp), 2)
  
  # Statistics that will be used throughout the calculations:
  # n = sample size of each group
  # groups = number of groups in data
  # Mean = means of each group sample
  # std = variance of each group sample
  n <- tapply(obs, grp, length)
  groups <- length(tapply(obs, grp, length))
  Mean <- tapply(obs, grp, mean)
  std <- tapply(obs, grp, var)
  
  statistics <- lapply(1:ncol(combs), function(x) {
    
    mean.diff <- Mean[combs[2,x]] - Mean[combs[1,x]]
    
    #t-values
    t <- abs(Mean[combs[1,x]] - Mean[combs[2,x]]) / sqrt((std[combs[1,x]] / n[combs[1,x]]) + (std[combs[2,x]] / n[combs[2,x]]))
    
    # Degrees of Freedom
    df <- (std[combs[1,x]] / n[combs[1,x]] + std[combs[2,x]] / n[combs[2,x]])^2 / # Numerator Degrees of Freedom
      ((std[combs[1,x]] / n[combs[1,x]])^2 / (n[combs[1,x]] - 1) + # Part 1 of Denominator Degrees of Freedom 
         (std[combs[2,x]] / n[combs[2,x]])^2 / (n[combs[2,x]] - 1)) # Part 2 of Denominator Degrees of Freedom
    
    #p-values
    p <- ptukey(t * sqrt(2), groups, df, lower.tail = FALSE)
    
    # Sigma standard error
    se <- sqrt(0.5 * (std[combs[1,x]] / n[combs[1,x]] + std[combs[2,x]] / n[combs[2,x]]))
    
    # Upper Confidence Limit
    upper.conf <- lapply(1:ncol(combs), function(x) {
      mean.diff + qtukey(p = 0.95, nmeans = groups, df = df) * se
    })[[1]]
    
    # Lower Confidence Limit
    lower.conf <- lapply(1:ncol(combs), function(x) {
      mean.diff - qtukey(p = 0.95, nmeans = groups, df = df) * se
    })[[1]]
    
    # Group Combinations
    grp.comb <- paste(combs[1,x], ':', combs[2,x])
    
    # Collect all statistics into list
    stats <- list(grp.comb, mean.diff, se, t, df, p, upper.conf, lower.conf)
  })
  
  # Unlist statistics collected earlier
  stats.unlisted <- lapply(statistics, function(x) {
    unlist(x)
  })
  
  # Create dataframe from flattened list
  results <- data.frame(matrix(unlist(stats.unlisted), nrow = length(stats.unlisted), byrow=TRUE))
  
  # Select columns set as factors that should be numeric and change with as.numeric
  results[c(2, 3:ncol(results))] <- round(as.numeric(as.matrix(results[c(2, 3:ncol(results))])), digits = 3)
  
  # Rename data frame columns
  colnames(results) <- c('groups', 'Mean Difference', 'Standard Error', 't', 'df', 'p', 'upper limit', 'lower limit')
  
  return(results)
} #code from userfriendlyscience package - package itself no longer supported, but function code available here: https://rpubs.com/aaronsc32/games-howell-test
oneway.test(meanratio~condition,data = alldata_summ, var.equal = F)
games.howell(alldata_summ$condition,alldata_summ$meanratio)

sink(file = "Rab5AreaLCV_Rabex5Rbpt5_stats.txt")
bartlett.test(meanratio~condition,data = alldata_summ)
print("Bartlett p<0.05, unequal variance")
oneway.test(meanratio~condition,data = alldata_summ, var.equal = F)
print("Games-Howell test- from userfriendlyscience package")
games.howell(alldata_summ$condition,alldata_summ$meanratio)
sink()

write.csv(alldata_summ,"LCVcountsRabex5Rabapatin5.csv", row.names = F)

# Figure 2B SidE fam Rab5 cloud - updated to include sidc/sdca KO DONE--------------------------------------------
setwd("/Volumes/AMS_ex/Paper 2 figures/Submission 2/Figure 2 Rab5 SidE fam submission 2/Area plot KO strains/")

image<- read.csv(file = "20250425_Rab5LCVarea_dsidCsdcAImage.csv")
legi<- read.csv(file = "20250425_Rab5LCVarea_dsidCsdcAmergedlegi.csv")
marker<-read.csv(file = "20250425_Rab5LCVarea_dsidCsdcAAF488.csv")
codes<-read.csv("codes.csv")

filenames<-data.frame("Filename"=image$FileName_lp, "ImageNumber"=image$ImageNumber)
codes$code<-toupper(substring(codes$code,1,4))
filenames$code<-substring(filenames$Filename,regexpr(paste(codes$code,collapse = "|"), filenames$Filename),regexpr(paste(codes$code,collapse = "|"),filenames$Filename)+3)
filenames<-merge(filenames, codes)
filenames$experiment<-substring(filenames$Filename, regexpr("MET|ERB",filenames$Filename),regexpr("MET|ERB",filenames$Filename)+5)

filenames[which(filenames$experiment=="ERB155"|filenames$experiment=="MET041"),"replicate"]<-"Rep 1"

filenames[which(filenames$experiment=="MET046"|filenames$experiment=="ERB134"|filenames$experiment=="ERB141"|filenames$experiment=="ERB180"),"replicate"]<-"Rep 2"

filenames[which(filenames$experiment=="MET062"|filenames$experiment=="ERB130"|filenames$experiment=="ERB185"|filenames$experiment=="ERB150"),"replicate"]<-"Rep 3"

legi<-data.frame("ImageNumber"=legi$ImageNumber,"LegiNumber"=legi$ObjectNumber,"LegiArea"=legi$AreaShape_Area)
marker<-data.frame("ImageNumber" = marker$ImageNumber,"LegiNumber" = marker$Parent_shrunkenlegi,"MarkerNumber" = marker$ObjectNumber, "MarkerArea" =marker$AreaShape_Area)

alldata<-merge(legi, marker)
alldata<-merge(alldata, filenames)

alldata<-alldata[which(alldata$strain!="LEG179"),]

#Find area ratio
alldata$ratio<-alldata$MarkerArea/alldata$LegiArea

#trim data
alldata<-alldata[which(alldata$ratio<quantile(alldata$ratio, probs = 0.995)),]

alldata$strain<-factor(alldata$strain,levels = c("WT","KO","pSdeB","pSdeB EE/AA", "LEG173", "LEG180","LEG181"))
alldata<-alldata[order(alldata$strain),]

#Summarize
alldata_summ<-summarise(group_by(alldata,replicate,strain), meanratio = mean(ratio), n = length(ratio))

alldata_xrep<-summarise(group_by(alldata,strain), meanratio = mean(ratio), n = length(ratio))

pdf(file = "Rab5Area_fullpaneljitter_wider.pdf",height = 8,width = 8)
ggplot() + geom_jitter(data = alldata, aes(x = strain, y = ratio, fill = replicate), pch = 21, color = "black", width = 0.3, size = 3) + 
  geom_jitter(data = alldata_summ, aes(x = strain, y = meanratio, fill = replicate), pch = 21, color = "black", width = 0.3, size = 8) +
  scale_fill_manual(values = c("slateblue","navy","lightblue"))+
  geom_errorbar(data = alldata_xrep,aes(x = strain, ymax = meanratio, ymin = meanratio), width = 0.7, linewidth = 2)+
  theme_bw(base_size = 30) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none")+
  labs(x = NULL, y = "Rab5 area/L.p. area") + ylim(0, 6.5)
dev.off()

bartlett.test(meanratio~strain,data = alldata_summ) #p>0.05, assume equal variance
model<-aov(meanratio~strain,data = alldata_summ)

sink(file = "FullPanelRab5area_ANOVA.txt")
bartlett.test(meanratio~strain,data = alldata_summ) 
print("Bartlett p>0.05, assume equal variance")
print("Summary - one way ANOVA")
summary(model)
print("Tukey Kramer post-hoc test results, confidence interval 0.95, alpha=0.05")
TukeyHSD(model, conf.level = 0.95)
sink()



# Figure 3B UI cell washout ------------------------------------------------
setwd("~/Desktop/Paper 2 figures/Figure 3 Rab5 detergent resistance/UI cell washout/")
met014<-read.csv("MET014_CTCF_alldata.csv")
met017<-read.csv("MET017_CTCF_alldata.csv")
met020<-read.csv("MET020_CTCF_alldata.csv")

met014$rep<-"Rep1"
met017$rep<-"Rep2"
met020$rep<-"Rep3"

alldata<-rbind(met014,met017,met020)
alldata$ctcf_scaled<-alldata$CTCF/10000

alldata_summ<-summarise(group_by(alldata,rep,perm), meanctcf= mean(ctcf_scaled), sdctcf = sd(ctcf_scaled), n = length(ctcf_scaled))
alldata_xrep<-summarise(group_by(alldata,perm),meanctcf= mean(ctcf_scaled), sdctcf = sd(ctcf_scaled), n = length(ctcf_scaled))

pdf("Rab5UIWashout.pdf",height = 6,width = 5)
ggplot() + geom_jitter(data = alldata, aes(x = perm, y = ctcf_scaled, fill = rep), pch = 21, color = "black", width = 0.3, size = 2) + 
  geom_jitter(data = alldata_summ, aes(x = perm, y = meanctcf, fill = rep), pch = 21, color = "black", width = 0.3, size = 8) +
  scale_fill_manual(values = c("white","lightgray","grey30"))+
  geom_errorbar(data = alldata_xrep,aes(x = perm, ymax = meanctcf, ymin = meanctcf), width = 0.8, linewidth = 2)+
  theme_bw(base_size = 30) + theme(legend.position = "none")+ ylim(0,10)+
  labs(y = "Norm. Rab5 intensity per cell (a.u.)")
dev.off()

sink(file = "Rab5_UIWashout_Ttest.txt")
t.test(meanctcf~perm, data = alldata_summ)
sink


# Figure 3D Rab5 detergent washout - updated log2fc DONE-----------------------------------------
setwd("~/Desktop/Paper 2 figures/Submission 2/Figure 3 Rab5 detergent resistance/")
met053<-read.csv("MET053_SDS_intensity_wholeFOV_final.csv")
met059<-read.csv("MET059_SDS_intensity_wholeFOV_final.csv")
met062<-read.csv("MET062_SDS_intensity_wholeFOV_final.csv")

alldata<-rbind(met053,met059, met062)

alldata$fc<-alldata$MeanLegiInt/alldata$BG
alldata$log2<-log2(alldata$fc)

alldata$condition<-factor(alldata$condition,levels = c("WT_saponin","ddotA_saponin","WT_SDS","ddotA_SDS"))
alldata<-alldata[order(alldata$condition),]

fc_summ<-alldata %>% group_by(experiment, condition) %>% summarise(mean = mean(log2), n = length(log2))
fc_xrep<-alldata %>% group_by(condition) %>% summarise(mean = mean(log2), perm = unique(perm))
# 
# pdf(file = "Rab5FC_violin_classic.pdf",height = 6,width = 5)
# ggplot() + geom_violin(data = alldata, aes(x = condition, y =fc, fill = perm), width = 1.3)+
#   geom_dotplot(data = fc_summ,aes(x = condition, y = meanfc), stackdir = "center",binaxis = "y",binwidth = 0.015)+
#   scale_fill_manual(values = c("deeppink1","turquoise3"))+
#   theme_classic(base_size = 25) + theme(legend.position = "none") +labs(x = NULL,y = "Norm. Rab5 intensity (a.u.)")+
#   ylim(0,2.5)
# dev.off()

pdf(file = "Rab5WashoutLog2.pdf", height = 6, width = 5)
ggplot()+ 
    geom_violin(data = alldata, aes(x = condition, y =log2, fill = perm), width = 1.1)+
    geom_dotplot(data = fc_summ,aes(x = condition, y = mean), stackdir = "center",binaxis = "y",binwidth = 0.015)+
    scale_fill_manual(values = c("deeppink1","turquoise3"))+
    theme_classic(base_size = 25) + theme(legend.position = "none") +labs(x = NULL,y = "log2(Rab5LCV/Rab5BG)")+
    scale_x_discrete(labels = c("WT_saponin" = "WT", "ddotA_saponin" = "dotA", "WT_SDS" = "WT", "ddotA_SDS" = "dotA"))+
    ylim(-0.4, 1.2)
dev.off()  

bartlett.test(mean~condition,data = fc_summ)
games.howell <- function(grp, obs) {
  
  #Create combinations
  combs <- combn(unique(grp), 2)
  
  # Statistics that will be used throughout the calculations:
  # n = sample size of each group
  # groups = number of groups in data
  # Mean = means of each group sample
  # std = variance of each group sample
  n <- tapply(obs, grp, length)
  groups <- length(tapply(obs, grp, length))
  Mean <- tapply(obs, grp, mean)
  std <- tapply(obs, grp, var)
  
  statistics <- lapply(1:ncol(combs), function(x) {
    
    mean.diff <- Mean[combs[2,x]] - Mean[combs[1,x]]
    
    #t-values
    t <- abs(Mean[combs[1,x]] - Mean[combs[2,x]]) / sqrt((std[combs[1,x]] / n[combs[1,x]]) + (std[combs[2,x]] / n[combs[2,x]]))
    
    # Degrees of Freedom
    df <- (std[combs[1,x]] / n[combs[1,x]] + std[combs[2,x]] / n[combs[2,x]])^2 / # Numerator Degrees of Freedom
      ((std[combs[1,x]] / n[combs[1,x]])^2 / (n[combs[1,x]] - 1) + # Part 1 of Denominator Degrees of Freedom 
         (std[combs[2,x]] / n[combs[2,x]])^2 / (n[combs[2,x]] - 1)) # Part 2 of Denominator Degrees of Freedom
    
    #p-values
    p <- ptukey(t * sqrt(2), groups, df, lower.tail = FALSE)
    
    # Sigma standard error
    se <- sqrt(0.5 * (std[combs[1,x]] / n[combs[1,x]] + std[combs[2,x]] / n[combs[2,x]]))
    
    # Upper Confidence Limit
    upper.conf <- lapply(1:ncol(combs), function(x) {
      mean.diff + qtukey(p = 0.95, nmeans = groups, df = df) * se
    })[[1]]
    
    # Lower Confidence Limit
    lower.conf <- lapply(1:ncol(combs), function(x) {
      mean.diff - qtukey(p = 0.95, nmeans = groups, df = df) * se
    })[[1]]
    
    # Group Combinations
    grp.comb <- paste(combs[1,x], ':', combs[2,x])
    
    # Collect all statistics into list
    stats <- list(grp.comb, mean.diff, se, t, df, p, upper.conf, lower.conf)
  })
  
  # Unlist statistics collected earlier
  stats.unlisted <- lapply(statistics, function(x) {
    unlist(x)
  })
  
  # Create dataframe from flattened list
  results <- data.frame(matrix(unlist(stats.unlisted), nrow = length(stats.unlisted), byrow=TRUE))
  
  # Select columns set as factors that should be numeric and change with as.numeric
  results[c(2, 3:ncol(results))] <- round(as.numeric(as.matrix(results[c(2, 3:ncol(results))])), digits = 3)
  
  # Rename data frame columns
  colnames(results) <- c('groups', 'Mean Difference', 'Standard Error', 't', 'df', 'p', 'upper limit', 'lower limit')
  
  return(results)
} #code from userfriendlyscience package - package itself no longer supported, but function code available here: https://rpubs.com/aaronsc32/games-howell-test
oneway.test(mean~condition,data = fc_summ, var.equal = F)
games.howell(fc_summ$condition,fc_summ$mean)

sink(file = "Rab5Washout_log2FC.txt")
bartlett.test(mean~condition,data = fc_summ)
print("p<0.05, unequal variance")
oneway.test(mean~condition,data = fc_summ, var.equal = F)
print("Games-Howell test, code from userfriendlyscience package")
games.howell(fc_summ$condition,fc_summ$mean)
sink()

# Figure 4A SidE fam Ub rec data DONE-------------------------------------------

setwd("/Volumes/AMS_ex/August 2024 trip analysis/Ub recruitment data/")
rec_data<-read.csv("Ubrec_compiled.csv")
rec_data$Strain<-factor(rec_data$Strain, levels = c("WT", "LEG151", "LEG151_SdeB", "LEG151_EEAA"))
rec_data<-rec_data[order(rec_data$Strain),]

rec_data_summ<-summarise(group_by(rec_data, Strain), mean = mean(Pct_Ub), sd = sd(Pct_Ub))

pdf(file = "UbRecSidEFam.pdf", height = 9, width = 5)
ggplot() + geom_errorbar(data = rec_data_summ, aes(y = mean, x = Strain, ymin = mean-1, ymax = mean + sd), width = 0.3, linewidth = 1)+ 
  geom_bar(data = rec_data_summ, aes(x =Strain, y = mean), stat = "identity",color="black",width = 0.8, linewidth=1, fill = "gray95") + 
  geom_dotplot(data = rec_data, aes(y = Pct_Ub, x = Strain), stackdir = "center", binaxis = "y", binwidth = 2) + 
  labs(y = "Percent Ub+ LCVs", x = "Strain") + 
  theme_classic(base_size = 35)+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))  + scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0,100))
dev.off()

# Figure 4C ub morph - updated how data was trimmed DONE------------------------------
setwd("~/Desktop/Paper 2 figures/Submission 2/Figure 4 ubiquitin morphology/")
met041<-read.csv("UbMorph_Rep1.csv")
met046<-read.csv("UbMorph_Rep2.csv")
met049<-read.csv("UbMorph_Rep3.csv")
met062<-read.csv("UbMorph_Rep4.csv")

#Trim top 0.5% of pooled data to remove extreme outliers 
alldata<-rbind(met041,met046,met049, met062)
alldata<-alldata[which(alldata$ratio<quantile(alldata$ratio, probs = 0.995)[[1]]),]

alldata$strain<-factor(alldata$strain, levels = c("WT","LEG151","LEG151_WT","LEG151_EEAA"))
alldata<-alldata[order(alldata$strain),]

alldata_summ<-summarise(group_by(alldata, experiment, strain), mean_ratio = mean(unique(ratio)), n = length(ratio))
alldata_xrep<-summarise(group_by(alldata, strain), mean_ratio = mean(ratio), n = length(ratio))

pdf("UbAreaSidEwider_retrimmed.pdf",width = 6,height = 8)
ggplot() + geom_jitter(data = alldata, aes(x = strain, y = ratio, fill = experiment), pch = 21, color = "black", width = 0.3, size = 2) + 
  geom_jitter(data = alldata_summ, aes(x = strain, y = mean_ratio, fill = experiment), pch = 21, color = "black", width = 0.3, size = 8) +
  scale_fill_manual(values = c("navy","slateblue","cornflowerblue","lightblue"))+
  geom_errorbar(data = alldata_xrep,aes(x = strain, ymax = mean_ratio, ymin = mean_ratio), width = 0.8, linewidth = 2)+
  theme_bw(base_size = 30) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none")+ ylim(0,5)+
  labs(x = NULL, y = "Ub area/L.p. area") + scale_x_discrete(labels = c("LEG151"="-",
                                                                            "LEG151_WT" = "+pSdeB WT",
                                                                            "LEG151_EEAA" = "+pSdeB EEAA"))
dev.off()


#stats
bartlett.test(mean_ratio~strain,data = alldata_summ) #p<0.05, need Games-Howell, bleh
games.howell <- function(grp, obs) {
  
  #Create combinations
  combs <- combn(unique(grp), 2)
  
  # Statistics that will be used throughout the calculations:
  # n = sample size of each group
  # groups = number of groups in data
  # Mean = means of each group sample
  # std = variance of each group sample
  n <- tapply(obs, grp, length)
  groups <- length(tapply(obs, grp, length))
  Mean <- tapply(obs, grp, mean)
  std <- tapply(obs, grp, var)
  
  statistics <- lapply(1:ncol(combs), function(x) {
    
    mean.diff <- Mean[combs[2,x]] - Mean[combs[1,x]]
    
    #t-values
    t <- abs(Mean[combs[1,x]] - Mean[combs[2,x]]) / sqrt((std[combs[1,x]] / n[combs[1,x]]) + (std[combs[2,x]] / n[combs[2,x]]))
    
    # Degrees of Freedom
    df <- (std[combs[1,x]] / n[combs[1,x]] + std[combs[2,x]] / n[combs[2,x]])^2 / # Numerator Degrees of Freedom
      ((std[combs[1,x]] / n[combs[1,x]])^2 / (n[combs[1,x]] - 1) + # Part 1 of Denominator Degrees of Freedom 
         (std[combs[2,x]] / n[combs[2,x]])^2 / (n[combs[2,x]] - 1)) # Part 2 of Denominator Degrees of Freedom
    
    #p-values
    p <- ptukey(t * sqrt(2), groups, df, lower.tail = FALSE)
    
    # Sigma standard error
    se <- sqrt(0.5 * (std[combs[1,x]] / n[combs[1,x]] + std[combs[2,x]] / n[combs[2,x]]))
    
    # Upper Confidence Limit
    upper.conf <- lapply(1:ncol(combs), function(x) {
      mean.diff + qtukey(p = 0.95, nmeans = groups, df = df) * se
    })[[1]]
    
    # Lower Confidence Limit
    lower.conf <- lapply(1:ncol(combs), function(x) {
      mean.diff - qtukey(p = 0.95, nmeans = groups, df = df) * se
    })[[1]]
    
    # Group Combinations
    grp.comb <- paste(combs[1,x], ':', combs[2,x])
    
    # Collect all statistics into list
    stats <- list(grp.comb, mean.diff, se, t, df, p, upper.conf, lower.conf)
  })
  
  # Unlist statistics collected earlier
  stats.unlisted <- lapply(statistics, function(x) {
    unlist(x)
  })
  
  # Create dataframe from flattened list
  results <- data.frame(matrix(unlist(stats.unlisted), nrow = length(stats.unlisted), byrow=TRUE))
  
  # Select columns set as factors that should be numeric and change with as.numeric
  results[c(2, 3:ncol(results))] <- round(as.numeric(as.matrix(results[c(2, 3:ncol(results))])), digits = 3)
  
  # Rename data frame columns
  colnames(results) <- c('groups', 'Mean Difference', 'Standard Error', 't', 'df', 'p', 'upper limit', 'lower limit')
  
  return(results)
} #code from userfriendlyscience package - package itself no longer supported, but function code available here: https://rpubs.com/aaronsc32/games-howell-test
oneway.test(mean_ratio~strain,data = alldata_summ, var.equal = F)
games.howell(alldata_summ$strain,alldata_summ$mean_ratio)

sink(file = "UbAreaLCV_SidEpanel_Retrimmed.txt")
bartlett.test(mean_ratio~strain,data = alldata_summ)
print("Bartlett p<0.05, unequal variance")
oneway.test(mean_ratio~strain,data = alldata_summ, var.equal = F)
print("Games-Howell test- from userfriendlyscience package")
games.howell(alldata_summ$strain,alldata_summ$mean_ratio)
sink()


# Figure 5B SDS washout Ub - updated log2fc DONE-------------------------------------------------
setwd("~/Desktop/Paper 2 figures/Submission 2/Figure 5 ubiquitin detergent/SDS_int/")
rep1<-read.csv("Rep1Ub_SDS_intensity_wholeFOV_final.csv")
rep2<-read.csv("Rep2Ub_SDS_intensity_wholeFOV_final.csv")
rep3<-read.csv("Rep3Ub_SDS_intensity_wholeFOV_final.csv")

alldata<-rbind(rep1,rep2,rep3)

alldata$fc<-(alldata$LegiInt/alldata$LegiArea)/alldata$BG
alldata$log2fc<-log2(alldata$fc)

alldata$condition<-factor(alldata$condition, levels = c("WT_saponin",
                                                        "LEG151_saponin",
                                                        "LEG151_WT_saponin",
                                                        "LEG151_EEAA_saponin",
                                                        "WT_SDS",
                                                        "LEG151_SDS",
                                                        "LEG151_WT_SDS",
                                                        "LEG151_EEAA_SDS"
))
alldata<-alldata[order(alldata$condition),]


fc_summ<-alldata %>% group_by(replicate,condition) %>% summarise(meanlog = mean(log2fc), n = length(log2fc))

pdf(file = "Ublog2FC_violin_NoTrim.pdf",width = 8,height = 6)
ggplot() + geom_violin(data=alldata, aes(x = condition, y = log2fc,fill = perm), trim = F)+
  geom_dotplot(data = fc_summ, aes(x = condition, y = meanlog), stackdir = "center",binaxis = "y",binwidth = 0.07)+
  theme_classic(base_size = 25) +  scale_fill_manual(values = c("deeppink1","turquoise3"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none")+
  labs(x = NULL, y = "log2(UbLCV/UbBG)") + ylim(-1.2, 5.5)+
  scale_x_discrete(labels = c("WT_saponin" = "WT",
                              "LEG151_saponin" = "-",
                              "LEG151_WT_saponin" = "+pSdeB WT",
                              "LEG151_EEAA_saponin" = "+pSdeB EEAA",
                              "WT_SDS" = "WT",
                              "LEG151_SDS" = "-",
                              "LEG151_WT_SDS" = "+pSdeB WT",
                              "LEG151_EEAA_SDS" = "+pSdeB EEAA"
                              ))
dev.off()



#Do stats
bartlett.test(meanlog~condition,fc_summ) #p>0.05, assume equal variance
model<-aov(meanlog~condition,fc_summ)

sink(file = "UbWashout_log2FC_ANOVA.txt")
bartlett.test(meanlog~condition,fc_summ)
print("Bartlett p>0.05, assume equal variance")
print("Summary - one way ANOVA")
summary(model)
print("Tukey Kramer post-hoc test results, confidence interval 0.95, alpha=0.05")
TukeyHSD(model, conf.level = 0.95)
sink()


# Figure 5C UI vs LP02 cyto Ub washout ------------------------------------
setwd("/Volumes/AMS_ex/Paper 2 figures/Submission 2/Figure 5 ubiquitin detergent/WT vs UI cyto ub SDS intensity/")

met046<-read.csv("MET046_WTUI_cyto.csv")
met049<-read.csv("MET049_WTUI_cyto.csv")
met062<-read.csv("MET062_WTUI_cyto.csv")

alldata<-rbind(met046,met049,met062)
alldata$rep<-substring(alldata$image,regexpr("MET",alldata$image), regexpr("MET",alldata$image)+5)

alldata$l2fc<-log2(alldata$fc_int)

alldata_summ<-alldata %>% group_by(rep,ID) %>% summarise(mean_l2fc = mean(l2fc))
alldata_xrep<-alldata %>% group_by(ID) %>% summarise(mean_l2fc = mean(l2fc))


pdf(file = "UI_WT_cyto_ub_SDS_violinl2.pdf", width = 4, height = 6)
ggplot() + geom_violin(data = alldata, aes(x = ID, y = l2fc), trim = F, fill = "lightgrey")+
  geom_dotplot(data = alldata_summ, aes(x = ID, y = mean_l2fc), binwidth = 0.02, binaxis = "y", stackdir = "center")+
  theme_classic(base_size = 30) + theme(legend.position = "none") + ylim(-0.1,1)+
  labs(x = NULL, y = "log2(Ubcyto/BG)")
dev.off()

sink(file = "UIvsWT_cyto_Ubwashout_ttestl2.txt")
t.test(mean_l2fc~ID,data = alldata_summ)
sink()


# Figure 5F sol insol blots SidE fam ------------------------------------------------

setwd("~/Desktop/Paper 2 figures/Submission 2/Figure 5 ubiquitin detergent/Sol insol - P4D1/")
#Had to re-normalize everything to dotA, because Chetan forgot to do UI for one of the SidC/SdcA ko reps lol
#Probably a better control anyways so that's fine

solfrac<-read.csv("SolInsol_SidE_dotAnorm.csv")
solfrac<-solfrac[which(solfrac!="HS"),] #don't actually care about quantifying this sample
solfrac<-na.omit(solfrac)
solfrac$condition<-paste(solfrac$Sample,solfrac$Fraction, sep = "_")
solfrac$condition<-factor(solfrac$condition,levels = c("WT_Sol","151_Sol","SdeB_Sol","EEAA_Sol",
                                                       "WT_Insol","151_Insol","SdeB_Insol","EEAA_Insol"))
solfrac<-solfrac[order(solfrac$condition),]

summ<-solfrac %>% group_by(condition) %>% summarise(avg = mean(FC_dotA),sd = sd(FC_dotA), fraction = unique(Fraction))

pdf("SidEfam_SolInsol.pdf",height = 8, width = 8)
ggplot() + geom_errorbar(data = summ, aes(y = avg, x = condition, ymin = avg-0.5, ymax = avg + sd), width = 0.3, linewidth = 1)+ 
  geom_bar(data = summ, aes(x =condition, y = avg, fill = fraction), stat = "identity",color="black",width = 0.8, linewidth=1) + 
  geom_dotplot(data = solfrac, aes(y = FC_dotA, x = condition), stackdir = "center", binaxis = "y", binwidth = .1) + 
  scale_fill_manual(values = c("lightblue", "lightgrey")) + labs(y = "FC HMW Ub over dotA", x = NULL) + 
  theme_classic(base_size = 35)+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))  + scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0,4))+
  scale_x_discrete(labels = c("WT","-","+pSdeB WT", "+pSdeB EEAA","WT","-","+pSdeB WT", "+pSdeB EEAA"))
dev.off()

#Stats

bartlett.test(FC_dotA~condition,data = solfrac)
#ok variance is unequal when pooling data, so just split into sol vs insol (since I don't care about comparisons between fractions - do this for all sol/insol experiments)

bartlett.test(FC_dotA~condition,data = solfrac[which(solfrac$Fraction=="Insol"),]) 
model<-aov(FC_dotA~condition,data = solfrac[which(solfrac$Fraction=="Insol"),])

sink(file = "InsolSol_ANOVA_dotA.txt")
bartlett.test(FC_dotA~condition,data = solfrac[which(solfrac$Fraction=="Insol"),])
print("Bartlett p>0.05, assume equal variance")
print("Summary - one way ANOVA")
summary(model)
print("Tukey Kramer post-hoc test results, confidence interval 0.95, alpha=0.05")
TukeyHSD(model, conf.level = 0.95)
sink()

bartlett.test(FC_dotA~Sample,data = solfrac[which(solfrac$Fraction=="Sol"),]) #p>0.05, assume equal variance
model<-aov(FC_dotA~Sample,data = solfrac[which(solfrac$Fraction=="Sol"),])

sink(file = "Sol_ANOVA_dotA.txt")
bartlett.test(FC_dotA~Sample,data = solfrac[which(solfrac$Fraction=="Sol"),])
print("Bartlett p>0.05, assume equal variance")
print("Summary - one way ANOVA")
summary(model)
print("Tukey Kramer post-hoc test results, confidence interval 0.95, alpha=0.05")
TukeyHSD(model, conf.level = 0.95)
sink()


# Figure 5I sol insol SidC SdcA --------------------------------------------
setwd("~/Desktop/Paper 2 figures/Submission 2/Figure 5 ubiquitin detergent/SidCSdcA SolInsol Infection Experiment/AMS analysis/")
solfrac<-read.csv("dSS_AMS_solinsol.csv")
solfrac<-solfrac[which(solfrac$Sample!="HS"),]
solfrac$condition<-paste(solfrac$Sample,solfrac$Fraction,sep = "_")
solfrac$condition<-factor(solfrac$condition, levels = c("WT_Sol","dSS_Sol","dSS_V_Sol","dSS_SdcA_Sol","dSS_SidC_Sol",
                                                        "WT_Insol","dSS_Insol","dSS_V_Insol","dSS_SdcA_Insol","dSS_SidC_Insol"))
summ<-solfrac %>% group_by(condition) %>% summarise(avg = mean(FC_dotA),sd = sd(FC_dotA), fraction = unique(Fraction))

pdf(file = "SolInsol_dSS.pdf",height = 8,width = 10)
ggplot() + geom_errorbar(data = summ, aes(y = avg, x = condition, ymin = avg-0.5, ymax = avg + sd), width = 0.3, linewidth = 1)+ 
  geom_bar(data = summ, aes(x =condition, y = avg, fill = fraction), stat = "identity",color="black",width = 0.8, linewidth=1) + 
  geom_dotplot(data = solfrac, aes(y = FC_dotA, x = condition), stackdir = "center", binaxis = "y", binwidth = .1) + 
  scale_fill_manual(values = c("lightblue", "lightgrey")) + labs(y = "FC HMW Ub over dotA", x = NULL) + 
  theme_classic(base_size = 35)+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))  + scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0,4))+
  scale_x_discrete(labels = c("WT","-","+vector", "+pSdcA","+pSidC","WT","-","+vector", "+pSdcA","+pSidC"))
dev.off()

sink(file = "dSSInsol_dotAnorm.txt")
bartlett.test(FC_dotA~condition, data = solfrac[which(solfrac$Fraction=="Insol"),])
print("p>0.05, assume equal variance")
model<-aov(FC_dotA~condition, data = solfrac[which(solfrac$Fraction=="Insol"),])
print("Summary - one way ANOVA")
summary(model)
print("Tukey Kramer post-hoc test results, confidence interval 0.95, alpha=0.05")
TukeyHSD(model, conf.level = 0.95)
sink()

sink(file = "dSSSol_dotAnorm.txt")
bartlett.test(FC_dotA~condition, data = solfrac[which(solfrac$Fraction=="Sol"),])
print("p>0.05, assume equal variance")
model<-aov(FC_dotA~condition, data = solfrac[which(solfrac$Fraction=="Sol"),])
print("Summary - one way ANOVA")
summary(model)
print("Not significant!")
sink()

# Figure 6B IF washout timecourse - updated log2fc DONE------------------------------------------

setwd("/Volumes/AMS_ex/Paper 2 figures/Submission 2/Figure 6 timecourse and DupA/IF timecourse intensity and area quant/")
met091<-read.csv("MET091.csv")
met094<-read.csv("MET094.csv")
met096<-read.csv("MET096.csv")

combo<-rbind(met091,met094,met096)
norm<-combo %>% filter(grepl("1HR",condition)) %>%group_by(exp,perm) %>% summarise(t1mean=mean(area_um))
combo<-merge(combo,norm)
combo$nbugs<-combo$area_um/combo$t1mean
combo$log2fc<-log2(combo$fc)
combo_summ<-summarise(group_by(combo, exp,condition), meanlog=mean(log2fc),max_fc = max(log2fc),n = length(log2fc),
                      mean_bugs=mean(nbugs), max_area = max(area_um),perm = unique(perm))
combo_xrep<-summarise(group_by(combo_summ, condition), avg_log=mean(meanlog),sd_log = sd(meanlog),avg_bugs=mean(mean_bugs), sd_bugs = sd(mean_bugs),perm = unique(perm))

#Combine area and log2FC intensity into same plot
pdf(file = "ComboUbAreaLog2intensity.pdf",height = 8,width = 6)
ggplot() + geom_errorbar(data = combo_xrep,aes(x = condition, ymin = avg_log-0.5, ymax = avg_log + sd_log), width = 0.3)+
  geom_bar(data = combo_xrep, aes(x = condition, y = avg_log, fill = perm), stat = "identity", color = "black") +
  geom_dotplot(data = combo_summ, aes(x = condition, y = meanlog),stackdir = "center",binaxis = "y",binwidth = 0.075)+
  geom_errorbar(data = combo_xrep,aes(x = condition, ymin = -avg_bugs-sd_bugs, ymax = -avg_bugs + 0.5), width = 0.3)+
  geom_bar(data = combo_xrep, aes(x = condition, y =-avg_bugs,fill = perm), stat = "identity", color = "black", linetype = "dashed")+
  geom_dotplot(data = combo_summ, aes(x = condition, y = -mean_bugs),stackdir = "center",binaxis = "y",binwidth = 0.075)+
  scale_fill_manual(values = c("deeppink1","turquoise3")) + theme_bw(base_size = 25)+ ylim(-3.5,2.5)+theme(legend.position = "none")+
  labs(x = "Hours post-infection",y=NULL)+
  scale_x_discrete(labels = c("saponin 1HR" = "1", 
                              "saponin 2HR" = "2",
                              "saponin 4HR" = "4",
                              "saponin 6HR" = "6",
                              "saponin 8HR" = "8",
                              "SDS 1HR" = "1", 
                              "SDS 2HR" = "2",
                              "SDS 4HR" = "4",
                              "SDS 6HR" = "6",
                              "SDS 8HR" = "8"))
dev.off()

#stats
bartlett.test(meanlog~condition,data = combo_summ) #p>0.05
model<-aov(meanlog~condition,data = combo_summ)

sink(file = "UbTimecourse_log2fc_ANOVA.txt")
bartlett.test(meanlog~condition,data = combo_summ)
print("Bartlett p>0.05, assume equal variance")
print("Summary - one way ANOVA")
summary(model)
print("Tukey Kramer post-hoc test results, confidence interval 0.95, alpha=0.05")
TukeyHSD(model, conf.level = 0.95)
sink()


# Figure 6E sol insol timecourse DONE -------------------------------------------
setwd("~/Desktop/Paper 2 figures/Submission 2/Figure 6 timecourse and DupA/Sol insol fractionation/AMS analysis/")
data<-read.csv("wtoverdota.csv")
data$condition<-paste(data$Sample,data$Fraction, sep = "_")
data$condition<-factor(data$condition, levels = c("WT1_sol",
                                                 "WT4_sol",
                                                 "WT8_sol",
                                                 "WT1_insol",
                                                 "WT4_insol",
                                                 "WT8_insol"))
data<-data[order(data$condition),]

data_summ<-summarise(group_by(data, condition), meanFC=mean(FC_dotA), sd = sd(FC_dotA), Fraction = unique(Fraction))


pdf(file = "TimecourseSolInsol.pdf",height = 8, width = 8)
ggplot() + geom_errorbar(data = data_summ, aes(y = meanFC, x = condition, ymin = meanFC-0.5, ymax = meanFC + sd), width = 0.3, linewidth = 1)+ 
  geom_bar(data = data_summ, aes(x =condition, y = meanFC, fill = Fraction), stat = "identity",color="black",width = 0.8, linewidth=1) + 
  geom_dotplot(data = data, aes(y = FC_dotA, x = condition), stackdir = "center", binaxis = "y", binwidth = .1) + 
  labs(y = "FC HMW Ub intensity (WT/dotA)", x = "Hours post-infection") + scale_fill_manual(values = c("lightgray","lightblue"))+
  theme_classic(base_size = 35)+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0,3.5))+
  scale_x_discrete(labels = c("WT1_sol" = "1", 
                              "WT4_sol" = "4",
                              "WT8_sol" = "8",
                              "WT1_insol" = "1",
                              "WT4_insol" = "4",
                              "WT8_insol" = "8"))
dev.off()


sink(file = "TimecourseInsolANOVA.txt")
bartlett.test(FC_dotA~condition,data = data[which(data$Fraction=="insol"),])
print("p>0.05, assume equal variance")
model<-aov(FC_dotA~condition,data=data[which(data$Fraction=="insol"),])
summary(model)
TukeyHSD(model)
sink()

sink(file = "TimecourseSolANOVA.txt")
bartlett.test(FC_dotA~condition,data = data[which(data$Fraction=="sol"),])
print("p>0.05, assume equal variance")
model<-aov(FC_dotA~condition,data=data[which(data$Fraction=="sol"),])
summary(model)
print("Not Significant!")
sink()

# Figure 6FG SDS washout DupA OE --------------------------------------------
setwd("~/Desktop/Paper 2 figures/Submission 2/Figure 6 timecourse and DupA/DupA OE FK2 intensity/")
met055<-read.csv("MET055_thresholded.csv")
met066<-read.csv("MET066_thresholded.csv")
met073<-read.csv("MET073_thresholded.csv")

met066<-met066 %>% select(!contains("405"))

alldata<-rbind(met055,met066,met073)
ext<-alldata[which(alldata$bug_pos=="ext"),]
alldata<-alldata[which(alldata$bug_pos=="int"),]

alldata$fc<-(alldata$LegiUbInt/alldata$LegiArea)/alldata$BGUb
alldata$log2fc<-log2(alldata$fc)

alldata$condition<-factor(alldata$condition, levels = c("saponin_EGFP","saponin_DupA","SDS_EGFP","SDS_DupA"))
alldata<-alldata[order(alldata$condition),]

medians<-alldata %>% group_by(experiment,perm) %>% summarise(median = median(log2fc))
alldata<-merge(alldata,medians)
alldata$score<-ifelse(alldata$log2fc>alldata$median, 1, 0)

threshold_summ<-alldata %>% group_by(experiment,condition) %>% summarise(n = length(log2fc), over = sum(score), under = length(which(score==0)),perm = unique(perm), construct = unique(construct))
forplot<-rbind(threshold_summ %>% select(-under) %>% rename(count = over) %>% mutate(category = "Over"),
               threshold_summ %>% select(-over) %>% rename(count = under) %>% mutate(category = "Under"))

forplot$pct<-forplot$count/forplot$n*100
forplot$label<-paste(forplot$category,forplot$condition,sep = "_")

forplot$label<-factor(forplot$label,levels = c("Under_saponin_EGFP","Under_saponin_DupA",
                                                       "Over_saponin_EGFP","Over_saponin_DupA",
                                                       "Under_SDS_EGFP","Under_SDS_DupA",
                                                       "Over_SDS_EGFP","Over_SDS_DupA"))

forplot<-forplot[order(forplot$label),]

forplot_x<-forplot %>% group_by(label) %>% summarise(mean = mean(pct), sd = sd(pct), perm = unique(perm), construct = unique(construct))

pdf("Simplified_bins_saponin_log2FC.pdf",height = 10,width = 7)
ggplot() +geom_errorbar(data = forplot_x[which(forplot_x$perm=="saponin"),],aes(x = label, ymin = mean -sd, ymax = mean + sd), width = 0.3)+
  geom_bar(data = forplot_x[which(forplot_x$perm=="saponin"),], aes(x = label, y = mean, fill = construct), stat = "identity", position = "dodge", color = "black")+
  geom_dotplot(data = forplot[which(forplot$perm=="saponin"),], aes(x = label, y = pct),binaxis = "y",stackdir = "center",binwidth = 2)+
  ylim(0,100) + theme_bw(base_size = 25) + scale_fill_manual(values = c("dodgerblue","grey95")) + labs(x = "Bin", y = "Percent LCVs per bin")
dev.off()

pdf("Simplified_bins_SDS_log2c.pdf",height =10,width = 7)
ggplot() +geom_errorbar(data = forplot_x[which(forplot_x$perm=="SDS"),],aes(x = label, ymin = mean -sd, ymax = mean + sd), width = 0.3)+
  geom_bar(data = forplot_x[which(forplot_x$perm=="SDS"),], aes(x = label, y = mean, fill = construct), stat = "identity", position = "dodge", color = "black")+
  geom_dotplot(data = forplot[which(forplot$perm=="SDS"),], aes(x = label, y = pct),binaxis = "y",stackdir = "center",binwidth = 2)+
  ylim(0,100) + theme_bw(base_size = 25) + scale_fill_manual(values = c("dodgerblue","grey95")) + labs(x = "Bin", y = "Percent LCVs per bin")
dev.off()

#plot paired means

fc_summ<-alldata %>% group_by(experiment,condition) %>% summarise(meanlog = mean(log2fc), n = length(log2fc), label = unique(paste(perm, experiment, sep = "_"))) 

pdf(file = "DupAmeansLog2FC.pdf",height = 8,width = 6)
ggplot(data = fc_summ, aes(x = condition, y = meanlog, group = label, color = experiment))+
  geom_line(linewidth = 2) + geom_point(size = 4) + scale_color_manual(values = c("black","grey","dodgerblue"))+
  theme_bw(base_size = 25) + theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  ylim(0,2.5) + labs(x = NULL, y = "Mean log2(UbLCV/UbBG)")+
  scale_x_discrete(labels = c("saponin_EGFP"="EGFP",
                              "saponin_DupA"="DupA",
                              "SDS_EGFP"="EGFP",
                              "SDS_DupA"="DupA"))
dev.off()

sink(file = "DupAWashoutMeans_ANOVA.txt")
bartlett.test(meanlog~condition,data = fc_summ)
print("Bartlett p>0.05, assume equal variance")
print("Summary - one way ANOVA")
model<-aov(meanlog~condition,data = fc_summ)
summary(model)
print("Not significant!")
sink()

# Figure 7BCD live imaging Ub only update colors and line weights---------------------------------------------------

setwd("~/Desktop/Paper 2 figures/Submission 2/Figure 7 live imaging/GFP Ub only/Area/")
live<-read.csv("all_manual_ub_live.csv")
live$bug<-as.character(live$bug)
live$t<-live$t-1
live$minUbpos<-live$t*7.5
live[which(live$experiment=="MET081"),"rep"]<-"R3"
live[which(live$experiment=="MET079"),"rep"]<-"R2"
live[which(live$experiment=="MET077"),"rep"]<-"R3"
live[which(live$experiment=="MET076"),"rep"]<-"R1"

norm<-live[which(live$t==0),c("bug","area.UB")]
norm<-rename(norm,"area_t0" = "area.UB")

live<-merge(live,norm)

live$ratiot0<-live$area.UB/live$area_t0
# 
# pdf(file = "WTUbTraces.pdf", height = 6,width = 4)
# ggplot(live[which(live$strain=="WT"),],aes(minUbpos,ratiot0,color = bug)) + geom_line() +
#   ylim(0,12) + theme_classic(base_size = 25) +theme(legend.position = "none")+
#   labs(x = "Time Ub positive (min)", y = "Ub area t0/t4") + geom_hline(yintercept = 1, linewidth = 2)
# dev.off()
# 
# pdf(file = "KOUbTraces.pdf", height = 6,width = 4)
# ggplot(live[which(live$strain=="KO"),],aes(minUbpos,ratiot0,color = bug)) + geom_line() +
#   ylim(0,12) + theme_classic(base_size = 25) +theme(legend.position = "none")+
#   labs(x = "Time Ub positive (min)", y = "Ub/LCV area ratio") + geom_hline(yintercept = 1, linewidth = 2)
# dev.off()

#Try plotting all timepoints

live$bin<-paste(live$strain,live$t,sep = "_")

live$bin<-factor(live$bin,levels = c("WT_0", "KO_0",
                                          "WT_1", "KO_1",
                                           "WT_2", "KO_2",
                                           "WT_3", "KO_3",
                                           "WT_4", "KO_4"))
live<-live[order(live$bin),]

bin_summ<-summarise(group_by(live,rep,bin), mean = mean(ratio), n = length(ratio), strain = unique(strain), t = unique(minUbpos))

pdf(file = "TimecourseAreaViolin.pdf",height = 6,width = 10)
ggplot() + geom_violin(data = live, aes(x = bin, y = ratio, fill = strain), width = 1.4) + 
  geom_dotplot(data = bin_summ, aes(x = bin, y = mean), binaxis = "y", stackdir = "center",binwidth = 0.2)+
  scale_fill_manual(values = c("grey90","coral"))+
  geom_hline(yintercept = 1, linetype = "dashed") + ylim(0,16) + theme_classic(base_size = 25) +
  labs(x = NULL, y = "Ub area/L.p. area")
dev.off()

model<-aov(mean~bin, bin_summ)

sink(file = "LiveUbArea_ANOVA.txt")
bartlett.test(mean~bin, bin_summ)
print("Bartlett p>0.05, assume equal variance")
print("Summary - one way ANOVA")
summary(model)
print("Tukey Kramer post-hoc test results, confidence interval 0.95, alpha=0.05")
TukeyHSD(model, conf.level = 0.95)
sink()

tfinal<-live[which(live$t==4),]
ggplot(tfinal,aes(x=strain,y=ratiot0)) + geom_jitter() + geom_hline(yintercept = c(1,2))

#ok make bins - just do over or under 1.5

tfinal$bin<-ifelse(tfinal$ratiot0<1.5,"bin3","bin1")

tfinal_summ<-summarise(group_by(tfinal,strain,bin),count = length(ratiot0))
N<-summarise(group_by(tfinal,strain),N = length(ratiot0))
tfinal_summ<-merge(tfinal_summ,N)
tfinal_summ$pct<-tfinal_summ$count/tfinal_summ$N*100

tfinal_summ$strain<-factor(tfinal_summ$strain,levels = c("WT","KO"))
tfinal_summ<-tfinal_summ[order(tfinal_summ$strain),]

pdf("Live_Bins.pdf",width = 5,height = 6)
ggplot(data = tfinal_summ,aes(x = strain,y=pct,fill=bin))+geom_bar(stat = "identity", color = "black")+
  scale_fill_manual(values = c("coral","grey95")) + theme_bw(base_size = 30) + labs(y = "% LCVs per category", x = "Strain")
dev.off()

#ok plot individual example bugs shown for inset - bug 23 for WT and 69 for KO
examples<-filter(live,bug ==23|bug==69)

pdf(file = "ExampleTraces.pdf",height = 6,width = 7)
ggplot(examples,aes(x = minUbpos,y=ratio,color = strain)) + geom_point(size = 4)+geom_line(linewidth = 2)+
  theme_bw(base_size = 30) + geom_hline(yintercept = 1, linetype = "dashed", linewidth = 1)+
  labs(y = "Ub area/L.p. area",x = "Time Ub positive (min)") + scale_y_continuous(breaks = seq(0,6,1))+
  scale_color_manual(values = c("black","coral"))
dev.off()




# Figure 7FI Rab5 Ub dual live intensity----------------------------------------------
#Import and format data
setwd("/Volumes/AMS_ex/Paper 2 figures/Submission 2/Figure 7 live imaging/Double positive/All manual/")

experiments<-c("MET069","MET070","MET078", "MET080", "MET081","MET089")

mch_lcv<-read.csv(paste(experiments[1], "_mCherry_foreground.csv", sep = ""))
egfp_lcv<-read.csv(paste(experiments[1], "_EGFP_foreground.csv", sep = ""))
mch_bg<-read.csv((paste(experiments[1], "_mCherry_background.csv", sep = "")))
egfp_bg<-read.csv((paste(experiments[1], "_mCherry_background.csv", sep = "")))

for (i in 2:length(experiments)){
  a<-read.csv(paste(experiments[i], "_mCherry_foreground.csv", sep = ""))
  mch_lcv<-rbind(mch_lcv,a)
  b<-read.csv(paste(experiments[i], "_EGFP_foreground.csv", sep = ""))
  egfp_lcv<-rbind(egfp_lcv,b)
  c<-read.csv(paste(experiments[i], "_mCherry_background.csv", sep = ""))
  mch_bg<-rbind(mch_bg,c)
  d<-read.csv(paste(experiments[i], "_EGFP_background.csv", sep = ""))
  egfp_bg<-rbind(egfp_bg,d)
}

mch_lcv$marker<-"Rab5"
mch_bg$marker<-"Rab5"
egfp_lcv$marker<-"Ub"
egfp_bg$marker<-"Ub"

lcv<-rbind(mch_lcv,egfp_lcv)
bg<-rbind(mch_bg,egfp_bg)

lcv<-rename(lcv, lp_area = Area)
lcv<-rename(lcv, lp_mean = Mean)
lcv<-rename(lcv, lp_intden=IntDen)

lcv<-subset(lcv,select = -c(X,RawIntDen))

bg<-rename(bg, bg_area = Area)
bg<-rename(bg, bg_mean = Mean)

bg<-subset(bg, select = -c(X,RawIntDen))

alldata<-merge(lcv,bg)

alldata$experiment<-substring(alldata$Label, regexpr("MET",alldata$Label), regexpr("MET", alldata$Label) + 5)
alldata$strain<-substring(alldata$Label,regexpr("LP|KO|LEG",alldata$Label), regexpr("LP|KO|LEG", alldata$Label) + 1)
alldata$strain<-gsub("LP","WT",alldata$strain)
alldata$strain<-gsub("LE","KO",alldata$strain)
alldata$timepoint<-substring(alldata$Label, regexpr("t:",alldata$Label)+2, regexpr("t:",alldata$Label)+3)
alldata$timepoint<-gsub("/","",alldata$timepoint)
alldata$timepoint<-as.numeric(alldata$timepoint)
alldata$series<-paste(alldata$experiment,alldata$ID, sep = "_")

alldata[which(alldata$experiment=="MET070"|alldata$experiment=="MET081"),"rep"]<-"Rep 1"
alldata[which(alldata$experiment=="MET080"),"rep"]<-"Rep 2"
alldata[which(alldata$experiment=="MET078"|alldata$experiment=="MET089"|alldata$experiment=="MET069"),"rep"]<-"Rep 3"

#Calculate fold change mean LCV region over mean BG
alldata$fc<-alldata$lp_mean/alldata$bg_mean

#Align start times in case I want to plot all traces, find max FC value
transformation<-summarise(group_by(alldata,series, marker), mint = min(timepoint), max_fc = max(fc))
alldata<-merge(alldata,transformation)
alldata$timepoint<-alldata$timepoint-alldata$mint-1 #want to start at -5

#Some of these timecourses are a little short - throw out any that are 10 or less.
trim<-summarise(group_by(alldata,series), maxt=length(timepoint))
trim[which(trim$maxt<=10),"drop"]<-"DROP"
alldata<-merge(alldata,trim)
alldata<-alldata[is.na(alldata$drop),]

#When does the max FC for each marker occur?
max_time_r5<-alldata[which(alldata$fc==alldata$max_fc & alldata$marker=="Rab5"),c("series","strain","timepoint", "rep", "experiment")]
max_time_ub<-alldata[which(alldata$fc==alldata$max_fc & alldata$marker=="Ub"),c("series","strain","timepoint", "rep", "experiment")]
max_time_r5<-rename(max_time_r5,r5tmax = timepoint)
max_time_ub<-rename(max_time_ub, ubtmax = timepoint)
max_time<-merge(max_time_r5,max_time_ub)
max_time$strain<-factor(max_time$strain, levels = c("WT","KO"))
max_time<-max_time[order(max_time$strain),]
#want in minutes
max_time$diff_min<-(max_time$r5tmax-max_time$ubtmax)*5
max_summ<-summarise(group_by(max_time,rep,strain), mean_diff = mean(diff_min), n = length(diff_min))
max_xrep<-summarise(group_by(max_time,strain), mean_diff = mean(diff_min))

pdf(file = "TmaxDiff.pdf",width = 5,height = 7)
ggplot() + geom_jitter(data = max_time, aes(x = strain, y = diff_min, fill = rep),pch=21, color = "black",width = 0.3, size = 3)+
  geom_jitter(data = max_summ, aes(x = strain, y = mean_diff, fill = rep), pch = 21, color = "black", width = 0.3, size = 8)+
  geom_errorbar(data=max_xrep, aes(x = strain, ymin = mean_diff,ymax = mean_diff), width = 0.6, linewidth = 2)+
  geom_hline(yintercept = 0, linetype = "dashed") + scale_fill_manual(values = c("gray30","gray","gray95"))+
  theme_bw(base_size = 25) + theme(legend.position = "none") + ylim(-100, 80) + labs(x = NULL, y = "TmaxRab5-TmaxUb (min)")
dev.off()

#Do stats on tmax difference
sink(file = "TmaxDiffTtest.txt")
t.test(mean_diff~strain,data = max_summ)
sink()

#Considered plotting all traces but it looks a mess
WT<-alldata[which(alldata$strain=="WT"),]
KO<-alldata[which(alldata$strain=="KO"),]

ggplot(WT[which(WT$marker=="Rab5"),],aes(x = timepoint, y = fc, color = series)) + geom_line() + theme(legend.position = "none")
ggplot(KO[which(KO$marker=="Rab5"),],aes(x = timepoint, y = fc, color = series)) + geom_line() #+ theme(legend.position = "none")

#Rep image/movie intensity plots

ko_plot<-alldata[which(alldata$series=="MET089_C"),]
ko_plot$tmin<-(ko_plot$timepoint)*5

wt_plot<-alldata[which(alldata$series=="MET080_I"),]
wt_plot$tmin<-(wt_plot$timepoint)*5

#Trim to equal length - WT series is shorter
ko_plot<-ko_plot[which(ko_plot$tmin<=max(wt_plot$tmin)),]

pdf(file = "KOMovie.pdf",width = 6, height = 6)
ggplot(ko_plot, aes(x=tmin,y = fc, color = marker)) + geom_line(linewidth = 1.5)+
  geom_point(size = 3)+
  theme_bw(base_size = 25) + geom_hline(yintercept = 1, linetype = "dashed")+
  scale_color_manual(values = c("goldenrod","deepskyblue")) + ylim(0.95,1.6)+
  labs(y = "Norm. intensity (LCV/BG)", x = "Time (minutes)")
dev.off()

pdf(file = "WTMovie.pdf",width = 6, height = 6)
ggplot(wt_plot, aes(x=tmin,y = fc, color = marker)) + geom_line(linewidth=1.5)+
  geom_point(size = 3)+
  theme_bw(base_size = 25) + geom_hline(yintercept = 1, linetype = "dashed")+
  scale_color_manual(values = c("goldenrod","deepskyblue")) +ylim(0.95,1.6)+
  labs(y = "Norm. intensity (LCV/BG)", x = "Time (minutes)")
dev.off()

# Figure 7G Rab5 Ub dual percentages ---------------------------------------
setwd("/Volumes/AMS_ex/Paper 2 figures/Submission 2/Figure 7 live imaging/Double positive/All manual/")

counts<-read.csv("manual_counts_live.csv")
counts$experiment<-substring(counts$Filename, regexpr("MET",counts$Filename), regexpr("MET",counts$Filename) + 5)
counts$strain<-substring(counts$Filename, regexpr("LP|LE|KO",counts$Filename), regexpr("LP|LE|KO",counts$Filename)+1)
counts$strain<-gsub("LP","WT",counts$strain)
counts$strain<-gsub("LE","KO",counts$strain)
counts$series<-c(1:length(counts$Filename))
counts$series<-paste(counts$experiment, counts$series,sep = "_")

count_summ<-summarise(group_by(counts,experiment,strain),total = sum(N_LCVs),RNUN = sum(RNUN), RPUN = sum(RPUN), RNUP = sum(RNUP), RPUP = sum(RPUP))
count_comb<-summarise(group_by(counts,strain),total = sum(N_LCVs),RNUN = sum(RNUN), RPUN = sum(RPUN), RNUP = sum(RNUP), RPUP = sum(RPUP))

category<-c("RNUN","RPUN","RNUP","RPUP")

pct_plot<-count_comb[,c("strain","total",category[1])]
pct_plot<-rename(pct_plot, count = sym(category[1]))
pct_plot$category<-category[1]

for (i in 2:length(category)) {
  x<-count_comb[,c("strain","total",category[i])]
  x<-rename(x, count = sym(category[i]))
  x$category<-category[i]
  pct_plot<-rbind(pct_plot,x)
}

pct_plot$pct<-pct_plot$count/pct_plot$total*100
pct_plot<-pct_plot[which(pct_plot$category!="RNUN"),]

pct_plot$strain<-factor(pct_plot$strain, levels = c("WT","KO"))
pct_plot<-pct_plot[order(pct_plot$strain),]

pct_plot$category<-factor(pct_plot$category, levels = c("RPUN","RNUP","RPUP"))
pct_plot<-pct_plot[order(pct_plot$category),]

pdf("ManualPercent.pdf",height = 8,width = 6)
ggplot(pct_plot, aes(x = category, y = pct, fill = strain)) + geom_bar(stat = "identity", position = "dodge", color = "black")+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0,35)) + theme_classic(base_size = 25) + 
  labs(y = "Percent LCVs per category", x = "Category") + scale_fill_manual(values = c("coral","gray95"))
dev.off()



# Figure 7H Rab5 dwell approx ----------------------------------------------
setwd("/Volumes/AMS_ex/Paper 2 figures/Submission 2/Figure 7 live imaging/Double positive/All manual/")
alldata<-read.csv("manual_rab5time.csv")

alldata$experiment<-substring(alldata$Label, regexpr("MET",alldata$Label), regexpr("MET", alldata$Label) + 5)
alldata$strain<-substr(alldata$Label,regexpr("LP|KO|LEG",alldata$Label), regexpr("LP|KO|LEG", alldata$Label) + 1)
alldata$strain<-gsub("LP","WT",alldata$strain)
alldata$strain<-gsub("LE","KO",alldata$strain)

alldata$dwell<-ifelse(alldata$t_rab5n=="X",5*(alldata$ttotal - alldata$t_rab5p), 5*(as.numeric(alldata$t_rab5n) - alldata$t_rab5p))
#idk this code made her mad but it works so ¿?

alldata$category<-ifelse(alldata$t_rab5n=="X","pos only","pos to neg")
alldata$strain<-factor(alldata$strain, levels = c("WT","KO"))

alldata %>% group_by(strain) %>% summarise(n=length(strain))

pdf(file = "DwellApprox_largerpoints.pdf",height = 8,width = 8)
ggplot(alldata, aes(x = strain, y = dwell, fill = category)) + geom_jitter(pch=21, size = 6, width = 0.3)+
  theme_bw(base_size = 25) + labs(y = "Rab5 LCV association time (min)", x = "Strain") +
  scale_fill_manual(values = c("grey95","purple"))
dev.off()



# Figure S2B dual stain Bap31 ----------------------------------------------
setwd("~/Desktop/Paper 2 figures/Figure S2 dual stain/Bap31 data/")
met046<-read.csv("MET046_Bap31_prepost_croppedImage.csv")
met059<-read.csv("20240912_MET059_Bap31_prepost_croppedImage.csv")

met046$experiment<-"MET046"
met059$experiment<-"MET059"

all_b31<-rbind(met046, met059)

all_b31$condition<-substring(all_b31$FileName_DAPI, regexpr("MET|OPT", all_b31$FileName_DAPI) + 7,  regexpr("MET|OPT", all_b31$FileName_DAPI) + 9)
all_b31$condition<-gsub("abf","AbFirst", all_b31$condition)
all_b31$condition<-gsub("CHA","AbFirst", all_b31$condition)
all_b31$condition<-gsub("per","PermFirst", all_b31$condition)
all_b31$condition<-gsub("CAB","PermFirst", all_b31$condition)


all_b31$condition<-factor(all_b31$condition, levels = c("PermFirst", "AbFirst"))
all_b31<-all_b31[order(all_b31$condition),]


b31_summ<-summarise(group_by(all_b31, experiment,condition), mean_pearson = mean(Correlation_Correlation_DAPI_FarRed), n = length(Correlation_Correlation_DAPI_FarRed))
b31_xrep<-summarise(group_by(all_b31, condition), mean_pearson = mean(Correlation_Correlation_DAPI_FarRed))

pdf(file = "Bap31_prepost.pdf",width = 6, height = 8)
ggplot() + geom_jitter(data = all_b31, aes(x = condition, y = Correlation_Correlation_DAPI_FarRed, fill = experiment), pch = 21, color = "black", width = 0.1, size = 4) + 
  geom_jitter(data = b31_summ, aes(x = condition, y = mean_pearson, fill = experiment), pch = 21, color = "black", width = 0.1, size = 8) +
  scale_fill_manual(values = c("white","dodgerblue"))+
  geom_errorbar(data = b31_xrep,aes(x = condition, ymax = mean_pearson, ymin = mean_pearson, width = 0.4, linewidth = 1))+
  theme_bw(base_size = 30) + theme(legend.position = "none")+ ylim(-0.1,1.1)+
  labs(y = "Pearson Correlation 405, 633", x = "GAR-405 staining step")
dev.off()

sink(file = "Bap31Pearson_ttest.txt")
t.test(mean_pearson~condition, data = b31_summ)
sink



# Figure S2DEFG dual stain Rab5 intensity- all WT Brodsky UPDATE log2fc-------------------------
#cropped either 405+ (ext) or 405- (int) LCVs

setwd("~/Desktop/Paper 2 figures/Submission 2/Figure S2 dual stain/intensity data/")
image<-read.csv("20241210_Rab5_intext_no405bgsImage.csv")
legi<-read.csv("20241210_Rab5_intext_no405bgsmergedlegi.csv")

filenames<-data.frame("ImageNumber" = image$ImageNumber, "filename" = image$FileName_GAR405, "BG_Rab5" = image$Intensity_MeanIntensity_MaskedBG, "BG_405" = image$Intensity_MeanIntensity_Masked405)
filenames$experiment<-substring(filenames$filename, regexpr("MET", filenames$filename), regexpr("MET", filenames$filename) +5)
filenames$bug_pos<-substring(filenames$filename, regexpr("int|ext", filenames$filename), regexpr("int|ext", filenames$filename) +2)

legi_min<-data.frame("ImageNumber" = legi$ImageNumber, "Rab5Integrated" = legi$Intensity_IntegratedIntensity_AF488, "GAR405Integrated"=legi$Intensity_IntegratedIntensity_GAR405, 
                     "legiArea" = legi$AreaShape_Area, "MeanRab5" = legi$Intensity_MeanIntensity_AF488, "Mean405" = legi$Intensity_MeanIntensity_GAR405)

alldata<-merge(filenames,legi_min, all = T)
alldata<-alldata[!is.na(alldata$legiArea),]

alldata$fcrab5<-log2(alldata$MeanRab5/alldata$BG_Rab5)
alldata$fc405<-log2(alldata$Mean405/alldata$BG_405)

thresholds<-summarise(group_by(alldata,experiment), q1 = quantile(fcrab5,probs = 0.25), q2 = median(fcrab5),q3 = quantile(fcrab5, probs = 0.75))
alldata<-merge(alldata,thresholds)

for (i in 1:length(alldata$fcrab5)){
  if (alldata$fcrab5[i]<alldata$q1[i]){
    alldata$quart[i]<-"Q1"
  }else if(alldata$fcrab5[i]<alldata$q2[i]){
    alldata$quart[i]<-"Q2"
  }else if(alldata$fcrab5[i]<alldata$q3[i]){
    alldata$quart[i]<-"Q3"
  }else{
    alldata$quart[i]<-"Q4"
  }
}

alldata_summ<-summarise(group_by(alldata,experiment,bug_pos), avg_r5=mean(fcrab5), avg_405 = mean(fc405), n = length(fc405))
alldata_xrep<-summarise(group_by(alldata,bug_pos), avg_r5=mean(fcrab5), avg_405 = mean(fc405))
count<-sum(alldata_summ$n)

totals<-summarise(group_by(alldata,experiment,bug_pos,quart), totals = length(fcrab5))
N<-summarise(group_by(alldata, experiment, bug_pos), N = length(fcrab5))
totals<-merge(totals, N)
totals$pct<-totals$totals/totals$N*100
totals$condition<-paste(totals$quart,totals$bug_pos,sep = "_")
totals_summ<-summarise(group_by(totals,condition), mean = mean(pct), sd = sd(pct), bug_pos = unique(bug_pos))

pdf(file = "QuartileFractionlog2FC.pdf", height = 8, width = 10)
ggplot() + geom_errorbar(data = totals_summ, aes(y = mean, x = condition, ymin = mean-1, ymax = mean + sd), width = 0.3, linewidth = 1)+ 
  geom_bar(data = totals_summ, aes(x =condition, y = mean, fill = bug_pos), stat = "identity",color="black",width = 0.8, linewidth=1) + 
  geom_dotplot(data = totals, aes(y = pct, x = condition), stackdir = "center", binaxis = "y", binwidth = 2) + 
  scale_fill_manual(values = c("steelblue","goldenrod"))+
  labs(y = "Percent per Rab5 intensity quartile", x = "Quartile") + 
  theme_classic(base_size = 35)+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(-2,75))
dev.off()

pdf("Rab5log2FCdual.pdf",height = 8, width=6)
ggplot() + geom_jitter(data = alldata, aes(x = bug_pos, y = fcrab5, fill = experiment), pch = 21, color = "black", width = 0.3, size = 3)+
  geom_jitter(data = alldata_summ, aes(x = bug_pos, y = avg_r5, fill = experiment), pch = 21, color = "black", width = 0.3, size = 8)+
  scale_fill_manual(values = c("blueviolet","white","magenta","grey80")) + geom_errorbar(data = alldata_xrep,aes(x = bug_pos, ymax = avg_r5, ymin = avg_r5),width = 0.7, linewidth = 2)+
  theme_bw(base_size = 30) + theme(legend.position = "none")+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  labs( y = "log2(Rab5Lp/Rab5BG)", x= NULL) + scale_x_discrete(labels = c("405+","405-")) + ylim(-0.3,1)
dev.off()

pdf("405log2FCdual.pdf",height = 8, width=6)
ggplot() + geom_jitter(data = alldata, aes(x = bug_pos, y = fc405, fill = experiment), pch = 21, color = "black", width = 0.3, size = 3)+
  geom_jitter(data = alldata_summ, aes(x = bug_pos, y = avg_405, fill = experiment), pch = 21, color = "black", width = 0.3, size = 8)+
  scale_fill_manual(values = c("blueviolet","white","magenta","grey80")) + geom_errorbar(data = alldata_xrep,aes(x = bug_pos, ymax = avg_405, ymin = avg_405),width = 0.7, linewidth = 2)+
  theme_bw(base_size = 30) + theme(legend.position = "none")+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  labs( y = "log2(405Lp/405BG)", x = NULL) + ylim(-0.1, 1.8) + scale_x_discrete(labels = c("405+","405-"))
dev.off()

sd_summ<-summarise(group_by(alldata_summ,bug_pos), sd_r5 = sd(avg_r5), sd_405 = sd(avg_405))
#ok sd is more than 2X in one group compared to other --> Welch's t test is appropriate

sink(file = "Rab5log2FCdualstain_ttest.txt")
t.test(avg_r5~bug_pos, data = alldata_summ)
sink

sink(file = "405log2FCdualstain_ttest.txt")
t.test(avg_405~bug_pos, data = alldata_summ)
sink

#Plot 405 vs Rab5 intensity

pdf(file = "log2Rab5405scatter.pdf",width = 6,height = 4)
ggplot(alldata,aes(x = fc405,y = fcrab5, fill = experiment)) + geom_point(pch=21,color = "black")+
  scale_fill_manual(values = c("black","lightgrey","blue","magenta"))+
  theme_bw(base_size = 25) + theme(legend.position = "NONE")+ labs(x = "log2(405Lp/405BG) ", y = "log2(Rab5Lp/Rab5BG)")
dev.off()


  