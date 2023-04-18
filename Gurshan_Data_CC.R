library(tidyverse)

data <- read.table(file="data/IrisLines_DLA_7296h2.csv", sep=",")

lineinfo = read.csv(file = "data/Iris_Crispr_Modeling.csv", header=TRUE)
lineinfo$Genotype <- as.factor(lineinfo$Genotype)

#### Adding the line info (gene coded as 0 or 1)

merge_data <- merge(data, lineinfo, by="Genotype", all.x = FALSE)

levels(as.factor(merge_data$Genotype))
levels(as.factor(lineinfo$Genotype))
levels(as.factor(data$Genotype))



##############################
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, draw_group = function(self, data, ..., draw_quantiles = NULL){
  data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
  grp <- data[1,'group']
  newdata <- plyr::arrange(transform(data, x = if(grp%%2==1) xminv else xmaxv), if(grp%%2==1) y else -y)
  newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
  newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x']) 
  if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
    stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 
                                              1))
    quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
    aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
    aesthetics$alpha <- rep(1, nrow(quantiles))
    both <- cbind(quantiles, aesthetics)
    quantile_grob <- GeomPath$draw_panel(both, ...)
    ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
  }
  else {
    ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
  }
})

geom_split_violin <- function (mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, position = position, show.legend = show.legend, inherit.aes = inherit.aes, params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

####################




data$Time <- as.factor(data$Time)

ggplot(merge_data_nocyp_nopad, aes(x=Genotype, y=Lesion.Size, fill=Time))+
  #geom_point()+
  #geom_line(size=1, aes(group= Isolate, color=Botrydial_HL)) +
  geom_split_violin()+
  #geom_boxplot(width=0.15)+ 
  theme_bw()


##### Filtering the positive controls

#data_nocyp_nopad <- data %>% 
#  filter(Genotype!="pad3", Genotype!="Cyp79B2/B3") 
#levels(as.factor(data_nocyp_nopad$Genotype))

merge_data_nocyp_nopad <- merge_data %>% 
  filter(Genotype!="pad3", Genotype!="Cyp79B2/B3") 
levels(as.factor(merge_data_nocyp_nopad$Genotype))

merge_data_nocyp_nopad$Genotime <- paste(merge_data_nocyp_nopad$Genotype, merge_data_nocyp_nopad$Time, sep="_")  

regression0 = lm(Leaf.Size ~ Genotype*Isolate + Tray + Time , data=merge_data_nocyp_nopad)
LM_Data<- as.data.frame(anova(regression0))
colnames(LM_Data) <- c("Df", "Sum_Sq","Mean_Sq", "F_value", "Pval")
Temp_Sum <- sum(LM_Data$Sum_Sq)
LM_Data$PercVar <- (LM_Data$Sum_Sq*100)/Temp_Sum
LM_Data$Labels <- paste(rownames(LM_Data), round(LM_Data$PercVar, digits=1), "%", sep="-" )

pdf("Pie_chart_LM.pdf")
pie(LM_Data$PercVar, labels=LM_Data$Labels, col=rainbow(length(LM_Data$Labels)))
dev.off()


Split <- split.data.frame(merge_data_nocyp_nopad, merge_data_nocyp_nopad$Time)
names(Split) <- c("NoCTRL72", "NoCTRL96")
list2env(Split, globalenv())

regression1 = lm(Leaf.Size ~ Genotype*Isolate + Tray, data=NoCTRL72)
anova(regression1)

regression2 = lm(Leaf.Size ~ Genotype*Isolate + Tray, data=NoCTRL96)
anova(regression2)

plot(merge_data_nocyp_nopad$Tray, merge_data_nocyp_nopad$Image)

regression_Gene = lm(Area.mm ~ Tray + Time + Isolate*AT1G22410_myb + Isolate*AT1G22410_evening + Isolate*AT4G33510, data=merge_data_nocyp_nopad)
RG <- as.data.frame(anova(regression_Gene))
colnames(RG) <- c("Df", "Sum_Sq","Mean_Sq", "F_value", "Pval")
RG_Sum <- sum(RG$Sum_Sq)
RG$PercVar <- (RG$Sum_Sq*100)/RG_Sum
RG$Labels <- paste(rownames(RG), round(RG$PercVar, digits=1), "%", sep="-" )

pdf("Pie_chart_LM_Gene.pdf")
pie(RG$PercVar, labels=RG$Labels, col=rainbow(length(RG$Labels)))
dev.off()


regression_Iso = lm(Area.mm ~ Tray + Time + Genotype*Botrydial_HL, data=merge_data_nocyp_nopad)
RI <- as.data.frame(anova(regression_Iso))
colnames(RI) <- c("Df", "Sum_Sq","Mean_Sq", "F_value", "Pval")
RI_Sum <- sum(RI$Sum_Sq)
RI$PercVar <- (RI$Sum_Sq*100)/RI_Sum
RI$Labels <- paste(rownames(RI), round(RI$PercVar, digits=1), "%", sep="-" )

pdf("Pie_chart_LM_Isolate.pdf")
pie(RI$PercVar, labels=RI$Labels, col=rainbow(length(RI$Labels)))
dev.off()


################################

data$key <- paste(data$Isolate, data$Genotype, data$Time, sep = "_")
Summary <- summarise(group_by(data, key),
                     mean_Area = mean(Area.mm))

Summary <- separate(Summary, col=key, into=c("Isolate", "Genotype", "Time"), sep = "_")

Summary_wide <- Summary %>% 
  pivot_wider(names_from = Genotype, values_from = mean_Area) %>% 
  rename("1-14-23" = "1.14.23", "1-14-1"= "1.14.1", "1-33-1" = "1.33.1", "5-10-2" = "5.10.2")

Summary_wide_diff <- Summary_wide %>% 
  mutate( "1-14-1-WT" = 1-14-1 - col0,
          "1-14-23-WT" = 1-14-23 - col0,
          "1-33-1-WT" = 1-33-1 - col0,
          "3-1-2-16-WT" = 3-1-2-16 - col0,
          "3-1-8-14-WT" =3-1-8-14- col0,
          "3-1-8-16-WT" = 3-1-8-16 - col0,
          "3-11-2-11-WT" =1-14-1 - col0,
          "5-10-2-WT" = 5-10-2 - col0,
          "CYP79B2/B3-WT" = get("Cyp79B2/B3") - col0,
          "pad3-WT" = pad3 - col0,
  )

Summary_wide_means <- Summary_wide_diff %>%
  select(-Isolate, -Time) %>% 
  colMeans()
Summary_wide_means <- as.data.frame(Summary_wide_means)



Diff_long <- Summary_wide_diff %>% 
  pivot_longer(
    cols = "1-14-1":"pad3-WT",
    values_to = "Lesion_Area"
  )
Diff_long %>%
  filter(name=="WT-1-14-1") %>% 
  ggplot(aes(x=Time, y=Lesion_Area, color=Isolate))+
  geom_line(aes(group=Isolate))+
  ggtitle("WT-1-14-1")

Diff_long %>%
  filter(name=="WT-1-14-23") %>% 
  ggplot(aes(x=Time, y=Lesion_Area, color=Isolate))+
  geom_line(aes(group=Isolate))+
  ggtitle("WT-1-14-23")

Diff_long %>%
  filter(name=="WT-1-33-1") %>% 
  ggplot(aes(x=Time, y=Lesion_Area, color=Isolate))+
  geom_line(aes(group=Isolate))+
  ggtitle("WT-1-33-1")

Diff_long %>%
  filter(name=="WT-3-1-2-16") %>% 
  ggplot(aes(x=Time, y=Lesion_Area, color=Isolate))+
  geom_line(aes(group=Isolate))+
  ggtitle("WT-3-1-2-16")

Diff_long %>%
  filter(name=="WT-3-1-8-14") %>% 
  ggplot(aes(x=Time, y=Lesion_Area, color=Isolate))+
  geom_line(aes(group=Isolate))+
  ggtitle("WT-3-1-8-14")

Diff_long %>%
  filter(name=="WT-3-1-8-16") %>% 
  ggplot(aes(x=Time, y=Lesion_Area, color=Isolate))+
  geom_line(aes(group=Isolate))+
  ggtitle("WT-3-1-8-16")

Diff_long %>%
  filter(name=="WT-3-11-2-11") %>% 
  ggplot(aes(x=Time, y=Lesion_Area, color=Isolate))+
  geom_line(aes(group=Isolate))+
  ggtitle("WT-3-11-2-11")

Diff_long %>%
  filter(name=="WT-5-10-2") %>% 
  ggplot(aes(x=Time, y=Lesion_Area, color=Isolate))+
  geom_line(aes(group=Isolate))+
  ggtitle("WT-5-10-2")

Diff_long %>%
  filter(name=="WT-Cyp79B2/B3") %>% 
  ggplot(aes(x=Time, y=Lesion_Area, color=Isolate))+
  geom_line(aes(group=Isolate))+
  ggtitle("WT-CYP79B2/B3")

Diff_long %>%
  filter(name=="WT-pad3") %>% 
  ggplot(aes(x=Time, y=Lesion_Area, color=Isolate))+
  geom_line(aes(group=Isolate))+
  ggtitle("WT-pad3")

Diff_long %>%
  filter(name=="1-14-1") %>% 
  ggplot(aes(x=Time, y=Lesion_Area, color=Isolate))+
  geom_line(aes(group=Isolate))+
  ggtitle("1-14-1")

Diff_long %>%
  filter(name=="1-14-23") %>% 
  ggplot(aes(x=Time, y=Lesion_Area, color=Isolate))+
  geom_line(aes(group=Isolate))+
  ggtitle("1-14-23")

Diff_long %>%
  filter(name=="1-33-1") %>% 
  ggplot(aes(x=Time, y=Lesion_Area, color=Isolate))+
  geom_line(aes(group=Isolate))+
  ggtitle("1-33-1")

Diff_long %>%
  filter(name=="3-1-2-16") %>% 
  ggplot(aes(x=Time, y=Lesion_Area, color=Isolate))+
  geom_line(aes(group=Isolate))+
  ggtitle("3-1-2-16")

Diff_long %>%
  filter(name=="3-1-8-14") %>% 
  ggplot(aes(x=Time, y=Lesion_Area, color=Isolate))+
  geom_line(aes(group=Isolate))+
  ggtitle("3-1-8-14")

Diff_long %>%
  filter(name=="3-1-8-16") %>% 
  ggplot(aes(x=Time, y=Lesion_Area, color=Isolate))+
  geom_line(aes(group=Isolate))+
  ggtitle("3-1-8-16")

Diff_long %>%
  filter(name=="3-11-2-11") %>% 
  ggplot(aes(x=Time, y=Lesion_Area, color=Isolate))+
  geom_line(aes(group=Isolate))+
  ggtitle("3-11-2-11")

Diff_long %>%
  filter(name=="5-10-2") %>% 
  ggplot(aes(x=Time, y=Lesion_Area, color=Isolate))+
  geom_line(aes(group=Isolate))+
  ggtitle("5-10-2")

Diff_long %>%
  filter(name=="CYP79B2/B3") %>% 
  ggplot(aes(x=Time, y=Lesion_Area, color=Isolate))+
  geom_line(aes(group=Isolate))+
  ggtitle("CYP79B2/B3")

Diff_long %>%
  filter(name=="pad3") %>% 
  ggplot(aes(x=Time, y=Lesion_Area, color=Isolate))+
  geom_line(aes(group=Isolate))+
  ggtitle("Pad3")
```
```{r}
Diff_long %>% 
  filter(!str_detect(name, 'WT-')) %>% 
  filter(!str_detect(name, 'CYP|pad3')) %>% 
  ggplot(aes(x=Time, y=Lesion_Area, color=name))+
  geom_line(aes(group=name))+
  facet_grid(.~Isolate)

Diff_long %>% 
  filter(!str_detect(name, '-WT')) %>% 
  #filter(!str_detect(name, 'Cyp|pad3')) %>% 
  ggplot(aes(x=Time, y=Lesion_Area, color=Isolate))+
  geom_boxplot()+
  geom_line(aes(group=Isolate))+
  facet_grid(.~name)
```
```{r}
new_merge <- merge_data_nocyp_nopad %>% 
  group_by(Genotype, Time) %>% 
  mutate(mean.genotype = mean(Lesion.Size))
```
```{r}
new_merge_pivot <- new_merge %>% 
  pivot_wider(names_from = Genotype,
              values_from = mean.genotype)

```
```{r}

Sum_line <- merge_data_nocyp_nopad %>% 
  group_by(Genotime) %>% 
  summarise(genotime = mean(Lesion.Size))

Sum_line <- separate(Sum_line, Genotime, into=c("Genotype", "Time"), sep="_")

sum_line_w <- Sum_line %>% 
  pivot_wider(
    names_from = "Genotype",
    values_from = "genotime"
  )
#write.csv(sum_line_w, file="data/Iris_lines_genotypeavgs.csv")
Mutant_WT_diffs = read.csv(file = "data/Iris_lines_genotypeavgs.csv", header=TRUE)
Mutant_WT_diff_long <-Mutant_WT_diffs %>%
  pivot_longer(
    cols = X1.14.1.WT:X5.10.2.WT,
    names_to = "Genotype",
    values_to = "mean.genotype"
  )

Mutant_WT_diff_long %>% 
  ggplot(aes(x=Time, y=mean.genotype, color=Genotype))+
  geom_line()

regression0 = lm(Leaf.Size ~ Genotype*Isolate + Tray + Time , data=merge_data_nocyp_nopad)
regression_Gene = lm(Area.mm ~ Tray + Time + Isolate*AT1G22410_myb + Isolate*AT1G22410_evening + Isolate*AT4G33510, data=merge_data_nocyp_nopad)
regression_Iso = lm(Area.mm ~ Tray + Time + Genotype*Botrydial_HL, data=merge_data_nocyp_nopad)