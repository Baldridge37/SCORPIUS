# SCORPIUS pseudotime analysis

```

#May need stringi - don't do compilation
library(SCORPIUS)
library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)
source("C:/Documents/SLM_wba_plotting/master_theme_ggplot.R")

setwd("C:/Documents/normalise_workflow/")

phn<- read.table("Combined_annotation2.txt", header=TRUE, sep="\t")
counts<-read.table("new_counts.txt", header = TRUE)
tpms<-read.table("protoplast_tpm.txt", header = TRUE)

#May not need feature data
##feat<-read.table("C:/Documents/SingleCell/gene_feature.txt", header=TRUE, sep="\t")
#fd <- new("AnnotatedDataFrame", data = feat)

#This works because I only have 1 and 20 cells. Wouldn't if I had any other number with a 1 in it
sc_index = grep(x = phn$cell_number, pattern ="1")
tpms<-tpms[,sc_index]
phn<-phn[sc_index,]

expression <-t(tpms)
group_name <-phn$Cell_type

dist <- correlation_distance(expression)

dim(dist)

plot(density(dist))

space <- reduce_dimensionality(dist)

draw_trajectory_plot(space)

draw_trajectory_plot(space, progression_group = group_name)
pre_filt_traj_plot<-draw_trajectory_plot(space, progression_group = group_name)
filt <- outlier_filter(dist)
expression <- expression[filt, ]
group_name <- group_name[filt]
dist <- dist[filt, filt]
space <- reduce_dimensionality(dist)

draw_trajectory_plot(space[, c(1, 2)])
draw_trajectory_plot(space[, c(1, 3)]) + labs(y = "Component 3")
draw_trajectory_plot(space[, c(2, 3)]) + labs(x = "Component 2", y = "Component 3")
traj <- infer_trajectory(space)

draw_trajectory_plot(space, progression_group = group_name, path = traj$path)
Trajectory_plot<-draw_trajectory_plot(space, progression_group = group_name, path = traj$path)

ggsave(filename="Trajectory_plot_prefilt_Scorpius.tiff", plot=pre_filt_traj_plot, dpi = 600, width=20, height=15, units="cm")
ggsave(filename="Trajectory_plot_postfilt_Scorpius.tiff", plot=Trajectory_plot, dpi = 600, width=20, height=15, units="cm")

#####Use the most informative genes to infer the trajectory

gimp <- gene_importances(expression, traj$time, num_permutations = 0, num_threads = 8)
gene_sel <- gimp[1:50,]
expr_sel <- expression[,gene_sel$gene]
#traj <- infer_trajectory(expr_sel)


#Check if necessary first
#traj<- reverse_trajectory(traj)
draw_trajectory_heatmap(expr_sel, traj$time, group_name)

draw_trajectory_plot(space, progression_group = group_name, path = traj$path)
#Is worse than the first one so I used the original trajectory

###Export expression
exprs<-t(expression)
write.csv(exprs, file="Expression.csv")
write.csv(traj, file="Scorpius_trajectory_final.csv")
options(stringsAsFactors = FALSE)

##########################################
############Expression vs time############
##########################################

#Manually added time from exported trajectory
ET<-read.table("Expression_time.txt", header=F, sep="\t")

ET<-t(ET)
colnames(ET)=ET[1,]
ET=ET[-1,]

ET<-as.data.frame(ET)
plot(ET$time, ET$AT1G01280)

ET$time<-as.numeric(as.character(ET$time))




ET_melt<-melt(ET, id.vars = c("Cell", "time"), measure.vars = )

gene<-"AT3G19150" #ICK4

ET_gene<-ET_melt[ET_melt$variable==gene,]


ps<-ggplot(ET_gene, aes(time, as.numeric(value)))

ICK4<-ps+geom_smooth(color="black", method="loess", size=2)+geom_point(color="red", size =1.5)+
  labs(title ="AT3G19150 ICK4", y = "Expression (TPM)", x = "Pseudotime(SCORPIUS)") + 
  theme_bw() + scale_x_continuous()



#exp_groups<-gimp[gimp$importance>0,]
exp_groups<-gimp[1:1000,]

exp_grp_sel <- expression[,exp_groups$gene]
exp_modules <- extract_modules(scale_quantile(exp_grp_sel), traj$time, verbose = F)

write.csv(exp_modules, file="Modules_1kgenes_final.csv")

##############Combined



#######Combined Modules 2 3 and 4 manually
Modules<-read.csv("Modules_1kgenes_edit.csv", header=T)
Module_1<-Modules[Modules$module=="1",]
Module_2<-Modules[Modules$module=="2",]


ET_Mod1<-ET_melt[ET_melt$variable%in%Module_1$feature,]
ET_Mod2<-ET_melt[ET_melt$variable%in%Module_2$feature,]


###Plot mean expression profiles for module one and 2
##Default span=0.75 seems to fit best. could make an argument for using 1 given variability
Fit1<-loess(formula= value~as.numeric(time), data=ET_Mod1, span=0.75)
Smooth1<-predict(Fit1)

M1<-ggplot(ET_Mod1, aes(x=as.numeric(time), y=as.numeric(value)))
Mod1_mean<-M1+geom_line(aes(y=Smooth1), size=2) +
  labs(title = "Module 1 mean expression", y = "Expression (TPM)", x = "Pseudotime(SCORPIUS)") + 
  theme_bw() + scale_x_continuous()

Fit2<-loess(formula= value~as.numeric(time), data=ET_Mod2, span=0.75)
Smooth2<-predict(Fit2)

M2<-ggplot(ET_Mod2, aes(x=as.numeric(time), y=as.numeric(value)))
Mod2_mean<-M2+geom_line(aes(y=Smooth2), size=2) +
  labs(title = "Module 2 mean expression", y = "Expression (TPM)", x = "Pseudotime(SCORPIUS)") + 
  theme_bw() + scale_x_continuous()

Comb_mean_plot<-grid.arrange(Mod1_mean, Mod2_mean)

ggsave(filename="Mean_expression_2_modules.tiff", plot=Comb_mean_plot, dpi=600, width = 15, height = 10, units = "cm")
ggsave(filename="Mean_expression_module1.tiff", plot=Mod1_mean, dpi=600, width=15, height=10, units="cm")
ggsave(filename="Mean_expression_module2.tiff", plot=Mod2_mean, dpi=600, width=15, height=10, units="cm")

####### Plotting marker genes from published data
#########################################################################


markers<-read.table("markers.txt", header=T, sep="\t")
mrk<-markers$gene
mrk<-as.character(mrk)
mark_sel <-subset(expression, select=mrk)

marker_heatmap<-draw_trajectory_heatmap(mark_sel, traj$time, group_name, show_labels_row = T)



############### Diagnostics

Pheno<-ET[,c(1,2)]
pass_cells<-phn[phn$sample_id%in%Pheno$Cell,]

Pheno_time<-cbind(Pheno, pass_cells)

gt<-ggplot(Pheno_time, aes(x=time, y=feature_number))
gene_time_plot<-gt+geom_smooth(color="black", method="loess", size=2, se= F)+geom_point(color="blue", size =1.5)+
  labs(title = "Change in genes expressed over pseudotime", y = "Genes expressed", x = "Pseudotime(SCORPIUS)") + 
  theme_bw() + scale_x_continuous()

ggsave(filename="gene_number_over_time.tiff", plot=gene_time_plot, dpi=600, width = 20, height = 15, units = "cm")

gr<-ggplot(Pheno_time, aes(x=time, y=read_count))
gene_read_plot<-gr+geom_smooth(color="black", method="loess", size=2, se= F)+geom_point(color="blue", size =1.5)+
  labs(title = "Read count variation with pseudotime", y = "Read count", x = "Pseudotime(SCORPIUS)") + 
  theme_bw() + scale_x_continuous()

ggsave(filename="read_count_over_time.tiff", plot=gene_read_plot, dpi=600, width = 20, height = 15, units = "cm")


```
