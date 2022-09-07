###DO YOU HAVE RRtools installed? 
###if not, unhash the bottom two lines,change the directory in second line and run both lines:
#library(devtools)
#install("C:/Users/yaps/Desktop/RRtools")

library(RRtools) 

load("C:/Users/yaps/Desktop/rrspecies/SYsaved_scripts.R")

###the 3 lines below are very important, as R look for data through the info you provide for these objects:
###RandRbase, Species, dataset

maindir <- "C:/Users/yaps/Desktop/"
setwd(maindir)

RandRbase <- "rrspecies/" #main directory 
species <- "PersHirs" #species name
dataset <- "DPers20-5061" #dart order

###CHECK your genotype file, note number of columns (nmetavar) and rows (topskip) until your genotype data starts , i.e,. the 0s,1s and 2s...

topskip   <- 6 
nmetavar  <- 18

d1        <- read.dart.xls.onerow(RandRbase,species,dataset,topskip, nmetavar, euchits=FALSE)
###this line will read like this:
###Reading data file: D:/test/dart_raw/Report_Dtest19-0000_SNP_singlerow_2.csv  
qc1       <- report.dart.qc.stats(d1, RandRbase, species, dataset, threshold_missing_loci = 0.8) # threshold_missing_loci: check which samples have 80% missing loci
###this line will generate reports in:  D:/test/qual_stat

d2        <- remove.poor.quality.snps(d1, min_repro=0.96, max_missing=0.2) 
###remove poor quality SNPs
###max missing: loci that 20% samples have missing info are removed. so the lower the score the more stringent..
###min_repro : minimum reproducibility score, a dart measure to test data quality, the better the data, the higher the score
qc2       <- report.dart.qc.stats(d2, RandRbase, species, dataset)

d3        <- sample.one.snp.per.locus.random(d2, seed=12345) 
###dart returns loci (sequences) that can have multiple SNPs, 
###this is because SNPs situated close together can occur due to linkage disequilibrium. our popgen studies dont rely on these assoiciations.. 
###so we only select one SNP to analyse if multiple SNPs occur on a locus
qc3       <- report.dart.qc.stats(d3, RandRbase, species, dataset)

# file <- write.dart.data(d3, RandRbase, species, dataset)
#this script writes the object, "d3" to a directory such as: rrspecies/PersHirs/dart_standard/raw_SNPFilt_1SNPperClone
#highlighted in yellow is the "treatment" given to the genotype data. this can be identified whenyou type "d3$treatment"
#you can save d1, d2, d3 in the same way
# load("rrspecies/PersHirs/dart_standard/raw_SNPFilt_1SNPperClone/PersHirs_DPers20-5061.rda")
#to  load d3

m1        <- read.meta.data(d3, RandRbase, species, dataset, fields=3)  
colnames(m1$analyses)
###import meta data from D:/test/meta/test_Dtest19-0000_meta.xlsx 
###make sure your meta data has at least 6 columns:
###sample, site, lat, long, and then 2 analysis columns with any name..
###the code looks for 2 analysis columns, if you want it to look for more columns, change:
###the "fields=2" to the number of your choice
dm        <- dart.meta.data.merge(d3, m1)
analysis <- "rrsites_SAM" #colnames(m1$analyses)[5]
fields    <- c(analysis) 
dms     <- data.by.meta.fields(dm, fields, RandRbase, species, dataset, object=analysis)
###dms contains the final cleaned data that you want to run popgen analysis with

#############generate a simple distance matrix of genotype data######################################
data <- as.matrix(dist(dms$gt))
#image(t(data)) #quick but ugly

library(reshape2)
melted_data <- melt(data)

library(ggplot2)
ggplot(data = melted_data, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()

#distance-based tree of data
library(ape)
plot(hclust(dist(dms$gt)))
par(mar=c(0,0,0,0)) #reset by: par(mar=c(2,2,2,2))
arbol <- nj(data); plot(arbol, cex=0.5)

library("igraph")
g1 = graph_from_incidence_matrix(data, directed=FALSE, weighted=TRUE)
hist(edge_attr(g1)$weight) #check histogram to find best cutoff, aim for the highest freq
mean(edge_attr(g1)$weight)
sd(edge_attr(g1)$weight)

cut.off <- 70
net.sp <- delete_edges(g1, E(g1)[weight<cut.off])
par(mai=c(0,0,0,0))
plot(g1,vertex.label=NA, layout=layout_with_fr,vertex.size=4,edge.width=0.1)
#plot(g1,vertex.label=NA, layout=layout_components,vertex.size=4,edge.width=0.1)#be vary of the results, it tries to find components
plot(net.sp,vertex.label=NA, layout=layout_as_star,vertex.size=4,edge.width=0.1) 

###downsample SNPs###########################################################################################################
dSVD <- downsample.snps(dms, number=50, seed=12345) #down sample SNPs

#######################PRINCIPAL COMPONENT ANALYSIS##########################################################################
library(adegenet) 
nd_gl  <- dart2gl(dms, RandRbase, species, dataset)
###converts the cleaned data to genlight format

nd_pca <- glPca(nd_gl, nf=5, parallel=FALSE) #nf indicates the number of principal components to be retained, a larger number (nf=300) is used for identifying clusters through k-means clustering / DAPC
scatter(nd_pca) #a quick way to check data


plot(nd_pca$scores[,1],nd_pca$scores[,2], xlab="PC1", ylab="PC2")
text(nd_pca$scores[,1],nd_pca$scores[,2], dms$meta$analyses[,analysis],pos=2,cex=0.3)

#########steps from here on makes pretty figs using the drawing program, "ggplot2"

sample_PC1_PC2 <- cbind(site=dms$meta$analyses[,analysis], #the cbind function combines different columns together
                        lat=dms$meta$lat,
                        long=dms$meta$long,
                        PC1=nd_pca$scores[,1], #PC1
                        PC2=nd_pca$scores[,2], #PC2
                        PC3=nd_pca$scores[,3]) #PC3

#write.table(sample_PC1_PC2, #this saves the table into a directory
#            PCAfile_directory, 
#            sep="\t",col.names=NA)

#remove rownames and add sample name to data, rownames do not get recognised by the drawing program
names <- rownames(sample_PC1_PC2);rownames(sample_PC1_PC2) <- NULL 
data <- cbind(names,sample_PC1_PC2)

df <- data.frame(data, stringsAsFactors = FALSE) #make data  from a matrix into dataframe, stringsAsFactors makes R view all columns as factors not numebrs

#make sure columns of numbers are viewed as numbers by R
df$PC1 <- as.numeric(df$PC1)
df$PC2 <- as.numeric(df$PC2)
df$PC3 <- as.numeric(df$PC3)
df$lat <- as.numeric(df$lat)
df$long <- as.numeric(df$long)

#make sure the data is ordered correctly, either by site or latitude
df2new <- df[order(-df$lat),] # this dataset is used in generating a summary of pops and their average latlongs for making maps/tables

df2new$site <- factor(df2new$site, levels = unique(df2new$site),ordered = TRUE)

library(plyr) #the data from df2new is summarised using ddply() in this R package
rep <- ddply(df2new, 
             .(site), 
             summarise,
             lat  = mean(lat),
             long = mean(long),
             PC1  = mean(PC1),
             PC2  = mean(PC2),
             PC3 = mean(PC3))

rep2 <- rep[order(-rep$lat),]
x <- length(unique(df2new$site)) #count how many pops there are
bgcols=rainbow(x) #generates X number of rainbow colours
rep2$nums_n_site <- paste0(1:x,": ",rep2$site)
rep2$num_site <- 1:x
rep222 <- rep2
rep222$nums_n_site <-factor(rep222$nums_n_site, levels=unique(rep222$nums_n_site)) #order data according to decreasing latitude

df3 <- merge(df,rep222[,c(1,7,8)], by="site") #note that merging reorders the data, so make sure both df and rep222 are in the same order

f <- nd_pca$eig[nd_pca$eig > sum(nd_pca$eig/length(nd_pca$eig))] #calculate variance explained
e <- round(f*100/sum(nd_pca$eig),1) #percent variance

library(ggplot2)
library(ggrepel)# learn ggplot, it is useful.
###http://www.sthda.com/english/wiki/ggplot2-scatter-plots-quick-start-guide-r-software-and-data-visualization

ggplot(data=df3, aes(x=PC1, y=PC2,col=nums_n_site))+ 
  geom_point(alpha=0.5, size=2) + #alpha is for transparency
  geom_text_repel(aes(label = num_site),col="black",size=2,data =df3,segment.color = NA)+
  theme(legend.text=element_text(size=8),
        legend.title = element_blank(),
        legend.position = "right",
        plot.title = element_text(face = "bold.italic", size=14))+
  #guides(col = guide_legend(ncol = 1)) + #how many columns in the legend
  ylab(paste0("PC2"," (",e[2],"%)"))+ 
  xlab(paste0("PC1"," (",e[1],"%)"))+
  scale_color_manual(labels=rep222$nums_n_site, 
                     values=c("blue","blue", "blue","red","grey","blue","blue"))
rainbow(length(unique(df2new$site)))
ggsave("D:/test/test_PC1_PC2_pops.tiff",width = 12, height = 8, dpi = 300, units = "in", device='tiff')
###this saves the ggplot to the directory of your choice

#########steps from here on to plot map  with dots and dot colour corresponds with PCA dot colour

###if you plan to save map, capture the image by running the tiff() below,
###tell it to stop capturing by running the dev.off() further down...
#tiff("E:/test/test_samp_distribPCA.tiff", units="in", width=12, height=8, res=480)

par(mai=c(0.1,0.1,0.2,0.1)) #par() sets where to plot figure, controls the amount of whitespace.

row_sub = apply(data.frame(dms$meta$long), 1, function(row) all(row !=0 )) #find any longitude in the data equals to zero 
newdmslong <- data.frame(dms$meta$long)[row_sub,]
row_sub = apply(data.frame(dms$meta$lat), 1, function(row) all(row !=0 )) #find any latitude in the data equals to zero 
newdmslat <- data.frame(dms$meta$lat)[row_sub,]

divxlims <- c(min(na.omit(newdmslong))-1,max(na.omit(newdmslong))+1) #find the min / max longitude
divylims <- c(min(na.omit(newdmslat))-1,max(na.omit(newdmslat))+1) #find the min / max latitude

library(oz) #draws australia coastlines and state boundaries
#try running:
#oz() or nsw()

oz(xlim=divxlims, ylim=divylims, lwd=0.5)
for (h in 1:x){
  points(rep222$long[h],rep222$lat[h],pch=21,lwd=0.2,bg=bgcols[h],cex=1) #put points on map
  }

rep23 <- rep222[complete.cases(rep222), ]#remove latlongs that are NAs

library(maptools) #run pointLabel, a function that moves around labels so that they dont overlap
pointLabel(rep23$long, rep23$lat, 
           labels = paste("  ", rep23$num_site, "  ", sep=""), 
           cex=1, font=4,offset=0.2)

dev.off()
#to revert back to default par(), type:
#par(mai=c(1,1,1,1))

############################################################################################
##Generate PCA plot with dots coloured according to latitude

PC1_PC2_plot <- ggplot(df, aes(PC1, PC2, col=lat))+ geom_point(alpha=0.5) + 
  theme(legend.position="none",aspect.ratio = 1,
        plot.title = element_text(face = "bold.italic", size=9),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 7),
        legend.title = element_text(size = 7))+
  ylab(paste0("PCA2"," (",e[2],"%)"))+ xlab(paste0("PCA1"," (",e[1],"%)"))+
  ggtitle(paste0(species, "\n"," PC1 vs PC2")) + scale_color_gradient(low="blue", high="red")

###PC1 vs PC3
PC1_PC3_plot <- ggplot(df, aes(PC1, PC3, col=lat))+ geom_point(alpha=0.5) + 
  theme(legend.position="none",aspect.ratio = 1,
        plot.title = element_text(face = "bold.italic", size=9),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 7))+
  ylab(paste0("PCA3"," (",e[3],"%)"))+ xlab(paste0("PCA1"," (",e[1],"%)"))+
  ggtitle(paste0("\n"," PC1 vs PC3")) + scale_color_gradient(low="blue", high="red")

###PC2 vs PC3
PC2_PC3_plot <- ggplot(df, aes(PC2, PC3, col=lat))+ geom_point(alpha=0.5) +
  ylab(paste0("PCA3"," (",e[3],"%)"))+ xlab(paste0("PCA2"," (",e[2],"%)"))+
  theme(legend.position="none",aspect.ratio = 1,
        plot.title = element_text(face = "bold.italic", size=9),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 7),
        legend.title = element_text(size = 7)) +    
  scale_color_gradient(low="blue", high="red")+
  ggtitle(paste0("\n"," PC2 vs PC3"))

PC_plot <- ggplot(df, aes(PC2, PC3, col=lat))+ geom_point() +
  ylab(paste0("PCA3"," (",e[3],"%)"))+ xlab(paste0("PCA2"," (",e[2],"%)"))+
  scale_color_gradient(low="blue", high="red")+
  ggtitle(paste0(" PC2 vs PC3"))
library(cowplot)
theme_set(theme_gray())

legend_b <- get_legend(PC_plot + theme(legend.position="bottom",
                                       legend.text=element_text(size=4))+ labs(col = "Latitude"))

# add the legend underneath the row we made earlier. Give it 10% of the height
# of one plot (via rel_heights).
prow <- plot_grid(PC1_PC2_plot,PC1_PC3_plot,PC2_PC3_plot, ncol=3)
p <- plot_grid( prow, legend_b, ncol = 1, rel_heights = c(1, .3))
p
ggsave(paste0("PCA/",species,"_PC1_PC2_PC3.tiff"),
       width = 8, height = 3, dpi = 300, units = "in", device='tiff')

###########################FST#############################################################################################

pFst      <- population.pw.Fst(dms, dms$meta$analyses[,analysis], RandRbase,species,dataset) #calculates genetic distance between populations
#fst is calculated using a function from R package SNPrelate, and Jason has written  to function into population.pw.Fst(), 
#calculates Fst using "W&H02" - relative beta estimator in Weir & Hill 2002, Jason's default
#an alternate estimator is "W&C84" - Fst estimator in Weir & Cockerham 1984
#to change the estimator, refer back to the population.pw.Fst() script

pS        <- population.pw.spatial.dist(dms, dms$meta$analyses[,analysis]) #calculates geographic distance between populations

####plot IBD plot

library(reshape2) #for melting data
library(vegan) #for mantel test

#tiff("E:/test/test fst plot.tiff", units="in", width=10, height=5, res=300)
Fst_sig <- cbind(melt(pS$S), unlist(as.list(pFst$Fst)))
colnames(Fst_sig)[3] <- "Geo_dist"
colnames(Fst_sig)[4] <- "Fst"
Fst_sig$Geo_dist2 <-Fst_sig$Geo_dist/1000 

plot(Fst_sig$Geo_dist2, Fst_sig$Fst, xlab="distance (km)", ylab="Fst",cex=1, 
     font=4, cex.main = 1) 

abline(lm(Fst_sig$Fst ~ Fst_sig$Geo_dist2), col = "blue")
title(main=paste0(species, " pairwise fst plots"),adj = 0.001,font.main=4) #from left (0) to right (1) with anything in between,
#line = positive values move title text up, negative - down)
man <- mantel(xdis = pS$S, ydis = pFst$Fst, permutations = 999, na.rm = TRUE) #mantel test for IBD
#if significant but R is low, maybe one population that is distant has low fst?
#if not significant but R is high, then sample size is low?
legend("bottomright", bty="n", cex=0.75,text.col="blue",
       legend=paste("Mantel statistic r is ", 
                    format(man$statistic, digits=4),
                    " P =",format(man$signif)))
#dev.off()

####plot geographic distance heatmap and fst heatmap 
library(ggplot2) #plots pretty figs
library(cowplot) #combines pretty figs

par(mfrow=c(2,1), oma=c(0,0,1,0))
geo_d <-pS$S #this is a square matrix
geo_d[upper.tri(geo_d)] <- NA #makes the upper triangular part of the matrix into nothing
rownames(geo_d) <- colnames(pS$S) #make sure rownames are the same as colnames

dimnames <- list (var1 = colnames(pS$S), var2 = colnames(pS$S)) 
mat <- matrix(geo_d, ncol=length(colnames(geo_d)), nrow=length(colnames(geo_d)), dimnames = dimnames)
df <- as.data.frame(as.table(mat))

p1 <- ggplot(df, aes(var1, var2)) + 
  geom_tile(aes(fill = Freq), colour = "white") +
  scale_fill_gradient(low = "white", high = "steelblue")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab("") + ylab("")+
  theme(legend.position="none",axis.text=element_text(size=5),plot.title = element_text(face = "bold.italic"))+
  ggtitle(paste0(species," heatmaps","\n","heatmap of pairwise geographic distance"))

genetic_d <-pFst$Fst
genetic_d[upper.tri(genetic_d)] <- NA
rownames(genetic_d) <- colnames(pFst$Fst)

dimnames2 <- list (var1 = colnames(pFst$Fst), var2 = colnames(pFst$Fst))
mat2 <- matrix(genetic_d, ncol=length(colnames(geo_d)), nrow=length(colnames(geo_d)), dimnames = dimnames)
df2 <- as.data.frame(as.table(mat2))

p2 <- ggplot(df2, aes(var1, var2, na.rm = TRUE)) + 
  geom_tile(aes(fill = Freq), colour = "white") +
  scale_fill_gradient(low = "white", high = "red")+
  geom_text(aes(label = round(Freq,3)),size=2,  df2)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("") + ylab("")+
  theme(legend.position="none",axis.text=element_text(size=7) ,plot.title = element_text(face = "bold.italic"))+
  ggtitle("heatmap of pairwise Fst")

plot_grid(p1,p2, ncol=1)

#ggsave("E:/test/test fst heatmaps.tiff"),
#       width = 8, height = 9, dpi = 300, units = "in", device='tiff')

###############################################LEA#######################################################################
#estimate individual allele frequencies through admixture proportions Q and population-specific
#allele frequencies F inferred from NMF (non-negative matrix factorisation)
#it is a dimension reduction and factor analysis method
#for finding a low-rank approximation of a matrix, which is similar to PCA

library(LEA)

kvalrange <- 1:16 

nd_lea <- dart2lea(dms, RandRbase, species, dataset)
snmf1=snmf(nd_lea, K=kvalrange, entropy = TRUE, repetitions = 1, project = "new") 
#snmf is a population structure method that uses non-negative matrix factorization algorithms
#if you wish to run this code a second time on your analysis, please go to your popgen folder
#and remove the entire lea folder. the code does not overwrite old outputs.

#tiff(entropy_file, units="in", width=8.6, height=5, res=300)
plot(snmf1, lwd = 6, col = "red", pch=1)#plot entropy graph
#dev.off()

###choose a k value to plot with colours chosen in the rainbow function
kval=2
kval_col = rainbow(kval) 

###choose the best run
snmf_file <- "C:\Users\yaps\Desktop\rrspecies\ToonAust\popgen\raw_SNPFilt_1SNPperClone_Field_rrsites\lea\ToonAust_DToo17-2792.snmf/K4/"
snmf_project <- load.snmfProject(snmf_file)
ce           <- cross.entropy(snmf_project, K = 2)
Rbest        <- which.min(ce) # with the lowest cross entropy criterion

#### open a Q file from LEA analysis for drawing barplot
##this Q file is in the popgen folder

leaQfile_dir <-paste0("C:/Users/yaps/Desktop/rrspecies/DorySass/popgen/raw_SNPFilt_1SNPperClone_Field_north/lea/DorySass_DDor16-2075.snmf/K2/run1/DorySass_DDor16-2075_r1.2.Q")
leaQfile <- read.table(leaQfile_dir)
par(mai=c(1.5,1,0,1))
barplot(t(as.matrix(leaQfile)), col=kval_col, ##draws barplot with labels underneath
        border = NA, space = 0, xlab = "", 
        ylab = "Admixture coefficients", 
        names.arg=dms$meta$analyses[, analysis], #this allows labels at the bottom of the barplots
        las=2,cex.names=0.5)

####draw LEA pies on map

lea_data <- data.frame(leaQfile)
lea_data2 <- cbind(dms$meta$analyses[, analysis], 
                   dms$meta$lat,dms$meta$long,lea_data
                   )
colnames(lea_data2)[1:3] <- c("site","lat","long")

coord <- data.frame(lea_data2$long, lea_data2$lat)
pop.factor <- factor(lea_data2$site)
pop = as.numeric(pop.factor)
npop <- length(unique(lea_data2$site))
qpop = matrix(NA, ncol = kval, nrow = npop)
coord.pop = matrix(NA, ncol = 2, nrow = npop)

#calculate aver age q values for each pop
for (c in unique(pop)){
  qpop[c,] = apply(leaQfile[pop == c,], 2, mean)#apply(m,2,mean), the "2" means get the mean of the columns
  coord.pop[c,] = apply(coord[pop == c,], 2, mean)}

#tiff(paste0(outputloc,"temp/", species," LEA pies_bar k=",kval,".tiff"), units="in", width=7, height=6, res=480)

library(oz)
library(mapplots)

par(mai=c(0.1,0.1,0.2,0.1))

#draw boundaries of australia, get min/max lat longs first
divxlims <- c(min(na.omit(dms$meta$long))-0.5,max(na.omit(dms$meta$long))+0.2)
divylims <- c(min(na.omit(dms$meta$lat))-0.5,max(na.omit(dms$meta$lat))+0.5)

oz(xlim=divxlims, ylim=divylims) 

for (n in 1:npop){
  add.pie(z = qpop[n,], x = coord.pop[n,1], y = coord.pop[n,2], 
          radius=(divylims[2]-divylims[1])/50, 
          labels = "",
          density = 100, # removing density fills pie slices
          col = kval_col)}

#dev.off()

##################################################Diversity######################################################

gp   <- dart2genepop(dms, RandRbase, species, dataset, dms$meta$analyses[,analysis])
#note that in dart2genepop, you can set your own min. allele frequency (maf_val). now it is set as maf_val=0.05
#the dart2genepop() also removes any loci that are missing for a population

#make sure that populations consisting of 1 individual have been removed as they can stuff up the scripts

library(diveRsity)
bs <- basicStats(infile = gp, outfile = NULL, #this step takes ages to run
                 fis_ci = FALSE, ar_ci = TRUE, 
                 ar_boots = 999, 
                 rarefaction = FALSE, ar_alpha = 0.05)

npop <- length(unique(dms$meta$analyses[,analysis]))
result <- mat.or.vec(npop,11)
measurement_names <- rownames(bs$main_tab[[1]])
population_names  <- names(bs$main_tab) #ls() rearranges the names 
rownames(result) <- population_names
colnames(result) <- measurement_names

for (r in 1:npop) {
  popstats <- bs$main_tab[[r]][,"overall"] ##extract from a list
  result[r,] <- popstats}

result <- as.data.frame(result)
result$sample_names <- rownames(result)

###getting latlongs into the diversity results
metnew2 <- cbind.data.frame(dms$meta$sample_names,dms$meta$analyses[,analysis], dms$meta$lat, dms$meta$long)
colnames(metnew2) <- c("sample_names", "site", "lat", "long")

data2 <- merge(result, metnew2, by= "sample_names", all.x=TRUE) #merge data

#write.table(data2, "diversity.csv",sep=",")

###finding number of Private alleles
gp2   <- dart2genepop0(dms, RandRbase, species, dataset, dms$meta$analyses[,analysis], maf_val = 0)

gp_genind <- read.genepop(gp2, ncode = 2)
gp_genind@other <- metnew2
strata(gp_genind) <- gp_genind@other
setPop(gp_genind) <- ~site
library(poppr)
p_allele <- data.frame(rowSums(private_alleles(gp_genind, locus ~ site, 
                                               count.alleles=F)))

p_allele$site <- rownames(p_allele)
rownames(p_allele) <- NULL
names(p_allele)[names(p_allele) == "rowSums.private_alleles.gp_genind..locus...site..count.alleles...F.."] <- "n_pa"
data <- merge(data2, p_allele, by.x="site",by.y="site")

metnew2_geo <- get_geodistinfo(metnew2)

data_last <- merge(data, metnew2_geo[,c(1,2,4)], by.x="site",by.y="site")

write.table(data_last, "C:/Users/yaps/Desktop/rrspecies/DorySass/diversity.csv",sep=",")

diversity_plots2(data)

source("C:/Users/yaps/Desktop/rrspecies/r scripts/diversity_plots_heat.R")
diversity_plots_heat(data,own_cols = c("yellow","purple")) #yellow is low, purple is high, make sure reset working directory

###############################Splitstree################################################################################
source("C:/Users/yaps/Desktop/RRtools/R/dart2splitstree.R") #save in species folder not popgen
snp <- dart2splitstree(dms, RandRbase, species, dataset,  dms$meta$analyses[,analysis], add_pop=TRUE)
#output nexus file will be in main directory after RandRbase

###############################SVDQuartet################################################################################
source("C:/Users/yaps/Desktop/RRtools/R/dart2svdquartetsSY.r") #save in species folder not popgen
dSVD <- dms
rownames(dSVD$gt) <- paste0(rownames(dms$gt),"_",dms$meta$analyses[,analysis]) #this allow changes to tree branches (need to copy the population info from previous file)
DorySVD <- dart2svdquartets(dSVD, RandRbase, species, dataset,  dSVD$meta$analyses[,analysis], add_pop=TRUE)
#the dart2svdquartets generates a nexus file of ATCGs and ambiguity codes

###############################SNAPP################################################################################
dSVD <- downsample.snps(dms, number=2000, seed=12345)
snp <- dart2snapp(dSVD, RandRbase, species, dataset,  dSVD$meta$analyses[,analysis], add_pop=TRUE)
#the dart2snapp generates a nexus file of 012?

###############################NeighbourJoining################################################################################

#NJ tree
library(ape)
gl_gm  <- dart2gl(dms, RandRbase, species, dataset)
tre <- nj(dist(as.matrix(gl_gm))) #generates  neighbor-joining tree estimate from euclidean matrix of the genotype data (genlight format)
plot(tre, type="fan", cex=0.5)
plot(tre, type="phylogram", cex=0.5)

#label tree tips according to sample and site name
#tre$tip.label <- paste0(dms$sample_name,"_",dms$meta$analyses[,analysis])
#plot(tre, type="fan", show.tip.label=TRUE, cex=0.3)
#plot(tre, type="phylogram", show.tip.label=TRUE, cex=0.3)

##########colour code according to PCA values
#PCA 
library(adegenet)
nd_gl  <- dart2gl(dms, RandRbase, species, dataset)
nd_pca <- glPca(nd_gl, nf=5, parallel=FALSE)
myCol <- colorplot(nd_pca$scores,nd_pca$scores, transp=TRUE, cex=2)
#scatter(nd_pca, xax=2,yax=3, posi="topleft", cex=0.1)
par(mfrow=c(1,2))
plot(nd_pca$scores[,1], nd_pca$scores[,2], xlab="PC1", ylab="PC2",col=myCol, pch=16)
text(nd_pca$scores[,1], nd_pca$scores[,2], labels=dms$meta$analyses[,analysis], pos=3, cex=0.5)

#NJ tree
plot(tre, typ="phylogram", show.tip=TRUE, no.margin=TRUE, cex=0.7)
tiplabels(pch=20, col=myCol, cex=2)

##########colour code according to latitude
dms_meta <- cbind.data.frame(dms$meta$sample_names, dms$meta$analyses[,analysis],dms$meta$lat, dms$meta$long)
colnames(dms_meta) <- c("sample_names","site", "lat", "long")

rbPal <- colorRampPalette(c('blue','red'))
dms_meta$col <- rbPal(50)[as.numeric(cut(dms_meta$lat,breaks = 50))]

#NJ tree
par(mfrow=c(1,2))
plot(tre, typ="phylogram", show.tip=TRUE, no.margin=TRUE, cex=0.6)
tiplabels(pch=20, col=dms_meta$col, cex=1)

#PCA
nd_gl  <- dart2gl(dms, RandRbase, species, dataset)
nd_pca <- glPca(nd_gl, nf=5, parallel=FALSE)
plot(nd_pca$scores[,1], nd_pca$scores[,2], xlab="PC1", ylab="PC2",col=dms_meta$col, pch=16)
text(nd_pca$scores[,1], nd_pca$scores[,2], labels=dms$meta$analyses[,analysis], pos=3, cex=0.5)

###############################TREEMIX########################################################################
source("C:/Users/yaps/Desktop/RRtools/R/dart2TreeMix.r")

tm_file <- dart2TreeMix(dms, RandRbase, species, dataset, pop=dms$meta$analyses[,analysis])

###############################Kinship analysis and UPGMA################################################################################

dC <- dms

#UPGMA
library(ape)
library(phangorn)

rownames(dC$gt) <- paste( rownames(dms$gt),dms$meta$analyses[,analysis],sep="_")
SNAPPC <- dist(dC$gt)
tCUNNINGHAMII <- upgma(SNAPPC)
tCUNNINGHAMII <- ladderize(tCUNNINGHAMII, right = FALSE)
is_tip <- tCUNNINGHAMII$edge[,2] <= length(tCUNNINGHAMII$tip.label)
ordered_tips <- tCUNNINGHAMII$edge[is_tip, 2]
tip_order <- tCUNNINGHAMII$tip.label[ordered_tips]

par(mar=c(1,1,1,1))
plot(tCUNNINGHAMII, cex=0.3)

#kin
#note: if your populations are highly differentiated, it is best to run kinship on individual populations not the entire species.
#expect that a clonal population might have an excess of heterozygotes because if the genet has 50% heterozygotes, 
#and the genet has multiple ramets, heterozygosity estimate of the population is more likely higher 
#as compared to a non-clonal population of genets with different number of heterozygotes 
#e.g. (50% x 3 vs 50%,30%,10%)

iIBD      <- individual.pw.IBD(dC,RandRbase,species,dataset)
kin       <-  iIBD$kinship #note if there is no kinship matrix, 
#it could be because there are no SNPs to compare
#or the snpgdsIBDMoM() in individual.pw.IBD script doesnt have: kinship=TRUE

#adds the row and column names corresponding to the tree to the spacial matrix
rownames(kin) <- paste0(dC$meta$sample_names,"_",dC$meta$analyses[,analysis])
colnames(kin) <- dC$meta$sample_names
# order kinship matrix using the tree
ik <- match(tip_order, rownames(kin))
ko <- kin[ik,ik]

#kinship heatmap
library(phytools)

png(file=paste0(RandRbase, species, "/",species,"_",analysis,"_kinship.png"), width = 6000, height = 6000,res=300)
phylo.heatmap(tCUNNINGHAMII,as.matrix(ko),fsize=0.3,ylim=c(0,1.1),split=c(0.25,0.75))#,split=c(0.25,0.75) ,,grid=T
dev.off()

write.table(as.matrix(ko), paste0(RandRbase, species, "/",species,"_",analysis,"_kinship.csv"),row.names = T,col.names = T)

####EXTRA KINSHIP ANALYSES IN CASE YOU WANT TO RUN KIN BY POP
b=6 #the number of the column with your analysis 
analysis <- colnames(m1$analyses)[b] 
colnames(m1$analyses)[b]
fields    <- c(analysis) 
dms     <- data.by.meta.fields(dm, fields, RandRbase, species, dataset, object=analysis)

####note this step is to calculate allele frequencies for each population separately and 
####then combine all data into one matrix
source("C:/Users/yaps/Desktop/rrspecies/r scripts/run_kinship.R")
###make sure there are no populations that only has one representative

big.mat <- run_kinship(d3_rm, species, RandRbase, b)
kin <- big.mat
image(kin)
write.csv(big.mat, paste0(RandRbase, species, "/",species,"_combinedkinship_",analysis,".csv"))

big.mat <- read.table(paste0(RandRbase, species, "/",species,"_combinedkinship_",analysis,".csv"),
                    header=T, sep=",", row.names = 1)
source("C:/Users/yaps/Desktop/rrspecies/r scripts/collect_ramet_sets.r")
source("C:/Users/yaps/Desktop/rrspecies/r scripts/ramets_to_prune.r")

rownames(kin) <- colnames(kin)
ls_ramets <- collect_ramet_sets(kin,0.45, lt=FALSE);print(ls_ramets)
write.csv(ls_ramets, paste0(RandRbase, species, "/",species,"_ls_ramets_",analysis,".csv"))

ramets_2p <- ramets_to_prune(ls_ramets, method="least_missing", gt=dms$gt)
write.csv(ramets_2p, paste0(RandRbase, species, "/",species,"_dontuse_ramets_",analysis,".csv"))

dg <- exclude.samples(dms, excluded_sample_names=ramets_2p, remove_fixed_loci=TRUE) 



###UPGMA
library(ape)
library(phangorn)
dC <- dms

rownames(dC$gt) <- paste( rownames(dms$gt),dms$meta$analyses[,analysis],sep="_")
SNAPPC <- dist(dC$gt)
tCUNNINGHAMII <- upgma(SNAPPC)
tCUNNINGHAMII <- ladderize(tCUNNINGHAMII, right = FALSE)
is_tip <- tCUNNINGHAMII$edge[,2] <= length(tCUNNINGHAMII$tip.label)
ordered_tips <- tCUNNINGHAMII$edge[is_tip, 2]
tip_order <- tCUNNINGHAMII$tip.label[ordered_tips]

par(mar=c(1,1,1,1))
plot(tCUNNINGHAMII, cex=0.7)

kin <- big.mat
ik <- match(tip_order, rownames(kin))
ko <- kin[ik,ik]
ik ###important to make sure there are no NAs, if so some thing has gone wrong

png(file=paste0(RandRbase, species, "/",species,"_combined_kinship.png"), width = 12000, height = 12000, res=300)
phylo.heatmap(tCUNNINGHAMII,as.matrix(ko),fsize=0.1,ylim=c(0,1.1),split=c(0.25,0.75))#,split=c(0.25,0.75) ,,grid=T
dev.off()

write.csv(ko, paste0(RandRbase, species, "/",species,"_combined_kinship_3.csv"))
###note that if you run this again, make sure you clear out your kinship folder or it will overwrite

dd <- read.table("rrspecies/RhodRube2/RhodRube2_pop_kinship.csv",
                      header=T, sep=",", row.names = 1)

dk <- match(colnames(big.mat), colnames(dd)); print(dk)
kd <- dd[dk,dk]
write.csv(kd, paste0("rrspecies/",species,"/",species,"_rrsites_SAM_kinship2.csv"))
write.csv(big.mat, paste0("rrspecies/",species,"/",species,"_pops_kinship.csv"))

diag(kd) <- NA
diag(big.mat) <- NA

kd_melt <- melt(as.matrix(kd))
big.mat_melt <- melt(big.mat)

par(mai=c(1,1,1,1))
plot(kd_melt$value,big.mat_melt$value, xlab="kd_melt", ylab="by pop")

all_data <- cbind(kd_melt, big.mat_melt)
write.csv(all_data, "rrspecies/RhodRube2/allkindata.csv")
###############################newhybrid
dms_less <- downsample.snps(dms, number=200, seed=12345) #down sample SNPs

source("RRtools/R/dart2newhy_SY.r") #script edited, so that output files all have an analysis label  
nhdir <- dart2newhy(dms_less, RandRbase, species, dataset, meta=dms_less$meta$analyses[,analysis])

#################################leaflet
library(ggplot2)
library(leaflet)
library(ggmap)
library(leaflet.minicharts)
library(htmlwidgets)

#############SATELLITE IMAGERY
dms_meta <- cbind(sample = dms$sample_names, site= dms$meta$analyses[,analysis], #the cbind function combines different columns together
                  lat=dms$meta$lat,long= dms$meta$long) 
write.table(dms_meta, paste0(RandRbase, species, "/meta/", species,"_",analysis,"_dms_meta.csv"))
dms_meta2 <- read.csv(paste0(RandRbase, species, "/meta/", species,"_",analysis,"_dms_meta.csv"))

# saveWidget(salty, file="E:/EUCALYPTUSCATTAI/saltwater_genet_map.html")
clumps <- unique(dms_meta2$genetic_clump)
dms_meta2_1 <- subset(dms_meta2, genetic_clump == clumps[1])
dms_meta2_2 <- subset(dms_meta2, genetic_clump == clumps[2])
dms_meta2_3 <- subset(dms_meta2, genetic_clump == clumps[3])
dms_meta2_4 <- subset(dms_meta2, genetic_clump == clumps[4])
dms_meta2_5 <- subset(dms_meta2, genetic_clump == clumps[5])
dms_meta2_6 <- subset(dms_meta2, genetic_clump == clumps[6])
dms_meta2_7 <- subset(dms_meta2, genetic_clump == clumps[7])
dms_meta2_8 <- subset(dms_meta2, genetic_clump == clumps[8])
dms_meta2_9 <- subset(dms_meta2, genetic_clump == clumps[9])

leaflet() %>% 
  addTiles() %>% 
  # setView(lng = min(na.omit(dms_meta2$long)), lat = min(dms_meta2$lat), zoom = 100) %>%
  # addProviderTiles("Esri.WorldImagery")%>%
  addCircles(lng= dms_meta2_1$long, lat=dms_meta2_1$lat, popup=dms_meta2_1$site, radius=rad, 
             stroke = TRUE, fillOpacity = 1,color= "red") %>% 
  addCircles(lng= dms_meta2_2$long, lat=dms_meta2_2$lat, popup=dms_meta2_2$site, radius=rad,
              stroke = TRUE, fillOpacity = 0.5,color= "orange") %>%
  addCircles(lng= dms_meta2_3$long, lat=dms_meta2_3$lat, popup=dms_meta2_3$site, radius=rad, 
             stroke = TRUE, fillOpacity = 1,color= "yellow") %>% 
  addCircles(lng= dms_meta2_4$long, lat=dms_meta2_4$lat, popup=dms_meta2_4$site, radius=rad,
             stroke = TRUE, fillOpacity = 0.5,color= "black") %>%
  addCircles(lng= dms_meta2_5$long, lat=dms_meta2_5$lat, popup=dms_meta2_5$site, radius=rad, 
             stroke = TRUE, fillOpacity = 1,color= "grey") %>% 
  addCircles(lng= dms_meta2_6$long, lat=dms_meta2_6$lat, popup=dms_meta2_6$site, radius=rad,
             stroke = TRUE, fillOpacity = 0.5,color= "green") %>%
  addCircles(lng= dms_meta2_7$long, lat=dms_meta2_7$lat, popup=dms_meta2_7$site, radius=rad, 
             stroke = TRUE, fillOpacity = 1,color= "darkgreen") %>% 
  addCircles(lng= dms_meta2_8$long, lat=dms_meta2_8$lat, popup=dms_meta2_8$site, radius=rad,
             stroke = TRUE, fillOpacity = 0.5,color= "blue") %>%
  addCircles(lng= dms_meta2_9$long, lat=dms_meta2_9$lat, popup=dms_meta2_9$site, radius=rad, 
             stroke = TRUE, fillOpacity = 1,color= "purple") %>% 
  addScaleBar(position = c("bottomright"), 
              options = scaleBarOptions(metric = TRUE, imperial = FALSE))

###note if you want to use other symbols make sure you start with e.g.,
# leaflet(mean_D_B) %>% 
#   addTiles() %>%
#   addProviderTiles("Esri.WorldImagery") %>% 
# to avoid this error:
#   Error in UseMethod("metaData") : 
#   no applicable method for 'metaData' applied to an object of class "NULL"

pchIcons = function(pch = 1, width = 30, height = 30, bg = "transparent", col = "black", ...) {
  n = length(pch)
  files = character(n)
  # create a sequence of png images
  for (i in seq_len(n)) {
    f = tempfile(fileext = '.png')
    png(f, width = width, height = height, bg = bg)
    par(mar = c(0, 0, 0, 0))
    plot.new()
    points(.5, .5, pch = pch[i], col = col[i], cex = min(width, height) / 8, ...)
    dev.off()
    files[i] = f
  }
  files
}

shapes = c(18) # base R plotting symbols (http://www.statmethods.net/advgraphs/parameters.html)
c1 <- c("#FF0000FF","#FFAA00FF")
#"yellow","#00FF00FF","grey","#00AAFFFF","#FF00AAFF","#AA00FFFF"
cols <- c(c1, "white")
iconFiles = pchIcons(shapes, 10, 10, 
                     col = c(cols), lwd = 1.5)
nums <- c(1,1,2,2,2,3,3,3,3,3,3,3)



leaflet(mean_D_B) %>% 
  addTiles() %>%
  addProviderTiles("Esri.WorldImagery") %>% 
  addMarkers(lng= mean_D_B$long, lat=mean_D_B$lat, popup=mean_D_B$gt,
             icon = ~ icons(iconUrl = iconFiles[nums]))%>% 
  addScaleBar(position = c("bottomright"), 
              options = scaleBarOptions(metric = TRUE, imperial = FALSE))

