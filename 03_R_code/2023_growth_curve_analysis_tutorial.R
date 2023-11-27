#### Code for calculating growth curve metrics from Bacterial OD data over time####
#### Authors: Jessica R. Bernardin ####
#### last update : November 27, 2023 ####

#### Load Required Packages ####
packages_to_load <- c(
  "ggplot2", "vegan", "lme4", "effects", "tidyverse", "growthcurver", "devtools",
  "plater", "janitor", "reshape", "readr"
)

# Load and install required packages
for (i in packages_to_load) { #Installs packages if not yet installed
  if (!require(i, character.only = TRUE)) install.packages(i)
}

#Set working directory
setwd(this.path::here())


#read in the parsed data and metadata
gc_data1 <- read.csv("../02_parsed_data/OD_11_10_23_long.csv", header=TRUE)
gc_meta1 <- read.csv("../05_metadata/plate_meta_11_16_2023.csv")

#remove time column
gc_data1 <- gc_data1[,-1]

#add a column for time
gc_data1$time<-0:192 #number of rows minus 1 bc 0 is the first time value
gc_data1$time <- 15*(gc_data1$time) # readings every 5 min starting at min 0

#Convert the "time" column from min to hours
gc_data1$time <- gc_data1$time / 60

#convert negative numbers to 0
gc_data1[gc_data1 < 0] <- 0

#put into long format
gc_data1_long <- melt(gc_data1, id=c("time"))
gc_data1_long$value <- as.numeric(gc_data1_long$value)

#remove wells that were empty
gc_data1_long_filt <- gc_data1_long %>%
  filter(variable %in% gc_meta1$Wells)

#visualize growth curves
ggplot(gc_data1_long_filt, aes(x=time, y=value)) + geom_point(alpha=0.5, size=.05) + theme_classic() +xlab("Time (hours)") + ylab("Bacterial Growth OD Plate 1 (Nov23)")

#filter data to look at negative controls
meta_blank <- gc_meta1 %>% filter(tube_id == "blank")
gc_data1_blanks <- gc_data1_long_filt %>%
  filter(variable %in% meta_blank$Wells)

ggplot(gc_data1_blanks, aes(x=time, y=value, color=variable)) + 
  geom_point(alpha=0.5, size=.5) + theme_classic() +xlab("Time (hours)") + 
  ylab("Bacterial Growth OD") + 
  ggtitle("Growth Curve for Blanks Plate 1 (Nov23)")

#looks like G8 was contaminated

### Let's try growthcurver to analyze the data
#Fits the data to a logistic function and calculates standard metrics like k(carrying capacity), N0(starting population), r(growthrate)

#filter meta data to remove negative controls
gcmeta_noblank <- gc_meta1 %>% filter(tube_id != "blank")

#filter wide data to match the samples above (no empty wells or blanks)
df_noblanks <- gc_data1 %>%
  select(intersect(names(.), gcmeta_noblank$Wells))

#convert data to numeric
df_noblanks <- as.data.frame(sapply(df_noblanks, as.numeric))

#isolate time column from original data
time <- gc_data1[,97]

#bind time column to our filtered dataframe in wide format
df_noblanks <- cbind(df_noblanks, time)

#uses data in wide format not long
gc_out1 <- SummarizeGrowthByPlate(df_noblanks, plot_fit = FALSE)

#Look at summary info and save as a txt
gc_file1 <- "../04_R_data_output/microbial_growth_metrics_Nov_2023.csv"
write.table(gc_out1, file = gc_file1, 
            quote = FALSE, sep = ",", row.names = FALSE)

#Look for normal distribution
hist(gc_out1$sigma, main = "Histogram of sigma values", xlab = "sigma")

# convert to dataframe for pca analysis
pca1_gc_out <- as_tibble(gc_out1)

# Prepare the gc_out data for the PCA
rownames(pca1_gc_out) <- pca1_gc_out$sample

# Do the PCA
pca1.res <- prcomp(pca1_gc_out %>% dplyr::select(k:sigma), center=TRUE, scale=TRUE)
summary(pca1.res)

# Plot the results
as_data_frame(list(PC1=pca1.res$x[,1],
                   PC2=pca1.res$x[,2],
                   samples = pca1_gc_out$sample)) %>% 
  ggplot(aes(x=PC1,y=PC2, label=samples)) + 
  geom_text(size = 3)


#Combine growth metric data with the metadata
gm1 <- merge(gc_out1, gc_meta1, by.x= "sample", by.y= "Wells")
gm1 <- gm1 %>% arrange(tube_id)
gm1$tube_id <- as.factor(gm1$tube_id)



######## STOP HERE
#there is a problem with the tube_id, it should be the same for the samples across rep 1, 2, 3


#ANOVA to look for differences in reps within tube_id
anova_growth_r_rep <- aov(r ~ tube_id*rep, data= gm1)
anova_growth_r_rep #no differences

dev.off()
ggplot(gm1, aes(x = rep, y = r, label = tube_id)) + geom_text() +
  geom_boxplot(outlier.shape = NA) +  geom_point() + 
  facet_grid(~tube_id) + 
  theme_classic() + ylab("Growth Rate") + xlab("Rep by Tube ID")


# aggregate and find mean for each bacterial isolate
gsum <- gm1 %>% group_by(tube_id) %>% 
  mutate(., r.mean = mean(r))

gsum <- gsum %>% group_by(tube_id) %>% 
  mutate(., n0.mean = mean(n0))

gsum <- gsum %>% group_by(tube_id) %>% 
  mutate(., k.mean = mean(k))

#remove columns that are not the averaged data
gsum_subset <- subset(gsum, select = -c(1:10))

#remove duplicates
gsum_subset <- gsum_subset %>%
  filter(!duplicated(r.mean))

#### Statistical Analysis ####
#read in metadata
strain_tax <- read.csv("05_metadata/Pitcher_plant_bacterial_isolates_metadata.csv", header=TRUE)
strain_enzyme <- read.csv("05_metadata/Pitcher_plant_bacterial_isolates_enzyme_rate.csv", header=TRUE)

#combine both metadata by strain id (we have more enzyme data than tax data so this will eliminate the extra)
list_df <- list(strain_tax,strain_enzyme)
df2 <- list_df %>% reduce(inner_join, by='Strain')

strain_meta_all <- merge(df2, gsum_subset, by.x = "ko_tube_number", by.y = "tube_id", all.x = TRUE)
strain_meta_all$genus <- as.factor(strain_meta_all$genus)
strain_meta_all$family <- as.factor(strain_meta_all$family)
strain_meta_all$order <- as.factor(strain_meta_all$order)

#look for differences in growth rates using taxonomy
aov_strain <- aov(r.mean~Strain, data=strain_meta_all)
summary(aov_strain)

aov_genus <- aov(r.mean~genus, data=strain_meta_all)
summary(aov_genus)

aov_fam <- aov(r.mean~family, data=strain_meta_all)
summary(aov_fam)

aov_order <- aov(r.mean~order, data=strain_meta_all)
summary(aov_order)

#look for differences in growth rates using enzyme activity
m1 <- glm(r.mean~Endochitinase + Lipase + Protease + β.N.acetylglucosaminidase + Chitobiosidase, data=strain_meta_all)
summary(m1)

## standardize the predictor variables using standardize package
strain_meta_all$lipase_scaled <- scale(strain_meta_all$Lipase)[, 1]
strain_meta_all$protease_scaled <- scale(strain_meta_all$Protease)[, 1]
strain_meta_all$bna_scaled <- scale(strain_meta_all$Endochitinase)[, 1]
strain_meta_all$endo_scaled <- scale(strain_meta_all$β.N.acetylglucosaminidase)[, 1]
strain_meta_all$cbio_scaled <- scale(strain_meta_all$Chitobiosidase)[, 1]

### model using brms and gamma dist.
#m2 <- brm(r.mean ~ endo_scaled + bna_scaled + cbio_scaled + lipase_scaled + protease_scaled,
#          data=strain_meta_all, family=Gamma(link="log"), iter = 10000, chains = 4, cores = 4) 

#saveRDS(m2, file = "~/Desktop/VIP/Bacteria_Growth_Curves/04_R_data_output/VIP2023_growthcurve_brm_model_enzymes.RDS")

m2 <- readRDS("~/Desktop/VIP/Bacteria_Growth_Curves/04_R_data_output/VIP2023_growthcurve_brm_model_enzymes.RDS") 

summary(m2)
#plot(m2)
#fixef(m2)
#plot(marginal_effects(m2), points = TRUE)
pp_check(m2)
dispersion <- function(x) {var(x)/mean(x)}
ppc_stat(y = strain_meta_all$r.mean,  
         # Compare the dispersion in the real data (y)
         yrep = posterior_predict(m2), 
         # to dispersion in predictions from our posterior (yrep)
         stat="dispersion")

posterior.m2 <- as.data.frame(m2)

## Visualize the credible intervals for our model's parameters:
mcmc_intervals(posterior.m2, 
               prob = 0.5, # Sets darker blue CI value
               prob_outer = 0.9, # Sets lighter blue CI value
               regex='^b') #Plot only the "betas" from our model
#density plots
color_scheme_set("teal")
#'arg' should be one of “blue”, “brightblue”, “darkgray”, “gray”, “green”, “orange”, “pink”, “purple”, “red”, “teal”, “yellow”, “viridis”, “viridisA”, “viridisB”, “viridisC”, “viridisD”, “viridisE”
plot_title <- ggtitle("Posterior distributions for Mean Bacterial Growth Rate",
                      "with medians and 80% intervals")

#Posterior distributions for Mean Bacterial Growth Rate with medians and 80% intervals
mcmc_areas(posterior.m2, pars=c("b_Intercept", "b_lipase_scaled",
                      "b_cbio_scaled", "b_protease_scaled", "b_endo_scaled", "b_bna_scaled"), prob=.8) +
  theme_classic() + xaxis_text(size=18, color="black")+ yaxis_text(size=18, color="black") +
  geom_vline(xintercept=0, linetype="dotted", color="darkgray")


#positive or negative effect
lipase_neg <- posterior.m2 %>% filter(b_lipase_scaled < 0)
#prob that b_lipase_scaled rate is pos correlated with growth rate
nrow(lipase_neg)/nrow(posterior.m2) #the probability of direction
#0.8604

#positive or negative effect
bna_pos <- posterior.m2 %>% filter(b_bna_scaled < 0)
#prob that b_bna_scaled rate is pos correlated with growth rate
nrow(bna_pos)/nrow(posterior.m2) #the probability of direction
#0.799
#evidence for a potential effect, don't have to ignore just because it crosses 0

# predictive distribution
plot(strain_meta_all$r.mean~strain_meta_all$lipase_scaled, 
     pch=21, bg=c("#f5bc42"),
     xlab="Lipase Rate (scaled)",
     ylab="Mean Bacterial Growth Rate")

## good news! Brms comes with its own function to predict from the posterior for you:
preds <- posterior_predict(m2, 
                           #Can set newdata = .... to any data you'd like to predict for
                           ndraws = 1000 # sets the number of posterior samples to use
)
#This "for loop" will plot all of the predictions from our posterior (for this set of values of our IVs), generated by each row (sample) from our posterior.
# The first part (for(i in...)) will repeat the subsequent action for each row [i] in our posterior (so each  MCMC sample):

x <- for(i in 1:nrow(preds)){ 
  curve(exp(      #Inverse Link function for gamma
    posterior.m2$b_Intercept[i]+ #Posterior for intercept
      posterior.m2$b_lipase_scaled[i]*0 + # holding prot density at mean
      posterior.m2$b_bna_scaled[i]*x),
    add=T, # Add to our existing plot
    col=rgb(0,0,0,0.01)) }

ggplot(posterior.m2, aes(x = b_lipase_scaled, fill = stat(x < 0))) + stat_halfeye() +
  scale_fill_manual(values=c( "grey50", "#20a198"))+
  geom_vline(aes(xintercept=0),  color="black", size=1, linetype="dashed")+
  ylab("density") + theme_classic() + guides(fill="none")

ggplot(posterior.m2, aes(x = b_bna_scaled, 
                         fill = stat(x > 0))) +stat_halfeye() +
  scale_fill_manual(values=c("#20a198", "grey50"))+
  geom_vline(aes(xintercept=0),  color="black", size=1, linetype="dashed")+
  ylab("density") +theme_classic() + guides(fill="none")

#visualize some of differences in taxonomy by r.mean
#reorder
strain_meta_all$r.mean <- as.numeric(strain_meta_all$r.mean)
strain_meta_all <- strain_meta_all[order(strain_meta_all$r.mean, decreasing = TRUE),]

ggplot(strain_meta_all, aes(x = fct_inorder(order), y = r.mean, fill = order)) + 
  geom_boxplot() +
  coord_flip() + ylab("Mean Growth Rate") + xlab("Order") +theme_classic()+ guides(color = guide_legend(title ="Order"))

ggplot(strain_meta_all, aes(x = fct_inorder(family), y = r.mean, fill = family)) + 
  geom_boxplot() +
  coord_flip() + ylab("Mean Growth Rate") + xlab("Family")+theme_classic()+ guides(color = guide_legend(title ="Family"))

ggplot(strain_meta_all, aes(x = fct_inorder(genus), y = r.mean, fill = genus)) + 
  geom_boxplot() +
  coord_flip() + ylab("Mean Growth Rate") + xlab("Genus")+theme_classic()+ guides(color = guide_legend(title ="Genus"))

#plots that look at rate by tax using tube number
ggplot(strain_meta_all, aes(x = fct_inorder(order), y = r.mean, label = ko_tube_number)) + geom_text(vjust = 0, nudge_x = 0.2, check_overlap = TRUE) +geom_point(aes(colour = order)) +
  theme_classic() + ylab("Mean Growth Rate") + xlab("Order") +
  coord_flip()+ guides(color = guide_legend(title ="Order"))

#plot for graph
ggplot(strain_meta_all, aes(x = fct_inorder(order), y = r.mean))+geom_point(aes(colour = order), show.legend = FALSE, cex=5) +
  theme_classic() + ylab("Mean Growth Rate") + xlab("Order") + scale_color_viridis(discrete=TRUE) +
  coord_flip()+ theme(text=element_text(size = 18, colour="black"), axis.text = element_text(color="black"))

ggplot(strain_meta_all, aes(x = fct_inorder(family), y = r.mean, label = ko_tube_number)) + geom_text(vjust = 0, nudge_x = 0.2, check_overlap = TRUE) +geom_point(aes(colour = family)) +
  theme_classic() + ylab("Mean Growth Rate") + xlab("Family") +
  coord_flip()+ guides(color = guide_legend(title ="Family"))

ggplot(strain_meta_all, aes(x = fct_inorder(genus), y = r.mean, label = ko_tube_number)) + geom_text(vjust = 0, nudge_x = 0.2, check_overlap = TRUE) +geom_point(aes(colour = genus)) +
  theme_classic() + ylab("Mean Growth Rate") + xlab("Genus") +
  coord_flip()+ guides(color = guide_legend(title ="Genus"))

strain_meta_all <- strain_meta_all[,-c(12:14)]
strain_meta_all <- strain_meta_all[,-17]

#all growth data from 6 reps are averaged, one weird rep from 16 filtered out
write.csv(strain_meta_all, "~/Desktop/VIP/Bacteria_Growth_Curves/04_R_data_output/2023_bacterial_isolates_mean_growth_data.csv", row.names=FALSE)

#### filter to only winners and losers (just 8 taxa) ####
win_los <- subset(strain_meta_all, !is.na(Winners_losers))
ggplot(win_los, aes(x = fct_inorder(genus), y = r.mean, label = Winners_losers)) + geom_text(vjust = 0, nudge_x = 0.2, check_overlap = TRUE) +geom_point(aes(colour = genus)) +
  theme_classic() + ylab("Mean Growth Rate") + xlab("Genus") +
  coord_flip()+ guides(color = guide_legend(title ="Genus"))

qqnorm(win_los$r.mean)
qqline(win_los$r.mean, col="green")
win_los$r.mean
win <- c(0.5695526, 0.5834786, 0.4926402, 0.8350893)
los <- c(0.3426037, 0.2607637, 0.4076297, 0.2342006)
var(win) #0.02212084
var(los) #0.006252523
var(win)/var(los) #3.537907 since less than 4 we can assume the variances between the two groups are appx =
df <- data.frame(win, los)
t.test(win, los) #t = 3.6676, df = 4.5705, p-value = 0.01697, true difference in means is not equal to 0

#### make PCA better
PC1 <- as_data_frame(list(PC1=pca1.res$x[,1],
                          PC2=pca1.res$x[,2],
                          samples = pca1_gc_out$sample))

PC1_merge <- merge(x=PC1, y=gc_meta1, by.x ="samples", by.y = "Wells")
PC1_merge <- merge(x=PC1_merge, y=strain_meta_all, by.x ="tube_id", by.y = "ko_tube_number")
PC1_merge$order <- as.factor(PC1_merge$order)
PC1_merge$family <- as.factor(PC1_merge$family)
PC1_merge$genus <- as.factor(PC1_merge$genus)

ggplot(PC1_merge, aes(x = PC1, y = PC2)) + 
  geom_point(size = 4, aes(colour = order))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "PC1", colour = "Order", y = "PC2")

ggplot(PC1_merge, aes(x = PC1, y = PC2)) + 
  geom_point(size = 4, aes(colour = family))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "PC1", colour = "Family", y = "PC2")

ggplot(PC1_merge, aes(x = PC1, y = PC2)) + 
  geom_point(size = 4, aes(colour = genus))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "PC1", colour = "Genus", y = "PC2")

#Enzymes
ggplot(PC1_merge, aes(x = PC1, y = PC2)) + 
  geom_point(size = 4, aes(colour = Endochitinase))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "Endochitinase", y = "NMDS2")+ scale_colour_gradientn(colours=rainbow(4))

ggplot(PC1_merge, aes(x = PC1, y = PC2)) + 
  geom_point(size = 4, aes(colour = β.N.acetylglucosaminidase))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "β.N.acetylglucosaminidase", y = "NMDS2")+ scale_colour_gradientn(colours=rainbow(4))

ggplot(PC1_merge, aes(x = PC1, y = PC2)) + 
  geom_point(size = 4, aes(colour = Lipase))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "Lipase", y = "NMDS2")+ scale_colour_gradientn(colours=rainbow(4))

ggplot(PC1_merge, aes(x = PC1, y = PC2)) + 
  geom_point(size = 4, aes(colour = Protease))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "Protease", y = "NMDS2")+ scale_colour_gradientn(colours=rainbow(4))

ggplot(PC1_merge, aes(x = PC1, y = PC2)) + 
  geom_point(size = 4, aes(colour = Chitobiosidase))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "Chitobiosidase", y = "NMDS2")+ scale_colour_gradientn(colours=rainbow(4))

### plotting growth curves with averages across all the reps ####
#put into long format
gc_data1.t_noblanks_long <- melt(gc_data1.t_noblanks, id=c("time"))
gc_data1.t_noblanks_long$value <- as.numeric(gc_data1.t_noblanks_long$value)

gc_data2.t_noblanks_long <- melt(gc_data2.t_noblanks, id=c("time"))
gc_data2.t_noblanks_long$value <- as.numeric(gc_data2.t_noblanks_long$value)

#visualize growth curves
ggplot(gc_data1.t_noblanks_long, aes(x=time, y=value)) + geom_point(alpha=0.5, size=.05) + theme_classic() +xlab("Time (hours)") + ylab("Bacterial Growth OD Plate 1 (2/1/23)")
ggplot(gc_data2.t_noblanks_long, aes(x=time, y=value)) + geom_point(alpha=0.5, size=.05) + theme_classic() +xlab("Time (hours)") + ylab("Bacterial Growth OD Plate 2 (2/8/23)")

metaplate_metastrain1 <- merge(gc_meta1, strain_tax, by.x = "tube_id", by.y = "ko_tube_number")
test <- merge(gc_data1.t_noblanks_long, metaplate_metastrain1, by.x = "variable", by.y = "Wells")

metaplate_metastrain2 <- merge(gc_meta2, strain_tax, by.x = "tube_id", by.y = "ko_tube_number")
metaplate_metastrain2$rep <- as.numeric(metaplate_metastrain2$rep)
metaplate_metastrain2$rep <- metaplate_metastrain2[, 3] + 3
test2 <- merge(gc_data2.t_noblanks_long, metaplate_metastrain2, by.x = "variable", by.y = "Wells")
#remove H03 from plate 2, not good readings
test2.filt <- test2 %>% filter(variable !="H03")
abs_all <- rbind(test, test2.filt)
abs_all_ave = abs_all %>%
  group_by(tube_id, time) %>%
  summarise(
    abs_mean = mean(value),
    abs_sd = sd(value))

abs_all_ave_meta <- merge(abs_all_ave, strain_meta_all, by.x = "tube_id", by.y = "ko_tube_number")

#visualize growth curves
ggplot(abs_all_ave_meta, aes(x=time, y=abs_mean, color=order)) + geom_point(alpha=0.5, size=0.7) + theme_classic() +
  xlab("Time (hours)") + ylab("Bacterial Growth Rate (r)")+ theme(text=element_text(size = 18)) + guides(color = guide_legend(title ="Order"))

#plot for poster
ggplot(abs_all_ave_meta, aes(x=time, y=abs_mean, color=order)) + geom_point(alpha=0.5, size=0.8) + theme_classic() +
  xlab("Time (hours)") + ylab("Mean Bacterial Growth Rate (r)")+ theme(text=element_text(size = 18)) + guides(color = guide_legend(title ="Order"))+
  scale_color_viridis(discrete=TRUE)

ggplot(abs_all_ave_meta, aes(x=time, y=abs_mean, color=order)) + geom_point(alpha=0.5, size=0.7) + theme_classic() +
  xlab("Time (hours)") + ylab("Bacterial Growth Rate (r)")+ theme(text=element_text(size = 18)) + guides(color = guide_legend(title ="Order"))+
  geom_dl(aes(label = Winners_losers), method = list(dl.combine("last.points")), cex = 0.7)

ggplot(abs_all_ave_meta, aes(x=time, y=abs_mean, color=family)) + geom_point(alpha=0.5, size=0.7) + theme_classic() +xlab("Time (hours)") + ylab("Bacterial Growth Rate (r)")+ theme(text=element_text(size = 18)) + guides(color = guide_legend(title ="Family"))
ggplot(abs_all_ave_meta, aes(x=time, y=abs_mean, color=genus)) + geom_point(alpha=0.5, size=0.7) + theme_classic() +xlab("Time (hours)") + ylab("Bacterial Growth Rate (r)")+ theme(text=element_text(size = 18)) + guides(color = guide_legend(title ="Genus"))

#rmarkdown::render("03_R_code/2023_growth_curve_analysis_34isolates.R", "html_document")
