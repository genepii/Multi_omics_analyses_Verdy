### FIGURE 1A

# Load necessary libraries
library(reshape2)
library(dplyr)
library(ggplot2)
library(ggsci)

# Load the relative abundance matrix
relab_matrix <- read.table("Relab_virome_fam_all.txt", header = TRUE, stringsAsFactors = FALSE)

# Melt the data for long-format processing
melted_data <- reshape2::melt(relab_matrix)
colnames(melted_data) <- c("Viral_family", "all", "Value")

# Create the pie chart
p <- ggplot(melted_data, aes(x = "", y = Value, fill = Viral_family)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  theme_void() +
  geom_text(aes(label = paste0(round(Value, 1), "%")),
            position = position_stack(vjust = 0.5),
            size = 6) +
  scale_fill_npg()
# Display the plot
print(p)

pdf(file = "fam_pie_chart.pdf", width = 10, height = 10)
p
dev.off()

### FIGURE 1B
##################################################### plot columns MRA and occurence

library(ggplot2)
library(reshape2)
library(readr)
library(RColorBrewer)
library(ggsci)

library(ggplot2)

otu <- read.table("votu.txt", header = TRUE, stringsAsFactors = FALSE)

otu_presence <- otu > 0

# Step 4: Summarize the presence of each OTU across samples as a percentage
otu_presence_count <- rowSums(otu_presence)  # Counts of presence across samples
total_samples <- ncol(otu)  # Total number of samples
otu_percentage <- (otu_presence_count / total_samples) * 100  # Convert to percentage

# Step 5: Filter OTUs based on minimum sample presence (optional)
min_samples <- 1  # Set your threshold here
otu_matrix_filtered <- otu[otu_percentage >= (min_samples / total_samples) * 100, ]
otu_percentage_filtered <- otu_percentage[otu_percentage >= (min_samples / total_samples) * 100]

# Step 6: Create a data frame for plotting, merging with Family annotation
otu_summary <- data.frame(
  OTU = rownames(otu_matrix_filtered),
  Occurrence_Percentage = otu_percentage_filtered
)

#write.table(otu_summary, file = "prevalence_vir.txt", sep = "\t")
data<-otu_summary
data  = read.delim(file="prevalence_vir.txt",header=T,check.names=F)

data$Occurrence_Percentage <- as.numeric(data$Occurrence_Percentage)

library(ggsci)

# Define custom colors for specific families
custom_colors <- c(
  "Viruses_Unclassified" = "yellow2",
  "Anelloviridae" = "#C4A484",
  "Autographiviridae" = "pink2",
  "Caudoviricetes_Unclassified" = "aquamarine4",
  "Potyviridae" ="aquamarine2",
  "Partitiviridae" ="slateblue1",
  "Papillomaviridae"="orange2",
  "Herelleviridae"="slateblue4"
)

sig <- ggplot(data, aes(x = MRA, y = Occurrence_Percentage, color = Family)) +
  geom_point(alpha = 0.4, size = 5) +
  theme_classic() +
  labs(title = "", x = "Log10MRA", y = "Occurrence Percentage") +
  theme(
    legend.position = "right",
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.title.x = element_text(size = 20)
  ) +
  scale_color_manual(
    values = c(custom_colors, ggsci::pal_npg()(length(unique(data$Family))))
  )

sig

pdf("All_otu_MRA_occurence2.pdf",width=8,height=5);
sig
dev.off()
ggsave("All_otu_MRA_occurence.jpeg", plot = sig, width = 7, height = 5, dpi = 300)


### FIGURE 1C

library(dplyr)
library(tidyr)
library(ggplot2)
library(PieGlyph)
library(ggiraph)
library(reshape2)
library(ggplot2)
library(ggforce)
library(patchwork)

df <- read.delim("Relab_fam.txt", sep = "\t", check.names = TRUE, row.names = 1)
meta <- read.delim("metadata.txt", sep = "\t", check.names = TRUE, row.names = 1)

# Convert df to a matrix and melt it
df <- as.matrix(df)
df <- melt(df)

# Rename columns
colnames(df) <- c("Viral_family", "Sample", "value")

# Ensure Sample names exist in metadata row names
df$minutes <- sapply(as.vector(df$Sample), function(x) if (x %in% rownames(meta)) meta[x, "min"] else NA)

# Display first few rows
head(df)

df$Viral_family <- factor(df$Viral_family)
color_map <- setNames(ifelse(levels(df$Viral_family) == "none", "white", pal_npg("nrc")(length(levels(df$Viral_family)))), levels(df$Viral_family))

piechart_fam <- ggplot(data = df, aes(x = minutes, y = reorder(Sample, minutes), fill = Viral_family)) +
  geom_pie_glyph(slices = 'Viral_family', values = 'value', radius = 0.5, colour = "black") +
  theme_linedraw() +
  scale_fill_manual(values = color_map) +
  scale_x_continuous(breaks = c(seq(0, 600, by = 50), seq(1200, 1300, by = 50)))+
  coord_cartesian(xlim = c(0, 1300)) +
  theme(panel.background = element_rect(fill = "white"),
        strip.text = element_text(size = 10),
        legend.title = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        plot.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.position = "right",
        legend.direction = "vertical") + 
  labs(
    y = "Samples"
  )

piechart_fam

############################
custom_colors <- c(
  "Viruses_Unclassified" = "yellow2",
  "Anelloviridae" = "#C4A484",
  "Autographiviridae" = "pink2",
  "Caudoviricetes_Unclassified" = "aquamarine4",
  "Potyviridae" = "aquamarine2",
  "Partitiviridae" = "slateblue1",
  "Papillomaviridae" = "orange2",
  "Herelleviridae" = "slateblue4",
  "none" = "white"
)


plot1 <- ggplot(data = df, aes(x = minutes, y = reorder(Sample, minutes))) +
  geom_pie_glyph(slices = 'Viral_family', values = 'value', radius = 0.5, colour = "black") +
  theme_linedraw() +
  scale_fill_manual(values = custom_colors) +
  scale_x_continuous(
    breaks = seq(0, 600, by = 50),  # Set the breaks for this range
    limits = c(0, 600)  # Limit the x-axis to 0 to 600
  ) +
  coord_cartesian(xlim = c(0, 600)) +  # Ensure the plot zooms into the 0 to 600 range
  theme(
    strip.text = element_text(size = 10),
    legend.title = element_blank(),
    axis.title.y = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    plot.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.position = "none",
    legend.direction = "vertical"
  ) +
  labs(y = "Samples")

# Second plot: showing the range from 1200 to 1300
plot2 <- ggplot(data = df, aes(x = minutes, y = reorder(Sample, minutes))) +
  geom_pie_glyph(slices = 'Viral_family', values = 'value', radius = 0.5, colour = "black") +
  theme_linedraw() +
  scale_fill_manual(values = custom_colors) +
  scale_x_continuous(
    breaks = seq(1250, 1300, by = 50),  # Set the breaks for this range
    limits = c(1250, 1300)  # Limit the x-axis to 1200 to 1300
  ) +
  coord_cartesian(xlim = c(1250, 1300)) +  # Ensure the plot zooms into the 1200 to 1300 range
  theme(
    panel.background = element_rect(fill = "white"),
    strip.text = element_text(size = 10),
    legend.title = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_blank(),
    plot.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.position = "right",
    legend.direction = "vertical"
  )

# Combine the two plots side by side using patchwork
library(gridExtra)
plot1 + plot2

pdf(file = "piechart_fam.pdf", width = 12, height = 10)  # Open PDF device
grid.arrange(plot1, plot2, ncol = 2, widths = c(0.7, 0.3))  # Place plotting code here
dev.off() 





### FIGURE 2

#### boxplot WILCOX NORMAL

library(ggplot2)
library(tidyverse)
library(rstatix)
library(ggplot2)
library(ggpubr)
library(dplyr)

dat <- read.delim("metadata - Copie.txt", stringsAsFactors = FALSE)
head(dat)

#custom_colors <- c("CS" = "PURPLE", "VD" = "PINK")
#custom_colors <- c("moderate_severe" = "red3", "No_mild" = "blue3")
custom_colors <- c("Oui" = "orange3", "Non" = "plum4")



#dat <- dat %>%
#  mutate(across(where(is.numeric), ~log10(.x+1)))

#dat <- dat %>%
#mutate(across(
#  .cols = where(is.numeric) & !all_of(c("DAY", "DAY2", "cluster")),
#  .fns = ~log10(.x + 1)
#))

dat %>% sample_n_by(ATB, size = 2)
dat %>%
  group_by(ATB) %>%
  get_summary_stats(VDR_total_vir, type = "median_iqr")

stat.test <- dat %>% 
  wilcox_test(VDR_total_vir ~ ATB) %>%
  add_significance()
stat.test
dat %>% wilcox_effsize(VDR_total_vir ~ ATB)
stat.test <- stat.test %>% add_xy_position(x = "ATB")

p <-ggplot(dat, aes(ATB, VDR_total_vir)) +  
  # Boxplot layer: will be drawn first, in the background
  geom_boxplot(aes(fill = ATB), width = 0.2, color = "black", outlier.shape = NA, alpha = 0.2) +  
  scale_color_manual(values = custom_colors) +  
  scale_fill_manual(values = custom_colors) +  
  
  # Add p-value annotation with stat_pvalue_manual()
  stat_pvalue_manual(stat.test, tip.length = 0, size = 10, bracket.size = 2) +
  
  # Customize labels and theme
  labs(subtitle = get_test_label(stat.test, detailed = TRUE)) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 20, colour = "black"),
    axis.text.y = element_text(size = 20),
    plot.subtitle = element_text(size = 20),
    legend.text = element_text(size = 20),
    legend.position = "none"
  ) +
  
  # Jittered points layer: will be drawn last, in the foreground
  geom_jitter(aes(color = ATB), width = 0.03, size = 10, alpha = 0.1) 
p

pdf("VDR_ATB_boxplot.pdf",width=5,height=5);
p
dev.off()
# kruskaL normal


library(tidyverse)
library(ggpubr)
library(rstatix)
library(ggsci)

data <- read.delim("metadata.txt", stringsAsFactors = FALSE)

head(data)

custom_colors <- c("24-25" = "gray", "26-27" = "firebrick", "28-29" = "blue")


data <- data %>%
  reorder_levels(GA_window, order = c("24-25", "26-27","28-29"))
head(data)
data %>% 
  group_by(GA_window) %>%
  get_summary_stats(Bacteriophages, type = "common")

ggboxplot(data, x = "GA_window", y = "Bacteriophages")

res.kruskal <- data %>% kruskal_test(Bacteriophages~GA_window )
res.kruskal

data %>% kruskal_effsize(Bacteriophages ~ GA_window )
pwc <- data %>% 
  dunn_test(Bacteriophages ~ GA_window, p.adjust.method = "bonferroni") 
pwc
pwc2 <- data %>% 
  wilcox_test(Bacteriophages ~ GA_window, p.adjust.method = "bonferroni")
pwc2

pwc <- pwc2 %>% add_xy_position(x = "GA_window")

# Plot
p <- ggboxplot(
  data, x = "GA_window", y = "Bacteriophages", fill = "GA_window",
  alpha = 0.5, width = 0.2
) +
  geom_jitter(aes(color = GA_window), width = 0.03, size = 10, alpha = 0.2) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  labs(
    subtitle = paste0("Kruskal-Wallis Test, Overall p-value: ", signif(res.kruskal$p, digits = 3))
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 20, colour = "black"),
    axis.title = element_text(size = 20, colour = "black"),
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 14)
  )

p

pdf(file = "vdr_Bacteriophages_GA_window.pdf", height = 6, width = 6)
p
dev.off()

##################################################### NMDS AND ANOSIM
library(vegan)
library(tidyverse)
set.seed(1)
data1<-read.delim("vdr_nmds.txt", row.names = 1)
dataTransposed1<-t(data1)
dist.1 <- vegdist(dataTransposed1, method = "bray")
metadata <- read.delim("metadata_nmds.txt")


ano = anosim(dataTransposed1, metadata$DM, distance = "bray", permutations = 9999)
ano
plot(ano)

nmds = metaMDS(dataTransposed1, distance = "bray")
nmds
plot(nmds)


ta.scores = as.data.frame(scores(nmds)$sites)
metadata$Samples = metadata$Samples
metadata$cluster = metadata$cluster
metadata$HAP_onset = metadata$HAP_onset

BPD_col <- c("moderate_severe" = "indianred", "No_mild" = "royalblue3")
BPD_col <- c("yes" = "orange3", "no" = "plum4")
DM_col <- c("CS" = "purple", "VD" = "pink4")

xx = ggplot(metadata, aes(x = ta.scores$NMDS1, y = ta.scores$NMDS2)) + 
  geom_point(aes(size = 10, colour = DM)) +  # Corrected color reference
  theme(
    axis.text.y = element_text(colour = "black", size = 20, face = "bold"), 
    axis.text.x = element_text(colour = "black", face = "bold", size = 20), 
    legend.text = element_text(size = 20, face = "bold", colour = "black"), 
    legend.position = "right", 
    axis.title.y = element_text(face = "bold", size = 20), 
    axis.title.x = element_text(face = "bold", size = 20, colour = "black"), 
    legend.title = element_text(size = 20, colour = "black", face = "bold"), 
    panel.background = element_blank(), 
    panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
    legend.key = element_blank()
  ) + 
  labs(x = "NMDS1", y = "NMDS2", colour = "DM", shape = "Type") + 
  scale_colour_manual(values = DM_col) +  # Added BPD colors here
  stat_ellipse(geom = "polygon", aes(group = DM, color = DM, fill = DM), alpha = 0.05)+
  annotate("text", x = -0.1, y = 0.5, label = paste0("stress: ", format(nmds$stress, digits = 4)), hjust = 3)+
  annotate("text", x = -0.1, y = 0.6, label = paste0("P=", format(ano$signif, digits = 4)), hjust = 3)

xx


pdf("nmds_virome_DM_OTU.pdf",width=14,height=10);
xx
dev.off()



### FIGURE 3

#### correlation

df <- read.delim("metadata.txt", stringsAsFactors = FALSE)

cor_result <- cor.test(df$VDR_total_vir, df$log_bdr, method = "s")
correlation_value <- round(cor_result$estimate, 2)
p_value <- round(cor_result$p.value, 4)

# Create the plot
p <- correlation_plot <- ggplot(df, aes(x = log_bdr, y = VDR_total_vir)) +
  geom_point(size = 3) +  # Scatter plot with NPG colors (red)
  geom_smooth(method = "lm", se = TRUE, size = 1, color="yellow") +  # Linear regression line with confidence interval
  labs(title = paste("r =", correlation_value, ", p =", p_value, ""),
       x = "log_bdr", 
       y = "VDR_total_vir"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Title styling
    axis.title.x = element_text(size = 20, face = "bold"),  # X-axis title styling
    axis.title.y = element_text(size = 20, face = "bold"),  # Y-axis title styling
    axis.text.x = element_text(size = 20),  # X-axis text styling
    axis.text.y = element_text(size = 20),  # Y-axis text styling
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Border around plot
    plot.margin = margin(10, 10, 10, 10),  # Adjust plot margins
    legend.position = "none"  # Remove legend
  ) +
  scale_color_npg()  # Define colors for elements

p

pdf("correl_vdr_bdr.pdf",width=6,height=6);
p
dev.off()



library(pheatmap)

data <- read.table("new_metarnaseq.txt", header = TRUE, sep = "\t", row.names = 1)

# Remove SAMPLEID from the data to set as rownames
data_matrix <- data[, -1]
rownames(data_matrix) <- data$SAMPLEID

# Create a heatmap
p <- pheatmap(t(log10(data+1)),
              cluster_rows = F,
              cluster_cols = F,
              cellheight = 60,
              color = colorRampPalette(c("white", "orange", "orange3"))(50),
              fontsize = 10,
              fontsize_row = 20,
              fontsize_col = 20,
              angle_col = 0# Font size for column names
)

pdf(file="Heatmap_density_newmetarnaseq.pdf", height = 5, width = 20)
p
dev.off()


### FIGURE 4


### cluster plot

library(ggplot2)
library(ggforce)
library(wesanderson)
library(ggrepel)
library(tidyverse)

meta = read.delim(file="metadata_r2.txt",header=T, check.names=F) 


cluster_colors <- c(
  "C1" = "blue2",
  "C2" = "purple3",
  "C3" = "yellow4",
  "C4" = "aquamarine4"
)

hulls <- meta %>%
  group_by(cluster) %>%
  slice(chull(PC1, PC2))

cluster_plot <- ggplot(meta, aes(x = PC1, y = PC2, color = cluster, shape = cluster)) +
  geom_point(size = 3, alpha = 0.4) +
  #geom_mark_ellipse(aes(fill = cluster, label = cluster), alpha = 0.1) +
  #xlim(-15,15)+ylim(-15,15)+
  #geom_mark_hull(aes(fill = cluster, label = cluster), alpha = 0.1) +
  geom_polygon(data = hulls, aes(fill = cluster), alpha = 0.25)+
  geom_text_repel(aes(label = UMAP.sample), size = 2.5) +
  #coord_cartesian(xlim = c(-10, 10), ylim = c(-10, 10))+
  scale_color_manual(values = cluster_colors) +
  scale_fill_manual(values = cluster_colors) +
  labs(
    title = "",
    x = "UMAP1",
    y = "UMAP2",
    color = "cluster",
    fill = "cluster"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 12),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 20),
    axis.line = element_line(color = "black", linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )+ coord_flip()

cluster_plot


pdf(file = "umapnew5.pdf", width = 12, height = 8)
cluster_plot
dev.off()


ggsave("umap.png", plot = cluster_plot, width = 10, height = 12, dpi = 300)


### bar plot

## cluster composition bpd/dm

library(ggplot2)
library(ggsci)
library(dplyr)
library(tidyr)

data <- read.delim("metadata.txt", stringsAsFactors = FALSE)

# Calculate the percentage of "yes" and "no" for each cluster
data_summary <- data %>%
  group_by(cluster, Delivrance_mode) %>%
  summarise(count = n()) %>%
  spread(key = Delivrance_mode, value = count, fill = 0) %>%
  mutate(total = CS + VD,
         CS = CS / total * 100,
         VD = VD / total * 100) %>%
  gather(key = "Delivrance_mode", value = "Percent", CS, VD)

# Generate color palette for Interferon
style_colors <- c("CS" = "BROWN", "VD" = "LIGHTGOLDENROD3")

cluster_colors <- c(
  "C1" = "blue2",
  "C2" = "brown",
  "C3" = "gold4",
  "C4" = "purple3"
)

# Create the plot
plot <- ggplot(data_summary, aes(x = cluster, y = Percent, fill = Delivrance_mode)) +
  geom_col(position = "stack", size = 1.2) +  # Border size
  scale_fill_manual(values = style_colors) +
  labs(x = "Cluster", y = "Delivrance mode (%)") +
  theme_classic() +
  geom_text(
    aes(label = round(Percent, 1)), 
    position = position_stack(vjust = 0.5), 
    size = 6, 
    color = "black"
  ) +
  theme(
    legend.text = element_text(size = 14),
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.title.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.text.y = element_text(face = "bold", size = 14),
    axis.text.x = element_text(face = "bold", size = 14, colour = "black")
  )

# Print the plot
print(plot)


pdf(file = "DM_clusters.pdf", height = 5, width = 6)
plot
dev.off()


ggsave("Undertermined_prevalence_Clusters.png", plot = plot, width = 12, height = 10, dpi = 300)

################## new colors budapset
library(ggplot2)
library(ggsci)
library(dplyr)
library(tidyr)
library(wesanderson)

data <- read.delim("metadata.txt", stringsAsFactors = FALSE)

# Calculate the percentage of "yes" and "no" for each cluster
data_summary <- data %>%
  group_by(cluster, BPD) %>%
  summarise(count = n()) %>%
  spread(key = BPD, value = count, fill = 0) %>%
  mutate(total = yes + no,
         yes = yes / total * 100,
         no = no / total * 100) %>%
  gather(key = "BPD", value = "Percent", yes, no)

# Generate color palette for Interferon
style_colors <- c("yes" = "indianred", "no" = "royalblue3")


# Create the plot
plot <- ggplot(data_summary, aes(x = cluster, y = Percent, fill = BPD)) +
  geom_col(position = "stack", size = 1.2) + 
  labs(x = "Cluster", y = "BPD (%)") +
  theme_classic() +
  geom_text(
    aes(label = round(Percent, 1)), 
    position = position_stack(vjust = 0.5), 
    size = 6, 
    color = "black"
  ) +
  theme(
    legend.text = element_text(size = 14),
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.title.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.text.y = element_text(face = "bold", size = 14),
    axis.text.x = element_text(face = "bold", size = 14, colour = "black")
  )

names(wes_palettes)

p <- plot + scale_fill_manual(values = wes_palette(n = 2, name = "Moonrise3"))

# Print the plot
print(p)


pdf(file = "BPD_clustersM.pdf", height = 5, width = 6)
plot
dev.off()

## effect size

####plot
# Load libraries
library(tidyverse)

df <- read.delim("results_all - Copie.txt")

# Create dataframe (copy your table here)


# Create significance label (based on adjusted p-value)
df <- df %>%
  mutate(stars = case_when(
    p_value < 0.0001 ~ "****",
    p_value < 0.001  ~ "***",
    p_value < 0.01   ~ "**",
    p_value < 0.05   ~ "*",
    TRUE ~ ""
  ))

# Sort within each cluster by decreasing effect size
df <- df %>%
  group_by(cluster) %>%
  arrange(desc(effect_size), .by_group = TRUE) %>%
  mutate(variable = factor(variable, levels = unique(variable))) %>%
  ungroup()

# Plot
p <- ggplot(df, aes(x = effect_size, 
                    y = reorder(variable, effect_size), 
                    fill = direction)) +
  geom_col(width = 0.5) +
  geom_text(aes(label = stars), 
            hjust = 0,
            size = 12) +
  facet_grid(~ cluster, scales = "free") +
  theme_classic() +
  labs(x = "Effect Size",
       y = "Variable",
       fill = "Direction") +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "none",
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 12, colour = "black")
  )

p

pdf("cluster_indcator2.pdf",width=10,height=7);
p
dev.off()







