## http://www.sthda.com/english/wiki/cluster-analysis-in-r-unsupervised-machine-learning
## http://www.sthda.com/english/wiki/static-and-interactive-heatmap-in-r-unsupervised-machine-learning

library("cluster")
library("factoextra") # install.packages("factoextra")
library("cluster")

setwd("~/Desktop/pharmacogenomics/scripts/")

my_data <- read.csv("../data/170824_zscores.csv", sep = ",", row.names = 1)
rownames(my_data)<-sub("_human", "", rownames(my_data))
# Remove any missing value (i.e, NA values for not available)
my_data <- na.omit(my_data)
# Scale variables
my_data <- scale(my_data)
row.names(my_data) <- toupper(row.names(my_data))

# View the firt 3 rows
# my_data <- t(my_data)
head(my_data, n = 3)

res.dist <- get_dist(my_data, stand = TRUE, method = "pearson")
fviz_dist(res.dist,
          gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

fviz_nbclust(my_data, kmeans, method = "gap_stat")

km.res <- kmeans(my_data, 5, nstart = 25)
# Visualize
fviz_cluster(km.res, data = my_data, frame.type = "convex")+
  theme_minimal()

# Heatmaps
if (!require("devtools")) install.packages("devtools")
devtools::install_github("rstudio/d3heatmap")
if (!require("devtools")) install.packages("devtools")
devtools::install_github("jokergoo/ComplexHeatmap")
library("gplots") # install.packages("gplots")
library("circlize")
library("d3heatmap") #
library("ComplexHeatmap") # install.packages("grid")
my_palette <- colorRampPalette(c("darkgreen", "white", "darkred"))(n = 1000)
hmap1 <- heatmap.2(my_data, scale = "column", col = my_palette,
          trace = "none", density.info = "none", margins = c(2, 8),
          key=FALSE, keysize=1.0, symkey=FALSE,
          colsep=1:10,
          sepcolor='white', sepwidth=0.01,
          cexRow=0.9,cexCol=2,
          lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lhei=c(0.25, 6, 0.25 ),
          srtCol=90)
dev.off()
hmap1
ggsave(file="~/Downloads/170829_clustering.eps", plot=hmap1, width=12, height=14)


# interactive
d3heatmap(scale(my_data), colors = "RdGr",
          k_row = 4, k_col = 2)

library(dendextend)
# order for rows
Rowv  <- my_data %>% scale %>% dist %>% hclust %>% as.dendrogram %>%
  set("branches_k_color", k = 3) %>% set("branches_lwd", 1.2) %>%
  ladderize
# Order for columns
# We must transpose the data
Colv  <- my_data %>% scale %>% t %>% dist %>% hclust %>% as.dendrogram %>%
  set("branches_k_color", k = 2, value = c("orange", "blue")) %>%
  set("branches_lwd", 1.2) %>%
  ladderize

heatmap.2(scale(my_data), scale = "none", col = bluered(100),
          Rowv = Rowv, Colv = Colv,
          trace = "none", density.info = "none")

if (!require("devtools")) install.packages("devtools")
devtools::install_github("jokergoo/ComplexHeatmap")
library("ComplexHeatmap")

library("RColorBrewer")
df <- scale(my_data)
Heatmap(df, name = "scores",
        col = colorRamp2(c(-2, 0, 2), brewer.pal(n=3, name="RdBu")))

library(dendextend)
row_dend = hclust(dist(my_data)) # row clustering
row_dend = color_branches(row_dend, k = 4)
Heatmap(my_data, name = "scores",
        cluster_rows = row_dend, split = 4)
