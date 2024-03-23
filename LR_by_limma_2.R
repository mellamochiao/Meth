####author: YAO

### Installation of GEOquery (if without GEOquery)

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Install remotes from CRAN
install.packages("remotes")

## Approach 1: Try installing from my GitHub repository
if (!requireNamespace("GEOquery", quietly = TRUE)) {
  tryCatch({
    library(remotes)
    install_github("curryhank08/GEOquery_with_modifiable_timeout_seconds", force = TRUE)
    library(GEOquery)
  }, error = function(e1) {
    ## Approach 2: If approach 1 fails, try installing from GEOquery's author's GitHub repository
    tryCatch({
      install_github("seandavi/GEOquery")
      library(GEOquery)
    }, error = function(e2) {
      stop("Failed to install GEOquery. Error messages:\nApproach 1:", e1, "\nApproach 2:", e2)
    })
  })
}

### Download dataset from NCBI GEO

# Load modified GEOquery or GEOquery 2.70 version
library(GEOquery)
# Setting the max timeout_seconds
options(timeout=100000)
# Check the input timeout_seconds
getOption("timeout")

# Download GSE30870 by a fuction getGEO() from GEOquery package.
gse30870 <- getGEO("GSE30870", GSEMatrix = TRUE, AnnotGPL = TRUE)
gse30870_matrix <- gse30870[[1]]
data <- exprs(gse30870_matrix)
head(data)


### Regressing age on methylation level for all CpGs by limma

library(limma)
## LR by limma
x1 <- gse30870_matrix$source_name_ch1
design_30870 <- model.matrix(~ x1)
colnames(design_30870)[colnames(design_30870) == "x1Nonagenarians"] <- "x1"
fit <- lmFit(data, design_30870)
fit <- eBayes(fit)
result_30870bylimma_x1 <- topTable(fit, coef="x1", number = Inf, adjust.method = "BH", sort.by = "P")
(result_30870bylimma_x1)

# 使用order函數按照P值由大到小排序
sorted_result <- result_30870bylimma_x1[order(-result_30870bylimma_x1$P.Value), ]
# 顯示排序後的前幾行
(sorted_result)

threshold <- 0.01# 設置閾值

# 篩選P值在閾值範圍內的資料
filtered_result <- sorted_result[sorted_result$P.Value >= (threshold - 0.05) & sorted_result$P.Value <= (threshold + 0.05), ]

# 顯示篩選後的結果
print(filtered_result)


significant_points <- subset(result_30870bylimma_x1, adj.P.Val < 0.05)
print(significant_points)
nrow(significant_points)
