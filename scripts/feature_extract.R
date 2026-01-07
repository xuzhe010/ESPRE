#!/usr/bin/env Rscript
# scripts/feature_extract.R

suppressMessages({
  library(tidyverse)
  library(optparse)
})

option_list <- list(
  make_option(c("-i", "--input_bed"), type = "character", help = "Input eccDNA BED"),
  make_option(c("-g", "--gc_file"), type = "character", help = "GC content file (from bedtools)"),
  make_option(c("-w", "--window"), type = "character", help = "Window BED file"),
  make_option(c("-o", "--output"), type = "character", help = "Output CSV file"),
  make_option(c("--min_length"), type = "integer", default = 100),
  make_option(c("--max_length"), type = "integer", default = 1000),
  make_option(c("--cutoff"), type = "integer", default = 500)
)

opt <- parse_args(OptionParser(option_list = option_list))

# 1. 读取数据
# 注意：bedtools nuc 输出的格式需要对应
gc_data <- read_tsv(opt$gc_file, col_names = FALSE, show_col_types = FALSE) %>%
    select(X1, X2, X3, X5) %>% # chr, start, end, GC
    rename(chr = X1, start = X2, end = X3, GC = X5)

bed_data <- read_tsv(opt$input_bed, col_names = FALSE, show_col_types = FALSE) %>%
    select(X1, X2, X3) %>% # chr, start, end
    rename(chr = X1, start = X2, end = X3)

# 给数据加上假的 sample ID，方便复用逻辑
sample_id <- "Target_Sample"
bed_data$sample <- sample_id
gc_data$sample <- sample_id

# 2. 读取窗口文件 (1Mb bins)
chr_list <- paste("chr", c(1:22, "X", "Y"), sep = "")
window <- read_tsv(opt$window, col_names = FALSE, show_col_types = FALSE) %>%
    rename(chr.Mb = X1, start.Mb = X2, end.Mb = X3) %>%
    filter(chr.Mb %in% chr_list)

# 3. 映射 eccDNA 到 1Mb 窗口 (这一步需要在 R 里做 overlap 或者假设输入已经是映射好的)
# *注意*：您原本的代码输入是 'bins_1Mbcompartments.csv'，说明在 R 之前已经做过 map。
# 为了简化，这里我们假设输入的 bed 已经是原始 eccDNA bed，我们需要在这里做 overlap。
# 使用 GenomicRanges 包会更严谨，但为了尽量复用您的 tidyverse 逻辑，我们这里做简化：
# *实际上，最高效的方法是在 Python 主程序里用 bedtools intersect -wa -wb 把 eccDNA 映射到 1Mb 窗口*
# 因此，我们假设输入的 --input_bed 已经是 bedtools intersect 之后的结果。
# 格式: [1Mb_chr, 1Mb_start, 1Mb_end, ecc_chr, ecc_start, ecc_end, GC_val]

# 读取整合后的数据 (由 Python 主程序生成)
raw_df <- read_tsv(opt$input_bed, col_names = FALSE, show_col_types = FALSE) %>%
    rename(chr.Mb = X1, start.Mb = X2, end.Mb = X3, 
           chr = X4, start = X5, end = X6, GC = X7) %>%
    mutate(sample = sample_id)

# 4. 长度过滤与分组
df <- raw_df %>%
    mutate(length = end - start) %>%
    filter(length >= opt$min_length & length <= opt$max_length) %>%
    mutate(group = ifelse(length <= opt$cutoff, "short", "long"))

# 5. 统计 Short, Long, GC
df.GC <- df %>%
    group_by(chr.Mb, start.Mb, end.Mb, sample) %>%
    summarise(GC = mean(GC, na.rm = TRUE), .groups = 'drop')

df.short <- df %>% filter(group == "short") %>%
    group_by(chr.Mb, start.Mb, end.Mb, sample) %>% summarise(short = n(), .groups = 'drop')

df.long <- df %>% filter(group == "long") %>%
    group_by(chr.Mb, start.Mb, end.Mb, sample) %>% summarise(long = n(), .groups = 'drop')

# 6. 合并
df.count <- df.GC %>%
    left_join(df.short, by = c("chr.Mb", "start.Mb", "end.Mb", "sample")) %>%
    left_join(df.long, by = c("chr.Mb", "start.Mb", "end.Mb", "sample")) %>%
    mutate(short = replace_na(short, 0), long = replace_na(long, 0)) %>%
    mutate(ratio = ifelse(short == 0 | long == 0, 0, short/long)) %>%
    filter(short > 0, long > 0)

# 7. GC 校正函数 (复用您的核心逻辑)
gc.correct <- function(coverage, bias) {
    if (length(coverage) < 10) return(coverage) # 样本太少不校正
    tryCatch({
        i <- seq(min(bias, na.rm = TRUE), max(bias, na.rm = TRUE), by = 0.001)
        coverage.trend <- loess(coverage ~ bias, span = 0.75)
        coverage.model <- loess(predict(coverage.trend, i) ~ i, span = 0.75)
        coverage.pred <- predict(coverage.model, bias)
        coverage.corrected <- coverage - coverage.pred + median(coverage, na.rm = TRUE)
        return(coverage.corrected)
    }, error = function(e) return(coverage))
}

# 8. 执行校正
if(nrow(df.count) > 0) {
    df.count$ratio.corrected <- gc.correct(df.count$ratio, df.count$GC)
} else {
    stop("Filtered data is empty. Cannot proceed.")
}

# 9. 格式化输出 (Wide format for Python)
# 输出两列: bin_id, value (ratio.corrected) 和 bin_id, value (GC)
output_df <- df.count %>%
    unite(bin_id, chr.Mb, start.Mb, end.Mb, sep = "_") %>%
    select(bin_id, ratio.corrected, GC)

write_csv(output_df, opt$output)
