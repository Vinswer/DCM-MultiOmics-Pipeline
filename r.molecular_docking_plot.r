# 清除当前环境所有变量并返回内存使⽤情况1
rm(list = ls()); gc()

#定义变量output和ORIGINAL_DIR
ORIGINAL_DIR <- "/data/nas1/yangyudong_OD/project/MD/844_1624/"
output <- file.path(ORIGINAL_DIR, "./")

#确保输出目录存在
if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}

#切换至项目目录
setwd(ORIGINAL_DIR) 


a <- file.path("./analysis_MAP2K1_MAP2K1_HBIN041139/")
# b <- file.path("./analysis_RPL14_Mesalamine/")
# c <- file.path("./analysis_PTGS2_Benzo_a_pyrene/")
# d <- file.path("./analysis_Prkcd_bassianins/")
# e <- file.path("./analysis_Spp1_quercetin/")

# 提取 analysis_ 后的内容
extract_name <- function(path) {
  sub("^\\.\\/analysis_[^_]+_([^/]+)/?$", "\\1", path)
}

an <- extract_name(a)

# 定义公共绘图主题------------
theme_pub <- function(legend_nrow = 1, x_expand = c(0, 0)) {       # legend_nrow: 图例行数，x_expand: 曲线与左右边框距离
  list(
    theme_bw(base_size = 12, base_family = "Helvetica"),            # 基础主题：白底黑框，字号12，Helvetica字体
    theme(
      panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.8),  # 绘图区黑色边框
      panel.grid.major  = element_blank(),                          # 去除主网格线
      panel.grid.minor  = element_blank(),                          # 去除次网格线
      axis.title        = element_text(size = 12, face = "bold"),   # 坐标轴标题：12号加粗
      axis.text         = element_text(size = 10, color = "black"), # 坐标轴刻度文字：10号黑色
      axis.ticks        = element_line(color = "black", linewidth = 0.4),  # 坐标轴刻度线：黑色细线
      axis.ticks.length = unit(0.15, "cm"),                         # 刻度线长度
      legend.position   = "top",                                    # 图例置于顶部
      legend.title      = element_text(size = 10, face = "bold"),   # 图例标题：10号加粗
      legend.text       = element_text(size = 9),                   # 图例文字：9号
      legend.key.width  = unit(1.2, "cm"),                          # 图例色块宽度
      legend.key.height = unit(0.4, "cm"),                          # 图例色块高度
      legend.margin     = margin(2, 6, 2, 6),                       # 图例内边距（上右下左）
      plot.margin       = margin(10, 15, 10, 10)                    # 图形外边距（上右下左）
    ),
    scale_x_continuous(expand = expansion(mult = x_expand)),        # 曲线与左右边框的距离，默认紧贴
    guides(color = guide_legend(nrow = legend_nrow, byrow = TRUE))  # 图例按指定行数排列，逐行填充
  )
}
# 公共绘图主题用法：
# p1 <- ggplot(...) + ... + theme_pub() 


# 自动配色脚本---------------
# 自动收集已定义的分组名
name_vars <- paste0(letters, "n")
name_vars[which(letters == "i")] <- "in1"   # 把 "in" 替换为 "in1"

names_list <- mget(name_vars, envir = environment(), ifnotfound = list(NULL))
group_names <- unlist(Filter(Negate(is.null), names_list))
cat("Groups:", paste(group_names, collapse = ", "), "\n")
cat("Count:", length(group_names), "\n")

library(RColorBrewer)
# 自动生成配色的函数
generate_pub_colors <- function(n, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)          # 可选：固定种子保证可复现
  
  if (n == 1) {
    # 随机选一个高饱和色
    hcl(h = runif(1, 0, 360), c = 50, l = 80) # 鲜艳明亮：130、65；淡雅清新：50、80；深沉高对比：90、55；莫兰迪：35、75（c饱和度，l亮度）
  } else {
    # 在色环上均匀取 n 个点，保证最大间隔（高区分度）
    # 然后加一个随机起始偏移（随机性）
    offset <- runif(1, 0, 360)                # 随机旋转起点
    hues <- (seq(0, 360 - 360/n, length.out = n) + offset) %% 360
    hues <- sample(hues)                      # 随机打乱顺序
    hcl(h = hues, c = 50, l = 80)             # 鲜艳明亮：130、65；淡雅清新：50、80；深沉高对比：90、55；莫兰迪：35、75（c饱和度，l亮度）
  }
}
# 每次运行颜色不同
# group_colors <- generate_pub_colors(length(group_names))
# 固定种子，保证结果可复现
group_colors <- generate_pub_colors(length(group_names), seed = 844)
# 映射分组和颜色
color_mapping <- setNames(group_colors, group_names)
# 查看配色
print(color_mapping)


library(dplyr)
library(ggplot2)



# 读取数据
read_xvg <- function(file_path) {
  # 读取文件所有行
  lines <- readLines(file_path)
  
  # 去除注释行(以#或@开头的行)
  data_lines <- lines[!grepl("^[#@]", lines)]
  
  # 从注释行中提取列名
  col_names <- c()
  
  # 提取x轴和y轴标签
  x_label <- lines[grepl("@\\s*xaxis\\s+label", lines)]
  if (length(x_label) > 0) {
    x_label <- gsub('.*label\\s*"(.*?)".*', "\\1", x_label[1])
    col_names <- c(col_names, x_label)
  }
  
  y_label <- lines[grepl("@\\s*yaxis\\s+label", lines)]
  if (length(y_label) > 0) {
    y_label <- gsub('.*label\\s*"(.*?)".*', "\\1", y_label[1])
    col_names <- c(col_names, y_label)
  }
  
  # 如果没有找到标签，使用默认列名
  if (length(col_names) == 0) {
    col_names <- c("X", "Y")
  }
  
  # 将数据读入数据框
  data <- read.table(text = data_lines, header = FALSE)
  names(data) <- col_names
  
  return(data)
}

# RMSD----------------
# 读取数据
a_rmsd <- read_xvg(file.path(a, "rmsd.xvg"))
# b_rmsd <- read_xvg(file.path(b, "rmsd.xvg"))
# c_rmsd <- read_xvg(file.path(c, "rmsd.xvg"))
# d_rmsd <- read_xvg(file.path(d, "rmsd.xvg"))
# e_rmsd <- read_xvg(file.path(e, "rmsd.xvg"))
# f_rmsd <- read_xvg("./analysis_STAT3_Regorafenib/rmsd.xvg")

# 为每个数据集添加基因名称标识
a_rmsd$Gene <- an
# b_rmsd$Gene <- bn
# c_rmsd$Gene <- cn
# d_rmsd$Gene <- dn
# e_rmsd$Gene <- en
# f_rmsd$Gene <- "STAT3_Regorafenib"

# 合并三个数据集
combined_rmsd <- rbind(a_rmsd)#, b_rmsd)#, c_rmsd)#, d_rmsd, e_rmsd)#, f_rmsd)

# 绘制折线图
p1 <- ggplot(combined_rmsd, aes(x = `Time (ns)`, y = `RMSD (nm)`, color = Gene)) +
  geom_line(linewidth = 0.6, alpha = 0.85) +
  labs(title = "",
       x = "Time (ns)",
       y = "RMSD (nm)",
       color = "Group") +
  scale_color_manual(values = color_mapping) +
  # 发表级主题
  theme_pub()
pdf(file = file.path(output, "01.RMSD_plot.pdf"), width = 8, height = 4)  # 打开PDF设备，设置宽度和高度（英寸）
print(p1)  # 绘制图形
dev.off()  # 关闭设备，保存文件

png(file = file.path(output, "01.RMSD_plot.png"), 
    width = 8, 
    height = 4, 
    units = "in",  # 单位：英寸
    res = 300)    # 分辨率：300 dpi
print(p1)  # 绘制图形
dev.off()  # 关闭设备，保存文件

# RMSF--------------------
# 读取数据
a_rmsf <- read_xvg(file.path(a, "rmsf.xvg"))
# b_rmsf <- read_xvg(file.path(b, "rmsf.xvg"))
# c_rmsf <- read_xvg(file.path(c, "rmsf.xvg"))
# d_rmsf <- read_xvg(file.path(d, "rmsf.xvg"))
# e_rmsf <- read_xvg(file.path(e, "rmsf.xvg"))
# f_rmsf <- read_xvg("./analysis_STAT3_Regorafenib/rmsf.xvg")

# 为每个数据集添加基因名称标识
a_rmsf$Gene <- an
# b_rmsf$Gene <- bn
# c_rmsf$Gene <- cn
# d_rmsf$Gene <- dn
# e_rmsf$Gene <- en
# f_rmsf$Gene <- "STAT3_Regorafenib"

# 氨基酸重新编号，从1开始
a_rmsf$Residue <- seq_len(nrow(a_rmsf))
# b_rmsf$Residue <- seq_len(nrow(b_rmsf))
# c_rmsf$Residue <- seq_len(nrow(c_rmsf))
# d_rmsf$Residue <- seq_len(nrow(d_rmsf))
# e_rmsf$Residue <- seq_len(nrow(e_rmsf))
# f_rmsf$Residue <- seq_len(nrow(f_rmsf))

# 合并三个数据集
combined_rmsf <- rbind(a_rmsf)#, b_rmsf)#, c_rmsf)#, d_rmsf, e_rmsf)#, f_rmsf)

# 绘制折线图
p2 <- ggplot(combined_rmsf, aes(x = Residue, y = `(nm)`, color = Gene)) +
  geom_line(linewidth = 0.8) +
  labs(title = "",
       x = "Residue",
       y = "RMSF (nm)",
       color = "Group") +
  scale_color_manual(values = color_mapping) +
  # 发表级主题
  theme_pub()
pdf(file = file.path(output, "02.RMSF_plot.pdf"), width = 8, height = 4)  # 打开PDF设备，设置宽度和高度（英寸）
print(p2)  # 绘制图形
dev.off()  # 关闭设备，保存文件

png(file = file.path(output, "02.RMSF_plot.png"), 
    width = 8, 
    height = 4, 
    units = "in",  # 单位：英寸
    res = 300)    # 分辨率：300 dpi
print(p2)  # 绘制图形
dev.off()  # 关闭设备，保存文件

# energy--------------------
read_multi_xvg <- function(file_path) {
  # 读取文件所有行
  lines <- readLines(file_path)
  
  # 去除注释行(以#或@开头的行)
  data_lines <- lines[!grepl("^[#@]", lines)]
  
  # 从注释行中提取列名
  col_names <- c()
  
  # 提取x轴标签
  x_label <- lines[grepl("@\\s*xaxis\\s+label", lines)]
  if (length(x_label) > 0) {
    x_label <- gsub('.*label\\s*"(.*?)".*', "\\1", x_label[1])
    col_names <- c(col_names, x_label)
  } else {
    col_names <- c(col_names, "Time")
  }
  
  # 提取y轴标签
  y_label <- lines[grepl("@\\s*yaxis\\s+label", lines)]
  if (length(y_label) > 0) {
    y_label <- gsub('.*label\\s*"(.*?)".*', "\\1", y_label[1])
  } else {
    y_label <- "Value"
  }
  
  # 提取图例信息作为列名
  legends <- lines[grepl("@\\s*s\\d+\\s+legend", lines)]
  if (length(legends) > 0) {
    legend_names <- gsub('.*legend\\s*"(.*?)".*', "\\1", legends)
    col_names <- c(col_names, legend_names)
  } else {
    # 如果没有图例信息，使用默认列名
    n_columns <- length(strsplit(data_lines[1], "\\s+")[[1]])
    if (n_columns > length(col_names)) {
      default_names <- paste0(y_label, "_", seq_len(n_columns - length(col_names)))
      col_names <- c(col_names, default_names)
    }
  }
  
  # 将数据读入数据框
  data <- read.table(text = data_lines, header = FALSE)
  
  # 确保列名数量与数据列数匹配
  if (ncol(data) > length(col_names)) {
    col_names <- c(col_names, paste0("V", seq(length(col_names)+1, ncol(data))))
  } else if (ncol(data) < length(col_names)) {
    col_names <- col_names[1:ncol(data)]
  }
  
  names(data) <- col_names
  
  return(data)
}
# 读取数据
a_en <- read_multi_xvg(file.path(a, "energy_total.xvg"))
# b_en <- read_multi_xvg(file.path(b, "energy_total.xvg"))
# c_en <- read_multi_xvg(file.path(c, "energy_total.xvg"))
# d_en <- read_multi_xvg(file.path(d, "energy_total.xvg"))
# e_en <- read_multi_xvg(file.path(e, "energy_total.xvg"))
# f_en <- read_multi_xvg("./analysis_STAT3_Capecitabine/energy_total.xvg")


# 为每个数据集添加基因名称标识
a_en$Gene <- an
# b_en$Gene <- bn
# c_en$Gene <- cn
# d_en$Gene <- dn
# e_en$Gene <- en
# f_en$Gene <- "STAT3_Capecitabine"

# 合并三个数据集
combined_en <- rbind(a_en)#, b_en)#, c_en)#, d_en, e_en)#, f_en)

# ps转ns
combined_en$`Time (ps)` <- combined_en$`Time (ps)` / 1000
# 重命名列
colnames(combined_en)[colnames(combined_en) == "Time (ps)"] <- "Time (ns)"

# 绘制折线图
library(ggrepel)
library(dplyr)

# 为每个分组选取标注位置
label_data <- combined_en %>%
  group_by(Gene) %>%
  slice(round(n() * 0.8)) %>%
  ungroup()

p3 <- ggplot(combined_en, aes(x = `Time (ns)`, y = `Total Energy`, color = Gene)) +
  geom_line(linewidth = 0.8, alpha = 0.7) +

  labs(title = "",
       x = "Time (ns)",
       y = "Energy (kcal/mol)",
       color = "Group") +
  scale_color_manual(values = color_mapping) +
  # 发表级主题
  theme_pub()
pdf(file = file.path(output, "03.energy_plot.pdf"), width = 8, height = 4)  # 打开PDF设备，设置宽度和高度（英寸）
print(p3)  # 绘制图形
dev.off()  # 关闭设备，保存文件

png(file = file.path(output, "03.energy_plot.png"), 
    width = 8, 
    height = 4, 
    units = "in",  # 单位：英寸
    res = 300)    # 分辨率：300 dpi
print(p3)  # 绘制图形
dev.off()  # 关闭设备，保存文件

# Hydrogen bond-------------
# 读取数据
a_hb <- read_xvg(file.path(a, "hbond_num.xvg"))
# b_hb <- read_xvg(file.path(b, "hbond_num.xvg"))
# c_hb <- read_xvg(file.path(c, "min_dist_protein_ligand.xvg"))
# d_hb <- read_xvg(file.path(d, "hbond_num.xvg"))
# e_hb <- read_xvg(file.path(e, "hbond_num.xvg"))
# f_hb <- read_xvg("./analysis_STAT3_Fruquintinib/hbond_num.xvg")

# 为每个数据集添加基因名称标识
a_hb$Gene <- an
# b_hb$Gene <- bn
# c_hb$Gene <- cn
# d_hb$Gene <- dn
# e_hb$Gene <- en
# f_hb$Gene <- "SORD_Aflatoxin_B1"

# 合并三个数据集
combined_hb <- rbind(a_hb)#, b_hb)#, c_hb)#, d_hb, e_hb)#, f_hb)

# ps转ns
combined_hb$`Time (ps)` <- combined_hb$`Time (ps)` / 1000
# 重命名列
colnames(combined_hb)[colnames(combined_hb) == "Time (ps)"] <- "Time (ns)"

# 绘制折线图
p4 <- ggplot(combined_hb, aes(x = `Time (ns)`, y = `Number`, color = Gene)) +
  geom_line(linewidth = 0.8) +  # 绘制折线
  labs(title = "",  # 标题
       x = "Time (ns)",          # x轴标签
       y = "Hydrogen bond numbers",          # y轴标签
       color = "Group") +         # 图例标题
  scale_color_manual(values = color_mapping) +
  # 发表级主题
  theme_pub()
pdf(file = file.path(output, "04.hbonds_plot.pdf"), width = 8, height = 4)  # 打开PDF设备，设置宽度和高度（英寸）
print(p4)  # 绘制图形
dev.off()  # 关闭设备，保存文件

png(file = file.path(output, "04.hbonds_plot.png"), 
    width = 8, 
    height = 4, 
    units = "in",  # 单位：英寸
    res = 300)    # 分辨率：300 dpi
print(p4)  # 绘制图形
dev.off()  # 关闭设备，保存文件

