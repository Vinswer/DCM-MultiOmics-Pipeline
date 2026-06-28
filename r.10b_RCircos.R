cat("[INFO] r.10b_RCircos.R 开始\n")

if (!requireNamespace("RCircos", quietly = TRUE)) {
  cat("[ERROR] RCircos 包未安装，请先安装\n")
  quit(status = 1)
}
library(RCircos)

WORK    <- ""
out_dir <- file.path(WORK, "09_RCircos")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

gene_positions <- data.frame(
  Chromosome = c("chr1",      "chr1",        "chr15"),
  chromStart = c(27341741L,   161158000L,    66427484L),
  chromEnd   = c(27360274L,   161162000L,    66568969L),
  Gene.Name       = c("FCN3",      "FCER1G",      "MAP2K1"),
  stringsAsFactors = FALSE
)
cat("[INFO] 基因位置:\n")
print(gene_positions)


data(UCSC.HG19.Human.CytoBandIdeogram)
cyto_info <- UCSC.HG19.Human.CytoBandIdeogram


RCircos.Set.Core.Components(
  cyto.info    = cyto_info,

  tracks.inside  = 2,
  tracks.outside = 0
)

rcircos.params <- RCircos.Get.Plot.Parameters()
rcircos.params$text.size <- 1.2
RCircos.Reset.Plot.Parameters(rcircos.params)


draw_rcircos <- function(out_file, width_px = 2400, height_px = 2400) {
  png(out_file, width = width_px, height = height_px, res = 200, bg = "white")
  RCircos.Set.Core.Components(
    cyto.info    = cyto_info,
    chr.exclude  = c("chrX", "chrY"),
    tracks.inside  = 2,
    tracks.outside = 0
  )
  rcircos.params <- RCircos.Get.Plot.Parameters()
  rcircos.params$text.size <- 1.2
  RCircos.Reset.Plot.Parameters(rcircos.params)


  RCircos.Set.Plot.Area()


  RCircos.Chromosome.Ideogram.Plot()


  tryCatch({
    RCircos.Gene.Connector.Plot(
      genomic.data = gene_positions,
      track.num    = 1,
      side         = "in"
    )
  }, error = function(e) {
    cat("[WARN] Connector plot 失败:", conditionMessage(e), "\n")
  })


  tryCatch({
    RCircos.Gene.Name.Plot(
      gene.data  = gene_positions,
      name.col   = 4,
      track.num  = 2,
      side       = "in"

    )
  }, error = function(e) {
    cat("[WARN] Gene name plot 失败:", conditionMessage(e), "\n")
  })


  title("Chromosomal Locations of Key Biomarkers\n(FCN3, FCER1G, MAP2K1)",
        cex.main = 1.3, font.main = 2, line = -1)
  dev.off()
}


png_file <- file.path(out_dir, "chromosome_location_RCircos.png")
draw_rcircos(png_file)
cat("[INFO] 已保存 chromosome_location_RCircos.png\n")


tryCatch({
  pdf_file <- file.path(out_dir, "chromosome_location_RCircos.pdf")
  pdf(pdf_file, width = 10, height = 10)
  RCircos.Set.Core.Components(
    cyto.info    = cyto_info,
    chr.exclude  = c("chrX", "chrY"),
    tracks.inside  = 2,
    tracks.outside = 0
  )
  rcircos.params <- RCircos.Get.Plot.Parameters()
  rcircos.params$text.size <- 1.2
  RCircos.Reset.Plot.Parameters(rcircos.params)
  RCircos.Set.Plot.Area()
  RCircos.Chromosome.Ideogram.Plot()
  tryCatch(RCircos.Gene.Connector.Plot(gene_positions, 1, "in"), error = function(e) NULL)
  tryCatch(RCircos.Gene.Name.Plot(gene_positions, 4, 2, "in"),   error = function(e) NULL)
  par(mar = c(2, 2, 2, 2))
  mtext("Chromosomal Locations of Key Biomarkers (FCN3, FCER1G, MAP2K1)",
        side = 3, line = 0.5, cex = 1.3, font = 2)
  dev.off()
  cat("[INFO] 已保存 chromosome_location_RCircos.pdf\n")
}, error = function(e) {
  cat("[WARN] PDF 保存失败:", conditionMessage(e), "\n")
})

cat("[DONE]\n")
