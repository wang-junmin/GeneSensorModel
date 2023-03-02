#############################
# <- Author: Junmin Wang -> #
#############################

## load libraries
library(ComplexHeatmap)
library(circlize)

## load model
params <- read.csv("best_params_min_err_tbl.csv")

## define functions
calc_drr_coef <- function(params, log10EBFP2, gRNA, uORF) {
  k_0 <- as.numeric(params["k0"])
  k_max <- as.numeric(
    params["c"] * gRNA / (gRNA + params["m"]) +
      params["d"] * uORF
  )
  H <- as.numeric(-log10(exp(1)) * (params["a1"] + 
                                      params["a12"] * log10EBFP2))
  K_D <- (as.numeric(exp(params["a0"] +
                          params["a2"] * log10EBFP2 + 
                          params["a3"] * gRNA +
                          params["a4"] * uORF
  )))
  res <- c(k_0, k_max, H, K_D)
  names(res) <- c("k_0", "k_max", "H", "K_D")
  return(res)
}

## initialize matrices to store results
k_0_mat_uORF_1 <- matrix(nrow = 20, ncol = 7)
k_max_mat_uORF_1 <- matrix(nrow = 20, ncol = 7)
H_mat_uORF_1 <- matrix(nrow = 20, ncol = 7)
K_D_mat_uORF_1 <- matrix(nrow = 20, ncol = 7)

k_0_mat_uORF_2 <- matrix(nrow = 20, ncol = 7)
k_max_mat_uORF_2 <- matrix(nrow = 20, ncol = 7)
H_mat_uORF_2 <- matrix(nrow = 20, ncol = 7)
K_D_mat_uORF_2 <- matrix(nrow = 20, ncol = 7)

#####################################
#### calculation of coefficients ####
#####################################
log10EBFP2_lst <- seq(2.1, 5.9, 0.2)
gRNA_lst <- c(1, 2, 3, 4, 8, 12, 16)
# iterate through all levels of EBFP2
for (i in 1:20) {
  log10EBFP2 <- log10EBFP2_lst[i]
  # iterate through all levels of gRNA
  for (j in 1:7) {
    gRNA <- gRNA_lst[j]
    # calculate for uORF = 1
    res1 <- calc_drr_coef(params, log10EBFP2, gRNA, 1)
    k_0_mat_uORF_1[i, j] <- res1[1]
    k_max_mat_uORF_1[i, j] <- res1[2]
    H_mat_uORF_1[i, j] <- res1[3]
    K_D_mat_uORF_1[i, j] <- res1[4]
    # calculate for uORF = 2
    res2 <- calc_drr_coef(params, log10EBFP2, gRNA, 2)
    k_0_mat_uORF_2[i, j] <- res2[1]
    k_max_mat_uORF_2[i, j] <- res2[2]
    H_mat_uORF_2[i, j] <- res2[3]
    K_D_mat_uORF_2[i, j] <- res2[4]
  }
}

######################
## make plots (K_D) ##
######################
colnames(K_D_mat_uORF_1) <- c(1, 2, 3, 4,
                              8, 12, 16)
rownames(K_D_mat_uORF_1) <- seq(2.1, 5.9, 0.2)

pdf("K_D_uORF_1.pdf", 
    width=7.2, height=5.4)
Heatmap(K_D_mat_uORF_1[nrow(K_D_mat_uORF_1):1, ], 
        cluster_rows = FALSE,
        row_title = expression(Log[10]*'EBFP2'),
        cluster_columns = FALSE, 
        column_title = "sgRNA",
        col = colorRamp2(c(0, 50), c("white", "blue")),
        column_names_rot = 0,
        heatmap_legend_param = list(title = 
                                      expression(K[D]),
                                    at=seq(0, 50, 10),
                                    legend_height=unit(4.5, "cm"),
                                    grid_width = unit(0.6, "cm"),
                                    title_gp=gpar(fontsize = 12),
                                    labels_gp=gpar(fontsize = 12)),
        cell_fun = function(j, i, x, y, width, height, fill) 
        {
          grid.text(sprintf("%.1f", 
                            K_D_mat_uORF_1[nrow(K_D_mat_uORF_1):1, ][i, j]), 
                    x, y, gp = gpar(fontsize = 10))
        })
dev.off()

colnames(K_D_mat_uORF_2) <- c(1, 2, 3, 4,
                              8, 12, 16)
rownames(K_D_mat_uORF_2) <- seq(2.1, 5.9, 0.2)

pdf("K_D_uORF_2.pdf", 
    width=7.2, height=5.4)
Heatmap(K_D_mat_uORF_2[nrow(K_D_mat_uORF_2):1, ], 
        cluster_rows = FALSE,
        row_title = expression(Log[10]*'EBFP2'),
        cluster_columns = FALSE, 
        column_title = "sgRNA",
        col = colorRamp2(c(0, 50), c("white", "blue")),
        column_names_rot = 0,
        heatmap_legend_param = list(title = 
                                      expression(K[D]),
                                    at=seq(0, 50, 10),
                                    legend_height=unit(4.5, "cm"),
                                    grid_width = unit(0.6, "cm"),
                                    title_gp=gpar(fontsize = 12),
                                    labels_gp=gpar(fontsize = 12)),
        cell_fun = function(j, i, x, y, width, height, fill) 
        {
          grid.text(sprintf("%.1f", 
                            K_D_mat_uORF_2[nrow(K_D_mat_uORF_2):1, ][i, j]), 
                    x, y, gp = gpar(fontsize = 10))
        })
dev.off()

####################
## make plots (H) ##
####################
colnames(H_mat_uORF_1) <- c(1, 2, 3, 4,
                            8, 12, 16)
rownames(H_mat_uORF_1) <- seq(2.1, 5.9, 0.2)

pdf("H_uORF_1.pdf", 
    width=7.2, 
    height=5.4)
Heatmap(H_mat_uORF_1[nrow(H_mat_uORF_1):1, ], 
        cluster_rows = FALSE,
        row_title = expression(Log[10]*"EBFP2"),
        cluster_columns = FALSE, 
        column_title = "sgRNA",
        col = colorRamp2(c(-0.35, 0, 0.35), 
                         c("red", "white", "blue")),
        column_names_rot = 0,
        heatmap_legend_param = list(title = 
                                      expression(H),
                                    at=seq(-0.36, 0.36,
                                           0.18),
                                    legend_height=unit(4.5, "cm"),
                                    grid_width = unit(0.6, "cm"),
                                    title_gp=gpar(fontsize = 12),
                                    labels_gp=gpar(fontsize = 12)),
        cell_fun = function(j, i, x, y, width, height, fill) 
        {
          grid.text(sprintf("%.2f", 
                            H_mat_uORF_1[nrow(H_mat_uORF_1):1, ][i, j]), 
                    x, y, gp = gpar(fontsize = 10))
        })
dev.off()

colnames(H_mat_uORF_2) <- c(1, 2, 3, 4,
                              8, 12, 16)
rownames(H_mat_uORF_2) <- seq(2.1, 5.9, 0.2)

pdf("H_uORF_2.pdf", 
    width=7.2, height=5.4)
Heatmap(H_mat_uORF_2[nrow(H_mat_uORF_2):1, ], 
        cluster_rows = FALSE,
        row_title = expression(Log[10]*"EBFP2"),
        cluster_columns = FALSE, 
        column_title = "sgRNA",
        col = colorRamp2(c(-0.35, 0, 0.35), 
                         c("red", "white", "blue")),
        column_names_rot = 0,
        heatmap_legend_param = list(title = 
                                      expression(H),
                                    at=seq(-0.36, 0.36,
                                           0.18),
                                    legend_height=unit(4.5, "cm"),
                                    grid_width = unit(0.6, "cm"),
                                    title_gp=gpar(fontsize = 12),
                                    labels_gp=gpar(fontsize = 12)),
        cell_fun = function(j, i, x, y, width, height, fill) 
        {
          grid.text(sprintf("%.2f", 
                            H_mat_uORF_2[nrow(H_mat_uORF_2):1, ][i, j]), 
                    x, y, gp = gpar(fontsize = 10))
        })
dev.off()

########################
## make plots (k_max) ##
########################
colnames(k_max_mat_uORF_1) <- c(1, 2, 3, 4,
                                8, 12, 16)
rownames(k_max_mat_uORF_1) <- seq(2.1, 5.9, 0.2)

pdf("k_max_uORF_1.pdf", 
    width=7.2, 
    height=5.4)
Heatmap(k_max_mat_uORF_1[nrow(k_max_mat_uORF_1):1, ], 
        cluster_rows = FALSE,
        row_title = expression(Log[10]*"EBFP2"),
        cluster_columns = FALSE, 
        column_title = "sgRNA",
        col = colorRamp2(c(1.6, 3.2), 
                         c("white", "blue")),
        column_names_rot = 0,
        heatmap_legend_param = list(title = 
                                      expression(k[max]),
                                    at=seq(1.6, 3.2, 
                                           0.4),
                                    legend_height=unit(4.5, "cm"),
                                    grid_width = unit(0.6, "cm"),
                                    title_gp=gpar(fontsize = 12),
                                    labels_gp=gpar(fontsize = 12)),
        cell_fun = function(j, i, x, y, width, height, fill) 
        {
          grid.text(sprintf("%.2f", 
                            k_max_mat_uORF_1[nrow(k_max_mat_uORF_1):1, ][i, j]), 
                    x, y, gp = gpar(fontsize = 10))
        })
dev.off()

colnames(k_max_mat_uORF_2) <- c(1, 2, 3, 4,
                                8, 12, 16)
rownames(k_max_mat_uORF_2) <- seq(2.1, 5.9, 0.2)

pdf("k_max_uORF_2.pdf", 
    width=7.2, 
    height=5.4)
Heatmap(k_max_mat_uORF_2[nrow(k_max_mat_uORF_2):1, ], 
        cluster_rows = FALSE,
        row_title = expression(Log[10]*"EBFP2"),
        cluster_columns = FALSE,
        column_title = "sgRNA",
        col = colorRamp2(c(1.6, 3.2), 
                         c("white", "blue")),
        column_names_rot = 0,
        heatmap_legend_param = list(title = 
                                      expression(k[max]),
                                    at=seq(1.6, 3.2,
                                           0.4),
                                    legend_height=unit(4.5, "cm"),
                                    grid_width = unit(0.6, "cm"),
                                    title_gp=gpar(fontsize = 12),
                                    labels_gp=gpar(fontsize = 12)),
        cell_fun = function(j, i, x, y, width, height, fill) 
        {
          grid.text(sprintf("%.2f", 
                            k_max_mat_uORF_2[nrow(k_max_mat_uORF_2):1, ][i, j]), 
                    x, y, gp = gpar(fontsize = 10))
        })
dev.off()

