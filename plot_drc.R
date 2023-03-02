#############################
# <- Author: Junmin Wang -> #
#############################

## load libraries
library(ggplot2)
library(scales)

## read data
train_data <- read.csv("train_res_tbl.csv")
test_data <- read.csv("test_res_tbl.csv")
best_params <- read.csv("best_params_min_err_tbl.csv")
  
## merge training and testing data
all_data <- rbind(train_data, test_data)

## plot mKate vs EYFP
for (gRNA_level in c(1, 2, 3, 4, 8, 12, 16)) {
  for (uORF_level in c(1, 2)) {
    p <- ggplot(subset(all_data, 
                       uORF==uORF_level & sgRNA==gRNA_level)) + 
      geom_line(aes(x=10^log10mKate, 
                    y=10^log10predEYFP, 
                    color=log10EBFP2,
                    group=log10EBFP2),
                linetype="dashed") +
      scale_color_gradientn(colours = rainbow(20),
                            breaks=c(3, 4, 5, 6),
                            labels=c(expression(10^3,
                                                10^4,
                                                10^5,
                                                10^6))) +
      labs(color="EBFP2\n[MEFL]", x="mKate [MEFL]", y="EYFP [MEFL]") +
      ggtitle(sprintf("[%sxuORF]-%sxsgRNA", uORF_level, 
                      gRNA_level)) +
      scale_y_continuous(trans="log10", 
                         limits=c(5e+4, 2e+8),
                         breaks=c(seq(5*10^4, 9*10^4, 10^4),
                                  seq(10^5, 9*10^5, 10^5),
                                  seq(10^6, 9*10^6, 10^6),
                                  seq(10^7, 9*10^7, 10^7),
                                  seq(10^8, 2*10^8, 10^8)),
                         labels=function(x) 
                         {ifelse(x %in% c(10^seq(5, 8, 1)),
                                 parse(text=gsub("[+]", "", 
                                                 gsub("1e", 
                                                      "10^", 
                                                      scientific_format()(x)
                                                 ))), ""
                         )}) +
      scale_x_continuous(trans="log10", 
                         limits=c(1e+3, 1e+6),
                         breaks=c(seq(10^3, 9*10^3, 10^3),
                                  seq(10^4, 9*10^4, 10^4),
                                  seq(10^5, 10^6, 10^5)
                         ),
                         labels=function(x) 
                         {ifelse(x %in% c(10^seq(3, 6, 1)),
                                 parse(text=gsub("[+]", "", 
                                                 gsub("1e", 
                                                      "10^", 
                                                      scientific_format()(x)
                                                 ))), ""
                         )}) +
      theme_bw() +
      theme(axis.line = element_line(color='black'),
            plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            plot.title=element_text(hjust=0.5),
            axis.text.x=element_text(color="black"),
            axis.text.y=element_text(color="black"),
            legend.position="right")
    ggsave(plot=p, 
           filename=sprintf("drc_gRNA_%s_uORF_%s.eps", 
                            gRNA_level,
                            uORF_level),
           width=4, 
           height=3)
  }
}

