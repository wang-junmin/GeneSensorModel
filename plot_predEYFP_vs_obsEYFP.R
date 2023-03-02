#############################
# <- Author: Junmin Wang -> #
#############################

## load libraries
library(ggplot2)
library(scales)

## read training data
train_data <- read.csv('train_res_tbl.csv')
test_data <- read.csv('test_res_tbl.csv')

## calculate R squared
train_RSS <- sum((train_data$log10predEYFP - train_data$log10obsEYFP)^2)
train_TSS <- sum((train_data$log10obsEYFP - mean(train_data$log10obsEYFP))^2)
train_Rsquared <- 1 - train_RSS / train_TSS

test_RSS <- sum((test_data$log10predEYFP - test_data$log10obsEYFP)^2)
test_TSS <- sum((test_data$log10obsEYFP - mean(test_data$log10obsEYFP))^2)
test_Rsquared <- 1 - test_RSS / test_TSS

## make plots
ggplot(data=train_data, aes(x=10^log10predEYFP, y=10^log10obsEYFP)) +
  geom_point(size=0.005) +
  labs(x="EYFP fitted [MEFL]", y="EYFP observed [MEFL]") +
  theme_bw() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title=element_text(hjust=0.5),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black")) +
  geom_abline(intercept = 0, slope = 1, linetype="dashed",
              color="mediumblue") +
  annotate("text", x=10^6, y=3*10^7,
           label=paste("R^{2} ==", round(train_Rsquared, 2)),
           parse=TRUE) +
  scale_x_continuous(trans="log10",
                     limits=c(7e+4, 2e+8),
                     breaks=c(seq(8*10^4, 9*10^4, 10^4),
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
  scale_y_continuous(trans="log10",
                     limits=c(7e+4, 2e+8),
                     breaks=c(seq(8*10^4, 9*10^4, 10^4),
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
                     )})
ggsave("yhat_vs_y_train.eps", width=4.5, height=4)

ggplot(data=test_data, aes(x=10^log10predEYFP, y=10^log10obsEYFP)) +
  geom_point(size=0.005) +
  labs(x="EYFP predicted [MEFL]", y="EYFP observed [MEFL]") +
  theme_bw() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title=element_text(hjust=0.5),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black")) +
  geom_abline(intercept = 0, slope = 1, linetype="dashed",
              color="mediumblue") +
  annotate("text", x=10^6, y=3*10^7,
           label=paste("R^{2} ==", round(test_Rsquared, 2)),
           parse=TRUE) +
  scale_x_continuous(trans="log10",
                     limits=c(5e+4, 2e+8),
                     breaks=c(seq(6*10^4, 9*10^4, 10^4),
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
  scale_y_continuous(trans="log10",
                     limits=c(5e+4, 2e+8),
                     breaks=c(seq(6*10^4, 9*10^4, 10^4),
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
                     )})
ggsave("yhat_vs_y_test.eps", width=4.5, height=4)

