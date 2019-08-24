
# Ben Stear
# Meta Analysis R code to reproduce Forest Plot
# 

# The R mirrors was not working correctly so I had to install from the cran website directly.
install.packages("metafor", dependencies=TRUE, repos="http://cran.rstudio.com/")
install.packages("meta", dependencies=TRUE, repos="http://cran.rstudio.com/")


setwd("~/Desktop/meta_analysis")

# Must load libraries after you've installed packages
library(readxl) # load in excel library reader
library(meta)
library(metafor)

data <- read_excel("data.xlsx")

# First step in meta-analysis is to calculate pooled effect sizes,
# meaning combine all effect sizes into one pooled effect size
# for this we need effect sizes and variances

raw <-metacont(N_post,
             mean_post,
             std_post,
             N_pre,
             mean_pre,
             std_pre,
             data=data,
             studlab=paste(Author),
             comb.fixed = FALSE,
             comb.random = TRUE,
             method.tau = "SJ",
             hakn = TRUE,
             prediction=TRUE,
             sm="SMD",
             title= 'The Effects of PA on QOL in the SCI Community',
             label.e='Physical Activity',
             label.c='Control')

# Create Forest Plot
forest(raw)

