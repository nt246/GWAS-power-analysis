## gwas power analysis
## for nsf postdoc proposal
## fall 2020

source("gwas-power-master/power_calc_functions.R")
library(tidyverse)

# power analysis
# from https://github.com/kaustubhad/gwas-power
# beta = effect size
pow <- power_beta_maf(beta = (5:40)/50, maf = (1:5)/10, n = 600, pval=5E-6) # iterating over beta and maf values
p1 <- power_plot(pow, "effect size", "maf")

pow2 <- power_beta_het(beta = (5:20)/50, het = (10:50)/100, n = 400, pval=5E-6) # iterating over beta and het values
p2 <- power_plot(pow2, "effect size", "het")


pow_var = power_n_hsq(n = (1:10)*100, qsq = (1:25)/100, pval=5E-8)
p3 <- power_plot(pow_var, "sample size", "% variance explained")



# transform into dataframes
pow_df <- rownames_to_column(as.data.frame(pow), var="effect_size") %>%
  as_tibble() %>% 
  pivot_longer(cols=c(2:ncol(pow)+1),names_to="maf") %>%
  mutate(effect_size=as.numeric(effect_size))

pow_df2 <- rownames_to_column(as.data.frame(pow2), var="effect_size") %>%
  as_tibble() %>% 
  pivot_longer(cols=c(2:ncol(pow2)),names_to="het") %>%
  mutate(effect_size=as.numeric(effect_size))

# heatmaps
heat1 <- ggplot(pow_df, aes(x=maf,y=effect_size,fill=value)) +
  geom_tile() +
  scale_fill_gradient(low="yellow",high="red",limits=c(0,1)) +
  xlab("Minor Allele Frequency") +
  ylab("Effect Size") +
  theme_bw() +
  theme(axis.title = element_text(size=16),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14,color="white"),
        panel.grid=element_blank(),
        panel.border=element_blank())

heat2 <- ggplot(pow_df2, aes(x=het,y=effect_size,fill=value)) +
  geom_tile() +
  scale_fill_gradient(low="yellow",high="red",limits=c(0,1)) +
  xlab("Heterozygous Genotype\nFrequency") +
  ylab("Effect Size") +
  theme_bw() +
  theme(axis.title = element_text(size=16),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14,color="white"),
        panel.grid=element_blank(),
        panel.border=element_blank())

lines <- ggplot(data=pow_df,aes(x=effect_size,y=value)) +
  geom_line(aes(col=maf), size=2) +
  geom_point(aes(col=maf), size=4) +
  xlab("Effect Size") +
  ylab("GWAS Power") +
  theme_bw() +
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=14))

# combine plots together
heat <- ggarrange(heat1,heat2,ncol=1,common.legend=TRUE, legend="bottom")
ggarrange(heat,lines,nrow=1,
          widths=c(1,2))
