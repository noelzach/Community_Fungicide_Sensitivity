``` r
rm(list = ls(all=TRUE)) # removes all variables in the global environment so you start fresh

Sys.time() # prints out the time and date you ran the code
```

    ## [1] "2016-04-28 12:27:00 EDT"

``` r
options(scipen = 999) # stops anything from being in scientific notation

op <- par(oma=c(5,7,1,1))
```

### Assay Development

Here we are analysing a dataset where we did different volumes of macerated agar in the wells of a 96 well plate to see which was going to be the closest to the actual poison plate data.

The following loop will generate EC50's for each isolate, agar volume, trial combination and save that output into a data frame called agardil, to use later for an ANOVA

``` r
library(drc)
agar.dil <- agar.dil[!agar.dil$is == "BLANK",]
agar.dil <- agar.dil[!agar.dil$is == "NTC",]
agar.dil$is <- factor(agar.dil$is)
agardil <- NULL
nm <- levels(agar.dil$is)
for (t in 1:2){
  for (i in seq_along(nm)){
agardil.drc <- drm(100*relgrow ~ conc, 
                   curveid = agarvol, 
                   data = agar.dil[agar.dil$is == nm[[i]] & agar.dil$trial == t,], 
                   fct = LL.4())

summary.mef.fit <- data.frame(summary(agardil.drc)[[3]])
#outputs the summary of just the EC50 data including the estimate, standard error, upper and lower bounds of the 95% confidence intervals around the EC50
print(nm[[i]])
print("RELATIVE EC50")
EC50.od.rel <- data.frame(ED(agardil.drc, 
                             respLev = c(50), 
                             type = "relative",
                             interval = "delta"),
                          level = 0.95)
rel.ec50 <- EC50.od.rel[1][[1]]
print("ABSOLUTE EC50")
EC50.od.abs <- data.frame(ED(agardil.drc, 
                             respLev = c(50), 
                             type = "absolute",
                             interval = "delta"),
                          level = 0.95)
abs.ec50 <- EC50.od.rel[1][[1]]
print(summary.mef.fit)
print("LACK FIT")
fit <- modelFit(agardil.drc)
lackfitpvalue <- fit[5]$`p value`[2]
print(fit)
print("COMP EC50")
SI(agardil.drc, c(50, 50), ci = "delta")
agar.vol <- unique(agar.dil$agarvol[agar.dil$is == nm[[i]]])
agardil_i <- data.frame(rep(nm[[i]], 4), c(agar.vol), rep(t, 4), c(abs.ec50))
agardil <- rbind.data.frame(agardil, agardil_i)
#species <- as.character(unique(agar.dil$species[agar.dil$is == nm[[i]]]))
#plot(agardil.drc, col = T, lty = 1, 
 #    ylim = c(0, 100), 
  #   ylab = expression(bold("% Growth, relative to control")), 
   #  xlab = expression(bold("Mefenoxam concentration (μg ml"^-1~")")))
#title(main = paste(species))
#plot(unique(agar.dil$agarvol), abs.ec50, 
 #    xlim = c(0, 100),
  #   ylab = expression(bold("Absolute EC"[50] ~ "(μg ml"^-1~")")),
   #  xlab = expression(bold("Agar Volume (μl)")))
  }
}
colnames(agardil) <- c("is", "agar.vol", "trial", "EC50")
```

Two-way ANOVA looking at the effect of isolate and agar volume added per well to the resulting EC50

``` r
agardil.aov <- lm(EC50 ~ agar.vol + is + is:agar.vol, data = agardil)
anova(agardil.aov)
```

    ## Analysis of Variance Table
    ## 
    ## Response: EC50
    ##             Df   Sum Sq  Mean Sq F value            Pr(>F)    
    ## agar.vol     1 0.052874 0.052874 159.504 0.000000000004306 ***
    ## is           3 0.073871 0.024624  74.282 0.000000000002743 ***
    ## agar.vol:is  3 0.019905 0.006635  20.016 0.000001016249563 ***
    ## Residuals   24 0.007956 0.000331                              
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Looks like the mean of everthing was significant. This means we a significant effect of agar volume within each isolate. This means we should look at just the main effects (agar volume) within isolate.

``` r
library(ggplot2)
library(Rmisc)
p1 <- ggplot(agardil, aes(x = agar.vol, y = EC50, color = factor(is))) +
  stat_summary(fun.y = mean, geom = "point") + 
  stat_summary(fun.y = mean, geom = "line") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 1.5) +
  theme_bw() +
  guides(fill = guide_legend(title = "Isolate", title.position = "top")) +
  theme(axis.text.x = element_text(size = 30, face = "bold.italic", 
                                   angle = 90, hjust = 1, vjust = 0.5),
        axis.ticks = element_line(colour = 'black', size = 2, linetype = 'dashed'),
        axis.text.y = element_text(size = 30, face = "bold"),
        axis.title.x = element_text(size = 40, face = "bold"),
        axis.title.y = element_text(size = 40, face = "bold"),
        legend.title = element_text(size  = 40, face = "bold"),
        legend.text = element_text(size  = 30, face = "bold.italic")) +
  labs(list(x = "Agar Volume (μl)",  
            y = expression(bold("Mean EC"[50] ~"(μg ml"^-1~")"))))
print(p1)
```

![](Assay_Development_files/figure-markdown_github/unnamed-chunk-5-1.png)<!-- -->
