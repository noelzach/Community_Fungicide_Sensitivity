``` r
rm(list = ls(all=TRUE)) # removes all variables in the global environment so you start fresh

Sys.time() # prints out the time and date you ran the code
```

    ## [1] "2016-06-25 13:51:54 EDT"

``` r
options(scipen = 999) # stops anything from being in scientific notation
```

### Assay Development

Here we are analysing a dataset where we inoculated different volumes of macerated agar in the wells of a 96 well plate to see how the fungicide and assay would perform under these conditions.

The following loop will generate EC50's for each isolate, agar volume, trial combination and save that output into a data frame called agardil, to use later for an ANOVA

``` r
agar.dil <- agar.dil[!agar.dil$is == "BLANK",]
agar.dil <- agar.dil[!agar.dil$is == "NTC",]
agar.dil$is <- factor(agar.dil$is)
agardil <- NULL
nm <- levels(agar.dil$is)
for (t in 1:2){
  for (r in 1:3){
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
   }
  }
}
colnames(agardil) <- c("is", "agar.vol", "trial", "EC50")
```

Looking at growth rates

``` r
agar.dil$gr <- agar.dil$od600meanblank/24
gr0 <- aov(gr ~ factor(agarvol) * is, data = agar.dil[agar.dil$conc == 0,])
gr1 <- lmer(gr ~ factor(agarvol) + (1|is), data = agar.dil[agar.dil$conc == 0,])
gr2 <- aov(gr ~ factor(agarvol), data = agar.dil[agar.dil$conc == 0,])
anova(gr1, gr2, gr0)
```

    ## refitting model(s) with ML (instead of REML)

    ## Data: agar.dil[agar.dil$conc == 0, ]
    ## Models:
    ## gr2: gr ~ factor(agarvol)
    ## gr1: gr ~ factor(agarvol) + (1 | is)
    ## gr0: gr ~ factor(agarvol) * is
    ##     Df     AIC     BIC logLik deviance  Chisq Chi Df            Pr(>Chisq)
    ## gr2  5 -741.01 -728.19 375.50  -751.01                                    
    ## gr1  6 -831.52 -816.13 421.76  -843.52  92.51      1 < 0.00000000000000022
    ## gr0 17 -916.03 -872.44 475.02  -950.03 106.52     11 < 0.00000000000000022
    ##        
    ## gr2    
    ## gr1 ***
    ## gr0 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(gr0)
```

![](Assay_Development_files/figure-markdown_github/unnamed-chunk-5-1.png)<!-- -->![](Assay_Development_files/figure-markdown_github/unnamed-chunk-5-2.png)<!-- -->![](Assay_Development_files/figure-markdown_github/unnamed-chunk-5-3.png)<!-- -->![](Assay_Development_files/figure-markdown_github/unnamed-chunk-5-4.png)<!-- -->

``` r
anova(gr1)
```

    ## Analysis of Variance Table
    ##                 Df     Sum Sq   Mean Sq F value
    ## factor(agarvol)  3 0.00031591 0.0001053  13.418

``` r
gr1_lsmeans <- lsmeans(gr1, c("agarvol"))
contrast(gr1_lsmeans, "pairwise")
```

    ##  contrast      estimate           SE df t.ratio p.value
    ##  25 - 50  -0.0033639757 0.0008086937 89  -4.160  0.0004
    ##  25 - 75  -0.0040700521 0.0008086937 89  -5.033  <.0001
    ##  25 - 100 -0.0047000000 0.0008086937 89  -5.812  <.0001
    ##  50 - 75  -0.0007060764 0.0008086937 89  -0.873  0.8187
    ##  50 - 100 -0.0013360243 0.0008086937 89  -1.652  0.3552
    ##  75 - 100 -0.0006299479 0.0008086937 89  -0.779  0.8638
    ## 
    ## P value adjustment: tukey method for comparing a family of 4 estimates

Specifying a linear models to describe these data:

model: EC50 as a function of agar volume and isolate as a random effect

We are treating isolate as a random effect because these four isolates were chosen at random from a larger population of isolates and we want to generalize over isolate to see the effect of agar volume on the EC50.

``` r
agardil$agar.vol <- as.factor(agardil$agar.vol)
class(agardil$agar.vol)
```

    ## [1] "factor"

``` r
agardil2 <- lmer(EC50 ~ agar.vol + (1|is), data = agardil)
summary(agardil2)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: EC50 ~ agar.vol + (1 | is)
    ##    Data: agardil
    ## 
    ## REML criterion at convergence: -357.4
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.8597 -0.6176  0.0147  0.6567  2.0545 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance  Std.Dev.
    ##  is       (Intercept) 0.0030401 0.05514 
    ##  Residual             0.0009084 0.03014 
    ## Number of obs: 96, groups:  is, 4
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept) 0.029382   0.028247   1.040
    ## agar.vol50  0.043585   0.008701   5.009
    ## agar.vol75  0.084341   0.008701   9.694
    ## agar.vol100 0.107606   0.008701  12.368
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) agr.50 agr.75
    ## agar.vol50  -0.154              
    ## agar.vol75  -0.154  0.500       
    ## agar.vol100 -0.154  0.500  0.500

Lets make sure we look at the regression diagnostics

``` r
plot(agardil2)
```

![](Assay_Development_files/figure-markdown_github/unnamed-chunk-7-1.png)<!-- -->

``` r
qqnorm(resid(agardil2)); qqline(resid(agardil2))
```

![](Assay_Development_files/figure-markdown_github/unnamed-chunk-7-2.png)<!-- -->

The residuals look normally destributed. One could argue the need for a log transformation for a "V" shaped residual plot, but I've already tried that and the contrasts do not change so I am keeping these data non-transformed.

``` r
# lsmeans for our linear model
lm_lsmeans <- lsmeans(agardil2, c("agar.vol"))
plot(lm_lsmeans)
```

![](Assay_Development_files/figure-markdown_github/unnamed-chunk-8-1.png)<!-- -->

``` r
lm_lsmeans
```

    ##  agar.vol     lsmean        SE   df    lower.CL  upper.CL
    ##  25       0.02938231 0.0282468 3.23 -0.05706177 0.1158264
    ##  50       0.07296744 0.0282468 3.23 -0.01347665 0.1594115
    ##  75       0.11372355 0.0282468 3.23  0.02727947 0.2001676
    ##  100      0.13698806 0.0282468 3.23  0.05054398 0.2234321
    ## 
    ## Confidence level used: 0.95

``` r
contrast(lm_lsmeans, "pairwise")
```

    ##  contrast    estimate          SE df t.ratio p.value
    ##  25 - 50  -0.04358513 0.008700589 89  -5.009  <.0001
    ##  25 - 75  -0.08434124 0.008700589 89  -9.694  <.0001
    ##  25 - 100 -0.10760575 0.008700589 89 -12.368  <.0001
    ##  50 - 75  -0.04075612 0.008700589 89  -4.684  0.0001
    ##  50 - 100 -0.06402062 0.008700589 89  -7.358  <.0001
    ##  75 - 100 -0.02326451 0.008700589 89  -2.674  0.0434
    ## 
    ## P value adjustment: tukey method for comparing a family of 4 estimates
