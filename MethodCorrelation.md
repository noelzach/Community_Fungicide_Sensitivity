A function to do a t-test to determine if varias values in a linear model is significantly different than others. Arguments: reg = the linear model coefnum = 1 = intercept, 2 = beta ... val = value you would like to test, the default lm tests if value is significantly different from zero.

This function will run a linear model of the percent relative growth, as well as correlations in both pearson and spearman and will also plot it in ggplot, if desired.

data = data.frame concentration = concentration plot = logical

``` r
RelGrowth.lm <- function(data, concentration, plot){
  data.string.y <- data$odrelgrowth[data$conc == concentration]
  data.string.x <- data$ppmeanrelgrowth[data$conc == concentration]
  lin.mod <- lm(data.string.y ~ data.string.x)
  sum.linmod <- summary(lin.mod)
  pearson <- cor.test(data.string.x, data.string.y, method = "pearson")
  spearman <- cor.test(data.string.x, data.string.y, method = "spearman")
  beta <- ttest(lin.mod, 2, 1) # t-test for Beta significantly different than one
  bias <- (sum.linmod[[4]][2] - 1)*100 #percent Bias
  coeff.var <- (sum.linmod[[4]][4]/sum.linmod[[4]][2])*100
  if(plot == TRUE){
  p <- ggplot(data[data$conc == concentration,], aes(y = odrelgrowth, x = ppmeanrelgrowth)) +
    geom_point(aes(colour = factor(species))) +
    guides(colour = guide_legend(title = "Species")) +
    scale_y_continuous(limits = c(0, 110), breaks = c(0, 25, 50, 75, 100)) +
    scale_x_continuous(limits = c(0, 110), breaks = c(0, 25, 50, 75, 100)) +
    geom_smooth(method = "lm", se = FALSE, fullrange = TRUE, col = "black") +
    xlab("% Poison Plate Growth, relative to control") + 
    ylab("% Optical Density growth, relative to control") + 
    theme_bw() +
    theme(axis.text.x = element_text(size = 10, face = "bold"),
          axis.text.y = element_text(size = 10, face = "bold"),
          axis.title.x = element_text(size = 15, face = "bold"),
          axis.title.y = element_text(size = 15, face = "bold"),
          legend.text = element_text(size = 10, face = "bold.italic"),
          legend.key = element_blank(),
          legend.title = element_text(size = 15, face="bold"))
  results <- list(sum.linmod, pearson, spearman, beta, bias, coeff.var, p)
  names(results) <- c("lm", "pearson", "spearman", "beta.t", "per.bias","coeff.var", "plot")
  return(results)
  } else {
    results <- list(sum.linmod, pearson, spearman, beta, bias, coeff.var)
    names(results) <- c("lm", "pearson", "spearman", "beta.t", "per.bias", "coeff.var")
    return(results)
  }
}
```

We want to look at the correlation of the two methods, poison plate and optical density, to determine if the optical density method is any good.

Ethaboxam relative growth correlations

``` r
eth.cor.pp <- ddply(cor[cor$chem == "ethaboxam" & cor$method == "poison_plate",], c("is", "species", "conc"), 
      summarize, 
      ppmeanrelgrowth = 100*mean(relgrowth, na.rm = TRUE))
eth.cor.od <- ddply(cor[cor$chem == "ethaboxam" & cor$method == "optical_density",], c("is", "species", "conc"), 
      summarize, 
      odmeanrelgrowth = 100*mean(relgrowth, na.rm = TRUE))
eth.cor <- cbind.data.frame(eth.cor.pp, eth.cor.od$odmeanrelgrowth)

colnames(eth.cor) <- c("is", "species", "conc", "ppmeanrelgrowth", "odrelgrowth")
eth.001 <- RelGrowth.lm(eth.cor, concentration = 0.01, plot = TRUE)
eth.001
```

    ## $lm
    ## 
    ## Call:
    ## lm(formula = data.string.y ~ data.string.x)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -15.0564  -1.5956   0.6452   1.9085   6.0361 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value       Pr(>|t|)    
    ## (Intercept)   14.87111    8.74019   1.701         0.0999 .  
    ## data.string.x  0.85958    0.09191   9.353 0.000000000413 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 3.953 on 28 degrees of freedom
    ## Multiple R-squared:  0.7575, Adjusted R-squared:  0.7489 
    ## F-statistic: 87.47 on 1 and 28 DF,  p-value: 0.0000000004127
    ## 
    ## 
    ## $pearson
    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  data.string.x and data.string.y
    ## t = 9.3528, df = 28, p-value = 0.0000000004127
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  0.7430963 0.9368618
    ## sample estimates:
    ##       cor 
    ## 0.8703579 
    ## 
    ## 
    ## $spearman
    ## 
    ##  Spearman's rank correlation rho
    ## 
    ## data:  data.string.x and data.string.y
    ## S = 1734, p-value = 0.0004029
    ## alternative hypothesis: true rho is not equal to 0
    ## sample estimates:
    ##      rho 
    ## 0.614238 
    ## 
    ## 
    ## $beta.t
    ## [1] 0.1377708
    ## 
    ## $per.bias
    ## [1] -14.04184
    ## 
    ## $coeff.var
    ## [1] 10.69199
    ## 
    ## $plot

![](MethodCorrelation_files/figure-markdown_github/Ethaboxam%20Relative%20Growth-1.png)

``` r
eth.01 <- RelGrowth.lm(eth.cor, concentration = 0.1, plot = TRUE)
eth.01
```

    ## $lm
    ## 
    ## Call:
    ## lm(formula = data.string.y ~ data.string.x)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -20.9532  -3.3288  -0.7175   4.8422  15.0876 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value           Pr(>|t|)    
    ## (Intercept)   29.26006    4.10788   7.123 0.0000000946325955 ***
    ## data.string.x  0.74368    0.05093  14.602 0.0000000000000128 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 7.187 on 28 degrees of freedom
    ## Multiple R-squared:  0.8839, Adjusted R-squared:  0.8798 
    ## F-statistic: 213.2 on 1 and 28 DF,  p-value: 0.00000000000001276
    ## 
    ## 
    ## $pearson
    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  data.string.x and data.string.y
    ## t = 14.602, df = 28, p-value = 0.00000000000001288
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  0.8769341 0.9714102
    ## sample estimates:
    ##       cor 
    ## 0.9401724 
    ## 
    ## 
    ## $spearman
    ## 
    ##  Spearman's rank correlation rho
    ## 
    ## data:  data.string.x and data.string.y
    ## S = 606, p-value = 0.0000005663
    ## alternative hypothesis: true rho is not equal to 0
    ## sample estimates:
    ##       rho 
    ## 0.8651835 
    ## 
    ## 
    ## $beta.t
    ## [1] 0.00002534679
    ## 
    ## $per.bias
    ## [1] -25.63232
    ## 
    ## $coeff.var
    ## [1] 6.848325
    ## 
    ## $plot

    ## Warning: Removed 2 rows containing missing values (geom_smooth).

![](MethodCorrelation_files/figure-markdown_github/Ethaboxam%20Relative%20Growth-2.png)

``` r
eth.05 <- RelGrowth.lm(eth.cor, concentration = 0.5, plot = TRUE)
eth.05
```

    ## $lm
    ## 
    ## Call:
    ## lm(formula = data.string.y ~ data.string.x)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -27.431  -4.319   1.050   5.400  21.930 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value         Pr(>|t|)    
    ## (Intercept)   17.56121    4.22120    4.16         0.000273 ***
    ## data.string.x  0.81890    0.06602   12.40 0.00000000000068 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 11.24 on 28 degrees of freedom
    ## Multiple R-squared:  0.846,  Adjusted R-squared:  0.8405 
    ## F-statistic: 153.9 on 1 and 28 DF,  p-value: 0.0000000000006796
    ## 
    ## 
    ## $pearson
    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  data.string.x and data.string.y
    ## t = 12.404, df = 28, p-value = 0.0000000000006795
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  0.8368378 0.9614634
    ## sample estimates:
    ##       cor 
    ## 0.9197995 
    ## 
    ## 
    ## $spearman
    ## 
    ##  Spearman's rank correlation rho
    ## 
    ## data:  data.string.x and data.string.y
    ## S = 586, p-value = 0.0000005246
    ## alternative hypothesis: true rho is not equal to 0
    ## sample estimates:
    ##       rho 
    ## 0.8696329 
    ## 
    ## 
    ## $beta.t
    ## [1] 0.01049333
    ## 
    ## $per.bias
    ## [1] -18.11014
    ## 
    ## $coeff.var
    ## [1] 8.062027
    ## 
    ## $plot

![](MethodCorrelation_files/figure-markdown_github/Ethaboxam%20Relative%20Growth-3.png)

``` r
eth.1 <- RelGrowth.lm(eth.cor, concentration = 1, plot = TRUE)
eth.1
```

    ## $lm
    ## 
    ## Call:
    ## lm(formula = data.string.y ~ data.string.x)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -35.804  -7.556   0.999  11.640  20.084 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value     Pr(>|t|)    
    ## (Intercept)    9.13243    5.28545   1.728        0.095 .  
    ## data.string.x  0.74812    0.09552   7.832 0.0000000157 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 16.69 on 28 degrees of freedom
    ## Multiple R-squared:  0.6866, Adjusted R-squared:  0.6754 
    ## F-statistic: 61.34 on 1 and 28 DF,  p-value: 0.00000001569
    ## 
    ## 
    ## $pearson
    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  data.string.x and data.string.y
    ## t = 7.8318, df = 28, p-value = 0.00000001569
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  0.6676330 0.9155582
    ## sample estimates:
    ##       cor 
    ## 0.8286008 
    ## 
    ## 
    ## $spearman
    ## 
    ##  Spearman's rank correlation rho
    ## 
    ## data:  data.string.x and data.string.y
    ## S = 1004, p-value = 0.000001908
    ## alternative hypothesis: true rho is not equal to 0
    ## sample estimates:
    ##       rho 
    ## 0.7766407 
    ## 
    ## 
    ## $beta.t
    ## [1] 0.01350023
    ## 
    ## $per.bias
    ## [1] -25.18779
    ## 
    ## $coeff.var
    ## [1] 12.76849
    ## 
    ## $plot

![](MethodCorrelation_files/figure-markdown_github/Ethaboxam%20Relative%20Growth-4.png)

``` r
eth.5 <- RelGrowth.lm(eth.cor, concentration = 5, plot = TRUE)
```

    ## Warning in cor.test.default(data.string.x, data.string.y, method =
    ## "spearman"): Cannot compute exact p-value with ties

``` r
eth.5
```

    ## $lm
    ## 
    ## Call:
    ## lm(formula = data.string.y ~ data.string.x)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -19.5394  -3.2108   0.2078   4.2866  14.9406 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value             Pr(>|t|)    
    ## (Intercept)   11.26559    1.45532   7.741         0.0000000197 ***
    ## data.string.x  1.01056    0.04476  22.579 < 0.0000000000000002 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 6.873 on 28 degrees of freedom
    ## Multiple R-squared:  0.9479, Adjusted R-squared:  0.9461 
    ## F-statistic: 509.8 on 1 and 28 DF,  p-value: < 0.00000000000000022
    ## 
    ## 
    ## $pearson
    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  data.string.x and data.string.y
    ## t = 22.579, df = 28, p-value < 0.00000000000000022
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  0.9447310 0.9875068
    ## sample estimates:
    ##      cor 
    ## 0.973621 
    ## 
    ## 
    ## $spearman
    ## 
    ##  Spearman's rank correlation rho
    ## 
    ## data:  data.string.x and data.string.y
    ## S = 2067.7, p-value = 0.00207
    ## alternative hypothesis: true rho is not equal to 0
    ## sample estimates:
    ##       rho 
    ## 0.5399933 
    ## 
    ## 
    ## $beta.t
    ## [1] 0.8152027
    ## 
    ## $per.bias
    ## [1] 1.055944
    ## 
    ## $coeff.var
    ## [1] 4.428861
    ## 
    ## $plot

    ## Warning: Removed 9 rows containing missing values (geom_smooth).

![](MethodCorrelation_files/figure-markdown_github/Ethaboxam%20Relative%20Growth-5.png)

``` r
eth.20 <- RelGrowth.lm(eth.cor, concentration = 20, plot = TRUE)
```

    ## Warning in cor.test.default(data.string.x, data.string.y, method =
    ## "spearman"): Cannot compute exact p-value with ties

``` r
eth.20
```

    ## $lm
    ## 
    ## Call:
    ## lm(formula = data.string.y ~ data.string.x)
    ## 
    ## Residuals:
    ##    Min     1Q Median     3Q    Max 
    ## -8.389 -3.986 -1.106  2.896 14.971 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value             Pr(>|t|)    
    ## (Intercept)   12.86069    1.18461   10.86      0.0000000000152 ***
    ## data.string.x  0.95357    0.04008   23.79 < 0.0000000000000002 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 5.816 on 28 degrees of freedom
    ## Multiple R-squared:  0.9529, Adjusted R-squared:  0.9512 
    ## F-statistic: 566.2 on 1 and 28 DF,  p-value: < 0.00000000000000022
    ## 
    ## 
    ## $pearson
    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  data.string.x and data.string.y
    ## t = 23.794, df = 28, p-value < 0.00000000000000022
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  0.9499662 0.9887137
    ## sample estimates:
    ##       cor 
    ## 0.9761532 
    ## 
    ## 
    ## $spearman
    ## 
    ##  Spearman's rank correlation rho
    ## 
    ## data:  data.string.x and data.string.y
    ## S = 2402.7, p-value = 0.009538
    ## alternative hypothesis: true rho is not equal to 0
    ## sample estimates:
    ##       rho 
    ## 0.4654791 
    ## 
    ## 
    ## $beta.t
    ## [1] 0.256445
    ## 
    ## $per.bias
    ## [1] -4.642783
    ## 
    ## $coeff.var
    ## [1] 4.202696
    ## 
    ## $plot

    ## Warning: Removed 6 rows containing missing values (geom_smooth).

![](MethodCorrelation_files/figure-markdown_github/Ethaboxam%20Relative%20Growth-6.png)

Mefenoxam relative growth correlations

``` r
mef.cor.pp <- ddply(cor[cor$chem == "mefenoxam" & cor$method == "poison_plate",], c("is", "species", "conc"), 
      summarize, 
      ppmeanrelgrowth = 100*mean(relgrowth, na.rm = TRUE))
mef.cor.od <- ddply(cor[cor$chem == "mefenoxam" & cor$method == "optical_density",], c("is", "species", "conc"), 
      summarize, 
      odmeanrelgrowth = 100*mean(relgrowth, na.rm = TRUE))
mef.cor.od <- mef.cor.od[!mef.cor.od$conc == 100,]
mef.cor.pp <- mef.cor.pp[!mef.cor.pp$conc == 100,]
mef.cor <- cbind.data.frame(mef.cor.pp, mef.cor.od$odmeanrelgrowth)
colnames(mef.cor) <- c("is", "species", "conc", "ppmeanrelgrowth", "odrelgrowth")
mef.001 <- RelGrowth.lm(mef.cor, concentration = 0.01, plot = TRUE)
mef.001
```

    ## $lm
    ## 
    ## Call:
    ## lm(formula = data.string.y ~ data.string.x)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -22.312  -7.813   1.125   4.147  43.193 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   -19.3623    24.0959  -0.804 0.430657    
    ## data.string.x   1.2286     0.2687   4.573 0.000165 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 12.51 on 21 degrees of freedom
    ## Multiple R-squared:  0.4989, Adjusted R-squared:  0.4751 
    ## F-statistic: 20.91 on 1 and 21 DF,  p-value: 0.0001652
    ## 
    ## 
    ## $pearson
    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  data.string.x and data.string.y
    ## t = 4.5729, df = 21, p-value = 0.0001652
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  0.4149891 0.8663221
    ## sample estimates:
    ##       cor 
    ## 0.7063615 
    ## 
    ## 
    ## $spearman
    ## 
    ##  Spearman's rank correlation rho
    ## 
    ## data:  data.string.x and data.string.y
    ## S = 628, p-value = 0.000388
    ## alternative hypothesis: true rho is not equal to 0
    ## sample estimates:
    ##       rho 
    ## 0.6897233 
    ## 
    ## 
    ## $beta.t
    ## [1] 0.4044911
    ## 
    ## $per.bias
    ## [1] 22.85705
    ## 
    ## $coeff.var
    ## [1] 21.86781
    ## 
    ## $plot

    ## Warning: Removed 2 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 2 rows containing missing values (geom_point).

    ## Warning: Removed 17 rows containing missing values (geom_smooth).

![](MethodCorrelation_files/figure-markdown_github/Mefenoxam%20relative%20growth-1.png)

``` r
mef.01 <- RelGrowth.lm(mef.cor, concentration = 0.1, plot = TRUE)
mef.01
```

    ## $lm
    ## 
    ## Call:
    ## lm(formula = data.string.y ~ data.string.x)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -32.564  -2.851   0.934   7.755  22.076 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value    Pr(>|t|)    
    ## (Intercept)    18.8134     5.3948   3.487      0.0022 ** 
    ## data.string.x   0.8039     0.1027   7.825 0.000000117 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 14.26 on 21 degrees of freedom
    ## Multiple R-squared:  0.7446, Adjusted R-squared:  0.7325 
    ## F-statistic: 61.23 on 1 and 21 DF,  p-value: 0.0000001171
    ## 
    ## 
    ## $pearson
    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  data.string.x and data.string.y
    ## t = 7.8253, df = 21, p-value = 0.0000001171
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  0.6995462 0.9405659
    ## sample estimates:
    ##       cor 
    ## 0.8629212 
    ## 
    ## 
    ## $spearman
    ## 
    ##  Spearman's rank correlation rho
    ## 
    ## data:  data.string.x and data.string.y
    ## S = 406, p-value = 0.000006679
    ## alternative hypothesis: true rho is not equal to 0
    ## sample estimates:
    ##       rho 
    ## 0.7994071 
    ## 
    ## 
    ## $beta.t
    ## [1] 0.0701064
    ## 
    ## $per.bias
    ## [1] -19.60604
    ## 
    ## $coeff.var
    ## [1] 12.77914
    ## 
    ## $plot

    ## Warning: Removed 1 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 1 rows containing missing values (geom_point).

![](MethodCorrelation_files/figure-markdown_github/Mefenoxam%20relative%20growth-2.png)

``` r
mef.05 <- RelGrowth.lm(mef.cor, concentration = 0.5, plot = TRUE)
```

    ## Warning in cor.test.default(data.string.x, data.string.y, method =
    ## "spearman"): Cannot compute exact p-value with ties

``` r
mef.05
```

    ## $lm
    ## 
    ## Call:
    ## lm(formula = data.string.y ~ data.string.x)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -22.754  -5.424   1.832   7.021  15.238 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value        Pr(>|t|)    
    ## (Intercept)    4.39083    2.81077   1.562           0.133    
    ## data.string.x  0.85091    0.06972  12.205 0.0000000000532 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 10.16 on 21 degrees of freedom
    ## Multiple R-squared:  0.8764, Adjusted R-squared:  0.8706 
    ## F-statistic:   149 on 1 and 21 DF,  p-value: 0.00000000005324
    ## 
    ## 
    ## $pearson
    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  data.string.x and data.string.y
    ## t = 12.205, df = 21, p-value = 0.00000000005324
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  0.8532511 0.9729345
    ## sample estimates:
    ##       cor 
    ## 0.9361856 
    ## 
    ## 
    ## $spearman
    ## 
    ##  Spearman's rank correlation rho
    ## 
    ## data:  data.string.x and data.string.y
    ## S = 654.61, p-value = 0.0003928
    ## alternative hypothesis: true rho is not equal to 0
    ## sample estimates:
    ##       rho 
    ## 0.6765746 
    ## 
    ## 
    ## $beta.t
    ## [1] 0.04439049
    ## 
    ## $per.bias
    ## [1] -14.90896
    ## 
    ## $coeff.var
    ## [1] 8.193344
    ## 
    ## $plot

![](MethodCorrelation_files/figure-markdown_github/Mefenoxam%20relative%20growth-3.png)

``` r
mef.1 <- RelGrowth.lm(mef.cor, concentration = 1, plot = TRUE)
```

    ## Warning in cor.test.default(data.string.x, data.string.y, method =
    ## "spearman"): Cannot compute exact p-value with ties

``` r
mef.1
```

    ## $lm
    ## 
    ## Call:
    ## lm(formula = data.string.y ~ data.string.x)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -14.9337  -4.1843  -0.1966   3.7998  17.1367 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value         Pr(>|t|)    
    ## (Intercept)    4.70494    2.25201   2.089            0.049 *  
    ## data.string.x  0.94986    0.06457  14.711 0.00000000000156 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 8.866 on 21 degrees of freedom
    ## Multiple R-squared:  0.9115, Adjusted R-squared:  0.9073 
    ## F-statistic: 216.4 on 1 and 21 DF,  p-value: 0.000000000001565
    ## 
    ## 
    ## $pearson
    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  data.string.x and data.string.y
    ## t = 14.711, df = 21, p-value = 0.000000000001565
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  0.8946287 0.9809135
    ## sample estimates:
    ##       cor 
    ## 0.9547498 
    ## 
    ## 
    ## $spearman
    ## 
    ##  Spearman's rank correlation rho
    ## 
    ## data:  data.string.x and data.string.y
    ## S = 769.79, p-value = 0.001613
    ## alternative hypothesis: true rho is not equal to 0
    ## sample estimates:
    ##     rho 
    ## 0.61967 
    ## 
    ## 
    ## $beta.t
    ## [1] 0.4460584
    ## 
    ## $per.bias
    ## [1] -5.01428
    ## 
    ## $coeff.var
    ## [1] 6.797615
    ## 
    ## $plot

![](MethodCorrelation_files/figure-markdown_github/Mefenoxam%20relative%20growth-4.png)

``` r
mef.10 <- RelGrowth.lm(mef.cor, concentration = 10, plot = TRUE)
```

    ## Warning in cor.test.default(data.string.x, data.string.y, method =
    ## "spearman"): Cannot compute exact p-value with ties

``` r
mef.10
```

    ## $lm
    ## 
    ## Call:
    ## lm(formula = data.string.y ~ data.string.x)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -12.865  -3.993  -1.948   2.909  15.851 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value          Pr(>|t|)    
    ## (Intercept)    4.83928    1.96862   2.458            0.0227 *  
    ## data.string.x  1.00122    0.06288  15.923 0.000000000000338 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 8.297 on 21 degrees of freedom
    ## Multiple R-squared:  0.9235, Adjusted R-squared:  0.9199 
    ## F-statistic: 253.6 on 1 and 21 DF,  p-value: 0.0000000000003381
    ## 
    ## 
    ## $pearson
    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  data.string.x and data.string.y
    ## t = 15.923, df = 21, p-value = 0.000000000000338
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  0.9087863 0.9835784
    ## sample estimates:
    ##       cor 
    ## 0.9609959 
    ## 
    ## 
    ## $spearman
    ## 
    ##  Spearman's rank correlation rho
    ## 
    ## data:  data.string.x and data.string.y
    ## S = 828.94, p-value = 0.003016
    ## alternative hypothesis: true rho is not equal to 0
    ## sample estimates:
    ##       rho 
    ## 0.5904454 
    ## 
    ## 
    ## $beta.t
    ## [1] 0.9847594
    ## 
    ## $per.bias
    ## [1] 0.1215481
    ## 
    ## $coeff.var
    ## [1] 6.280042
    ## 
    ## $plot

    ## Warning: Removed 1 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 1 rows containing missing values (geom_point).

![](MethodCorrelation_files/figure-markdown_github/Mefenoxam%20relative%20growth-5.png)

From the analysis above it looks like there are some isolates with percent relative growth above 50% at the highest concentration tested. Therefore those isolates will not have an EC50 and cannot be calculated. We will express these isolates in terms of their relative growth.

``` r
# briefly we are going to split up the data frame into each concentration to make it easier to take out individual isolates
cor_eth <- cor[cor$chem == "ethaboxam",]
cor_mef <- cor[cor$chem == "mefenoxam",]

# taking out insensitive isolates, otherwise convergence would not occur with the ll.4 model
insens_iso_eth <- c("MISO_8-29.1", "C-MNSO2_2-21", "NDSO_L-8-6", "ILSO_6-15C")
insens_iso_mef <- c("V-MISO2_6-46", "23.4B", "1.18A")

# creating a dataframe with only those isolates
insenscor_eth <- cor_eth[cor_eth$is %in% insens_iso_eth,]
insenscor_mef <- cor_mef[cor_mef$is %in% insens_iso_mef,]

# taking out those isolates
cor_eth <- cor_eth[!cor_eth$is %in% insens_iso_eth,]
cor_mef <- cor_mef[!cor_mef$is %in% insens_iso_mef,]
cor_eth$is <- factor(cor_eth$is)
cor_mef$is <- factor(cor_mef$is)

cor <- rbind.data.frame(cor_eth, cor_mef)
```

This code iterates through every isolate and generates a relative and absolute EC50 using the LL.4 model of drc and saves the output. I already went through and picked out the isolates for both mefenoxam and ethaboxam that have an EC50 beyond the concentration range tested in this study.

``` r
relative <- function(data){
  EC50.pp.rel <- data.frame(ED(data, 
                             respLev = c(50), 
                             type = "relative",
                             interval = "delta"),
                          level = 0.95)
rel.ec50.pp <- EC50.pp.rel[1][[1]]
return(rel.ec50.pp)
}
absolute <- function(data){
  EC50.pp.abs <- data.frame(ED(data, 
                             respLev = c(50), 
                             type = "absolute",
                             interval = "delta"),
                          level = 0.95)
abs.ec50.pp <- EC50.pp.abs[1][[1]]
return(abs.ec50.pp)
}
plug_EC <- function(chemistry){
  nm <- unique(cor$is[cor$chem == as.character(chemistry)])
dataframe_names <- c("is", "species", "method", "trial", "chem", "absolute", "relative")
meth.cor <- NULL
for (i in seq_along(nm)){
mefcor.pp.drc <- drm(100*relgrowth ~ conc, data = cor[cor$is == nm[[i]] & cor$method == "poison_plate" & cor$chem == as.character(chemistry),], curveid = trial, fct = LL.4(), na.action = na.omit)
mefcor.od.drc <- drm(100*relgrowth ~ conc, data = cor[cor$is == nm[[i]] & cor$method == "optical_density" & cor$chem == as.character(chemistry),], curveid = trial, fct = LL.4(), na.action = na.omit)
paste(print(nm[[i]]))
# RELATIVE 
rel.pp.mef <- relative(mefcor.pp.drc)
rel.od.mef <- relative(mefcor.od.drc)

# ABSOLUTE 
abs.pp.mef <- absolute(mefcor.pp.drc)
abs.od.mef <- absolute(mefcor.od.drc)

mef.cor_od <- data.frame(cbind.data.frame( 
                  nm[[i]], 
                  as.character(unique(cor$species[cor$is == nm[[i]]])),
                  as.character(unique(cor$method[cor$is == nm[[i]] & cor$method == "optical_density"])),
                  unique(cor$trial[cor$is == nm[[i]] & cor$method == "optical_density"]),
                  as.character(chemistry),
                  as.numeric(abs.od.mef),
                  as.numeric(rel.od.mef)))
colnames(mef.cor_od) <- dataframe_names
mef.cor_pp <- data.frame(cbind.data.frame( 
                  nm[[i]], 
                  as.character(unique(cor$species[cor$is == nm[[i]]])),
                  as.character(unique(cor$method[cor$is == nm[[i]] & cor$method == "poison_plate"])),
                  unique(cor$trial[cor$is == nm[[i]] & cor$method == "poison_plate"]),
                  as.character(chemistry),
                  as.numeric(abs.pp.mef),
                  as.numeric(rel.pp.mef)))
colnames(mef.cor_pp) <- dataframe_names

meth.cor <- rbind.data.frame(meth.cor, mef.cor_od, mef.cor_pp)
  }
return(meth.cor)
}

mefenoxam_ec50 <- plug_EC("mefenoxam")
ethaboxam_ec50 <- plug_EC("ethaboxam")

ec50_cor <- rbind.data.frame(ethaboxam_ec50, mefenoxam_ec50)

# not sure why one of the isolates had duplicates but they did, so this code will find duplicates in the absolute ec50s and delete the entire row. 
ec50_cor <- ec50_cor[-which(duplicated(ec50_cor$absolute)),] 
```

Table of EC50 values for every isolate tested for both mefenoxam and ethaboxam on both poison plate and optical density method.

| is               | species                             | chem      | method           |   mean.abs|    std.abs|   mean.rel|    std.rel|
|:-----------------|:------------------------------------|:----------|:-----------------|----------:|----------:|----------:|----------:|
| AR\_127.S.2.3.A  | Pythium irregulare                  | ethaboxam | optical\_density |  0.6017397|  0.4608413|  0.3258476|  0.2118367|
| AR\_127.S.2.3.A  | Pythium irregulare                  | ethaboxam | poison\_plate    |  0.6474189|  0.0242181|  0.7462384|  0.1600082|
| AR\_127.S.2.3.A  | Pythium irregulare                  | mefenoxam | optical\_density |  0.1831292|  0.0722302|  0.1036954|  0.0178995|
| AR\_127.S.2.3.A  | Pythium irregulare                  | mefenoxam | poison\_plate    |  0.1259046|  0.0503855|  0.1228201|  0.0557412|
| AR\_262.S.1.6.A  | Pythium spinosum                    | ethaboxam | optical\_density |  0.4890281|  0.0532135|  0.4472674|  0.0672784|
| AR\_262.S.1.6.A  | Pythium spinosum                    | ethaboxam | poison\_plate    |  0.2093399|  0.0147487|  0.1991737|  0.0172343|
| ARS\_284.S1.4.A  | Pythium lutarium                    | ethaboxam | optical\_density |  0.2669444|  0.0983561|  0.2488793|  0.0977931|
| ARS\_284.S1.4.A  | Pythium lutarium                    | ethaboxam | poison\_plate    |  0.2260918|  0.0293079|  0.1949081|  0.0054857|
| C-NDSO2\_1-11    | Pythium perplexum                   | ethaboxam | optical\_density |  0.1673536|  0.0045666|  0.1324687|  0.0026915|
| C-NDSO2\_1-11    | Pythium perplexum                   | ethaboxam | poison\_plate    |  0.0983410|  0.0112713|  0.0977276|  0.0100928|
| IASO\_10-37.8RT  | Pythium aff. dissotocum             | ethaboxam | optical\_density |  1.1413192|  0.3019554|  1.0115718|  0.2570090|
| IASO\_10-37.8RT  | Pythium aff. dissotocum             | ethaboxam | poison\_plate    |  1.0220161|  0.6068654|  1.0072584|  0.5928676|
| IASO\_10-38.14RT | Pythium aff. dissotocum             | ethaboxam | optical\_density |  1.0125252|  0.2669542|  0.9121016|  0.1950198|
| IASO\_10-38.14RT | Pythium aff. dissotocum             | ethaboxam | poison\_plate    |  1.7250786|  0.2861206|  1.8877662|  0.3800344|
| IASO\_10-39.16rt | Pythium aff. dissotocum             | ethaboxam | optical\_density |  0.6248467|  0.0759800|  0.5562301|  0.0206426|
| IASO\_10-39.16rt | Pythium aff. dissotocum             | ethaboxam | poison\_plate    |  1.6056028|  0.1531086|  1.7370847|  0.1606027|
| IASO\_6-10.15h   | Pythium oopapillum                  | ethaboxam | optical\_density |  0.5136203|  0.0286944|  0.4836774|  0.0102900|
| IASO\_6-10.15h   | Pythium oopapillum                  | ethaboxam | poison\_plate    |  1.1130479|  0.2173577|  1.0919763|  0.2239991|
| ILSO\_1-31       | Pythium irregulare                  | ethaboxam | optical\_density |  0.8616940|  0.1247902|  0.8372427|  0.0556342|
| ILSO\_1-31       | Pythium irregulare                  | ethaboxam | poison\_plate    |  0.8058171|  0.0387544|  1.0671747|  0.0765512|
| ILSO\_3-48C      | Pythium irregulare                  | ethaboxam | optical\_density |  0.9968918|  0.2276545|  0.5450930|  0.1574642|
| ILSO\_3-48C      | Pythium irregulare                  | ethaboxam | poison\_plate    |  0.7286362|  0.0211646|  0.9552630|  0.0460969|
| ILSO\_3-48C      | Pythium irregulare                  | mefenoxam | optical\_density |  0.2218693|  0.0780914|  0.1509015|  0.0575959|
| ILSO\_3-48C      | Pythium irregulare                  | mefenoxam | poison\_plate    |  0.1263038|  0.0133000|  0.1136575|  0.0248296|
| ILSO\_5-42C      | Pythium oopapillum                  | ethaboxam | optical\_density |  1.0026325|  0.4214322|  0.9374770|  0.3795077|
| ILSO\_5-42C      | Pythium oopapillum                  | ethaboxam | poison\_plate    |  1.3388329|  0.0166538|  1.5284298|  0.1501364|
| ILSO\_5-42C      | Pythium oopapillum                  | mefenoxam | optical\_density |  0.0934891|  0.0727647|  0.0828173|  0.0621251|
| ILSO\_5-42C      | Pythium oopapillum                  | mefenoxam | poison\_plate    |  0.0176811|  0.0013969|  0.0190659|  0.0028424|
| ILSO\_6-2B       | Pythium oopapillum                  | ethaboxam | optical\_density |  0.7590296|  0.0573760|  0.6391139|  0.0124021|
| ILSO\_6-2B       | Pythium oopapillum                  | ethaboxam | poison\_plate    |  0.9836553|  0.0599357|  0.9708861|  0.0039751|
| ILSO\_6-2B       | Pythium oopapillum                  | mefenoxam | optical\_density |  0.1616391|  0.0430944|  0.1512861|  0.0385509|
| ILSO\_6-2B       | Pythium oopapillum                  | mefenoxam | poison\_plate    |  0.1374047|  0.0218105|  0.1122750|  0.0353251|
| INSO\_1-10C      | Pythium sylvaticum                  | ethaboxam | optical\_density |  0.6021485|  0.2261376|  0.4896181|  0.2288380|
| INSO\_1-10C      | Pythium sylvaticum                  | ethaboxam | poison\_plate    |  0.1023615|  0.0310532|  0.1085730|  0.0178761|
| INSO\_1-10C      | Pythium sylvaticum                  | mefenoxam | optical\_density |  0.0556356|  0.0142788|  0.0540463|  0.0146426|
| INSO\_1-10C      | Pythium sylvaticum                  | mefenoxam | poison\_plate    |  0.0240773|  0.0109447|  0.0240587|  0.0109250|
| INSO\_1-8A       | Pythium lutarium                    | ethaboxam | optical\_density |  0.6975005|  0.1292836|  0.6548995|  0.0999101|
| INSO\_1-8A       | Pythium lutarium                    | ethaboxam | poison\_plate    |  1.4568919|  0.1234575|  1.5417715|  0.0704994|
| KSSO\_6-1        | Pythium ultimum var. ultimum        | ethaboxam | optical\_density |  2.1738744|  0.1268949|  1.8935903|  0.0802024|
| KSSO\_6-1        | Pythium ultimum var. ultimum        | ethaboxam | poison\_plate    |  1.4760269|  0.1444716|  1.6207549|  0.0955815|
| KSSO\_6-47       | Phytophthora sansomeana             | ethaboxam | optical\_density |  0.0634017|  0.0315296|  0.0526138|  0.0286371|
| KSSO\_6-47       | Phytophthora sansomeana             | ethaboxam | poison\_plate    |  0.0324276|  0.0043185|  0.0324150|  0.0039692|
| MISO\_1-4        | Pythium lutarium                    | ethaboxam | optical\_density |  0.9515964|  0.1695581|  0.8880759|  0.1894814|
| MISO\_1-4        | Pythium lutarium                    | ethaboxam | poison\_plate    |  0.6058080|  0.0382837|  0.6033125|  0.0375648|
| MISO\_5-19H      | Pythium glomeratum                  | ethaboxam | optical\_density |  1.3033651|  0.3326300|  0.9839949|  0.2931130|
| MISO\_5-19H      | Pythium glomeratum                  | ethaboxam | poison\_plate    |  1.1133215|  0.0441768|  1.1914970|  0.0877672|
| MISO\_8-10       | Pythium ultimum var. ultimum        | ethaboxam | optical\_density |  0.7494623|  0.2384325|  0.6596120|  0.2148583|
| MISO\_8-10       | Pythium ultimum var. ultimum        | ethaboxam | poison\_plate    |  1.2705559|  0.0327570|  1.3791302|  0.0529770|
| NDSO\_1-42       | Pythium sylvaticum                  | ethaboxam | optical\_density |  0.5083377|  0.0455087|  0.4342847|  0.0397514|
| NDSO\_1-42       | Pythium sylvaticum                  | ethaboxam | poison\_plate    |  0.5029008|  0.0315540|  0.5163838|  0.0454677|
| NDSO\_1-42       | Pythium sylvaticum                  | mefenoxam | optical\_density |  0.0843794|  0.0395403|  0.0776105|  0.0309785|
| NDSO\_1-42       | Pythium sylvaticum                  | mefenoxam | poison\_plate    |  0.0584407|  0.0271372|  0.0592125|  0.0288098|
| NESO\_2-13       | Pythium sylvaticum                  | ethaboxam | optical\_density |  0.5351099|  0.1587134|  0.4859459|  0.1842366|
| NESO\_2-13       | Pythium sylvaticum                  | ethaboxam | poison\_plate    |  0.2755640|  0.0195702|  0.3037442|  0.0022381|
| NESO\_2-13       | Pythium sylvaticum                  | mefenoxam | optical\_density |  0.1626739|  0.0513231|  0.1500946|  0.0446556|
| NESO\_2-13       | Pythium sylvaticum                  | mefenoxam | poison\_plate    |  0.0659475|  0.0136608|  0.0657198|  0.0179300|
| NESO\_4-29       | Pythium perplexum                   | ethaboxam | optical\_density |  0.0997645|  0.0070532|  0.0880328|  0.0143253|
| NESO\_4-29       | Pythium perplexum                   | ethaboxam | poison\_plate    |  0.1238281|  0.0031861|  0.1199439|  0.0006862|
| V-KSSO2\_1-7     | Phytophthora sansomeana             | ethaboxam | optical\_density |  0.0477238|  0.0204292|  0.0408790|  0.0145521|
| V-KSSO2\_1-7     | Phytophthora sansomeana             | ethaboxam | poison\_plate    |  0.0164446|  0.0004557|  0.0161274|  0.0007975|
| V-KSSO2\_3-6     | Phytophthora sansomeana             | ethaboxam | optical\_density |  0.0527088|  0.0075370|  0.0435138|  0.0123592|
| V-KSSO2\_3-6     | Phytophthora sansomeana             | ethaboxam | poison\_plate    |  0.0331940|  0.0092453|  0.0336671|  0.0089418|
| V-MISO2\_2-57    | Pythium perplexum                   | ethaboxam | optical\_density |  0.1575673|  0.0323676|  0.1342038|  0.0319532|
| V-MISO2\_2-57    | Pythium perplexum                   | ethaboxam | poison\_plate    |  0.0746280|  0.0000147|  0.0735610|  0.0014948|
| WISO\_4-13       | Pythium ultimum var. ultimum        | ethaboxam | optical\_density |  1.4114181|  0.3025662|  1.2496867|  0.2475802|
| WISO\_4-13       | Pythium ultimum var. ultimum        | ethaboxam | poison\_plate    |  1.6811566|  0.0761205|  1.7426576|  0.0073123|
| AR\_96.S.2.1.A   | Pythium spinosum                    | mefenoxam | optical\_density |  0.0640711|  0.0155281|  0.0552741|  0.0148170|
| AR\_96.S.2.1.A   | Pythium spinosum                    | mefenoxam | poison\_plate    |  0.0408055|  0.0210795|  0.0416514|  0.0219245|
| C-KSSO2\_1-25    | Phytopythium litorale               | mefenoxam | optical\_density |  1.1289415|  0.2611411|  0.6216587|  0.0899566|
| C-KSSO2\_1-25    | Phytopythium litorale               | mefenoxam | poison\_plate    |  0.9949350|  0.1711380|  0.6115304|  0.1309216|
| C-SDSO2\_5-35    | Pythium intermedium                 | mefenoxam | optical\_density |  0.1038246|  0.0015915|  0.1042345|  0.0085674|
| C-SDSO2\_5-35    | Pythium intermedium                 | mefenoxam | poison\_plate    |  0.2039187|  0.0848682|  0.1791958|  0.0740958|
| IASO\_3-41.17    | Phytophthora sojae                  | mefenoxam | optical\_density |  0.0564770|  0.0357976|  0.0646445|  0.0374231|
| IASO\_3-41.17    | Phytophthora sojae                  | mefenoxam | poison\_plate    |  0.0193404|  0.0062906|  0.0193835|  0.0063338|
| IASO\_6-10.15H   | Pythium oopapillum                  | mefenoxam | optical\_density |  0.2521354|  0.0684946|  0.2051422|  0.0498267|
| IASO\_6-10.15H   | Pythium oopapillum                  | mefenoxam | poison\_plate    |  0.1369200|  0.0624607|  0.1212678|  0.0535314|
| ILSO\_3-21A      | Pythium ultimum var. sporangiiferum | mefenoxam | optical\_density |  0.0171647|  0.0110661|  0.0164300|  0.0103643|
| ILSO\_3-21A      | Pythium ultimum var. sporangiiferum | mefenoxam | poison\_plate    |  0.0444185|  0.0319651|  0.0486010|  0.0361461|
| INSO\_3-10       | Pythium pleroticum                  | mefenoxam | optical\_density |  0.0354960|  0.0029223|  0.0307489|  0.0010692|
| INSO\_3-10       | Pythium pleroticum                  | mefenoxam | poison\_plate    |  0.0199470|  0.0063061|  0.0201349|  0.0064933|
| INSO\_3-43       | Pythium spinosum                    | mefenoxam | optical\_density |  0.2043519|  0.0807916|  0.1498540|  0.0461110|
| INSO\_3-43       | Pythium spinosum                    | mefenoxam | poison\_plate    |  0.0983442|  0.0584485|  0.0984697|  0.0614379|
| INSO\_4-40       | Pythium spinosum                    | mefenoxam | optical\_density |  0.0905818|  0.0403220|  0.0878968|  0.0396769|
| INSO\_4-40       | Pythium spinosum                    | mefenoxam | poison\_plate    |  0.0502867|  0.0279427|  0.0467290|  0.0252301|
| INSO\_5-50       | Pythium ultimum var. sporangiiferum | mefenoxam | optical\_density |  0.0193989|  0.0011695|  0.0195310|  0.0012364|
| INSO\_5-50       | Pythium ultimum var. sporangiiferum | mefenoxam | poison\_plate    |  0.0394601|  0.0130432|  0.0409852|  0.0183301|
| KSSO\_6-30       | Pythium ultimum var. sporangiiferum | mefenoxam | optical\_density |  0.0155855|  0.0055027|  0.0140713|  0.0044018|
| KSSO\_6-30       | Pythium ultimum var. sporangiiferum | mefenoxam | poison\_plate    |  0.0729989|  0.0409378|  0.0843618|  0.0553740|
| V-IASO2\_6-55\_1 | Phytopythium litorale               | mefenoxam | optical\_density |  0.9960691|  0.2530901|  0.4135866|  0.1211882|
| V-IASO2\_6-55\_1 | Phytopythium litorale               | mefenoxam | poison\_plate    |  2.4032623|  0.9107058|  0.9779876|  0.5341688|
| V-SDSO2\_1-53    | Phytophthora sojae                  | mefenoxam | optical\_density |  0.0846261|  0.0036287|  0.0807847|  0.0015759|
| V-SDSO2\_1-53    | Phytophthora sojae                  | mefenoxam | poison\_plate    |  0.0525604|  0.0387617|  0.0529119|  0.0391127|

Lets do an ANOVA for each chemistry and the effect of method. We have log transformed these data for homogeneity of variance.

model 1: method as fixed effect and isolate as random effect

We are treating isolate as a random effect because we we sampled these isolates from a larger possible population of isolates and we want to generalize over all isolates

``` r
lm_mef <- lmer(absolute ~ method + (1|is), data = ec50_cor[ec50_cor$chem == "mefenoxam",])
hist(residuals(lm_mef)) # not normally distributed residuals 
```

![](MethodCorrelation_files/figure-markdown_github/unnamed-chunk-4-1.png)

``` r
qqnorm(resid(lm_mef), main = "not log transformed"); qqline(resid(lm_mef))
```

![](MethodCorrelation_files/figure-markdown_github/unnamed-chunk-4-2.png)

``` r
lm_mef <- lmer(log(absolute) ~ method + (1|is), data = ec50_cor[ec50_cor$chem == "mefenoxam",])
hist(residuals(lm_mef)) # log transformation is good
```

![](MethodCorrelation_files/figure-markdown_github/unnamed-chunk-4-3.png)

``` r
qqnorm(resid(lm_mef), main = "log transformed"); qqline(resid(lm_mef))
```

![](MethodCorrelation_files/figure-markdown_github/unnamed-chunk-4-4.png)

``` r
lmerTest::anova(lm_mef, test.statistic="F", type = 2) # using type II ANOVA for unbalanced data. Some isolates have more technical replicates than others. So the mean over all isolates is different.  
```

    ## Analysis of Variance Table of type II  with  Satterthwaite 
    ## approximation for degrees of freedom
    ##         Sum Sq Mean Sq NumDF  DenDF F.value Pr(>F)
    ## method 0.71795 0.71795     1 84.159  1.1889 0.2787

``` r
plot(lm_mef, type = c("p", "smooth"), id = 0.05) # regression diagnostics
```

![](MethodCorrelation_files/figure-markdown_github/unnamed-chunk-4-5.png)

``` r
lsmeans_mef <- lsmeans::lsmeans(lm_mef, "method")
plot(lsmeans_mef)
```

![](MethodCorrelation_files/figure-markdown_github/unnamed-chunk-4-6.png) This is what we expected, no significant differnces for the method.

Lets do the same for ethaboxam.

``` r
lm_eth <- lmer(absolute ~ method + (1|is), data = ec50_cor[ec50_cor$chem == "ethaboxam",])
hist(residuals(lm_eth)) # not normally distributed residuals 
```

![](MethodCorrelation_files/figure-markdown_github/unnamed-chunk-5-1.png)

``` r
qqnorm(resid(lm_eth), main = "not log transformed"); qqline(resid(lm_eth))
```

![](MethodCorrelation_files/figure-markdown_github/unnamed-chunk-5-2.png)

``` r
lmerTest::anova(lm_eth, test.statistic="F", type = 2) # using type II ANOVA for unbalanced data. Some isolates have more technical replicates than others. So the mean over all isolates is different.  
```

    ## Analysis of Variance Table of type II  with  Satterthwaite 
    ## approximation for degrees of freedom
    ##         Sum Sq Mean Sq NumDF  DenDF F.value Pr(>F)
    ## method 0.10793 0.10793     1 82.088  1.1123 0.2947

``` r
plot(lm_eth, type = c("p", "smooth"), id = 0.05)# regression diagnostics
```

![](MethodCorrelation_files/figure-markdown_github/unnamed-chunk-5-3.png)

``` r
lsmeans_eth <- lsmeans::lsmeans(lm_eth, "method")
plot(lsmeans_eth)
```

![](MethodCorrelation_files/figure-markdown_github/unnamed-chunk-5-4.png)

Now lets put chemistry in as a fixed effect and fit the interaction bewtween chemistry and method.

We have log-transformed these data for homogienity of variance

``` r
lm3 <- lm(log(absolute) ~ is * chem * method, data = ec50_cor)
hist(residuals(lm3)) # log transformation is good
```

![](MethodCorrelation_files/figure-markdown_github/unnamed-chunk-6-1.png)

``` r
qqnorm(resid(lm3), main = "log transformed"); qqline(resid(lm3))
```

![](MethodCorrelation_files/figure-markdown_github/unnamed-chunk-6-2.png)

``` r
anova(lm3)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(absolute)
    ##                 Df Sum Sq Mean Sq  F value                Pr(>F)    
    ## is              38 365.52   9.619  17.6018 < 0.00000000000000022 ***
    ## chem             1  64.77  64.770 118.5218 < 0.00000000000000022 ***
    ## method           1   0.76   0.758   1.3869              0.241245    
    ## is:chem          6  11.28   1.880   3.4411              0.003579 ** 
    ## is:method       38  20.44   0.538   0.9843              0.505606    
    ## chem:method      1   0.24   0.241   0.4403              0.508217    
    ## is:chem:method   6   0.91   0.152   0.2778              0.946444    
    ## Residuals      121  66.12   0.546                                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(lm3)
```

![](MethodCorrelation_files/figure-markdown_github/unnamed-chunk-6-3.png)![](MethodCorrelation_files/figure-markdown_github/unnamed-chunk-6-4.png)![](MethodCorrelation_files/figure-markdown_github/unnamed-chunk-6-5.png)![](MethodCorrelation_files/figure-markdown_github/unnamed-chunk-6-6.png)

``` r
lsmeans3 <- lsmeans::lsmeans(lm3, c("is", "chem", "method"))
plot(lsmeans3)
```

![](MethodCorrelation_files/figure-markdown_github/unnamed-chunk-6-7.png)

Lets do correlation analysis between the two methods separated by chemistry.

We are testing the correlation of the absolute EC50s between the two methods. We are going to use spearman's correlation coeffiecient since it is rank based it can handle outliers with high leverage.

``` r
cor_mef <- lm(mean.abs.pp ~ mean.abs.od, data = EC50[EC50$chem == "mefenoxam",])
summary(cor_mef)
```

    ## 
    ## Call:
    ## lm(formula = mean.abs.pp ~ mean.abs.od, data = EC50[EC50$chem == 
    ##     "mefenoxam", ])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.69793 -0.08873  0.00269  0.05529  0.91905 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value    Pr(>|t|)    
    ## (Intercept) -0.07988    0.07755  -1.030       0.317    
    ## mean.abs.od  1.57027    0.21667   7.247 0.000000972 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2866 on 18 degrees of freedom
    ## Multiple R-squared:  0.7448, Adjusted R-squared:  0.7306 
    ## F-statistic: 52.52 on 1 and 18 DF,  p-value: 0.000000972

``` r
par(mfrow = c(2,2))
plot(cor_mef)
```

![](MethodCorrelation_files/figure-markdown_github/unnamed-chunk-8-1.png)

``` r
ttest(cor_mef, 1, 0) # tests if intercept is significantly different than 0
```

    ## [1] 0.3166105

``` r
ttest(cor_mef, 2, 1) # tests if slope (beta) is significantly different than 1
```

    ## [1] 0.01692273

``` r
(summary(cor_mef)[[4]][2] - 1)*100 #percent Bias
```

    ## [1] 57.02711

There is a significant linear relationship between the mean absolute EC50 using either method. Since it looked like there were some points with a bit of leverage we will use spearman's correlation to test the significance of the correlation.We will also look at the spearman:pearson correlation ratio to see if the correlation is more monotonic or linear.

``` r
spear.cor <- cor.test(EC50$mean.abs.pp[EC50$chem == "mefenoxam"], 
         EC50$mean.abs.od[EC50$chem == "mefenoxam"], 
         method = "spearman")
pear.cor <- cor.test(EC50$mean.abs.pp[EC50$chem == "mefenoxam"], 
         EC50$mean.abs.od[EC50$chem == "mefenoxam"], 
         method = "pearson")
spear.cor
```

    ## 
    ##  Spearman's rank correlation rho
    ## 
    ## data:  EC50$mean.abs.pp[EC50$chem == "mefenoxam"] and EC50$mean.abs.od[EC50$chem == "mefenoxam"]
    ## S = 358, p-value = 0.0003757
    ## alternative hypothesis: true rho is not equal to 0
    ## sample estimates:
    ##       rho 
    ## 0.7308271

``` r
pear.cor
```

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  EC50$mean.abs.pp[EC50$chem == "mefenoxam"] and EC50$mean.abs.od[EC50$chem == "mefenoxam"]
    ## t = 7.2471, df = 18, p-value = 0.000000972
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  0.6802565 0.9447289
    ## sample estimates:
    ##       cor 
    ## 0.8629927

Since the spearman correlation coeficient is lower than the pearson coefficient, this indicates we have more of a linear relationship than a monotonic one. This is a good thing because we would expect a perfect linear relationship between the methods.

``` r
cor_eth <- lm(mean.abs.pp ~ mean.abs.od, data = EC50[EC50$chem == "ethaboxam",])
summary(cor_eth)
```

    ## 
    ## Call:
    ## lm(formula = mean.abs.pp ~ mean.abs.od, data = EC50[EC50$chem == 
    ##     "ethaboxam", ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -0.5862 -0.1956 -0.1347  0.2654  0.9172 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value  Pr(>|t|)    
    ## (Intercept)   0.1342     0.1364   0.984     0.335    
    ## mean.abs.od   0.8869     0.1625   5.457 0.0000131 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.4029 on 24 degrees of freedom
    ## Multiple R-squared:  0.5538, Adjusted R-squared:  0.5352 
    ## F-statistic: 29.78 on 1 and 24 DF,  p-value: 0.00001312

``` r
par(mfrow = c(2,2))
plot(cor_eth)
```

![](MethodCorrelation_files/figure-markdown_github/unnamed-chunk-11-1.png)

``` r
ttest(cor_eth, 1, 0) # tests if intercept is significantly different than 0
```

    ## [1] 0.3350091

``` r
ttest(cor_eth, 2, 1) # tests if slope (beta) is significantly different than 1
```

    ## [1] 0.4931334

``` r
(summary(cor_eth)[[4]][2] - 1)*100 #percent Bias
```

    ## [1] -11.31059

There is a significant linear relationship between the mean absolute EC50 using either method. Since it looked like there were some points with a bit of leverage we will use spearman's correlation to test the significance of the correlation.

``` r
spear.cor <- cor.test(EC50$mean.abs.pp[EC50$chem == "ethaboxam"], 
         EC50$mean.abs.od[EC50$chem == "ethaboxam"], 
         method = "spearman")
pear.cor <- cor.test(EC50$mean.abs.pp[EC50$chem == "ethaboxam"], 
         EC50$mean.abs.od[EC50$chem == "ethaboxam"], 
         method = "pearson")
spear.cor
```

    ## 
    ##  Spearman's rank correlation rho
    ## 
    ## data:  EC50$mean.abs.pp[EC50$chem == "ethaboxam"] and EC50$mean.abs.od[EC50$chem == "ethaboxam"]
    ## S = 482, p-value = 0.000001932
    ## alternative hypothesis: true rho is not equal to 0
    ## sample estimates:
    ##       rho 
    ## 0.8352137

``` r
pear.cor
```

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  EC50$mean.abs.pp[EC50$chem == "ethaboxam"] and EC50$mean.abs.od[EC50$chem == "ethaboxam"]
    ## t = 5.4573, df = 24, p-value = 0.00001312
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  0.5012888 0.8783242
    ## sample estimates:
    ##       cor 
    ## 0.7441461

Since the spearman correlation coeficient is higher than the pearson coefficient, this indicates we have more of a monotonic relationship than a linear one. This is ok, but we would like to see it flipped around.

``` r
EC50_spec <- ddply(EC50, c("species_od","chem", "method"), 
      summarize, 
      mean.abs.pp = mean(mean.abs.pp, na.rm = TRUE),
      std.abs.pp = std.error(std.abs.pp, na.rm = TRUE),
      mean.abs.od = mean(mean.abs.od, na.rm = TRUE),
      std.abs.od = std.error(std.abs.od, na.rm = TRUE))

ggplot(EC50_spec, aes(mean.abs.pp, mean.abs.od)) + 
  geom_point(aes(colour = species_od, shape = chem)) +
  #geom_errorbarh(aes(xmax = mean.abs.pp + std.abs.pp, xmin = mean.abs.pp - std.abs.pp, height = .01)) +
  #geom_errorbar(aes(ymax = mean.abs.od + std.abs.od, ymin = mean.abs.od - std.abs.od, width = .01)) +
  geom_smooth(method = "lm", se = FALSE) + 
  theme_bw() +
  guides(colour = guide_legend(title = "Species"),
         shape = guide_legend(title = "Chemistry")) + 
  xlab(expression(bold("Poison Plate EC"[50]))) + 
  ylab(expression(bold("Optical Density EC"[50]))) + 
  scale_y_continuous(limits = c(0, 2), breaks = c(0, 0.5, 1, 1.5, 2)) +
  scale_x_continuous(limits = c(0, 2), breaks = c(0, 0.5, 1, 1.5, 2)) +
  theme(axis.text.x = element_text(family = "Times New Roman", size = 10, face = "bold"),
          axis.text.y = element_text(family = "Times New Roman", size = 10, face = "bold"),
          axis.title.x = element_text(family = "Times New Roman", size = 10, face = "bold"),
          axis.title.y = element_text(family = "Times New Roman", size = 10, face = "bold"),
          legend.text = element_text(family = "Times New Roman", size = 10, face = "bold.italic"),
          legend.key = element_blank(),
          legend.title = element_text(family = "Times New Roman", size = 15, face="bold"),
   strip.text.x = element_text(family = "Times New Roman",size = 15, face = "bold"))
```

![](MethodCorrelation_files/figure-markdown_github/unnamed-chunk-14-1.png)

``` r
summary(lm(mean.abs.pp ~ mean.abs.od, data = EC50_spec))
```

    ## 
    ## Call:
    ## lm(formula = mean.abs.pp ~ mean.abs.od, data = EC50_spec)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.36128 -0.11452 -0.02151  0.06983  0.50524 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value     Pr(>|t|)    
    ## (Intercept) -0.04456    0.07842  -0.568        0.577    
    ## mean.abs.od  1.16557    0.12007   9.707 0.0000000239 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2351 on 17 degrees of freedom
    ## Multiple R-squared:  0.8472, Adjusted R-squared:  0.8382 
    ## F-statistic: 94.23 on 1 and 17 DF,  p-value: 0.00000002389

Based on the above analysis, we will select the relative growth at

``` r
mef.cor$chem <- "Mefenoxam"
eth.cor$chem <- "Ethaboxam"
cor.plot.relgrowth <- rbind.data.frame(eth.cor[eth.cor$conc == 0.5,], mef.cor[mef.cor$conc == 0.5,])
EC50$chem2 <- factor(EC50$chem, labels = c("Ethaboxam", "Mefenoxam"))
length(levels(cor.plot.relgrowth$is))
```

    ## [1] 46

``` r
 p <- ggplot(cor.plot.relgrowth, aes(y = odrelgrowth, x = ppmeanrelgrowth)) +
    geom_point(aes(color = factor(species))) +
    guides(colour = guide_legend(title = "Species")) +
    scale_y_continuous(limits = c(0, 125), breaks = c(0, 25, 50, 75, 100, 125)) +
    scale_x_continuous(limits = c(0, 125), breaks = c(0, 25, 50, 75, 100, 125)) +
     geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x, fullrange = TRUE) +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(..eq.label.., ..adj.rr.label.., sep = "~~~")), 
               parse = TRUE) +
    xlab("% Relative Growth, Poison Plate") + 
    ylab("% Relative Growth, Optical Density") + 
    theme_classic() +
    theme(axis.text.x = element_text(size = 10, face = "bold"),
          axis.text.y = element_text(size = 10, face = "bold"),
          axis.title.x = element_text(size = 10, face = "bold"),
          axis.title.y = element_text(size = 10, face = "bold"),
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          legend.text = element_text(size = 6, face = "bold.italic"),
          legend.key = element_blank(),
          legend.title = element_text(size = 10, face="bold"),
          legend.position = "bottom", 
          strip.text.x = element_text(size = 15, face = "bold")) + 
   facet_wrap(~chem)
p1 <- ggplot(EC50, aes(mean.abs.pp, mean.abs.od)) + 
  geom_point(aes(color = factor(species_od))) +
 geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x, fullrange = TRUE) +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(..eq.label.., ..adj.rr.label.., sep = "~~~")), 
               parse = TRUE) +
  guides(colour = guide_legend(title = "Species")) +
  xlab(expression(bold("Poison plate EC"[50]))) + 
  ylab(expression(bold("Optical density EC"[50]))) + 
  theme_classic() + 
  theme(axis.text.x = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 10, face = "bold"),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        legend.text = element_text(size = 9, face = "bold.italic"),
        legend.key = element_blank(),
        legend.title = element_text(size = 10, face="bold"),
        strip.text.x = element_text(size = 15, face = "bold")) +
  facet_wrap(~chem2)
```

``` r
#extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(p)
```

``` r
p3 <- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                         p + theme(legend.position="none"),
                         nrow=2),
             mylegend, nrow=2,heights=c(5, 1))
```

![](MethodCorrelation_files/figure-markdown_github/unnamed-chunk-18-1.png)

``` r
print(p3)
```

    ## TableGrob (2 x 1) "arrange": 2 grobs
    ##   z     cells    name              grob
    ## 1 1 (1-1,1-1) arrange   gtable[arrange]
    ## 2 2 (2-2,1-1) arrange gtable[guide-box]

``` r
ggsave(file="FigureCorrelation.pdf", plot = p3, width = 9, height = 9, dpi = 300)
```
