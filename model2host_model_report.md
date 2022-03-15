Risk mapping by climate, host trait and distribution
================
Yanjie Xu
15/03/2022

# Salmonella example

This is a report for the risk modelling using *Salmonella* as an
example. We only show the model training process here, as the model
prediction/extrapolation takes long (\~950,000 data points - 2 days).
The raw code including the full process is available at:
<https://github.com/yanjie-xu/pathogen.riskmap>.

## 1. Check the colinearity between predictors

use `vif` fuction:

``` r
vif(lmer(Npositive ~ scale(Temp)+scale(Prec)+
      scale(Clutch_MEAN)+scale(Maximum.longevity)+
      scale(Mass)+Habitat+Migration+Ntested+
        (1|animal), data=preva_trait, REML=FALSE))
```

    ##                              GVIF Df GVIF^(1/(2*Df))
    ## scale(Temp)              1.088479  1        1.043302
    ## scale(Prec)              1.075523  1        1.037074
    ## scale(Clutch_MEAN)       1.242252  1        1.114563
    ## scale(Maximum.longevity) 1.942145  1        1.393609
    ## scale(Mass)              1.515324  1        1.230985
    ## Habitat                  1.862775  5        1.064182
    ## Migration                1.373422  2        1.082558
    ## Ntested                  1.013381  1        1.006668

## 2.MCMCglmm

To enable multiple entries for single species, we define `species` as
`animal`. And define the subset birdtree as the `pedigree`.

The `animal` column should be the `Tiplabel` in birdtree instead of
scientific name.

``` r
df = preva_trait[preva_trait$Pathogen == pathogen[12],]
pathogen[12]
```

    ## [1] "Salmonella"

The model:

``` r
model = MCMCglmm(cbind(Npositive,Ntested-Npositive) ~ scale(Temp)+scale(Prec)+
                   scale(Clutch_MEAN)+scale(Maximum.longevity)+
                   scale(Mass)+Migration+Habitat, 
                 random=~animal, 
                 pedigree = phy,
                 family =  "multinomial2",
                 data=df,
                 nitt=133000, 
                 burnin=3000, 
                 thin=100)
```

``` r
summary(model)
```

    ## 
    ##  Iterations = 3001:132901
    ##  Thinning interval  = 100
    ##  Sample size  = 1300 
    ## 
    ##  DIC: 7145.854 
    ## 
    ##  G-structure:  ~animal
    ## 
    ##        post.mean l-95% CI u-95% CI eff.samp
    ## animal     3.107   0.1906    6.345      353
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units     2.837    1.694    4.009     1136
    ## 
    ##  Location effects: cbind(Npositive, Ntested - Npositive) ~ scale(Temp) + scale(Prec) + scale(Clutch_MEAN) + scale(Maximum.longevity) + scale(Mass) + Migration + Habitat 
    ## 
    ##                           post.mean l-95% CI u-95% CI eff.samp   pMCMC    
    ## (Intercept)                 -5.3532  -8.4239  -1.7519   1300.0 0.00462 ** 
    ## scale(Temp)                  0.8982   0.3940   1.3656   1300.0 < 8e-04 ***
    ## scale(Prec)                  0.2030  -0.1462   0.5148   1160.5 0.21077    
    ## scale(Clutch_MEAN)          -0.3935  -1.1038   0.3402   1300.0 0.25077    
    ## scale(Maximum.longevity)    -0.1011  -0.7719   0.5082   1300.0 0.78615    
    ## scale(Mass)                  0.5975   0.1164   1.0364   1136.9 0.01385 *  
    ## MigrationPartialMigratory    1.5243   0.1569   2.7929   1157.3 0.01538 *  
    ## MigrationSedentary           1.5320   0.1011   2.9129    965.3 0.02615 *  
    ## HabitatForest               -0.2588  -2.9297   2.4442   1300.0 0.82615    
    ## HabitatFreshwater           -0.5233  -3.7693   2.3591   1300.0 0.72615    
    ## HabitatGrassland            -0.7130  -3.3601   1.8777   1082.2 0.59231    
    ## HabitatHuman Modified       -0.6695  -3.2809   2.3043   1300.0 0.60308    
    ## HabitatMarine               -1.3539  -4.6529   1.7546   1300.0 0.38769    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Effective sample size \> 1000, OK: Your effective sample size should be
quite high aiming at 1000-2000. If not, add iterations to achieve a
comparable effective sample size. The `R structure` is the residual
structure. The `G structure` is the random effects structure.

## 3. Check model convergence

### 3.1 Check for model convergence for fixed effects

To make sure your model has converged, the trace plot should look like a
fuzzy caterpillar. Well mixed in our case, OK!

![](model2host_model_report_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->![](model2host_model_report_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->![](model2host_model_report_files/figure-gfm/unnamed-chunk-5-3.png)<!-- -->![](model2host_model_report_files/figure-gfm/unnamed-chunk-5-4.png)<!-- -->

### 3.2 Check for model convergence for random effects

![](model2host_model_report_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

Again mixed well - a little bit pattern in animal (species), OK?

### 3.3 Calculate R squared

``` r
R2(model)
```

    ## [1] 0.637
