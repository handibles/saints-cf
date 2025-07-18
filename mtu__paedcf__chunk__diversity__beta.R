
##   beta diversity & co.
      
      rm(list=ls())
      
      source("analysis/mtu__paedcf__spreadsheeting__shortcuts.R")
      source("analysis/background_code/R__fns_jfg/fn_definitions.R")
      
      set.seed(2205)
      
      library(ggplot2)
      library(ggpubr)
      library(vegan)
      library(parallel)
      
      
## =======================================
      
      print("applying theme_update in spreadsheeting_shortcuts")
      theme_set( theme_minimal())
      theme_update(
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 14, colour = "grey20"),
        axis.title.x = element_text(size = 14, colour = "grey20"),
        axis.line = element_line(colour = "grey80", size = 0.2),
        #
        legend.box = "horizontal",
        legend.position = "right",
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.75, "cm"),
        #
        panel.grid.major.x = element_line(colour = "lightskyblue3", size = 0.1),
        panel.grid.major.y = element_line(colour = "lightskyblue3", size = 0.1),
        panel.grid.minor.x = element_line(colour = "lightskyblue3", size = 0.1),
        panel.grid.minor.y = element_line(colour = "lightskyblue3", size = 0.1),
        panel.spacing.x = unit(1, "lines"),
        #
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),  # if plotting class using super
        plot.tag = element_text(size = 25, face = "bold", margin = margin()),  #  margin = margin(l = 10)),  # if plotting class using super
        plot.title = element_text(size = 18, face = "bold"),  #  margin = margin(l = 10)),  # if plotting class using super
        plot.subtitle = element_text(size=14),
        plot.background = element_rect(fill = "white", colour = "white"),
        #
        strip.background = element_blank(),
        strip.text.x = element_text(angle = 0, size = 16),   # Top
        strip.text.y = element_text(angle = 0, size =14),
        NULL
      )
    
    
      samp_type <- c("BAL", "NS", "TS")
      
      head(mgdat) ; dim(mgdat)
      mgdat <- dplyr::filter(mgdat, !grepl("LIBNEG", rownames(mgdat)))
      mgfeat <- mgfeat[ , rownames(mgdat) ]
      mgfeat_ra <- mgfeat_ra[ , rownames(mgdat) ]
      
      dim( mgfeat_raka <- k_A( mgfeat_ra, k = 0.001, A = 0.1))
      bc_dist_raka <- vegdist( t(mgfeat_raka), method = "bray")
      ja_dist_raka <- vegdist( t(mgfeat_raka), method = "jaccard")

      
##  removing H. sapiens DNA and retest...   ====================================

      dim(nhsfeat_ra)
      nhsfeat_raka  <- k_A( nhsfeat_ra, k = 0.001, A = 0.1)
      dim(nhsfeat_raka <- nhsfeat_raka[ , colSums(nhsfeat_raka) > 0 ])

      ja_dist_raka__nHsap <- vegdist( t(nhsfeat_ra[ , colSums(nhsfeat_ra) > 0 ]), method = "jaccard")

      euc_dist_raka__nHsap <- vegdist( t(nhsfeat_clr[ ,  ]), method = "euclidean")

            
##  heteroskedasticity  ---  noHsap   ==========================================

      # try it with euclidean / CLR
        ## total
        permutest( betadisper(
          euc_dist_raka__nHsap,
          mgdat$substrate ),
          permutations = 100 )
      
            # >  Permutation test for homogeneity of multivariate dispersions
            # >  Permutation: free
            # >  Number of permutations: 100
            # >  
            # >  Response: Distances
            # >            Df Sum Sq Mean Sq      F N.Perm   Pr(>F)   
            # >  Groups     2 9742.2  4871.1 118.82    100 0.009901 **
            # >  Residuals 85 3484.5    41.0                          
            # >  ---
            # >  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
      
      # BAL-NS homoskedastic
      lapply( samp_type, function(aa){
        bb_dat <- dplyr::filter( mgdat, substrate != aa )
        cc_feat <- nhsfeat_clr[ , rownames(bb_dat)]
        dd_dist <- vegdist( t(cc_feat), method = "euclidean")
        permutest( betadisper( dd_dist, bb_dat$substrate ), permutations = 100 )
      })
      
            # >  [[1]]
            # >  
            # >  Permutation test for homogeneity of multivariate dispersions
            # >  Permutation: free
            # >  Number of permutations: 100
            # >  
            # >  Response: Distances
            # >            Df Sum Sq Mean Sq      F N.Perm   Pr(>F)   
            # >  Groups     1 7978.0    7978 242.02    100 0.009901 **
            # >  Residuals 54 1780.1      33                          
            # >  ---
            # >  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
            # >  
            # >  [[2]]
            # >  
            # >  Permutation test for homogeneity of multivariate dispersions
            # >  Permutation: free
            # >  Number of permutations: 100
            # >  
            # >  Response: Distances
            # >            Df Sum Sq Mean Sq      F N.Perm   Pr(>F)   
            # >  Groups     1 7600.9  7600.9 170.98    100 0.009901 **
            # >  Residuals 52 2311.7    44.5                          
            # >  ---
            # >  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
            # >  
            # >  [[3]]
            # >  
            # >  Permutation test for homogeneity of multivariate dispersions
            # >  Permutation: free
            # >  Number of permutations: 100
            # >  
            # >  Response: Distances
            # >            Df  Sum Sq Mean Sq      F N.Perm Pr(>F)
            # >  Groups     1    1.42   1.421 0.0316    100 0.8218
            # >  Residuals 64 2877.24  44.957                     
      
      
##  adonis by substrate   =================================================================
      
      # BAL-NS no blocking structure at this high level
      lapply( samp_type, function(aa){
        bb_dat <- dplyr::filter( mgdat[ labels(euc_dist_raka__nHsap) , ], substrate != aa )
        cc_dist <- as.matrix(euc_dist_raka__nHsap)[ rownames(bb_dat) , rownames(bb_dat)]
        adonis2( cc_dist ~ substrate, data = bb_dat )
      })    
            # >  [[1]]
            # >  Permutation test for adonis under reduced model
            # >  Permutation: free
            # >  Number of permutations: 999
            # >  
            # >  adonis2(formula = cc_dist ~ substrate, data = bb_dat)
            # >           Df SumOfSqs      R2    F Pr(>F)    
            # >  Model     1    11096 0.31472 24.8  0.001 ***
            # >  Residual 54    24160 0.68528                
            # >  Total    55    35256 1.00000                
            # >  ---
            # >  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
            # >  
            # >  [[2]]
            # >  Permutation test for adonis under reduced model
            # >  Permutation: free
            # >  Number of permutations: 999
            # >  
            # >  adonis2(formula = cc_dist ~ substrate, data = bb_dat)
            # >           Df SumOfSqs      R2      F Pr(>F)    
            # >  Model     1     9824 0.28533 20.761  0.001 ***
            # >  Residual 52    24605 0.71467                  
            # >  Total    53    34429 1.00000                  
            # >  ---
            # >  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
            # >  
            # >  [[3]]
            # >  Permutation test for adonis under reduced model
            # >  Permutation: free
            # >  Number of permutations: 999
            # >  
            # >  adonis2(formula = cc_dist ~ substrate, data = bb_dat)
            # >           Df SumOfSqs      R2      F Pr(>F)  
            # >  Model     1    144.2 0.02791 1.8374  0.041 *
            # >  Residual 64   5023.4 0.97209                
            # >  Total    65   5167.6 1.00000                
            # >  ---
            # >  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
      
      
##  heteroskedasticity     ================================================================

## total
permutest( betadisper( 
  vegdist( t(mgfeat_ra), method = "jaccard"), 
  mgdat$substrate ), 
  permutations = 100 )


# BAL-NS homoskedastic
lapply( samp_type, function(aa){
  bb_dat <- dplyr::filter( mgdat, substrate != aa )
  cc_feat <- mgfeat_raka[ , rownames(bb_dat)]
  dd_dist <- vegdist( t(cc_feat), method = "jaccard")
  permutest( betadisper( dd_dist, bb_dat$substrate ), permutations = 100 )
})    


  ## MTU rewrite 

    # >  [[1]]     no BAL
    # >  
    # >  Permutation test for homogeneity of multivariate dispersions
    # >  Permutation: free
    # >  Number of permutations: 100
    # >  
    # >  Response: Distances
    # >  Df Sum Sq Mean Sq      F N.Perm Pr(>F)  
    # >  Groups     1 0.3500 0.34999 3.8979    100 0.0396 *
    # >    Residuals 54 4.8487 0.08979                       
    # >  ---
    # >    Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    # >  
    # >  
    # >  
    # >  [[2]]     no NS
    # >  
    # >  Permutation test for homogeneity of multivariate dispersions
    # >  Permutation: free
    # >  Number of permutations: 100
    # >  
    # >  Response: Distances
    # >  Df Sum Sq Mean Sq      F N.Perm Pr(>F)  
    # >  Groups     1 0.6400 0.63996 8.7329    100 0.0297 *
    # >    Residuals 52 3.8106 0.07328                       
    # >  ---
    # >    Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    # >  
    # >  
    # >  
    # >  [[3]]    no TS
    # >  
    # >  Permutation test for homogeneity of multivariate dispersions
    # >  Permutation: free
    # >  Number of permutations: 100
    # >  
    # >  Response: Distances
    # >  Df Sum Sq  Mean Sq      F N.Perm Pr(>F)
    # >  Groups     1 0.0588 0.058767 0.4495    100  0.505
    # >  Residuals 64 8.3674 0.130740                     
    
      

##  adonis by substrate   =================================================================

# BAL-NS no blocking structure at this high level
lapply( samp_type, function(aa){
  bb_dat <- dplyr::filter( mgdat, substrate != aa )
  cc_dist <- as.matrix(ja_dist)[ rownames(bb_dat) , rownames(bb_dat)]
  adonis2( cc_dist ~ substrate, data = bb_dat )
})    

    ## MTU re-write

      # >  
      # >  [[1]]       ## no BAL
      # >  
      # >  Permutation test for adonis under reduced model
      # >  Permutation: free
      # >  Number of permutations: 999
      # >  
      # >  adonis2(formula = cc_dist ~ substrate, data = bb_dat)
      # >  Df SumOfSqs      R2      F Pr(>F)    
      # >  Model     1   6.4282 0.34288 28.176  0.001 ***
      # >  Residual 54  12.3196 0.65712                  
      # >  Total    55  18.7478 1.00000                  
      # >  ---
      # >  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
      # >  
      # >  
      # >  [[2]]       ## no NS
      # >  
      # >  Permutation test for adonis under reduced model
      # >  Permutation: free
      # >  Number of permutations: 999
      # >  
      # >  adonis2(formula = cc_dist ~ substrate, data = bb_dat)
      # >  Df SumOfSqs      R2      F Pr(>F)    
      # >  Model     1    6.543 0.37641 31.388  0.001 ***
      # >  Residual 52   10.840 0.62359                  
      # >  Total    53   17.383 1.00000                  
      # >  ---
      # >  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
      # >  
      # >  
      # >  [[3]]    ## no TS
      # >  
      # >  Permutation test for adonis under reduced model
      # >  Permutation: free
      # >  Number of permutations: 999
      # >  
      # >  adonis2(formula = cc_dist ~ substrate, data = bb_dat)
      # >  Df SumOfSqs      R2      F Pr(>F)
      # >  Model     1   0.2967 0.02548 1.6736  0.141
      # >  Residual 64  11.3466 0.97452              
      # >  Total    65  11.6433 1.00000              
      # >  


