## 

    rm(list=ls())

    source("analysis/mtu__paedcf__spreadsheeting__shortcuts.R")
    set.seed(2205)

        
    library(rtk)
    library(ggplot2)
    library(vegan)
    library(parallel)
    
    head(mgdat) ; dim(mgdat)
    mgdat <- dplyr::filter(mgdat, !grepl("LIBNEG", rownames(mgdat)))
    mgfeat <- mgfeat[ , rownames(mgdat) ]
    mgfeat_ra <- mgfeat_ra[ , rownames(mgdat) ]
    
    
  ## v 0.3.5.    no more humans, hot lizards only   -----------------
    
    ##   - done in spreadsheets

    dim(nhsfeat)
    dim(nhsfeat_ra)
    
    any( rowSums(nhsfeat) == 0)   # WELL THATS GOOD
    
    
## rarefaction   ========================================================

    ## repeat for non-Homo sapiens microbiome estimates. A lot more patchy, but probably a fairer estiamte.
    
      mgdat$substrate <- gsub("[a,b]$", "", mgdat$substrate, perl = TRUE)
      # ggplot(data=mgdat, aes(x=reorder(sample, countTotal), y=log10(countTotal), color=substrate)) + geom_point(size = 4) + labs(title = "sample depth")
    
      rtk_usable_depth <- 1000
      sum( rtk_usable_depth > sort(colSums(mgfeat)))   # two. not bad.
      sum( rtk_usable_depth > sort(colSums(nhsfeat)))  # 43. way worse.
      
      plot( log10(colSums(mgfeat)), log10(colSums(nhsfeat)), pch = 19, 
            col = c("dodgerblue", "orange", "chartreuse")[ factor(mgdat[ colnames(mgfeat), "substrate"])],
            xlab = "all reads", ylab = "non-H.sap reads")

      rtk_its <- 250
      # not locally manageable, unstable
      mgrar <- rtk(mgfeat,
                   margin = 2,
                   repeats = rtk_its,
                   threads = 7,
                   #              depth = seq(100, 7000, by = 1000),
                   depth = seq(1000, 2000, by = 100),   # from the min to the mean - no real need to go all the way out?
                   ReturnMatrix = rtk_its)
      ## PITILESS they called it
      saveRDS(mgrar, "output/mtu__paedcf__krak2-hostile-88__rar_1K-2K-1C__noHsap.RDS")
      mgrar <- readRDS("output/mtu__paedcf__krak2-hostile-88__rar_1K-2K-1C__noHsap.RDS")

      ## bad cuve
      # par( mfrow = c(1,3))
      # plot(mgrar, div = "richness", groups = mgdat$Gender,  )
      # plot(mgrar, div = "shannon", groups = mgdat$substrate)
      # plot(mgrar, div = "simpson", groups = mgdat$substrate)

      ## lines do not stabilise
      d1k <- t(sapply(mgrar[[ as.character(rtk_usable_depth) ]][["divvs"]], function(aa){   # aa <- mgrar[["6100"]][["divvs"]][[15]]
        bb_dmed <- lapply( aa[2:7], function(aaa){ median(aaa)})
        unlist(c( aa[[1]], bb_dmed))
      } ))
      rownames(d1k) <- d1k[,1] ; d1k <- t(apply(d1k[,-1], 1, as.numeric))
      colnames(d1k) <- c("richness", "shannon", "simpson", "invsimpson", "chao1", "evenness")
      head(d1k)

      miss_div <- rownames(d1k)[apply(d1k, 1, function(aa) all(is.na(aa)) )]
      miss_div_fill <- cbind(
        "richness" = apply( mgfeat[ , miss_div], 2, function(aa) sum(aa > 0)),
        "shannon" = vegan::diversity( t(mgfeat[ , miss_div]), index = "shannon"),
        "simpson" = vegan::diversity( t(mgfeat[ , miss_div]), index = "simpson"),
        "invsimpson" = vegan::diversity( t(mgfeat[ , miss_div]), index = "invsimpson"),
        "chao1" = vegan::estimateR( t(mgfeat[ , miss_div]))[2,],
        "evenness" = NA
      )
      head(miss_div_fill)
      # # Hsap
      # miss_div_fill <- rep(NA, 6)

    
## non pareil ::  sequence entropy   ================================================

    ## flipping great programme
    ## run in conda env geno
    # install.packages('Nonpareil');
    library(Nonpareil)
    files <- 'output/nonpareil_hostile'
    samps <- list.files(files)

    # NP rules
    Nonpareil.read_data( Nonpareil.curve( file = paste0("output/nonpareil_hostile/",samps[3]) ), correction.factor = TRUE )
    Nonpareil.curve( paste0("output/nonpareil_hostile/",samps[3]) )$diversity
    predict.Nonpareil.Curve(Nonpareil.curve( paste0("output/nonpareil_hostile/",samps[3]) ))   # the estimated coverage...
    # for a multiple sample set of curves, see...
    np_dat <- summary.Nonpareil.Set( Nonpareil.set( paste0("output/nonpareil_hostile/",samps) ))

    View(np_dat)

    ## summary.Nonpareil.Curve {Nonpareil}
      # >   -    kappa: "Redundancy" value of the entire dataset.
      # >   -    C: Average coverage of the entire dataset.
      # >   -    LRstar: Estimated sequencing effort required to reach the objective average coverage (star, 95
      # >   -    LR: Actual sequencing effort of the dataset.
      # >   -    modelR: Pearson's R coefficient betweeen the rarefied data and the projected model.
      # >   -    diversity: Nonpareil sequence-diversity index (Nd). This value's units are the natural logarithm
      # >   - - -  of the units of sequencing effort (log-bp), and indicates the inflection point of the fitted model
      # >   - - -  for the Nonpareil curve. If the fit doesn't converge, or the model is not estimated, the value is zero (0).

    saveRDS(np_dat, "output/mtu__paedcf__nonpar-hostile-88.RDS")

    ggplot( np_dat) +
        geom_point( shape = 21, alpha = 0.6, aes(x = log10(LR), y = log10(LRstar), size = C, fill = diversity)) +
        scale_fill_viridis_c() +
        labs( x = "LR - actual seq", y = "LRstar - amount of further seq required", title = "most are grossly underCHARACTERISED, with C(overage) < 0.2") +
    ggplot( np_dat) +
        geom_point( shape = 21, alpha = 0.6, aes(size = log10(LR), fill = log10(LRstar), x = C, y = diversity)) +
        labs( x = "Coverage", y = "diversity", title = "coverage and diversity inversely related (possibly?)")
    
    np_dat <- readRDS("output/mtu__paedcf__nonpar-hostile-88.RDS")

      	
## Combine  =================================================================================
    
    # bring together in correct order - first the d1k, minus the NA samples
    
    head(d1k)
    head(miss_div_fill)
    head(np_dat)

    dim(  d1k_miss <- rbind(
                      d1k[ !apply(d1k, 1, function(aa){all(is.na(aa))}) , ],
                      miss_div_fill
      ) )
    rownames( np_dat ) <- gsub( "-",  ".", gsub("_S.*", "", rownames( np_dat )))

    dim(
      mgdd <- data.frame( mgdat[ intersect( rownames(mgdat), rownames(d1k_miss) ) , ],
                        d1k_miss[  intersect( rownames(mgdat), rownames(d1k_miss) ), ],    # intersect( rownames(mgdat), rownames(d1k_miss) )
                        apply( np_dat[ intersect( rownames(mgdat), rownames(d1k_miss)) , ], 2, as.numeric)
      ) )

    head(mgdd)
    saveRDS(mgdd, "output/mtu__paedcf__metadata_88-43_A-Div__noHsap.RDS")
    mgdd <- readRDS("output/mtu__paedcf__metadata_88-43_A-Div__noHsap.RDS")


## reads per substrate tested, and don';'t differ
    
                        
## alpha testing   -------------------------------------
    
    # can use models for this if residuals are norm, n'es pas?
    head(mgdd)
    mgdd$BMI.z.Score <- as.numeric(mgdd$BMI.z.Score)
    
    colnames(mgdd)
    
    ## NP previously referred to , "np.a" and "coverage"  - coverage is now C, but np.a is.. Diversity, based on the labels previoulsy applied. 
    divs <- c("richness", "shannon", "invsimpson", "chao1", "C", "diversity")   # #"faith_pd", 
    mdivs <- mgdd[ , c("substrate", "Genotype", "participant", "Gender", "age.at.sampling", "countTotal", "proph_bin", "modul_bin", "BMI.z.Score", divs) ] 
    head(mdivs)
    
    # # not particularly parametric. and probably wrong given the differences in ID
    md_m <- reshape2::melt(mdivs)
    ggplot(md_m, aes(x = value)) +
      facet_wrap(variable~., scale = "free") +
      geom_density(aes(fill = substrate), alpha = 0.5)
    
    
    library(car)
    library(e1071)
    # normal?
    t(sapply( divs, function(aa){ # aa <- divs[3]
      # hist(mdivs[ , aa])
      bb_lm <- lm( mdivs[,aa] ~ mgdat[ rownames(mdivs), "substrate"] )
      cc_v <- unlist(car::leveneTest(bb_lm))[ c(3,5)]
      dd_v <- unlist(shapiro.test(residuals(bb_lm)))[ c(1,2)]
      c(cc_v, dd_v)
    }))
    
    ## non-normal, non-normal residuals - use KW + post hoc
    adiv_kw <- t(sapply( divs, function(aa){ # aa <- divs[3]
      # hist(mdivs[ , aa])
      bb_kw <- kruskal.test( mdivs[,aa] ~ mgdat[ rownames(mdivs), "substrate"] )
      unlist(bb_kw)
    }))
    # flatten - all indices show significant differences in FDR
    adiv_kwf <- apply(adiv_kw, 2, function(aa){aa})[, c(1,3)]
    adiv_kwff <- cbind( adiv_kwf, "FDR" = p.adjust(adiv_kwf[,"p.value"] , method="fdr"))
    
    ## all tests, post hoc
    u1 <- unique(t(apply( expand.grid( unique(mgdat$substrate), unique(mgdat$substrate)), 1, function(aa)sort(aa))))
    (u2 <- u1[ apply(u1, 1, function(aa){length(unique(aa)) == 2}),])
    adiv_wt <- lapply( divs, function(aa){ # aa <- divs[2]
      bb_test <- t(apply( u2, 1, function(aaa){ # aaa <- u2[1,]
        c(
          "comparison" = paste0(aaa, collapse = "-"),
          unlist(wilcox.test(
            dplyr::filter(mdivs, substrate == aaa[1])[ , aa],
            dplyr::filter(mdivs, substrate == aaa[2])[ , aa]))[c(1,2)]
          )
      }))
      data.frame(
        "index" = aa, bb_test
      )
    })
    adiv_wtm <- do.call("rbind", adiv_wt)
    adiv_wtmf <- cbind(adiv_wtm, "FDR" = p.adjust( adiv_wtm[,"p.value"] , method = "fdr"))
    knitr::kable( adiv_out <- dplyr::filter( data.frame(adiv_wtmf,stringsAsFactors = FALSE), FDR < 0.05) )
    
      # >  |index      |comparison |statistic.W |p.value              |       FDR|
      # >  |:----------|:----------|:-----------|:--------------------|---------:|
      # >  |richness   |BAL-TS     |14.5        |1.90349621639817e-09 | 0.0000000|
      # >  |richness   |NS-TS      |17          |1.43725066036356e-09 | 0.0000000|
      # >  |shannon    |BAL-TS     |12          |1.52645146530306e-09 | 0.0000000|
      # >  |shannon    |NS-TS      |11          |8.82450280343067e-10 | 0.0000000|
      # >  |invsimpson |BAL-TS     |8           |9.79428594886615e-10 | 0.0000000|
      # >  |invsimpson |NS-TS      |6           |5.16915453448855e-10 | 0.0000000|
      # >  |chao1      |BAL-TS     |14.5        |1.90349621639817e-09 | 0.0000000|
      # >  |chao1      |NS-TS      |17          |1.43725066036356e-09 | 0.0000000|
      # >  |C          |BAL-NS     |344         |0.00983028158008985  | 0.0147454|
      # >  |C          |BAL-TS     |121         |2.07599374113239e-05 | 0.0000374|
      # >  |C          |NS-TS      |230         |0.0151792254288567   | 0.0210174|
      # >  |diversity  |BAL-TS     |133         |6.14062613586712e-05 | 0.0001005|
      # >  |diversity  |NS-TS      |89          |2.66036789416476e-07 | 0.0000005|
        
    # View(adiv_out)


## alpha div plots ================================================================
    
    ## evaluate, across parameters above
    # 
    # ggplot(mgdd, aes(x = substrate, y  = shannon)) +
    #   theme_bw() +
    #   facet_grid(.~modul_bin, space = "free") +
    #   # scale_fill_manual(values = var_cols) +
    #   geom_boxplot(aes(fill = substrate), outlier.shape = NA) +
    #   geom_point(aes(fill = substrate), position = position_jitterdodge(jitter.width = 0.25), shape = 21, size = 3.5, alpha = 0.8) + 
    #   # ggbeeswarm::geom_beeswarm(aes(fill = BMI.z.Score)) + 
    #   labs(y = "Shannon Diversity", x = "") +
    #   theme( axis.text.x = element_text(angle = 0))
    # 
    # 
    # ggplot(mgdd, aes(x = BMI.z.Score, y  = shannon)) +
    #   theme_bw() +
    #   facet_grid(.~proph_bin, space = "free") +
    #   # scale_fill_manual(values = var_cols) +
    #   # geom_boxplot(aes(fill = substrate), outlier.shape = NA) +
    #   geom_point(aes(fill = substrate), shape = 21, size = 3.5, alpha = 0.8) + 
    #   # ggbeeswarm::geom_beeswarm(aes(fill = BMI.z.Score)) + 
    #   labs(y = "Shannon Diversity", x = "") +
    #   theme( axis.text.x = element_text(angle = 0))
    #     
    #     
    
    # bees? not violins, sadly..
    # col_gend <-nr_col[c(13,14)]   # show_col(nr_col)
    
    mgdd_nona <- dplyr::filter( mgdd, !is.na(age.at.sampling))

    theme_update(
      plot.title = element_text(size = 16, face = "plain", hjust = 0.5),
      NULL
    )
        
    ad1 <-  ggplot(mgdd_nona, aes(x = substrate, y  = richness)) +
      geom_boxplot(aes(fill = substrate), outlier.shape = NA) +
      geom_point(aes(fill = substrate), position = position_jitterdodge(jitter.width = 0.25), shape = 21, size = 2, alpha = 0.8) + 
      labs(y = "n species observed", title = "species\nrichness", x = "") +
      theme( axis.text.x = element_text(angle = 90),
             axis.title.y = element_text(size = 10),
             legend.position = "none")
    
    ad2 <-  ggplot(mgdd_nona, aes(x = substrate, y  = shannon)) +
      geom_boxplot(aes(fill = substrate), outlier.shape = NA) +
      geom_point(aes(fill = substrate), position = position_jitterdodge(jitter.width = 0.25), shape = 21, size = 2, alpha = 0.8) + 
      labs(y = "", title = "Shannon\nDiversity", x = "") +
      theme( axis.text.x = element_text(angle = 90), legend.position = "none")
    
    ad3 <-  ggplot(mgdd_nona, aes(x = substrate, y  = invsimpson)) +
      geom_boxplot(aes(fill = substrate), outlier.shape = NA) +
      geom_point(aes(fill = substrate), position = position_jitterdodge(jitter.width = 0.25), shape = 21, size = 2, alpha = 0.8) + 
      labs(y = "", title = "Inverse\nSimpsons", x = "") +
      theme( axis.text.x = element_text(angle = 90), legend.position = "none")
    
    ad5 <-  ggplot(mgdd_nona, aes(x = substrate, y  = log10(countTotal))) +
      geom_boxplot(aes(fill = substrate), outlier.shape = NA) +
      geom_point(aes(fill = substrate), position = position_jitterdodge(jitter.width = 0.25), shape = 21, size = 2, alpha = 0.8) +
      labs(y = "", title = "log10 total\ncounts / sample", x = "") +
      theme( axis.text.x = element_text(angle = 90), legend.position = "none")
    
    
    ad6 <-  ggplot(mgdd_nona, aes(x = substrate, y  = diversity)) +
      geom_boxplot(aes(fill = substrate), outlier.shape = NA) +
      geom_point(aes(fill = substrate), position = position_jitterdodge(jitter.width = 0.25), shape = 21, size = 2, alpha = 0.8) +
      labs(y = "", title = "k-mer\ndiversity", x = "") +
      theme( axis.text.x = element_text(angle = 90), legend.position = "none")
    
    ad7 <-  ggplot(mgdd_nona, aes(x = substrate, y  = C)) +
      geom_boxplot(aes(fill = substrate), outlier.shape = NA) +
      geom_point(aes(fill = substrate), position = position_jitterdodge(jitter.width = 0.25), shape = 21, size = 2, alpha = 0.8) +
      labs(y = "", title = "k-mer\ncoverage", x = "") +
      theme( axis.text.x = element_text(angle = 90), legend.position = "none")
    
    
    # ad_plot <- ggpubr::ggarrange(ad1, ad2, ad3, ad5,  ad6, ad7, common.legend = TRUE, legend = "right", ncol = 3, nrow = 2 )
    library("patchwork")
    ad_plot <- ad1 + ad2 + ad3 + ad6 + plot_layout(guides = "collect", ncol = 4, nrow = 1 )
    ad_plot
    # ggsave(ad_plot, filename = "vis/mtu__paedcf_alpha-divx4.png", device = "png", width = 14, height = 6)
    # ggsave(ad_plot, filename = "vis/mtu__paedcf-Hsap_alpha-divx4.png", device = "png", width = 14, height = 6)
    # ggsave(ad_plot, filename = "vis/mtu__paedcf_alpha-divx6.png", device = "png", width = 12, height = 10)

      
## a g e   a t   s a m p l i n g   ## ====================================
    
    summary(mgdd$age.at.sampling)
    
    # needs a mean value, and boxplots at each
    
    divs.a <- c("richness", "shannon", "invsimpson", "diversity", "C")   #"faith_pd", , "coverage" 
    
    mgdd.l <- mgdd ; mgdd.l[ , "richness"] <- log10(mgdd.l[ , "richness"])
    head(mg <- reshape2::melt( mgdd[ , c(divs.a[ ], "age.at.sampling", "substrate")], id.vars = c("age.at.sampling", "substrate")))
    
    (ad1 <- ggplot(mg, aes(x = age.at.sampling, y  = value)) +
      facet_wrap(.~variable, scale = "free") +
      geom_point(aes(fill = substrate), position = position_jitterdodge(jitter.width = 0.1), shape = 21, size = 3.5, alpha = 0.6) +
      geom_smooth(method = lm, aes(colour = substrate, group = substrate), size = 1.5, alpha = 0.2) +
      scale_colour_manual(values = those_relic_cols) +
      scale_fill_manual(values = those_relic_cols) +
      labs(y = "diversity unit", x = "age at sampling") +
      theme( axis.text.x = element_text(angle = 90),
             legend.position = "bottom") )
    # ggsave( ad1, filename = "vis/mtu__paedcf__alpha_div_wrt_Age.png", width = 6, height = 3)

  
  # it does look like a duck, but shows no signal
    
    lapply( divs.a, function(aa){   # aa <- "richness"
      bb_df <- mgdd[ , c( "age.at.sampling", "substrate")]
      bb_df$dep_var <- mgdd[ , aa]
      # anova( lm(dep_var ~ age.at.sampling:substrate + age.at.sampling, bb_df ))
      summary( lm(dep_var ~ age.at.sampling*substrate + age.at.sampling, bb_df ))
      
      # emmeans::pwpm( emmeans::emmeans( lm(dep_var ~ age.at.sampling*substrate, bb_df ), specs = c( "substrate", "age.at.sampling")))
      # plot( emmeans::emmeans( lm(dep_var ~ age.at.sampling*substrate, bb_df ), specs = c( "substrate", "age.at.sampling")))
      
      })

    
  # consider subsetting and cor() --  A: still nothing
    
    sapply( divs.a, function(aa){   # aa <- "richness"   aa <- "invsimpson"
      
      sapply( c("BAL", "NS", "TS"), function(aaa){   # aaa <- "TS" aaa <- "NS"
        
        bb_df <- dplyr::filter( mgdd[ , c( "age.at.sampling", "substrate", "simpson")], substrate == aaa, simpson != 0 )
        bb_df$dep_var <- mgdd[ rownames(bb_df), aa]
        bb_df <- bb_df[ !apply(bb_df, 1, function(bbb) any(is.na(bbb))) , ]
        cc_t <- cor.test( bb_df[ , 1], bb_df[ , 3], method = "spearman" ) 
        # ifelse( cc_t$p.value < 0.1, cc_t$estimate, NA )
        cc_t$estimate
        })

      })
    
      ## spearman having removed the NAs, and 0-diversity samples
        # >        richness    shannon invsimpson  diversity          C
        # >  BAL  0.1604503  0.1604503  0.1604503  0.1604503  0.1604503
        # >  NS  -0.2097436 -0.2097436 -0.2097436 -0.2097436 -0.2097436
        # >  TS   0.1384467  0.1384467  0.1384467  0.1384467  0.1384467
      
    
  # consider subsetting and proper testing() --  A: still nothing still! 
    
    cf.model <- lm(dep_var ~ age.at.sampling * substrate, data = bb_df)
    summary(cf.model)
    anova(cf.model)

    library(emmeans)
    emtrends(cf.model, ~ substrate, var = "age.at.sampling")
    emtrends(cf.model, pairwise ~ substrate, var = "age.at.sampling")
    
    
## mgdd versus antibiotics   ---   a bit but not much
    
    unique(mgdd$Prohphylaxis.ABs)
    sum(mgdd$Prohphylaxis.ABs == "No", na.rm = TRUE)
    
    adiv_abx <- lapply( divs, function(aa){
      # bb_df <- mgdd_nona
      # bb_df <- dplyr::filter(mgdd_nona, substrate != "TS")
      bb_df <- dplyr::filter(mgdd_nona, !is.na(Prohphylaxis.ABs))
      bb_df$ifso <- bb_df[ , aa]
      bb_df$richness <- log10(bb_df$richness)
      ggplot(bb_df, aes(x = substrate, y  = ifso)) +
        facet_grid(.~Prohphylaxis.ABs == "No", space = "free") +
        geom_boxplot(aes(fill = substrate), outlier.shape = NA) +
        geom_point(aes(fill = substrate), position = position_jitterdodge(jitter.width = 0.25), shape = 21, size = 2, alpha = 0.8) + 
        # ggbeeswarm::geom_beeswarm(aes(fill = BMI.z.Score)) + 
        labs(y = "", title = aa, x = "") +
        theme( axis.text.x = element_text(angle = 90), legend.position = "none")
    })
    ggsave(
      adiv_abx[[1]] + adiv_abx[[2]] + adiv_abx[[3]] + adiv_abx[[4]] + adiv_abx[[5]] + adiv_abx[[6]] + plot_layout(ncol = 3),
      file = "vis/mtu__paedcf__adiv_abx.png",device = "png", 
      width = 15, height = 8)

    