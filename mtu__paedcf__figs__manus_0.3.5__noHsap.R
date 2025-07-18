## mtu__paedcf__manus_fig


    ## v.0.3.5  ---   remove the human read sfrom kraken2 diversity estimates. 


      # fig1   # PPV/NPV versus culture, and versus mBAL
      # fig2   # alpha, beta, relative abundance
      # fig3   # DA, pathogen rel ab.
      #        # missing the red arrows ...
      # fig4   # culture NMDS + tresh curves


      rm(list=ls())
      set.seed(2205)
      
      
      library(vegan)
      library(ggplot2)
      library(patchwork)
      
      
      ab_col_pal <- ( c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a", "#FF1493",
                    "#b15928", "#737f33", "#8B008B", "#32fbd8", "#fdbf6f",
                    # RColorBrewer::brewer.pal(5, "Spectral")[1],
                    "#b2df8a", "#fb9a99", "#d9e627", "#EE82EE", "#DEB887",
                    "#a6cee3",
                    RColorBrewer::brewer.pal(4,'Accent')  ))

            
      source("../../grp__project/analysis/background_code/R__fns_jfg/fn_definitions.R")
      # does this in microbiology below  ~~~source("analysis/mtu__paedcf__spreadsheeting__shortcuts.R")
      source("analysis/mtu__paedcf__chunk__microbiology.R")
      
      substrate_substitute <- function(aa_df){ # aa_df <- mgdat
        bb_v <- colnames(aa_df)[ apply( aa_df, 2, function(aaa){ any(aaa %in% c("BAL", "NS", "TS"))  }) ]
        cc_df <- aa_df
        cc_df[ , bb_v] <- sapply( bb_v, function(aaaa){ # aaaa <- bb_v[1]
                                dplyr::recode( 
                                "NS" = "MMS",     # middle meatal swabs
                                "TS" = "OPS",    # oropharyngeal swabs
                                "BAL" = "BAL",   # i don't know much about the gold standard
                                cc_df[ , aaaa]) })
        cc_df
      }
      
      mgdat <- substrate_substitute(mgdat)
      naughtiers <- micro_uniq
      

### figure 1:   =======================================================================
  
    ##   PPV / NPV   -------------------------------------------------------------
      ## compare only to BAL culture results
      mgdat$micro_unif_BAL <- unlist(sapply( as.character(mgdat$participant), function(aaaa){
        bb <- dplyr::filter( mgdat, participant == aaaa, substrate == "BAL")
        ifelse( is.null(unlist(bb$micro_unif)), NA, unlist(bb$micro_unif) )
      }))
      
      # we know that this measure is different for each bug
      micro_uniq
      # and that the different methods perform differently
      data_stacks <- c( "micro.BAL", "micro.NS", "micro.TS")
      data_stacks2 <- c( "micro_unif_BAL", "micro_unif_NS", "micro_unif_TS")
      
      ## we also need to standardise TS and NS
      mgdat$micro_unif_NS <- unlist(sapply( as.character(mgdat$participant), function(aaaa){   # aaaa <- "SC114C5"
        bb <- dplyr::filter( mgdat, participant == aaaa, substrate == "MMS")
        ifelse( is.null(unlist(bb$micro_unif)), NA, unlist(bb$micro_unif) )
      }))
      
      mgdat$micro_unif_TS <- unlist(sapply( as.character(mgdat$participant), function(aaaa){
        bb <- dplyr::filter( mgdat, participant == aaaa, substrate == "OPS")
        ifelse( is.null(unlist(bb$micro_unif)), NA, unlist(bb$micro_unif) )
      }))
      
      n_events <- length(unique(mgdat$participant))
      
      cult_cult <- as.data.frame( do.call("rbind", 
                                          lapply( data_stacks2[-1], function(aa){   # aa <- "micro_unif_BAL"   aa <- "micro_unif_NS"   
                                            
                                            bb_df <- do.call("rbind", lapply( micro_uniq[-c(1,2,6)], function(aaa){   # aaa <- micro_uniq[8]
                                              
                                              
                                              # for all samples where pathogen WAS in reference:
                                              bbb_df <- unique(dplyr::filter( mgdat, grepl( aaa, micro_unif_BAL))[ , c( "participant", aa, "micro_unif_BAL")])
                                              true_pos <- sum( grepl( aaa, bbb_df[ , aa ]) )/n_events
                                              false_neg <- sum( !grepl( aaa, bbb_df[ , aa ]) )/n_events
                                              
                                              # for all samples where pathogen was NOT in reference:
                                              ccc_df <- unique(dplyr::filter( mgdat, !grepl( aaa, micro_unif_BAL))[ , c( "participant", aa, "micro_unif_BAL")])
                                              true_neg <- sum( !grepl( aaa, ccc_df[ , aa ]) )/ n_events
                                              false_pos <- sum( grepl( aaa, ccc_df[ , aa ]) )/ n_events
                                              
                                              # if( (true_neg + true_pos + false_neg + false_pos) == 1 ){
                                              c( 
                                                "true_pos" = true_pos,
                                                "true_neg" = true_neg,
                                                "false_pos" = false_pos,
                                                "false_neg" = false_neg,
                                                "ppv" = true_pos/(true_pos + false_pos + 0.000001),
                                                "npv" = true_neg/(true_neg + false_neg + 0.000001)
                                              )
                                              
                                              
                                            }))
                                            rownames(bb_df) <- paste0( gsub("micro_unif_", "", aa), "_", micro_uniq[-c(1,2,6)] )
                                            bb_df
                                            
                                          }) ), stringsAsFactors = FALSE)
      
      cult_cult[ , "micro"] <- gsub(".*_", "", rownames(cult_cult))
      cult_cult[ , "substrate"] <- gsub("_.*", "", rownames(cult_cult))
      cult_cult[ , "freq"] <- sapply( cult_cult$micro, function(aa){ length(grep(aa, mgdat$micro_unif_BAL)) })
      
      cult_cult$substrate <- gsub("NS", "MMS", gsub("TS", "OPS", cult_cult$substrate))
      
      
      (fig1a <- ggplot( cult_cult, aes( x = ppv, y = npv, fill = micro, size = freq )) +
          theme_minimal() +
          coord_fixed( ratio = 1, xlim = c(0,1), ylim = c(0,1)) + 
          # facet_wrap(substrate ~ . ,ncol = 1) +
          facet_wrap(substrate ~ ., nrow = 1) +
          scale_x_continuous( breaks = c(0,0.25,0.5,0.75,1,10) ) +
          scale_y_continuous( breaks = c(0,0.25,0.5,0.75,1,10) ) +
          scale_fill_manual(values = nr_col, NULL) +
          geom_abline(intercept = c(0,0), slope = 1, colour = "grey65", lty = 2) +
          geom_text_repel(data = dplyr::filter(cult_cult, ppv < 0.3), aes( label = gsub("(.).* (.*)", "\\1. \\2", micro, perl = TRUE) ),
                          direction = "both", max.time = 2, max.iter = Inf,
                          nudge_y = -0.5, nudge_x = 0.2,
                          size = 4, fontface = "italic", segment.color = "grey80", max.overlaps = 100) +
          geom_text_repel(data = dplyr::filter(cult_cult, ppv > 0.3), aes( label = gsub("(.).* (.*)", "\\1. \\2", micro, perl = TRUE) ),
                          direction = "both", max.time = 2, max.iter = Inf,
                          # nudge_y = -0.4, nudge_x = 0.1,
                          size = 4, fontface = "italic", segment.color = "grey80", max.overlaps = 100) +
          geom_point(shape = 22, alpha = 0.7) + # , position = position_jitter(width = 0.02, height = 0.02))
          labs(
            x = "positive predictive value", #\n(true culture positive / all culture positive in BAL)",
            y = "negative predictive value", #\n(true culture negative / all culture negative in BAL)",
            title="MMS and OPS culture as predictors of BAL culture:",
            # tag = "B",
            NULL) +
          theme(
            axis.text.x = element_text(size = 11),
            axis.text.y = element_text(size = 11),
            axis.title.x = element_text(size = 13, face = "plain"),
            axis.title.y = element_text(size = 13, face = "plain"),
            legend.direction = "vertical",
            legend.position = c(0.95, 0.22),
            legend.spacing = unit(units = "cm", 0.001),
            legend.title = element_text( size = 8 ),
            legend.text = element_text( size = 7 ),
            panel.border = element_rect(colour = "grey30", fill= NA, size = 0.6),
            panel.spacing.x = unit(1, "cm"),
            plot.title = element_text(size = 18, face = "plain", hjust = 0.5),
            strip.text.x = element_text(size = 14, face = "bold"),
          ) + 
          guides(
            # tag = "A",
            fill = "none",
            size = guide_legend( "n BAL\nsamples:", title.position = "top")
          ))
      
      
      
    ## comparing efficiency with BAL ----------------------------
      
      # ppv / npv
      
      sB_acc_dat <- substrate_substitute(sB_acc_dat)
      
      (fig1b <- ggplot( sB_acc_dat, aes( ppv, npv, fill = micro )) +
          coord_fixed(xlim = c(0,1.2), ylim = c(0,1.2)) +
          geom_abline(intercept = c(0,0), slope = 1, colour = "grey65", lty = 2) +
          facet_grid( . ~substrate ) +
          geom_point(shape = 21, aes(size = mean_ra), alpha = 0.7 ) +
          geom_text_repel(aes( label = gsub("(.).* (.*)", "\\1. \\2", micro, perl = TRUE) ),
                          direction = "both", max.time = 2, max.iter = Inf,
                          size = 4, fontface = "italic", segment.color = "grey80", max.overlaps = 100) +
          scale_x_continuous( breaks = c(0,0.25,0.5,0.75,1,10) ) +
          scale_y_continuous( breaks = c(0,0.25,0.5,0.75,1,10) ) +
          scale_fill_manual(values = ab_col_pal, NULL) +
          labs(
            x = "positive predictive value", #\n(true culture positive / all culture positive in BAL)",
            y = "negative predictive value", #\n(true culture negative / all culture negative in BAL)",
            title="\n\n\nBAL, MMS, and OPS culture as predictors of metagenomic BAL:", subtitle = NULL,
            NULL) +
          theme(
            axis.text.x = element_text(size = 8),
            axis.text.y = element_text(size = 8),
            axis.title.x = element_text(size = 12, face = "plain"),
            axis.title.y = element_text(size = 12, face = "plain"),
            # legend.position = "none",
            legend.position = c(0.95, 0.22),
            legend.title = element_text( size = 8 ),
            legend.text = element_text( size = 7 ), 
            panel.border = element_rect(colour = "grey30", fill= NA, size = 0.6),
            plot.title = element_text(size = 18, face = "plain", hjust = 0.5),
            # plot.title.position = element_text ),                                                   ###    can we shift the plots/titles closer together?
            strip.text.x = element_text(size = 14, face = "bold"),
          )+ 
          guides(
            # tag = "B",
            fill = "none",
            size = guide_legend( "average\nabundance (%):", title.position = "top")
          ))
      
      
      ## abandoned plot with some heatmaps in it...
      
      fig1 <- (fig1a / fig1b  ) & plot_annotation(tag_levels = "A") & theme( plot.tag = element_text(size =  12, face = "plain"))
      ggsave( fig1, filename = "vis/mtu__paedcf__figs_manus__fig1AB-cul-mgx-predict.png", device = "png", height = 11, width = 11, dpi = 200)
      
      write.table(sB_acc_dat, file = "output/mtu__SAINTS-CF__PPV-NPV__cult-to-BALmgx.tsv", sep = "\t")

      
### figure 2:   =======================================================================
      
    ## ---   alpha   --------------------------------------------------
      
      head(mgdd <- substrate_substitute(readRDS("output/mtu__paedcf__metadata_88-43_A-Div__noHsap.RDS")) )   #   mtu__paedcf__metadata_79-38_A-Div.RDS"))
      
      # can use models for this if residuals are norm, n'es pas?
      mgdd$BMI.z.Score <- as.numeric(mgdd$BMI.z.Score)
      
      colnames(mgdd)
      
    ## in this version, we re-estimate the microbial community absent the human reads. 
      
      divs <- c("richness", "shannon", "invsimpson", "chao1", "kappa", "C", "diversity")   #"faith_pd", 
      mdivs <- mgdd[ , c("substrate", "Genotype", "participant", "Gender", "age.at.sampling", "countTotal", "proph_bin", "modul_bin", "BMI.z.Score", divs) ]   # [1:4] 
      head(mdivs)
      
      mgdd <- dplyr::filter( mgdd, !is.na(age.at.sampling))
      
      ad1 <- ggplot(mgdd, aes(x = substrate, y  = chao1, fill = substrate, shape = substrate)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.55, size = 0.8) +
        geom_point(position = position_jitterdodge(jitter.width = 0.25), size = 2, stroke = 1.2, alpha = 0.8, aes(colour = substrate)) + 
        scale_shape_manual(values = c(4,15,17)) +
        scale_colour_manual(values = c("#FBA90A","#9A0794","#4D854C")) +
        scale_fill_manual(values = c("#FBA90A","#9A0794","#4D854C")) +
        labs(y = "estimated n species ", title = "species richness\n(Chao1 index)", x = "") +
        theme( axis.text.x = element_text(angle = 0),
               axis.title.y = element_text(size = 10),
               legend.position = "none")
      
      ad2 <- ggplot(mgdd, aes(x = substrate, y  = shannon, fill = substrate, shape = substrate)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.55, size = 0.8) +
        geom_point(position = position_jitterdodge(jitter.width = 0.25), size = 2, stroke = 1.2, alpha = 0.8, aes(colour = substrate)) + 
        scale_shape_manual(values = c(4,15,17)) +
        scale_colour_manual(values = c("#FBA90A","#9A0794","#4D854C")) +
        scale_fill_manual(values = c("#FBA90A","#9A0794","#4D854C")) +
        labs(y = "", title = "Shannon's H\ndiversity", x = "") 
      
      ad3 <- ggplot(mgdd, aes(x = substrate, y  = invsimpson, fill = substrate, shape = substrate)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.55, size = 0.8) +
        geom_point(position = position_jitterdodge(jitter.width = 0.25), size = 2, stroke = 1.2, alpha = 0.8, aes(colour = substrate)) + 
        scale_shape_manual(values = c(4,15,17)) +
        scale_colour_manual(values = c("#FBA90A","#9A0794","#4D854C")) +
        scale_fill_manual(values = c("#FBA90A","#9A0794","#4D854C")) +
        labs(y = "", title = "Inverse Simpson's\ndiversity", x = "") 
      
      ad4 <- ggplot(mgdd, aes(x = substrate, y  = kappa, fill = substrate, shape = substrate)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.55, size = 0.8) +
        geom_point(position = position_jitterdodge(jitter.width = 0.25), size = 2, stroke = 1.2, alpha = 0.8, aes(colour = substrate)) + 
        scale_shape_manual(values = c(4,15,17)) +
        scale_colour_manual(values = c("#FBA90A","#9A0794","#4D854C")) +
        scale_fill_manual(values = c("#FBA90A","#9A0794","#4D854C")) +
        labs(y = "", title = "sequence diversity\n(NonPareil k-mer index)", x = "")  # 
      
      ad5 <- ggplot(mgdd, aes(x = substrate, y  = C, fill = substrate, shape = substrate)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.55, size = 0.8) +
        geom_point(position = position_jitterdodge(jitter.width = 0.25), size = 2, stroke = 1.2, alpha = 0.8, aes(colour = substrate)) + 
        scale_shape_manual(values = c(4,15,17)) +
        scale_colour_manual(values = c("#FBA90A","#9A0794","#4D854C")) +
        scale_fill_manual(values = c("#FBA90A","#9A0794","#4D854C")) +
        labs(y = "", title = "kmer coverage\n(NonPareil completeness, C)", x = "") 
      
      ad6 <- ggplot(mgdd, aes(x = substrate, y  = diversity, fill = substrate, shape = substrate)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.55, size = 0.8) +
        geom_point(position = position_jitterdodge(jitter.width = 0.25), size = 2, stroke = 1.2, alpha = 0.8, aes(colour = substrate)) + 
        scale_shape_manual(values = c(4,15,17)) +
        scale_colour_manual(values = c("#FBA90A","#9A0794","#4D854C")) +
        scale_fill_manual(values = c("#FBA90A","#9A0794","#4D854C")) +
        labs(y = "", title = "sequence (kmer)\ndiversity", x = "")    # NonPareil sequence Kmer diversity, Nd
      
      (ad1 + ad2 + ad3 + ad4 + ad6) +    #  + ad5
        plot_layout(guides = "collect" ) +
        plot_annotation(title = "here we go again") +
        guide_area()
      

      ad1.a <- branch_droop(the_plot = ad1, branch_y = 110, droop_y = 110*0.8, droops_x = c(1.5,3),
                            droop_size = 0.6, b.lab_y = 120, b.lab_size = 3, b.lab_label = "FDR < 0.00001")
      ad1.b <- branch_droop(the_plot = ad1.a, branch_y = c(110*0.95), droops_x = c(1,2), droop_size = 0.6, b.lab_size = 0, droop_factor = 0.00, b.lab_label = "")


      ad2.a <- branch_droop(the_plot = ad2, branch_y = 4.5, droop_y = 4.5*0.8, droops_x = c(1.5,3),
                            droop_size = 0.6, b.lab_y = 5, b.lab_size = 3, b.lab_label = "FDR <0.00001 ", droop_factor = 0.025)
      ad2.b <- branch_droop(the_plot = ad2.a, branch_y = c(4.5*0.97), droops_x = c(1,2), droop_size = 0.6, b.lab_size = 0, droop_factor = 0.00, b.lab_label = "" )
      

      ad3.a <- branch_droop(the_plot = ad3, branch_y = 30, droop_y = 30*0.8, droops_x = c(1.5,3),
                            droop_size = 0.6, b.lab_y = 33, b.lab_size = 3, b.lab_label = "FDR < 0.00001 ", droop_factor = 0.05)
      ad3.b <- branch_droop(the_plot = ad3.a, branch_y = c(30*0.95), droops_x = c(1,2), droop_size = 0.6, b.lab_size = 0, droop_factor = 0.00, b.lab_label = "" )


      ad6.a <- branch_droop(the_plot = ad6, branch_y = 20, droop_y = 20*0.8, droops_x = c(1.5,3),
                            droop_size = 0.6, b.lab_y = 21, b.lab_size = 3, b.lab_label = "FDR < 0.0002 ", droop_factor = 0.05)
      ad6.b <- branch_droop(the_plot = ad6.a, branch_y = c(20*0.95), droops_x = c(1,2), droop_size = 0.6, b.lab_size = 0, droop_factor = 0.00, b.lab_label = "" )
      
            
      (ad_plot.annot <- ad1.b + ad2.b + ad3.b + ad6.b + plot_layout(guides = "collect", ncol = 4, nrow = 1 ) &
          theme(
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        legend.position = "none",
        NULL) )
      
      
     ###  ---   bars    -------------------------------------------------------------------- 
      
        colnames(mgtax) <- c("K", "P", "C", "O", "F", "G", "S", "ID")
        plot_by <- "G"
        plot_n <- 15    # 20 to remove betapolyomavirus
        facet_by <- "P"
        facet_n <- 4       # currently doing nothing
        colour_by <- "O"
        feats <- mgfeat_ra
        dim(feats)
        dim(feats)
        hier <- mgtax[ !apply(mgtax, 1, function(aa) any(is.na(aa))), ]
        moniker <- "other"
        na_term <- "unkn."
        

      ## crush to plottable ------------------------------------------
        
        top_Nplot <- names( sort(sapply( unique(hier[ , plot_by]), function(aa){  ## aa <- "g__Hydromonas"
          sum( feats[ rownames(hier)[grep(aa, hier[,plot_by]) ] , ])
        }), decreasing = TRUE)[ 1:plot_n])
        
      # fiddle-on-the-fly
        top_Nplot <- top_Nplot[-c(10)]
        
        top_Nfacet <- names( sort(sapply( unique(hier[ , facet_by]), function(aa){
          sum( feats[ rownames(hier)[grep(aa, hier[,facet_by]) ] , ])
        }), decreasing = TRUE)[ 1:facet_n])
        
        glomf <- t(do.call("rbind", lapply( top_Nplot, function(aa){  # aa <- "g__Rothia"
          if( sum( hier[ , plot_by] == aa) == 1 ){
            feats[ hier[,plot_by]==aa ,  ]
          }else{
            colSums(feats[ grep(aa, hier[,plot_by]) ,  ])
          }
        }))  )
        ## need fn here to glom by facet also, but see below too
        
        glomfo <- cbind( glomf, colSums(feats) - rowSums(glomf))
        colnames(glomfo) <- c( top_Nplot, moniker)
        dim(glomfo) ; glomfo[1:5, 1:5]
        
        
      ## plot  
        
        mg_m <- reshape2::melt( as.matrix(glomfo[,]))
        mg_m$Genotype <- mgdat[ mg_m$Var1 , "Genotype"]
        mg_m$substrate <- mgdat[ mg_m$Var1 , "substrate"]
        mg_m$Gender <- mgdat[ mg_m$Var1 , "Gender"]
        mg_m$sample <- mgdat[ mg_m$Var1 , "sample"]
        mg_m$participant <- mgdat[ mg_m$Var1 , "participant"]
        mg_m$hcluster <- factor(mgdat[ mg_m$Var1 , "hclust"], levels = mgdat$hclust)
        head(mg_m)
        mg_m$plot_var <- hier[ match( mg_m$Var2, hier[ , plot_by]) , plot_by]
        mg_m$facet_var <- hier[ match( mg_m$Var2, hier[ , plot_by]) , facet_by]
        mg_m[ is.na(mg_m[ , "plot_var"]) , "plot_var" ] <- moniker
        mg_m[ is.na(mg_m[ , "facet_var"]) , "facet_var" ] <- moniker
        #
        mg_m$plot_var <- gsub("g__", "", mg_m$plot_var)
        # mg_m$facet_var <- gsub("p__", "", mg_m$facet_var)
        mg_m$facet_var <- factor( gsub("p__", "", mg_m$facet_var), levels = c("Actinomycetota", "Bacillota", "Bacteroidota", "Campylobacterota", 
                                                                              "Fusobacteriota", "Pseudomonadota",
                                                                              "Cossaviricota", "Chordata", "other"))
      
      ## HOW FUNCTIONABLE
        mg_m$plot_var <- factor(mg_m$plot_var,
                                levels = gsub("g__", "", as.character( aggregate(value ~ Var2, mg_m, FUN = mean)[  order( aggregate(value ~ Var2, mg_m, FUN = mean)$value, decreasing = TRUE ), "Var2"] ))
        )
        
        ab_col_pal2 <- ab_col_pal
        names(ab_col_pal2) <- levels(mg_m$plot_var)

                
        # mg_m <- dplyr::filter( mg_m, !(facet_var %in% c("Chordata")) )
        (abund <- ggplot( mg_m, aes(fill = plot_var, x = sample, y = value*100)) +
            theme_minimal() + 
            facet_grid(facet_var ~ substrate, space = "free_x", scale = "free") +
            geom_col(color = "black") +
            # labs( x = "samples, ordered by similarity (Ward's D2 clustering on Jaccard dissim.)", ylab = "relative abundance") +
            labs( 
              # title = "", 
              y = "relative abundance (%)", 
              x = "", #  Kaiju Taxonomic output",
              # tag = "B",
              NULL ) +
            scale_fill_manual(values = ab_col_pal2, "15 most abundant genera") +
            theme(
              axis.text.x = element_blank(), #element_text(angle = 90),
              panel.border = element_rect(size = 0.5, linetype = 1, colour = "grey50", fill = NA),
              legend.text = element_text(face = "italic", size =10),
              legend.title = element_text(size =14, hjust = 0.5),
              # legend.position = c(0.55, 1.1), #"bottom",
              legend.position = "top",
              # axis.line = element_line(), 
              panel.grid.minor = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_line(size = 0.5, linetype = 3, colour = "grey70"),
              panel.spacing.y = unit(0.4, "cm"),
              plot.margin = unit(units = "cm", c(0.2, 0.2, 0.2, 1)),
              strip.text.x = element_text(size = 16, face = "bold"),
              strip.text.y = element_text(angle = 0, size=14, face = "italic", hjust = -0.05)
            ) + 
            guides( fill=guide_legend(nrow = 3, title.position = "top")))
        
      
      ##  beta :: CHIRes / Manus version  =============================
        
        # head(mgdat) ; dim(mgdat)
        
        mgdat <- dplyr::filter(mgdat, !grepl("LIBNEG", rownames(mgdat)))
        mgfeat <- mgfeat[ , rownames(mgdat) ]
        mgfeat_ra <- mgfeat_ra[ , rownames(mgdat) ]
        
      ## n counts in the NMDS are the number of species corresponding to that genus, not direct abund/incidence
        dim( mgfeat_raka <- k_A( mgfeat_ra, k = 0.0005, A = 0.05))
        ja_nmds_ra <- metaMDS( t(mgfeat_raka[ , colSums(mgfeat_raka) > 0]), distance = "jaccard", try = 20, trymax = 100)
      ## NOTE ::: 
        print("  + + +    < ! >    v 0.3.5  ---    death to all humans' reads :   nhsfeat / _ra   + + +")
        dim( nhsfeat_raka <- k_A( nhsfeat_ra, k = 0.0005, A = 0.05))
        nhs_ja_nmds_ra <- metaMDS( t(nhsfeat_raka[ , colSums(nhsfeat_raka) > 0]), distance = "jaccard", try = 20, trymax = 100)
        
        nhs_eu_nmds_clr <- metaMDS( t(nhsfeat_clr[ ,]), distance = "euclid", try = 20, trymax = 100)

        
      ## hacking of the ord_plot function - much anciliary
        
        # vegout <-  nhs_ja_nmds_ra
        vegout <-  nhs_eu_nmds_clr
        pdatf <- mgdat
        ptax <- apply( mgtax, 2, function(aa){ gsub( "\\w__", "", aa )})
        pfeat <- mgfeat
        
        samples_are_rows = FALSE
        
        ## testers
        alt_title <- NULL
        ax_lab <- "axis"
        alt_title <- "paed_cf mgfeat_ra Jac"
        OMIT = TRUE
        DO_VARS = FALSE
        DO_SPEC = TRUE
        tax_cent_thresh = 0
        taxon_level = "G"              # "g"           ## why did this change ???
        filt_asv = FALSE
        
        require(ggrepel)
        require(ggvegan)
        if(!samples_are_rows){ pfeat <- t(pfeat) }    
        ptax <- as.matrix( ptax )
        
        if( any ("matrix" %in% class(vegout) | grepl("MDS", class(vegout)) ) ){
          RD1 <- paste0( ax_lab, " 1")
          RD2 <- paste0( ax_lab, " 2")
          ord_df = data.frame(
        ##      
        ##      
        ## ---- < ! > ------   version conflict here
            # data.frame(pdatf[rownames(scores(vegout)$sites),], stringsAsFactors = FALSE),
            # "axis1" = scores(vegout)$sites[,1],
            # "axis2" = scores(vegout)$sites[,2])
        ## ---- < ! > ------   version conflict here
            
      ## v.0.3.5 euc-CLR
          data.frame(pdatf[rownames(scores(vegout)),], stringsAsFactors = FALSE),
          "axis1" = scores(vegout)[,1],
          "axis2" = scores(vegout)[,2])
        }
        ## ---- < ! > ------   version conflict here
        ##      
        ##      
        
      ## v.0.3.5 euc-CLR
        # asv_cent <- data.frame(vegout$species)
        asv_cent <- as.data.frame(wascores( scores(vegout), t(nhsfeat_ra)))
        
        colnames(asv_cent) <- c("asvAxis1", "asvAxis2")
        asv_cent$fill_var <- sapply( rownames(asv_cent), function(aa){
          if( aa %in% rownames(ptax) && !grepl("NA", ptax[ aa, taxon_level ])){
            ptax[ aa, taxon_level ] 
          }else{
            "undefined taxa"
          } })   # aa <- rownames(asv_cent)[21]
        asv_cent$size_var <- colMeans(pfeat[ , rownames(asv_cent)])
        ## better still to use conf ints, these v. small ad limited
        taxon_cent <- data.frame(
          aggregate( asvAxis1 ~ fill_var, FUN = mean, data = asv_cent)[ , ],
          "asvAxis2" = aggregate( asvAxis2 ~ fill_var, FUN = mean, data = asv_cent)[ , 2],
          "label_var" = apply(aggregate( . ~ fill_var, FUN = length, data = asv_cent)[ , 1:2], 1, function(aa){paste0(aa[[1]], " (sp=", aa[[2]], ")" )}),
          "size_var" = aggregate( size_var ~ fill_var, FUN = sum, data = asv_cent)[ , 2]
        )
        
        tax_count <- sapply( unique(asv_cent$fill_var), function(aa){ sum(asv_cent$fill_var == aa)}) 
        taxon_cent <- taxon_cent[ taxon_cent$fill_var %in% names( tax_count[ tax_count > tax_cent_thresh]) , ]
        tax_removed <- names( tax_count[ tax_count <= tax_cent_thresh])
        tax_removed_sub <- paste("omitted", length(tax_removed), "taxa from labels as incidence <", tax_cent_thresh,":\n", paste(abbreviate(tax_removed, minlength = 7,strict = TRUE, dot = TRUE), collapse = ", "))
        print(tax_removed_sub)
        
        if( filt_asv ){
          asv_cent <- dplyr::filter( asv_cent, !(fill_var %in% tax_removed))   # remove altogether
          print( paste0("   + + +    hiding all features with incidence below ", tax_cent_thresh))
        }else{
          print("   + + +    showing all features, even if not in centroids - will clutter up legend\n   + + +    consider \\'filt_asv = TRUE\\'")
        }
          
        
        if(is.null(alt_title)){plot_title <-  deparse(substitute(vegout)) }else{plot_title <- alt_title}
        plot_y_lim <- c( -0.65, 0.5)
        plot_x_lim <- c( -0.45, 1.75)
        
        
      ## --    B E T A   P L O T S   -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -   ##
        
        ## ---- A  
        
        ## i think take out undefined taxa - simply confuses the message, and largest to boot
        taxon_cent_no.unk <- dplyr::filter(taxon_cent, fill_var %in% levels(mg_m$plot_var) )
        taxon_cent_no.unk$label_var <- gsub( "\\w__", "", taxon_cent_no.unk$label_var )

        (bd_plot <- ggplot(ord_df, aes(x = axis1, y = axis2)) +  
            geom_hline(yintercept = 0, colour = "grey80")  +   # weight is the wrong arg
            geom_vline(xintercept = 0, colour = "grey80")  +
            
            ## add samples
            geom_point(aes(colour=substrate, shape = substrate), size = 4, stroke=1.2, alpha = 0.8) +    
            annotate("rect", alpha = .4, fill = "white", ) +
            scale_shape_manual(values = c(4,15,17)) +
            scale_colour_manual(values = c("#FBA90A","#9A0794","#4D854C")) +
            scale_fill_manual(values = ab_col_pal2 ) +
            ggrepel::geom_text_repel(data = taxon_cent_no.unk,
                                     aes(x = asvAxis1, y = asvAxis2, fill = NULL), 
                                     label = taxon_cent_no.unk$label_var,
                                     fontface = "italic",
                                     size = 4,
                                     force =0.9,
                                     force_pull =0,
                                     box.padding = 0.75,
                                     direction = "both",
                                     max.overlaps = Inf,
                                     min.segment.length = 0.5,
                                     segment.size = 0.2,
                                     hjust = 0.7,
                                     alpha = 0.85,
                                     max.time = 300, max.iter = 1e6) +
            # label.r = unit(0.2, "cm")) +
            geom_point(data = taxon_cent_no.unk, aes(x = asvAxis1, y = asvAxis2, fill = fill_var, size = size_var), alpha = 1, shape = 21, stroke = 0.55, colour = "black") +
            theme(
              axis.text.x=element_text(size=10),
              axis.text.y=element_text(size=10),
              axis.title.x=element_text(size=10),
              axis.title.y=element_text(size=10),
              legend.position = c(0.22, 0.95),
              legend.text=element_text(size=12),
              legend.title=element_text(size=12),
              panel.border = element_rect(colour = "grey30", fill= NA, size = 0.8),
              panel.grid.minor.x= element_blank(),
              panel.grid.minor.y= element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_blank(),
              # plot.margin = unit(c(1,-3,-10,0), "cm"),  # this doesnt work - spacing is overwritten
              plot.tag = element_text(size = 14, face = "plain", hjust = -2), 
              # plot.title = element_text(size = 14, face = "plain", hjust = 0.5), 
              plot.subtitle = element_text(size = 14, face = "plain", hjust = 0.5),
            ) +
            guides(
              fill = "none", 
              colour = guide_legend(ncol = 3,direction = "horizontal", "sample type", override.aes = list(size = 3, stroke = 3) ),
              shape = guide_legend(ncol = 3,direction = "horizontal", "sample type"), # , override.aes = list(size = 3, stroke = 3)
              size = "none"
            ) +
            labs(
              # tag = "A",
              # title = "metagenomic NMDS of matched BAL, MMS, and OPS samples\nshows strong similarity (overlap) between BAL and MMS",
              subtitle = "NMDS plot of metagenomic beta diversity (Jaccard dissimilaity)", #BAL and MMS are similar, while OPS is highly variable",  #balances better than the title
              x = "NMDS axis 1",
              y = "NMDS axis 2"
            )
        )
        
        
      ### fig 2 combo    --------------------------------------------------------------------
        
        # (ad_plot_quad <- plot_spacer() + ad1.b + ad2.c + ad3.b + ad4.c + plot_layout(guides = "collect", ncol = 1, heights = c(-0.2, 0.25, 0.25, 0.25, 0.25)) & 
        (ad_plot_quad <- plot_spacer() + ad1.b + ad2.b + ad3.b + ad6.b + plot_layout(guides = "collect", ncol = 1, heights = c(-0.2, 0.25, 0.25, 0.25, 0.25)) & 
           theme_void() +
           theme(
             axis.text = element_text(angle = 90, vjust = 1),
             axis.line.y.left = element_line(colour = "grey10", size = 0.5),
             axis.line.x.bottom = element_line(colour = "grey10", size = 0.5),
             plot.margin = unit(units = "cm", c(0, 0, 0.15, 0.15)),
             plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
             panel.grid.major.x = element_blank(),
             panel.grid.major.y = element_blank(),
             panel.grid.minor.x = element_blank(),
             panel.grid.minor.y = element_blank(),
             legend.position = "none") )
        
        fig2.a <- ad_plot_quad #& theme(plot.margin = unit(units = "cm", c(0.2, 0.2, -5, 0.2)))
        fig2.b <- ( abund / (bd_plot + plot_spacer() + plot_layout(widths = c(1, -0.15))) ) + plot_layout(heights = c(0.6, 0.4))
        fig2 <- (fig2.a | fig2.b) + plot_layout(widths = c(0.22,0.78), tag_level = "new") + plot_annotation(tag_level = c("A"))
        ggsave( fig2, filename = "vis/mtu__paedcf__figs_manus__fig2-NMDS-adiv-relab.png", device = "png", height = 12, width = 12, dpi = 200, bg = "white")

          
### figure 3:   ===========================================================================
      
    # 	- pathogen abundance
    # 	- DA in pathogens

        
    ## new abund fn ----  
        
        colnames(mgtax) <- c("D", "P", "C", "O", "F", "G", "S", "k_n", "ID")
        plot_by <- "G"
        facet_by <- "G"
        facet_n <- 4       # currently doing nothing
        colour_by <- "S"
        
        ## subset to pathogens only, grabbing Achro on the fly
        naughtiers %in% rownames(mgfeat)
        ## this ISNT what you want as it won't index the feat_ra
        naut_actual <- naughtiers[ unlist(lapply( naughtiers, function(aa){ any( grepl( gsub(" " , "_", aa), rownames(mgfeat)) ) })) ]    # aa <- naughtiers[4]
        ## instead, this is the matching naut terms
        naut <- unlist(lapply( naughtiers, function(aa){ grep( gsub(" " , "_", aa), rownames(mgfeat), value = TRUE)  }))    # aa <- naughtiers[4]
        names(naut) <- naut_actual
                        
        mg_p <- reshape2::melt(
          rbind(
            mgfeat_ra[ naut , ]
          ## these dont feature (Burk) or area already caught (Achro)
            # "Achromobacter misc." = colSums(mgfeat_ra[ grep("Achromo", rownames(mgfeat), value = TRUE) , ])
            # "Burkholderia misc." = colSums(mgfeat_ra[ grep("Burkholderia", rownames(mgfeat), value = TRUE) , ]) )     # not found
        ))
        mg_p$Var1 <- as.character(mg_p$Var1)
        mg_p$Var2 <- as.character(mg_p$Var2)
        mg_p$value <- as.numeric(mg_p$value)
        str(mg_p)
        
        
      ## assemble vars  
        
        str( mg_p <- cbind(mg_p, mgdat[ mg_p$Var2 , c( "Genotype", "substrate", "Gender", "sample", "participant", "hclust")]))
        mg_p$sample <- as.character(mg_p$sample)
        mg_p$participant <- as.character(mg_p$participant)
        mg_p$hclust <- factor(mgdat[ mg_p$Var2 , "hclust"], levels = mgdat$hclust)
        mg_p$plot_var <- gsub( " .*", "", mg_p$Var1)
        mg_p$facet_var <- factor( 
          gsub( " .*", "", mg_p$Var1),
          levels = c("Achromobacter misc.", sort(unique(gsub( " .*", "", mg_p$Var1))), "pathogen total")
        )

        ## again, in the revision this is not necessa s A.xylo alreadty caught...        
        mg_p[ is.na(mg_p[ , "plot_var"]) , "plot_var" ] <-   "Achromobacter misc." # moniker
        mg_p[ is.na(mg_p[ , "facet_var"]) , "facet_var" ] <- "Achromobacter misc."
        mg_p$value <- ifelse(mg_p$value == 0, NA, mg_p$value*100)
        str(mg_p)
        head(mg_p)
        
        ## highlight detected taxa  - could just have appended the micro_unif to mg_pb, but...
        ## change :: make the   -----------------------------------------------------------
        head(mg_p)
        mg_p$det <- as.numeric(apply( mg_p, 1, function(aa){  # aa <- mg_p[1,]      aa <- mg_p[141,]      aa <- mg_p[640,]      
          if( grepl( 
              # this hack to catch Achromo cases....
              gsub( "(Achromobacter) .*", "\\1", 
                gsub("s__(.*)_(.*)", "\\1 \\2", aa["Var1"], perl = TRUE)), mgdat[ unlist(aa["sample"]) , "micro_unif" ] ) ){
            # manage cases where dectected in culture but not in MGX
            ifelse( is.na(aa["value"]), 0.000000000000001, aa["value"])
            }else
            if( !grepl(aa["Var1"], mgdat[ unlist(aa["sample"]) , "micro_unif" ] ) ){NA}else  # not ideal
            {NA}
        }))
        
        
      ## -----------------------------------------------------------------------
        
        ## as factor here with accounting for additional groups      
        mg_p$plot_var <- factor( mg_p$plot_var, levels = c( sort(unique(mg_p$plot_var)), "pathogen total", "detected in culture"))
        mg_p$facet_var <- gsub( "Haemophi.*", "H. haemolyticus /\nH. influenzae", gsub("Strepto.*", "Streptococcus pneumoniae", mg_p$Var1))
        mg_p$facet_var <- factor( mg_p$facet_var, levels = c( sort(unique(mg_p$facet_var)), "pathogen total"))
        
        ## total pathogen burden?
        mg_pb <- mg_p
        mg_pb$facet_var <- "pathogen total"
        mg_pb$det <- NA
        mg_pb <- rbind(
          mg_p, 
          mg_pb
        )
        
        ## get a legend for the highlights
        mg_pb$plot_type <- mg_pb$plot_var
        mg_pb_leg <- mg_pb[ 1, ]
        mg_pb_leg$value <- 0.000000000001
        mg_pb_leg$plot_var <- "detected in culture"
        mg_pb_leg$plot_type <- "detected in culture"
        mg_pb <- rbind(
          mg_pb,
          mg_pb_leg
        )
        
        thresh <- 0
        
        
    ##   plot fig 3   ## =================
        
        ## centralise
        
        theme_update(
          plot.margin = unit(c(0,0,0,0), "cm"),  # if plotting class using super
          legend.position = "bottom",
          panel.border = element_rect(size = 0.5, linetype = 1, colour = "grey50", fill = NA),
        )  
        
        mg_pbs <- substrate_substitute(mg_pb)
        mg_pbs$Var1 <- gsub("s__(.*)_(.*)", "\\1 \\2", mg_pbs[,"Var1"], perl = TRUE)
        mg_pbs$facet_var <- factor( gsub("_"," ", gsub("s__", "", mg_pbs$facet_var)), levels = c(  "Achromobacter xylosoxidans", "Escherichia coli", "H. haemolyticus /\nH. influenzae", 
                                                                                            "Moraxella catarrhalis","Staphylococcus aureus", "Streptococcus pneumoniae", "pathogen total"))
        
        ##  in plot call, d.i.c. is an additional call via plot_type from above, with allowances made in scale_manual_xxx
        
        (path_abund_s_dic <- ggplot( mg_pbs, aes(fill = Var1, x = sample, y = value)) +
            theme_minimal() + 
            # facet_grid( facet_var ~ substrate, space = "free_x", scale = "free") +
            facet_grid( facet_var ~ substrate, space = "free_x", scale = "free") +
            geom_col(colour = "black", size = 0.25, alpha = 0.9) + # aes(color = plot_type)
            geom_text(aes(y = det, colour = is.na(value), label = ifelse( is.na(value), "", "*")), size = 12, na.rm = TRUE, alpha = 1) +
            geom_text(aes(y = det, colour = is.na(value), label = ifelse( is.na(value), "↓", "")), size = 8, na.rm = TRUE, alpha = 1, vjust = -0.25, face = "bold") +
            labs(
              x = NULL, 
              y = "relative abundance (%)", 
              NULL
            ) +
            scale_colour_manual(values = c("black", "red3")) + #c(ab_col[1:10], "yellow"), NULL) +
            scale_fill_manual(values = ab_col_pal, NULL) + #c(ab_col[1:10], "yellow"), NULL) +
            scale_y_continuous(minor_breaks = c(-10*100, thresh*100, 10*100), breaks = scales::breaks_extended(n = 3, only.loose = FALSE)) +
            theme(
              axis.text.x = element_blank(), #element_text(angle = 90),
              axis.text.y = element_text(size = 8), #element_text(angle = 90),
              axis.ticks.y = element_line(size = 0.5 ),
              axis.title = element_text(size = 12, face = "plain"),
              legend.position = "none",
              panel.border = element_rect(size = 0.5, linetype = 1, colour = "grey50", fill = NA),
              panel.grid.minor = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_blank(), #element_line(size = 0.5, linetype = 3, colour = "grey85"),
              panel.spacing.y = unit(0.25, "cm"),
              panel.spacing.x = unit(0.7, "cm"),
              #
              panel.ontop=TRUE, 
              panel.background = element_rect(fill = NA),
              #
              plot.margin = unit(c(0,0,0.2,0), "cm"),
              plot.title = element_text(size = 18, face = "plain", hjust = 0.5), 
              plot.subtitle = element_text(size = 18, face = "plain", hjust = 0.5), 
              strip.text.x = element_text(size = 16, face = "bold"),
              strip.text.y = element_text(angle = 0, size=14, face = "italic", hjust = -0.04)
            ) +
            guides( 
              fill=guide_legend(ncol = 10, nrow = 2, direction = "horizontal", title.position = "top", title.hjust = 0.5),
              colour=guide_legend(ncol = 10, nrow = 2, direction = "horizontal", title.position = "top", title.hjust = 0.5)  # , override.aes = list(shape = c(rep(22,10),21))
            ))
        
                
    ## test data      ---------------------------------------------------------
        
        test_df <- readRDS("output/mtu__paedcf__DA_LMM-EMM__kA-0.001-0.1.RDS")
        test_df$feat <- as.vector(test_df$feat)
        test_df$id <- as.vector(test_df$id)
        test_df$facet_by <- mgtax[ test_df$id , "P" ]
        test_df$plot_by <- mgtax[ test_df$id , "G" ]
        test_df$contrast <- gsub("TS", "OPS", 
                                 gsub("NS", "MMS",
                                      test_df$contrast))
        test_df$stable <-  "unit"
        test_out <- test_df[ , c( "feat", "id", "plot_by", "facet_by", "contrast", "stable", "estimate", "t.ratio", "p.value", "FDR")]
        
        # get abundance values
        ms_ids <- unique(test_out$id)
        head(ms_idm <- reshape2::melt( as.matrix( mgfeat_clr[ ms_ids , ])))
        colnames(ms_idm) <- c("id_melt", "sample", "clr")
        
        dim(test_out_0.05 <- dplyr::filter(test_out, FDR < 0.05))
        test_out_0.05$contrast <- dplyr::recode(
          "BAL - MMS" = "BAL w.r.t. MMS",
          "BAL - OPS" = "BAL w.r.t. OPS",
          "MMS - OPS" = "MMS w.r.t. OPS",
          test_out_0.05$contrast
        )
        # View(test_out)
        head(test_out)
        
        
    ##   P A T H O G E N   D A   P L O T   ## ----------------------------------------------------------------------------
 
        # dim(test_out_0.05path <- dplyr::filter( test_out_0.05, feat %in% micro_uniq_named ))
        dim(test_out_0.05path <- dplyr::filter( test_out_0.05, id %in% naut ))
        test_out_0.05path$plot_by <- mgtax[ test_out_0.05path$id , "S" ]
        
        length(unique(test_out_0.05path$id))
        
        test_out_0.05path$plot_by <- gsub("_"," ", gsub("s__", "", test_out_0.05path$plot_by))
        test_out_0.05path$facet_by <- gsub("_"," ", gsub("p__", "", test_out_0.05path$facet_by))

        # vector late to the game        
        test_out_0.05path$which <- ifelse( test_out_0.05path$estimate > 0, gsub("(\\w*) .*", "\\1", test_out_0.05path$contrast), gsub(".* (\\w*)", "\\1", test_out_0.05path$contrast))
        
        (lollypop_plot <- ggplot( test_out_0.05path, aes(fill = which, x = reorder(plot_by, estimate), y= estimate)) +    # , alpha = FDR < 0.05
            facet_grid(facet_by~contrast, scale = "free_y", space = "free_y") +
            coord_flip() +
            geom_hline(yintercept = 0, lty = 1, colour = "grey10", size = 0.25) +
            geom_point(shape = 21, aes(size = -log10(FDR)), alpha = 0.9) +
            geom_col(width = 0.045) +
            scale_fill_manual(values = those_effing_cols, "more abundant in:") + #c(ab_col[1:10], "yellow"), NULL) +
            theme_minimal() +
            theme(
              axis.text.x = element_text(size = 14),
              axis.text.y = element_text(face = "italic", size = 14),
              legend.position = "top",
              legend.direction = "horizontal",
              plot.title = element_text(size = 18, face = "plain", hjust = 0.5), 
              panel.background = element_rect( colour = "grey15", size = 0.30, fill = NA), #grey15
              panel.grid.major.y = element_line(colour = "grey40", size = 0.10),
              strip.text.y = element_text(angle = 0, size = 14, face = "italic", hjust = 0),
              strip.text.x = element_text(angle = 0, size=14, face = "bold"),
              NULL
            ) +
            labs(y = " difference in log-abundance to log-mean (CLR)", #"log-mean abundance (CLR)", #log-ratio to mean sample abundance" ,  #"ratio of modelled abundances",
                 x = "",
                 NULL) +
            guides(
              # fill = "none",
              size = guide_legend("significance: -log10(FDR)", title.position = "left", nrow = 1)
            ))
        
        
    ## make fig 3   =============
        
        fig3 <- lollypop_plot +         
            (plot_spacer() +
            path_abund_s_dic + 
            plot_layout(widths = c(-0.3,0.99)) ) +
          plot_layout(nrow = 2, heights = c(0.32, 0.68)) +
          plot_annotation(tag_levels = "A")
        # fig3
        ggsave( fig3, filename = "vis/mtu__paedcf__figs_manus__fig3-PathAb-PathDA.png", device = "png", height = 10, width = 12, dpi = 200)
        
        
### figure 4:   ============================================================================
    
      ## NMDS of Jaccardian culture  ----------------------------------------------- 
        
        mgdat <- substrate_substitute( mgdat )
        bug_dat <- substrate_substitute( bug_dat )
        
        # use bug dat and micro_uniq to create a matrix
        head(bug_dat)
        micro_uniq
        
        micro_uniq %in% bug_dat$micro
        unique(bug_dat$micro) %in% micro_uniq
        
      ## out of interest (for the ms) how many samples are each of the pathogens seen in across the study?
              #   sapply( micro_uniq[-c(1,2,6, 13)], function(aa){
              #     sum( mgfeat[ aa , ] > 0)
              #   })/86
              # ## %
              #   sapply( micro_uniq[-c(1,2,6, 13)], function(aa){
              #     sum( mgfeat_ra[ aa , ] > 0.01)
              #   })/86
        sapply( naut, function(aa){
          sum( mgfeat[ aa , ] > 0)
        })/86   ## mtu__paedcf
                # >  Staphylococcus aureus   Haemophilus influenzae Streptococcus pneumoniae Haemophilus haemolyticus         Escherichia coli    Moraxella catarrhalis  Streptococcus anginosus            Achromobacter 
                # >  0.01162791               0.27906977               0.32558140               0.18604651               0.03488372               0.26744186               0.02325581               0.01162791 
      ## %
        sapply( naut, function(aa){
          sum( mgfeat_ra[ aa , ] > 0.01)
        })/86   ## mtu__paedcf
                # >  Staphylococcus aureus   Haemophilus influenzae Streptococcus pneumoniae Haemophilus haemolyticus         Escherichia coli    Moraxella catarrhalis  Streptococcus anginosus            Achromobacter 
                # >  0.01162791               0.22093023               0.29069767               0.09302326               0.03488372               0.23255814               0.00000000               0.01162791 
        
        
        # removing NA only is interesting, but for comparability remove Norm, misc, NA
        path_mat <- sapply( unique(micro_uniq[-c(2)]), function(aa){   # aa <- micro_uniq[3]
          bb_df <- dplyr::filter(bug_dat, micro == aa)
          rownames(mgdat) %in% bb_df$sample
          # *1 to make it numeric!! genius.
        })*1
        rownames(path_mat) <- rownames(mgdat)
        path_mat <- path_mat[ rowSums(path_mat) != 0 , colSums(path_mat) != 0 ]
        
        path_mds_jac <- vegan::metaMDS( path_mat, distance = "jacc") # makes little diff : try = 1e4, trymax = 1e4, )
        p_m <- as.data.frame(vegan::scores( path_mds_jac )$sites)
        p_m$sample <- rownames(p_m)
        p_m$sample_proc <- gsub("SC(.*)C\\d-.*", "\\1", rownames(p_m))
        p_m$substrate <- mgdat[ p_m$sample , "substrate"]
        p_m$parti <- as.character(mgdat[ p_m$sample , "participant"])
        
        p_s <- as.data.frame(vegan::scores( path_mds_jac )$species)
        p_s$count <- colSums( path_mat[  , rownames(p_s)]) 
        rownames(p_s) <- gsub("Other Species Unidentified", "misc.unkn", rownames(p_s))
        p_s$nom <- gsub("^(\\w)\\w* ", "\\1. ", rownames(p_s))
        

        (mb_bd_plot <- ggplot(p_m, aes(x = NMDS1, y = NMDS2)) +
            geom_hline(yintercept = 0, colour = "grey80")  +   # weight is the wrong arg
            geom_vline(xintercept = 0, colour = "grey80")  +
            # coord_fixed(ratio = 1) + #, ylim = plot_y_lim, xlim = plot_x_lim  ) +
            
            ## add samples
            geom_jitter(aes(colour=substrate, shape = substrate), size = 4, stroke=3, alpha = 0.8, position = position_jitter(width = 0.05, height = 0.05)) +
            annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = .4, fill = "white", ) +
            scale_shape_manual(values = c(4,15,17)) +
            scale_colour_manual(values = c("#FBA90A","#9A0794","#4D854C")) +
            scale_fill_manual(values = c(ab_col_pal, "dodgerblue") ) +
            ggrepel::geom_text_repel(data = p_s,
                                     aes(x = NMDS1, y = NMDS2, fill = NULL),
                                     label = p_s$nom,
                                     fontface = "italic",
                                     size = 4.5,
                                     force =0.9,
                                     force_pull =0,
                                     box.padding = 0.75,
                                     direction = "both",
                                     max.overlaps = Inf,
                                     min.segment.length = 0,
                                     segment.size = 0.4,
                                     hjust = 0.7,
                                     alpha = 0.85,
                                     max.time = 300, max.iter = 1e6) +
            # label.r = unit(0.2, "cm")) +
            geom_point(data = p_s, aes(x = NMDS1, y = NMDS2, fill = nom, size = count), alpha = 1, shape = 21, stroke = 0.55, colour = "black") +
            theme(
              axis.text.x=element_text(size=10),
              axis.text.y=element_text(size=10),
              axis.title.x=element_text(size=10),
              axis.title.y=element_text(size=10),
              legend.position = c(0.08, 0.15),
              # legend.direction = "vertical",
              legend.text=element_text(size=12),
              legend.title=element_text(size=12),
              panel.border = element_rect(colour = "grey30", fill= NA, size = 0.8),
              panel.grid.minor.x= element_blank(),
              panel.grid.minor.y= element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_blank(),
              # plot.margin = unit(c(0,0,0,-2), "cm"),
              plot.tag = element_text(size = 14, face = "plain"),
              plot.title = element_text(size = 12, face = "plain", hjust = 0.5),
              plot.caption = element_text(size = 9, face = "italic", hjust = 0, colour = "grey20"),
            ) +
            guides(
              fill = "none",
              colour = guide_legend(nrow = 3,direction = "vertical", "sample type", title.position = "top", override.aes = list(size = 3, stroke = 3) ),
              shape = guide_legend(nrow = 3,direction = "vertical", "sample type", title.position = "top"), # , override.aes = list(size = 3, stroke = 3)
              size = "none"
            ) +
            labs(
              # tag = "B",
              # title = "NMDS of pathogens detected through clinical culture surveillance,\nshowing separation of possible opportunists (left) and likely commensals (right)",
              caption = "note: NMDS samples shifted slightly to reduce overplotting",
              x = "NMDS axis 1",
              y = "NMDS axis 2"
            )
          
        )

        # ggsave( mb_bd_plot, filename="vis/mtu__paedcf__micro_culture__jacc.png", device = "png", width = 8, height = 6)
        
        
  ##   A C C U R A C Y   R A N G E    =====================================================================
        
      micro_gu <- unique( gsub(" .*", "", micro_uniq) )
      thresh_vec <- seq(0.000, 0.05, 0.0005)
      
      ##  check range       
      # thresh_range <- parallel::mclapply( thresh_vec, function(aa_thresh){  # aa_thresh <- thresh_vec[34]
      thresh_range <- lapply( thresh_vec, function(aa_thresh){  # aa_thresh <- thresh_vec[34]
        
        thresh <- aa_thresh

        ##  G E N U S   x   S U B S T R A T E   level agglomerate   -------------------------------------

        g_acc_list <- lapply( c("BAL", "MMS", "OPS"), function(mm){                 # mm <- "BAL"

          nn_df <- dplyr::filter( mgdat, substrate == mm)
          oo_feat <- mgfeat_ra[ , rownames(nn_df)]

          pp_dat <- t( sapply( micro_gu[ -c(1,2,6)], function(aa){         # aa <- micro_gu[9]   aa <- micro_gu[7]   aa <- micro_gu[4]

            if( !any(grepl(aa, rownames(oo_feat))) ){
              c( "micro" = aa, "sens" = NA, "spec" = NA, "ppv" = NA, "npv" = NA, "mean_ra" = NA)
            }else{
              if( is.null(dim(mgfeat_ra[ grepl(aa, rownames(mgfeat_ra)), ])) ){
                bb_df <- mgfeat_ra[ grepl(aa, rownames(mgfeat_ra)), ]
              }else{
                bb_df <- colSums( mgfeat_ra[ grepl(aa, rownames(mgfeat_ra)), ] )
              }
              
              true_pos <- sum( sapply( rownames(nn_df), function(bbb){  #  bbb <- rownames(nn_df)[30]
                bb_df[ bbb ] > thresh & grepl( aa, nn_df[ bbb, "micro_unif" ])
              }))
              true_neg <- sum( sapply( rownames(nn_df), function(aaa){   # aaa <- rownames(nn_df)[30]
                bb_df[ aaa ] < thresh & !(grepl( aa, nn_df[ aaa, "micro_unif" ]))
              }))

              false_neg <- sum( sapply( rownames(nn_df), function(aaa){   # aaa <- rownames(nn_df)[30]
                bb_df[ aaa ] > thresh & !(grepl( aa, nn_df[ aaa, "micro_unif" ]))
              }))
              false_pos <- sum( sapply( rownames(nn_df), function(bbb){  #  bbb <- rownames(nn_df)[30]
                bb_df[ bbb ] < thresh & grepl( aa, nn_df[ bbb, "micro_unif" ])
              }))


              mean_ra <- sum(bb_df)

              c(
                "micro" = aa ,
                "sens" = (true_pos/ (true_pos + false_neg)) ,
                "spec" = (true_neg/ (true_neg + false_pos)) ,
                "ppv" = (true_pos/ (true_pos + false_pos)) ,
                "npv" = (true_neg/ (true_neg + false_neg)),
                "mean_ra" = mean_ra
                # (true_neg/ all_pres) ,
                # (true_pos/ all_pres) ,
                # (( true_neg + true_pos) / nrow(nn_df)))
              )
            }
          }) )

          cbind( "substrate" = mm, pp_dat)
        })

        g_acc_datish <-as.data.frame( do.call("rbind", g_acc_list), stringsAsFactors = FALSE )
        g_acc_dat <- data.frame( t( apply( g_acc_datish, 1, function(aa){ gsub("NaN", "0", aa) })), stringsAsFactors = FALSE )
        g_acc_dat[ , 2:6] <- apply( g_acc_dat[ , 2:6], 2, function(aa){ as.numeric(as.character(aa)) })
        # head(g_acc_dat)
        
        
    ##  G E N U S, No Substrate   level agglomerate   -------------------------------------
        
        g_acc_list__NoSub <- t( sapply( micro_gu[ -c(1,2,6)], function(aa){         #  aa <- micro_gu[9]   aa <- micro_gu[10]   aa <- micro_gu[7]
          
          if( !any(grepl(aa, rownames(mgfeat_ra))) ){
            c( "micro" = aa, "sens" = NA, "spec" = NA, "ppv" = NA, "npv" = NA, "mean_ra" = NA)
          }else{
            if( is.null(dim(mgfeat_ra[ grepl(aa, rownames(mgfeat_ra)), ])) ){
              bb_df <- mgfeat_ra[ grepl(aa, rownames(mgfeat_ra)), ]
            }else{
              bb_df <- colSums( mgfeat_ra[ grepl(aa, rownames(mgfeat_ra)), ] )
            }
            true_pos <- sum( sapply( rownames(mgdat), function(bbb){  #  bbb <- rownames(mgdat)[30] 
              bb_df[ bbb ] > thresh & grepl( aa, mgdat[ bbb, "micro_unif" ]) }))
            true_neg <- sum( sapply( rownames(mgdat), function(aaa){   # aaa <- rownames(mgdat)[30]
              bb_df[ aaa ] < thresh & !(grepl( aa, mgdat[ aaa, "micro_unif" ])) }))
            false_neg <- sum( sapply( rownames(mgdat), function(aaa){   # aaa <- rownames(mgdat)[30]
              bb_df[ aaa ] > thresh & !(grepl( aa, mgdat[ aaa, "micro_unif" ])) }))
            false_pos <- sum( sapply( rownames(mgdat), function(bbb){  #  bbb <- rownames(mgdat)[30] 
              bb_df[ bbb ] < thresh & grepl( aa, mgdat[ bbb, "micro_unif" ]) }))
            mean_ra <- sum(bb_df)
            c(
              "micro" = aa , 
              "sens" = (true_pos/ (true_pos + false_neg)) ,
              "spec" = (true_neg/ (true_neg + false_pos)) ,
              "ppv" = (true_pos/ (true_pos + false_pos)) ,
              "npv" = (true_neg/ (true_neg + false_neg)),
              "mean_ra" = mean_ra
            )
              }
            
          }) )
            
        
        g_acc_NoSub_datish <-as.data.frame( g_acc_list__NoSub, stringsAsFactors = FALSE )
        g_acc_NoSub_dat <- data.frame( t( apply( g_acc_NoSub_datish, 1, function(aa){ gsub("NaN", "0", aa) })), stringsAsFactors = FALSE )
        g_acc_NoSub_dat[ , 2:6] <- apply( g_acc_NoSub_dat[ , 2:6], 2, function(aa){ as.numeric(as.character(aa)) })
        
        
        list(
          cbind( aa_thresh, g_acc_dat),
          cbind( aa_thresh, g_acc_NoSub_dat)
        )
        
      # not any more we dont    }, mc.cores = 6, mc.preschedule = TRUE)
      })
      
      
      # View(thresh_range)
      
      tax_thresh_range <- do.call("rbind", lapply(micro_gu[-c(1,2,6)], function(aa){  # aa <- micro_gu[4]
        reshape2::melt(
          do.call( "rbind", lapply( thresh_range, function(aaa){  
            bbb <- aaa[[1]][ grep( aa, rownames(aaa[[1]] )) , ]
            bbb[ , "micro"] <- aa
            bbb
          }) )[, -8],   # aaa <- thresh_range[[25]]
          id.vars = c("substrate" , "micro", "aa_thresh") )
      }) )
      
      tax_thresh_range_NoSub <- do.call("rbind", lapply(micro_gu, function(aa){  # aa <- micro_gu[4]
        reshape2::melt(
          do.call( "rbind", lapply( thresh_range, function(aaa){  aaa[[2]][ grep( aa, rownames(aaa[[1]] )) , ] }) )[, -7],   # aaa <- thresh_range[[25]]
          id.vars = c("micro", "aa_thresh")
        ) }) )
      
      
      head(tax_thresh_range)
      str(tax_thresh_range)
      #
      head(tax_thresh_range_NoSub)
      str(tax_thresh_range_NoSub)
      tax_thresh_range_NoSub <- tax_thresh_range_NoSub[ !apply( tax_thresh_range_NoSub, 1, function(aa) any(is.na(aa)) ) , ]
      
      
      ## plot accuracy parameters  -  substrate  ----------------
      
      ttr_pv <- dplyr::filter(tax_thresh_range, grepl("pv", variable))
      ttr_pv$variable2 <- dplyr::recode( ttr_pv$variable, 
                                         "ppv" = "PPV\n(% correct pos)",
                                         "npv" = "NPV\n(% correct neg)",
      )
      
      ttr_pv$aa_thresh_100 <- ttr_pv$aa_thresh*100
      (acc_range <- ggplot(
            dplyr::filter(ttr_pv, micro != "Chryseobacterium", micro != "Pseudomonas", ), 
            aes(lty = substrate, colour = substrate, x = aa_thresh_100)) + 
          geom_line(aes(y = value), size = 0.8) +
          facet_grid(micro  ~ variable2) + 
          scale_y_continuous(breaks = c(0, 0.5, 1.0)) +
          scale_colour_manual(values = c("#FBA90A","#9A0794","#4D854C")) +
          labs(
            x = "% relative abundance cut-off",
            y = "predictive value (0-1)\n",
            NULL) +
        theme_void() +
        theme(
            axis.text.x = element_text(angle = 0, size = 9, colour = "grey20", hjust = 0),
            axis.text.y = element_text(angle = 0, size = 9, colour = "grey20", hjust = 0),
            axis.title.y = element_text(angle = 90, size = 14, colour = "grey20", hjust = 0.5),
            axis.title.x = element_text(size = 14, colour = "grey20", hjust = 0.5),
            axis.line.y = element_line(colour = "grey10", size = 0.5),
            axis.line.x = element_line(colour = "grey10", size = 0.5),
            plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
            panel.spacing.y = unit(0.5, "cm"),  
            panel.grid.major.x = element_line(size = 0.5, linetype = 3, colour = "grey85"),
            strip.text.x = element_text(face = "bold", hjust = 0.5, size = 12),
            strip.text.y = element_text(face = "italic", hjust = 0.0, size = 14),
            legend.position = "none") + 
        guides( colour = guide_legend(override.aes = list(size = 5.0)),
                  linetype = guide_legend(override.aes = list(size = 5.0)),
                  NULL)
      )
      
      # ggsave( acc_range, filename = "vis/mtu__som__figs_CHIRes__acc_range.png", device = "png", height = 8, width = 6, dpi = 200)
      
      
      (fig4 <- mb_bd_plot + acc_range + plot_annotation(tag_level = "A") + plot_layout(widths = c(0.7, 0.3)))
      ggsave( fig4, filename = "vis/mtu__paedcf__figs_manus__fig4-acc_range.png", device = "png", height = 6, width = 12, dpi = 200)
      
        
##   wild statements   - Supp Mat   ================================================================
      
      ##
      ## question here as to whether rates of detection are negatively correlated with the amount of host DNA (presumably so)
      ##
      ## could give each sample a total score to show how well they performed (e.g. 4/7, 0/1 = 0, 0/0 = NA) and then spearman with host abundance?
      ##
      ## that would give you an aggregated score, and rare taxa will inflate that? plot human v. pathogens and lm?
      ##
      ## More likely: logit regression between correctly identified (1,0) and CLR/RA human abundance?
      ##
      
      
      # more MGX pathogns than Cult pathogens
      sum(apply(micro_pc, 2, function(aa){ sum( aa > 0 )})[-c(1:2)])
      sum(aggregate(sample ~ micro,  data = bug_dat, FUN = function(aa){ (length(aa)) })[ -c(1:2), 2])
      
      # put a value on that
      bug_dat_sums <- aggregate( micro ~ sample, FUN = length, data = bug_dat)
      rownames(bug_dat_sums) <- bug_dat_sums[ , 1]
      # presume non param, but
      hist( apply( micro_pc[ rownames(bug_dat_sums) , ], 1, function(aa){ sum( aa > 0 )}) )
      hist( bug_dat_sums[ , 2] )
      wilcox.test( paired = TRUE, 
                    bug_dat_sums[  , 2],
                    apply( micro_pc[ rownames(bug_dat_sums) , ], 1, function(aa){ sum( aa > 0 )})
                    )
      
      ## correlation depth and n feats
      cor( log10(colSums(mgfeat)), apply( micro_pc[ , -c(1,2,6)], 1, function(aa){sum(aa>0)}) , method = "spearman" )
      plot( log10(colSums(mgfeat)), apply( micro_pc[ , -c(1,2,6)], 1, function(aa){sum(aa>0)})  )
      
      
  ## filtering effects
      # cat $WRK/2__filt/*log | grep 'Input Read Pairs' | sed -E "s/.*Pairs: ([0-9]*) Both Surviving: ([0-9]*) .*Forward Only Surviving: ([0-9]*) .* Reverse Only Surviving: ([0-9]*) .* Dropped: ([0-9]*) .*/\1,\2,\3,\4,\5/g" > $MAT/mtu__paedcf_bt2logcounts.csv
      # scp -i $key -o ProxyJump=$jgary $jdaed:/mnt/workspace2/jamie/mtu__paedcf/Materials/mtu__paedcf_bt2logcounts.csv outpout/
      head(trimmfilt <- read.csv( file = "output/mtu__paedcf_bt2logcounts.csv", header = FALSE))
      colnames(trimmfilt) <- c("raw", "f_surv", "r_surv", "both_surv", "dropped")
      summary(trimmfilt)
            # >        raw               f_surv             r_surv          both_surv         dropped       
            # >  Min.   :  520607   Min.   :  181518   Min.   :  14869   Min.   : 13469   Min.   : 162343  
            # >  1st Qu.: 3129192   1st Qu.: 1679139   1st Qu.: 114632   1st Qu.: 94983   1st Qu.: 474815  
            # >  Median : 4448924   Median : 2949450   Median : 192754   Median :148081   Median :1110450  
            # >  Mean   : 5176767   Mean   : 3317408   Mean   : 224856   Mean   :174004   Mean   :1460498  
            # >  3rd Qu.: 5996288   3rd Qu.: 4259176   3rd Qu.: 268396   3rd Qu.:210368   3rd Qu.:1661820  
            # >  Max.   :22237272   Max.   :14106997   Max.   :1181526   Max.   :816343   Max.   :7663006              
      ggplot( reshape2::melt(trimmfilt), aes(x = variable, y = log10(value))) + 
        coord_flip() + 
        theme_minimal() +
        geom_boxplot(shape = 21, outlier.shape = NA) + 
        ggbeeswarm::geom_beeswarm( shape = 21, fill = rep(c("yellow","black"), 220))
      
      
      head(bt2filt <- read.table( file = "output/mtu__paedcf_filt_multi__multiqc_general_export.tsv", header = TRUE, sep = "\t"))
      summary(bt2filt$Seqs)
            # >  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
            # >  0.1815  1.6791  2.9495  3.3174  4.2592 14.1070 

      ggplot( bt2filt, aes(x = "all", y = Seqs)) + 
        coord_flip() + 
        theme_minimal() +
        geom_boxplot(shape = 21, outlier.shape = NA) + 
        ggbeeswarm::geom_beeswarm( shape = 21, fill = sample( c("yellow", "black"), size = nrow(bt2filt), replace = TRUE)) + 
        labs( title = "bt2 outputs")
      
      
  ## inter versus intra distances (euc-CLR)   ----------------------------------
      
      clr_dist <- as.dist( vegdist( t(nhsfeat_clr[ , ]), method = "euclidean"))
      dist_mat <- as.matrix(clr_dist)
      dist_mat[ dist_mat == 0 ] <- NA
      
      parties <- unique( mgdat$participant )
      
      dmmelt <- reshape2::melt( dist_mat )
      dmmelt$isself <- apply( dmmelt, 1, function(aa){ ifelse( gsub("\\..*", "", aa["Var1"]) == gsub("\\..*", "", aa["Var2"]), "self", "other" ) 
      })
      dmmelt$compar <- apply( dmmelt, 1, function(aa){ gsub( "a|b", "", 
                                                             paste( sort( c(gsub(".*\\.", "", unlist(aa["Var1"])), gsub(".*\\.", "", unlist(aa["Var2"]))) ), collapse = "_") )
      })
      
      dmmelt_o <- dplyr::filter( dmmelt, !is.na(value), compar %in% c("BAL_NS", "BAL_TS", "NS_TS") )    # isself != "self", 
      
      lapply( unique(dmmelt_o$compar), function(aa){   # aa <- "BAL_NS"
        bb_df <- dplyr::filter( dmmelt_o, compar == aa )
        cc_t <- wilcox.test( dplyr::filter( bb_df, isself == "self")$value,
                             dplyr::filter( bb_df, isself == "other")$value
        )
        
        dd_p <- ggplot( bb_df, aes( x = isself, y = value, fill = isself)) +
          geom_violin( draw_quantiles = c(0.50) , ) + 
          # ggbeeswarm::geom_beeswarm(shape = 21) + 
          scale_fill_manual(values = c("grey40", "yellow"))
        
        list(cc_t, dd_p)
        
      })
      
      
      ## euclidean-CLR on nhsfeat   ----------------------------------------------
      
      ###  BAL_NS                                                                                          
      # >  statistic   58700                                                                                           
      # >  parameter   NULL                                                                                            
      # >  p.value     0.02967165                    * thats a single star!!                                                                  
      # >  null.value  0                                                                                               
      # >  alternative "two.sided"                                                                                     
      # >  method      "Wilcoxon rank sum test with continuity correction"                                             
      # >  data.name   "dplyr::filter(bb_df, isself == "self")$value and dplyr::filter(bb_df, isself == "other")$value"
      ###  BAL_TS                                                                                          
      # >  statistic   29024                                                                                           
      # >  parameter   NULL                                                                                            
      # >  p.value     0.7110175                                                                                       
      # >  null.value  0                                                                                               
      # >  alternative "two.sided"                                                                                     
      # >  method      "Wilcoxon rank sum test with continuity correction"                                             
      # >  data.name   "dplyr::filter(bb_df, isself == "self")$value and dplyr::filter(bb_df, isself == "other")$value"
      ###  NS_TS                                                                                           
      # >  statistic   34084                                                                                           
      # >  parameter   NULL                                                                                            
      # >  p.value     0.799277                                                                                        
      # >  null.value  0                                                                                               
      # >  alternative "two.sided"                                                                                     
      # >  method      "Wilcoxon rank sum test with continuity correction"                                             
      # >  data.name   "dplyr::filter(bb_df, isself == "self")$value and dplyr::filter(bb_df, isself == "other")$value"
      
            
            ## jaccard on mgfeat_ra_ka   -----------------------------------------------
            
            # ja_nmds_ra_dist <- as.dist( vegdist( t(mgfeat_raka[ , colSums(mgfeat_raka) > 0]), method = "jaccard"))
            # dist_mat <- as.matrix(ja_nmds_ra_dist)
            
            ###  BAL_NS                                                                                          
            # >  statistic   8522.5                                                                                          
            # >  parameter   NULL                                                                                            
            # >  p.value     0.2751722                                                                                       
            # >  null.value  0                                                                                               
            # >  alternative "two.sided"                                                                                     
            # >  method      "Wilcoxon rank sum test with continuity correction"                                             
            # >  data.name   "dplyr::filter(bb_df, isself == "self")$value and dplyr::filter(bb_df, isself == "other")$value"
            ###  BAL_TS                                                                                          
            # >  statistic   6661.5                                                                                          
            # >  parameter   NULL                                                                                            
            # >  p.value     0.370851                                                                                        
            # >  null.value  0                                                                                               
            # >  alternative "two.sided"                                                                                     
            # >  method      "Wilcoxon rank sum test with continuity correction"                                             
            # >  data.name   "dplyr::filter(bb_df, isself == "self")$value and dplyr::filter(bb_df, isself == "other")$value"
            ###  NS_TS                                                                                           
            # >  statistic   8617                                                                                            
            # >  parameter   NULL                                                                                            
            # >  p.value     0.784462                                                                                        
            # >  null.value  0                                                                                               
            # >  alternative "two.sided"                                                                                     
            # >  method      "Wilcoxon rank sum test with continuity correction"                                             
            # >  data.name   "dplyr::filter(bb_df, isself == "self")$value and dplyr::filter(bb_df, isself == "other")$value"
            
            