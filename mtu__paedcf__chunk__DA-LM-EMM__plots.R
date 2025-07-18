
## mtu__som__DA__LMM-EMM-plots_feats_boxplot

    rm(list=ls())
    
    library(ggplot2)
    library(vegan)
    source("analysis/background_code/R__fns_jfg/fn_definitions.R")
    source("analysis/mtu__paedcf__spreadsheeting__shortcuts.R")

    ## add in pathogens here to avoid overwriting theme etc.
    get_mun <- function(){ 
      source("analysis/mtu__paedcf__chunk__microbiology.R")
      return( micro_uniq_named)
    }
    micro_uniq_named <- get_mun()
    
    ## test data  
    test_df <- readRDS("output/mtu__paedcf__DA_LMM-EMM__kA-0.001-0.1.RDS")

    
## ============================================================

    test_df$feat <- as.vector(test_df$feat)
    test_df$id <- as.vector(test_df$id)
    test_df$facet_by <- mgtax[ test_df$id , "p" ]
    test_df$plot_by <- mgtax[ test_df$id , "g" ]
    
    head(test_df)

    # test_df$stable <-  gsub("^ ", "", gsub(" $", "", gsub("\\.", "", test_df$contrast, perl = TRUE ), perl = TRUE ), perl = TRUE )
    test_df$stable <-  "unit"
    test_out <- test_df[ , c( "feat", "id", "plot_by", "facet_by", "contrast", "stable", "estimate", "t.ratio", "p.value", "FDR")]

        
    ## LOOK AT THE OUTPUT WITH YOUR EYES
    head(test_out)
    head(dplyr::filter(test_out, FDR < 0.05))
    # # View(test_out)
    # View(dplyr::filter(test_out, FDR < 0.05))
    
    # out of interest
    sum(test_out$FDR < 0.01)  # 241
    sum(test_out$FDR < 0.05)  # 270
    
    # get abundance values
    ms_ids <- unique(test_out$id)
    head(ms_idm <- reshape2::melt( as.matrix( mgfeat_clr[ ms_ids , ])))
    colnames(ms_idm) <- c("id_melt", "sample", "clr")
    
    
  ## complicate - three groups
      
      head(test_out)
      dim(test_out_0.05 <- dplyr::filter(test_out, FDR < 0.05))
    
      # View(test_out)
      test_out_0.05$contrast <- dplyr::recode(
        "BAL - NS" = "BAL w.r.t. MMS",
        "BAL - TS" = "BAL w.r.t. OPS",
        "NS - TS" = "MMS w.r.t. OPS",
        test_out_0.05$contrast
        )


##  P L O T S   =============================================================================================

      
  ##   theme update for common plots   ## ----------------------------------------------------------------------------

      # could probably assemble the entire plot less the data, and apply that way...
      
      theme_update(
        axis.text.x = element_text(angle = 0, size = 9),
        axis.text.y = element_text(face = "italic", size = 10),
        legend.position = "bottom",
        plot.title = element_text(face = "plain"),
        plot.subtitle = element_text(size = 8),
        panel.background = element_rect( colour = "grey15", fill = NA), #grey15
        panel.grid.major.y = element_line(colour = "grey70", size = 0.10),
        panel.spacing = unit(0.25, "cm"),
        strip.background = element_rect(fill = "grey90", colour = NA),
        strip.text.y = element_text(angle = 0, size = 14), #element_blank() #,
        strip.text.x = element_text(size = 14),
        text=element_text(size = 13)
      )

      ## boohoo advice boohoo
      # geom_point(aes(size = -log10(FDR)), shape = 21, alpha = 0.6, position = position_jitterdodge(jitter.width =0.1) ) +

               
  ##   all   ## ----------------------------------------------------------------------------
    
      # View(test_out_0.05[, c(1,5,7,8,10)])
      length(unique(test_out_0.05$id))
      length(unique(dplyr::filter(test_out_0.05, contrast!="BAL w.r.t. NS", estimate >0)$id))
      length(unique(dplyr::filter(test_out_0.05, contrast!="BAL w.r.t. NS", estimate >0)$plot_by))
      length(unique(dplyr::filter(test_out_0.05, contrast!="BAL w.r.t. NS", estimate <0)$id))

        those_effing_cols <- c("#FBA90A","#9A0794","#4D854C")
        names(those_effing_cols) <- c("BAL", "MMS", "OPS")
      
      (boxplot_BNT <- ggplot( test_out_0.05, aes(fill = facet_by, x = reorder(plot_by, estimate), y= estimate)) +    # , alpha = FDR < 0.05
        facet_grid(facet_by~contrast, scale = "free_y", space = "free_y") +
        theme_minimal() +
        coord_flip() +
        geom_hline(yintercept = 0, lty = 1, colour = "grey35", size = 0.25) +
        geom_boxplot(aes(), shape = 21, alpha = 0.6, outlier.shape = NA ) +
        # geom_point(shape = 21, aes(size = -log10(FDR), alpha = -log10(FDR)), position = position_jitterdodge(jitter.width =0.1) ) +
        geom_point(shape = 21, aes(size = -log10(FDR)), alpha=0.6, position = position_jitterdodge(jitter.width =0.1) ) +
        theme(
          legend.position = "bottom",
        ) +
        labs(y = "modelled abundance (CLR)", #log-ratio to mean sample abundance" ,  #"ratio of modelled abundances",
             x = "",
             title = "differentially abundant taxa between sample types",
             subtitle = "FDR threshold of 0.5%, prevalence threshold of 0.1% in 10% of samples (105 taxa)",
             NULL) +
        guides(
          # fill = guide_legend(nrow = 5, "______ at D13", override.aes = c( fill.size = 20, alpha = 1)),
          fill = "none"
          # size = guide_legend(nrow = 1, "FDR value (-log10)")
        ))
      
            
      ## no longer descends to such resolutions   (mtu__paedcf)       
        # # 1% thresh

              those_effing_cols <- c("#FBA90A","#9A0794","#4D854C")
              names(those_effing_cols) <- c("BAL", "MMS", "OPS")
              dim( ka_fae <- k_A( (mgfeat), k = 0.01, A = 0.1))
              test_out_0.05$which <- ifelse( test_out_0.05$estimate < 0, gsub("(\\w*) .*", "\\1", test_out_0.05$contrast), gsub(".* (\\w*)", "\\1", test_out_0.05$contrast))
              (boxplot_BNT2 <- ggplot(
                # test_out_0.05[ test_out_0.05$feat %in% rownames(ka_fae) , ],
                test_out_0.05,
                  aes(fill = which, x = reorder( gsub("g__", "", plot_by), estimate), y= estimate)) +    # , alpha = FDR < 0.05
                facet_grid( gsub("p__", "", facet_by)~contrast, scale = "free_y", space = "free_y") +
                theme_minimal() +
                coord_flip() +
                geom_hline(yintercept = 0, lty = 1, colour = "grey35", size = 0.25) +
                # geom_point(shape = 21, aes(size = -log10(FDR), alpha = -log10(FDR)), position = position_jitterdodge(jitter.width =0.1) ) +
                # geom_col(colour = "grey35", width = 0.001, fill = NA, position = "dodge") + # aes(color = plot_type)
                geom_boxplot(aes(), shape = 21, alpha = 0.6, outlier.shape = NA, key_glyph = "point" ) +
                geom_point(shape = 21, aes(size = -log10(FDR)), alpha=0.6, position = "dodge", key_glyph = "point" ) +
                # geom_point(shape = 21, aes(size = -log10(FDR)), alpha=0.6, position = position_jitterdodge(jitter.width =0.1), key_glyph = "point" ) +
                scale_fill_manual(values = those_effing_cols, "higher in: ") +
                theme( panel.background = element_rect( colour = "grey35"),
                  strip.text.x = element_text(angle = 0, size = 14, face = "bold"),
                  axis.text.y = element_text(angle = 0, size = 13, face = "italic"),
                  strip.text.y = element_text(angle = 0, size = 14),
                  legend.position = "bottom",
                ) +
                labs(y = "modelled abundance (CLR)", #log-ratio to mean sample abundance" ,  #"ratio of modelled abundances",
                     x = "",
                     # title = "differentially abundant taxa between sample types",
                     # subtitle = "FDR threshold of 0.5%, prevalence threshold of 0.1% in 10% of samples (105 taxa)",
                     NULL) +
                guides(
                  fill = guide_legend( override.aes = c( size = 7, alpha = 1, colour = NA, shape = 22)),
                  # fill = "none"
                  size = guide_legend( "FDR value (-log10)")
                ))
            
                
            
      
      
  ##   P A T H O G E N S   O N L Y   ## ----------------------------------------------------------------------------
      
      dim(test_out_0.05path <- dplyr::filter( test_out_0.05, feat %in% micro_uniq_named ))
      test_out_0.05path$plot_by <- mgtax[ test_out_0.05path$id , "s" ]
      
      length(unique(test_out_0.05path$id))
      
      (boxplot_path <- ggplot( test_out_0.05path, aes(fill = plot_by, x = reorder(plot_by, estimate), y= estimate)) +    # , alpha = FDR < 0.05
          facet_grid(facet_by~contrast, scale = "free_y", space = "free_y") +
          coord_flip() +
          geom_hline(yintercept = 0, lty = 1, colour = "grey35", size = 0.25) +
          # geom_boxplot(aes(), shape = 21, alpha = 0.6, outlier.shape = NA ) +
          geom_point(shape = 21, aes(size = -log10(FDR)), alpha = 0.6, position = position_jitterdodge(jitter.width =0.1) ) +
          theme(
          ) +
          labs(y = "modelled abundance (CLR)", #log-ratio to mean sample abundance" ,  #"ratio of modelled abundances",
               x = "",
               title = "differences in pathogenic taxa between sample types",
               subtitle = "FDR threshold of 0.5%, prevalence threshold of 0.1% in 10% of samples (8 taxa)",
               NULL) +
          guides(
            # fill = guide_legend(nrow = 5, "______ at D13", override.aes = c( fill.size = 20, alpha = 1)),
            fill = "none",
            size = guide_legend("significance of difference: -log10(FDR)")
            # size = guide_legend(nrow = 1, "FDR value (-log10)")
          ))
      
      
##   S A V E    O U T   ================================================================

      # ggsave(boxplot_BNT, filename = "vis/mtu__paedcf_boxplDA-BAL-NS-TS_0.001-0.1.png", device = "png", width = 8, height = 8)
      # ggsave(boxplot_path, filename = "vis/mtu__paedcf_boxplDA-BAL-NS-TS_0.001-0.1__path.png", device = "png", width = 8, height = 3.5)

      # write.table( cbind( mgtax_mgmt[ test_df$id, c(7, 1:6)], test_df[ , c(2,3,4,5,6,7,9,15)] ), file = "output/mtu__van__DA-LMM-EMM_total.tsv", sep = "\t", row.names = FALSE)


##   C H E C K    E M   ===============================================================

            # check_feat <- "asv__0012"  # Akkermansia
            # check_feat <- "asv__0738"  # Listeria
            check_feat <- "Actinomyces odontolyticus"
            check_feat <- "s__Homo_sapiens"
            check_feat <- "s__Rothia_aeria"
            View( cbind(
              mgdat[ colnames(mgfeat_ra) , "substrate"] ,
              mgfeat_ra[ check_feat , colnames(mgfeat_ra) ]
              ))
            # png(filename = paste0("vis/mtu__van__abund_plot__", check_feat, "__", mgtax[ check_feat, "g"], ".png"), width = 500, height = 300)
            plot(mgfeat_clr[ check_feat , ],
                 col = c(1:5)[factor( mgdat[ colnames(mgfeat_ra), "substrate"])],
                 # pch = c(1:5)[factor( mgdat[ rownames(mgfeat_ra), "______"])],
                 ylab = "% abundance",
                 xlab ="sample",
                 main = paste0("abundance of ", mgtax[ check_feat , "g" ], " - ", check_feat)
            )
            legend("topleft",                    # Add legend to plot
                   legend = c("BAL", "NS", "TS"),
                   col = c(1:3),
                   pch = c(19,19,19))
            # dev.off()


      # ##   B L A S T    E M   ===============================================================
      # 
      #       uknk_set <- dplyr::filter(test_out, facet_by == "Bacteroidia", feat == "unkn.", estimate  < -5)$id
      #       # all Muribaculaceae
      #       mgtax[ unique(uknk_set) , 1:7 ]
      # 
      #       unkn_df <- data.frame(mgtax[ unique(uknk_set) , 8 ])
      #       rownames(unkn_df) <- unique(uknk_set)
      #       # View(unkn_df)
      #       # write.table(unkn_df, file="output/mtu__van__DA_unkn_bacteroidia.tsv", sep = "\t" )
      # 
      #     ## BLAST!
      #       # - all outputs are most closely related to Muribaculaceae (<97%),
      #       # - except for asv__0033, which is Parabacteroides goldsteinii (cov 100%/ID 100%)


      