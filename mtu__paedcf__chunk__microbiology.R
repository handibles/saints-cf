## category count off

      # rm(list=ls())
      source("analysis/mtu__paedcf__spreadsheeting__shortcuts.R")
  
      library(ggplot2)
      library(ggrepel)
      
      ## i mean, why not
      # load("output/mtu__paedcf__microbiology_data.RData")
      
          
    ## unify, i.e. spellcheck  -----------------------
      
      # where did $micro come from??? it came from kaiju_spreadsheeting
      
      mgdat$micro <- gsub( ", ", ",", mgdat$micro)
      mgdat$micro_unif <- sapply( mgdat$micro, function(aa){
          gsub("Achromobacter Species", "Achromobacter",
          gsub("Chryseobacterium Indologenese", "Chryseobacterium indologenes",
          gsub("E Coli", "Escherichia coli",
          gsub("Ecoli", "Escherichia coli",
          gsub("Haemophillus Haemolytical", "Haemophilus haemolyticus",
          gsub("Haemophillus Haemolyticus", "Haemophilus haemolyticus",
          gsub("Haemophilus influenzaee", "Haemophilus influenzae",
          gsub("Haemophillius Influenza", "Haemophilus influenzae",
          gsub("Haemophillus Infleunzae", "Haemophilus influenzae",
          gsub("Haemophillus Influenzae", "Haemophilus influenzae",
          gsub("Haemophillus Influenza", "Haemophilus influenzae",
          gsub("Haemophilus Influenzae", "Haemophilus influenzae",
          gsub("Haemphillus Influenza", "Haemophilus influenzae",
          gsub("Moroxeela Cattarhalis", "Moraxella catarrhalis",
          gsub("Moroxella Catarrhalis", "Moraxella catarrhalis",
          gsub("MRSA", "Staphylococcus aureus",
          gsub("NA", "no data", perl = TRUE,
          gsub("Normal", "normal/misc", perl = TRUE,
          gsub("Normal Flora", "normal/misc", perl = TRUE,
          gsub("Other Organisms unspecified", "Other Species Unidentified", perl = TRUE,
          gsub("Other Species$", "Other Species Unidentified", perl = TRUE,
          gsub("Pseudomonas Aerug.", "Pseudomonas aeruginosa",
          gsub("Staph Aureus", "Staphylococcus aureus",
          gsub("Step Anginosis", "Streptococcus anginosus",
          gsub("Strep Pnemoniae", "Streptococcus pneumoniae",
          gsub("Strep Pneumoniae", "Streptococcus pneumoniae",
          gsub("Strep.Pneumoniae", "Streptococcus pneumoniae",  aa)))))))))))))))))))))))))))
      })

      
      
  ## enumerate  ---------------------------
      
    # atomise      
      micro_uniq <- unique( unlist(sapply( mgdat$micro_unif, function(aa){ strsplit( aa, ",") })) )
      micro_uniq[ micro_uniq %in% mgtax[ , 6] ]
      micro_uniq[ !(micro_uniq %in% mgtax[ , 6]) ]
      # exclude classes with no rational representation in NGS
      micro_uniq_named <- micro_uniq[ -c(1,2,6)]

    # abundance in sampling
      rownames(mgtax) == rownames(mgfeat)
      micro_pc <- sapply( micro_uniq, function(aa){    # aa <- micro_uniq[3]
        bb_df <- mgfeat_ra[ grep( gsub( "(.*) .*", "\\1", aa), mgtax[ rownames(mgfeat_ra) , "g"] ) , ]
        if( is.null( dim(bb_df))){ bb_df }else{ colSums(bb_df)}
      })
      mp_m <- reshape2::melt(as.matrix(micro_pc) )
      mp_m$value <- as.numeric( mp_m$value)
      mp_m$substrate <- mgdat[ mp_m$Var1, "substrate" ]

      
    ## Do standard NMDS plot of the relevant pathogens
      micro_pc0 <- micro_pc[ , ( !colSums(micro_pc) == 0 & colnames(micro_pc) != "Achromobacter") ]
      # micro_MJ <- vegan::metaMDS( comm = micro_pc0, distance = "jaccard", try = 200)
      # plot3b <- ord_plot( micro_MJ, alt_title = "NMDS of surveillance samples")
      
      colnames(mgtax) <- c("k", "p", "c", "o", "f", "g", "s", "ID", "kid")
      # plot3a <- ord_plot( micro_MJ,
      #           alt_title = "NMDS of surveyed pathogen genera",
      #           pfeat = micro_pc0,
      #           taxon_level = "g",
      #           DO_SPEC = TRUE,
      #           samples_are_rows = TRUE,
      #           filt_asv = FALSE, 
      #           tax_cent_thresh = 0)
      # 
      # plot3 <- ggpubr::ggarrange( plot3b, plot3a, ncol = 2, nrow = 1, common.legend = FALSE)      
      # ggsave( plot3, filename="vis/mtu__paedcf__micro_culture_NMDS-JA.png", device = "png", width = 18, height = 6, bg = "white")
      
      
    ## 
                  
    # incidence in sampling  
      micro_binary <- sapply( micro_uniq, function(aa){
        bb_df <- mgfeat_ra[ grep( gsub( "(.*) .*", "\\1", aa), mgtax[ rownames(mgfeat_ra) , "g"] ) , ]
        if( is.null( dim(bb_df))){ bb_df > 0 }else{ colSums(bb_df) > 0 }
      })
      mp_bin_m <- reshape2::melt(as.matrix(micro_binary) )
      mp_bin_m$value <- as.numeric( mp_bin_m$value)
      mp_bin_m$substrate <- mgdat[ mp_bin_m$Var1, "substrate" ]
      

    # n times in metadata    
      metad_melt <- do.call("rbind",
            apply( mgdat, 1, function(aa){ # aa <- mgdat[30,]
                    bb_df <- as.data.frame( strsplit( aa["micro_unif"], ","))
                    bb_df[,"substrate"] <- aa[ "substrate"]
                    colnames(bb_df) <- c("micro", "substrate")
                    bb_df
                  })
        )
      # metad_melt$micro <- ifelse( grepl("NA", metad_melt$micro), NA, metad_melt$micro)
      # metad_melt <- dplyr::filter( metad_melt, !( micro %in% micro_uniq[ c(8)] ))
      
      
    # sorting vector  -  explicitly use as sorted
      levels( mp_m$Var2 )
      levels(mp_bin_m$Var2)
      metad_melt$micro <- factor( metad_melt$micro, levels = levels(mp_bin_m$Var2) )

      
  # does cult match metag? contingency table 
      
      mp_m$decis <- sapply( as.character(mp_m$Var1), function(aa){   # aa <-  as.character(mp_m$Var1)[39]
        if( grepl("^no data$", mgdat[ aa, "micro_unif"], perl = TRUE)){"no data"}else if(grepl("^normal/misc$", mgdat[ aa, "micro_unif"], perl = TRUE)){"normal/misc"}else{"pot. pathogen"}
        })
      mp_bin_m$decis <- sapply( as.character(mp_bin_m$Var1), function(aa){
        if( grepl("^no data$", mgdat[ aa, "micro_unif"], perl = TRUE)){"no data"}else if(grepl("^normal/misc$", mgdat[ aa, "micro_unif"], perl = TRUE)){"normal/misc"}else{"pot. pathogen"} 
        })
      metad_melt$decis <- sapply( gsub("\\.\\d$", "", rownames(metad_melt), perl = TRUE), function(aa){
        if( grepl("^no data$", mgdat[ aa, "micro_unif"], perl = TRUE)){"no data"}else if(grepl("^normal/misc$", mgdat[ aa, "micro_unif"], perl = TRUE)){"normal/misc"}else{"pot. pathogen"} 
        })
        
      # mp_m$decis <- sapply( as.character(mp_m$Var1), function(aa){ if( grepl("no data", mgdat[ aa, "micro"])){"no data"}else if(grepl("normal/misc", aa)){"normal/misc"}else{"pot. pathogen"} })
      # mp_bin_m$decis <- sapply( as.character(mp_bin_m$Var2), function(aa){ if( grepl("no data", aa)){"no data"}else if(grepl("normal/misc", aa)){"normal/misc"}else{"pot. pathogen"} })
      # metad_melt$decis <- sapply( as.character(metad_melt$micro), function(aa){ if( grepl("no data", aa)){"no data"}else if(grepl("normal/misc", aa)){"normal/misc"}else{"pot. pathogen"} })
      
      
      bug_dat <- as.data.frame( do.call("rbind", 
              apply( mgdat, 1, function(aa){                          # aa <- mgdat[20,]
                bb_samp <- as.character(aa[ "sample"][[1]])
                cc_micro <- unlist( strsplit( aa[ "micro_unif"][[1]], ","))
                do.call("rbind", lapply( cc_micro, function(aaa){                      # aaa <- "Escherichia coli"
                    if( all(!( grepl( gsub( " ", "_", aaa), mgtax[ rownames(mgfeat_ra) , "s"] ) )) ){
                      c( aaa, bb_samp, NA )
                    }else{  c( aaa, bb_samp, sum(mgfeat_ra[ grep( gsub(" ","_",aaa), mgtax[ rownames(mgfeat_ra) , "s"]) , bb_samp] )) }
                }) )
              }) ) , stringsAsFactors = FALSE)
      colnames(bug_dat) <- c( "micro", "sample", "percent")
      bug_dat$percent <- as.numeric( bug_dat$percent)
      bug_dat$incidence <- ifelse( bug_dat$percent > 0, 1, 0)
      bug_dat$substrate <- mgdat[ bug_dat$sample, "substrate"]
      bug_dat$micro <- factor( bug_dat$micro, levels = levels(mp_bin_m$Var2) )
            
            
  ##   #############################################################  ##      
  ##   ################   november 2023   ##########################  ##      
  ##   #############################################################  ##      
      
       microbiology_incidence <- t(
        sapply( 
          unique( micro_uniq), function(aa){ 
            sapply( unique( mgdat$substrate), function(aaa){
              sum( grepl(aa, dplyr::filter( mgdat, substrate == aaa)$micro_unif)) 
                }) 
              })
      )
      # write.table(microbiology_incidence, file = "output/mtu__paedcf__micro_uniq_incdence__Nov2023-Mar2025.tsv", sep = "\t")

      
      # there IS micro data for each substrate!      
      microbiology_binary <- t(
        sapply( 
          unique( micro_uniq), function(aa){ 
            
            sapply( rownames( mgdat ), function(aaa){
              grepl(aa, mgdat[ aaa, "micro_unif"]) 
                }) 
              })
      )

       library(ComplexHeatmap)
       library(circlize)
       library(RColorBrewer)
      
      proph_v <- sapply( mgdat$Prohphylaxis.ABs, function(aa){
        if( grepl( "No|NA", aa)){
          "none"
        }else if( grepl( "zithro", aa)){
          "azithromycin"
        }else if( grepl( "lucloxaci|loxpen", aa, perl = TRUE)){
          "flucloxacillin"
        }else{"none"}
      })
      modul_v <- sapply( mgdat$Modulator...Start.Date, function(aa){
        if( grepl( "alydeco", aa)){
          "kalydeco"
        }else if( grepl( "rkambi", aa)){
          "orkambi"
        }else{"none"}
      })
      
       sevr <- apply( mgdat, 1, function(aa){ sum( micro_uniq %in% unlist(strsplit( gsub("no data", "", aa[ "micro_unif"]), split=",")))})/ (mgdat$age.at.sampling*mgdat$BMI.z.Score)
                                                                                                                       
       # gender_names <- RColorBrewer::brewer.pal(2, "Dark2")[1:2] # terhor[1:2]
       gender_names <- c("#9CB4D3", "#86DEB7")
       names(gender_names) <- c("Female", "Male")
       P.aero_names <- c(RColorBrewer::brewer.pal(2, "Set1")[1], "white", "grey90") # terhor[3:5]
       names(P.aero_names) <- c("Positive", "Free", "NA")
       modul_names <- c(RColorBrewer::brewer.pal(2, "Set2")[1:2], "white") # terhor[6:8]
       names(modul_names) <- c("kalydeco", "orkambi", "none")
       proph_names <- c(RColorBrewer::brewer.pal(4, "Set3")[c(1,4)], "white") # terhor[9:11]
       names(proph_names) <- c("azithromycin", "flucloxacillin", "none")
       location_names <- c("#FBA90A","#9A0794","#4D854C") # terhor[12:14]
       names(location_names) <- c("BAL", "MMS", "OPS")

     # ramp off
       bmi_col_fun = colorRamp2(c(-1, 0, 2), c( "white", brewer.pal(2, "YlGn")[2:3]))
       age_col_fun = colorRamp2(c(0, 2, 6), c( "white", brewer.pal(2, "RdPu")[2:3]))
       load_col_fun = colorRamp2(c(0, 2, 4, 8), c( "white", brewer.pal(3, "Blues")))

       ## why are the colours so wrong?
       
       # https://jokergoo.github.io/ComplexHeatmap-reference/book/legends.html#the-side-of-legends
                       
       column_ha <- HeatmapAnnotation(
         na_col = "grey90",
         gender = mgdat$Gender,
         age = mgdat$age.at.sampling,
         BMI = mgdat$BMI.z.Score,
         modul. = modul_v,
         proph. = proph_v,
         P.aero = mgdat$Pseudomonas.Status,
         location = gsub("NS", "MMS", gsub("TS", "OPS", mgdat$substrate)),
         `n cultured` = apply( mgdat, 1, function(aa){ sum( micro_uniq %in% unlist(strsplit( gsub("no data", "", aa[ "micro_unif"]), split=",")))}),
         
         annotation_legend_param = list(direction = "horizontal",
                                        # foo1 = list(direction = "horizontal"),
                                        # bar1 = list(nrow = 1),
                                        age = list(direction = "horizontal"),
                                        BMI = list(direction = "horizontal"),
                                        `n cultured` = list(direction = "horizontal"),
                                        gender = list(nrow = 1),
                                        modul. = list(nrow = 1),
                                        proph. = list(nrow = 1),
                                        P.aero = list(nrow = 1),
                                        location = list(nrow = 1) ),
         
         col = list(
           gender = gender_names,
           age = age_col_fun,
           age = age_col_fun,
           BMI = bmi_col_fun,
           P.aero = P.aero_names,
           modul. = modul_names,
           proph. = proph_names,
           location = location_names,
           `n cultured` = load_col_fun
         ))

       # row_ha = rowAnnotation(foo2 = runif(10), bar2 = anno_barplot(runif(10)))
       
       hm1 <- ComplexHeatmap::Heatmap( microbiology_binary+0,
                                       show_heatmap_legend = FALSE,
                                       column_title = "anthropometrics:",
                                col = c("white", "black"),
                                column_labels = rep("", nrow(mgdat)),
                                top_annotation = column_ha, 
                                row_names_gp = gpar(fontface = "italic"),
                                heatmap_legend_param = list(direction = "horizontal"))
       draw( hm1, annotation_legend_side = "bottom")
       
       
       # png(filename = "vis/mtu__paedcf__Anthropometrics_heatmap.png", width = 500, height = 500)
       draw(hm1, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
       # dev.off()
       

## SPECIES level Sensitivity / Specificity ,  PPV/NPV   ==============================================================
    
      # you were very close! 
      # https://uk.cochrane.org/news/sensitivity-and-specificity-explained-cochrane-uk-trainees-blog
        
          ## Bug[i]_accuracy = (true_neg + true_pos)/nsamp 
      
      ## PPV = TP / TP + FP  (true_pos below)
      ## NPV = TN / TN + FN  (true_neg below)1
      ## sensitivity = TP / (TP + FN) , i.e. diagnosed diseased / total actual diseased
      ## specificity = TN / (TN + FP),  i.e. diagnosed healthy / total actual healthy
    
        
      # abundances below thresh don't count as being present, which is why it boosts specificity so much
        ##  ---  low values give a very hihg NPV and nothing else   -  i.e, super IN-sensitive MGX, so never there and all Negatives map as true
        ##  ---  high values give a very high PPV and nothing else  -  i.e, super sensitive MGX, so always there and all Positives map as true
        ##  ---  sens/spec for both remain lpoor
      # thresh <- 0.025
      # thresh <- 0.001
      thresh <- 0.0001
      # thresh <- 0.00    ## issue here is that nothing ends nothing ever ends

      acc_list <- lapply( c("BAL", "NS", "TS"), function(mm){                 # mm <- "BAL"
        
        nn_df <- dplyr::filter( mgdat, substrate == mm)
        oo_feat <- mgfeat_ra[ , rownames(nn_df)]
        
        # pp_dat <- do.call( "rbind" , lapply( micro_uniq[ -c(1,2,6,13)], function(aa){         # aa <- micro_uniq[10]
        pp_dat <- t( sapply( micro_uniq[ -c(1,2,6,13)], function(aa){         # aa <- micro_uniq[3]
          # aa_kid <- mgtax[ grep( gsub(" ", "_", aa), mgtax[ rownames(oo_feat) , "s" ]) , "kid" ] 
          aa_kid <- mgtax[ grep( gsub(" ", "_", aa), mgtax[ rownames(oo_feat) , "s" ]) , "s" ]          ## changed for Kraken2..
          if( length(aa_kid) == 0 ){ 
            c( "micro" = aa, "sens" = NA, "spec" = NA, "ppv" = NA, "npv" = NA, "mean_ra" = NA)
          }else{
            bb_df <- mgfeat_ra[ aa_kid, ]
            # all_abs <- sum( !grepl( aa, nn_df$micro_unif ))
            # all_pres <- sum( grepl( aa, nn_df$micro_unif ))
  
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
            
            mean_ra <- mean(bb_df)
            
            c( "micro" = aa , 
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
      
      acc_datish <-as.data.frame( do.call("rbind", acc_list), stringsAsFactors = FALSE )
      acc_dat <- data.frame( t( apply( acc_datish, 1, function(aa){ gsub("NaN", "0", aa) })), stringsAsFactors = FALSE )
      acc_dat$spec <- as.numeric( as.character(acc_dat$spec))
      acc_dat$sens <- as.numeric( as.character(acc_dat$sens))
      acc_dat$ppv <- as.numeric( as.character(acc_dat$ppv))
      acc_dat$npv <- as.numeric( as.character(acc_dat$npv))
      acc_dat$mean_ra <- as.numeric( as.character(acc_dat$mean_ra))
      head(acc_dat)
      
      
      head( acc_dat_m <- reshape2::melt( acc_dat, id.vars = c("micro", "substrate")) )
      acc_dat_m$value <- as.numeric(acc_dat_m$value)
      
      
##  G E N U S   level - simply need ot agglomerate abundances matching glom_rank   =============================
      
      micro_gu <- unique( gsub(" .*", "", micro_uniq) )
      g_acc_list <- lapply( c("BAL", "NS", "TS"), function(mm){                 # mm <- "BAL"
        
        nn_df <- dplyr::filter( mgdat, substrate == mm)
        oo_feat <- mgfeat_ra[ , rownames(nn_df)]
        
        pp_dat <- t( sapply( micro_gu[ -c(1,2,6)], function(aa){         # aa <- micro_gu[ -c(1,2,6)][6]
          
          if( !any( grepl( aa, mgtax[ rownames(oo_feat), "g"] )) ){         
            c( "micro" = aa, "sens" = NA, "spec" = NA, "ppv" = NA, "npv" = NA, "mean_ra" = NA)
          }else{ 
            if( sum(grepl( aa, mgtax[ rownames(oo_feat), "g"] )) == 1 ){                                ## changed for kraken2 - manage case of 1-row matrices
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
            
            c( "micro" = aa , 
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
      g_acc_dat$spec <- as.numeric( as.character(g_acc_dat$spec))
      g_acc_dat$sens <- as.numeric( as.character(g_acc_dat$sens))
      g_acc_dat$ppv <- as.numeric( as.character(g_acc_dat$ppv))
      g_acc_dat$npv <- as.numeric( as.character(g_acc_dat$npv))
      g_acc_dat$mean_ra <- as.numeric( as.character(g_acc_dat$mean_ra))
      head(g_acc_dat)
      
      head( g_acc_dat_m <- reshape2::melt( g_acc_dat, id.vars = c("micro", "substrate")) )
      g_acc_dat_m$value <- as.numeric(g_acc_dat_m$value)
      
      
##  v e r s u s    B A L   ===============================================================
      
      micro_gu <- unique( gsub(" .*", "", micro_uniq) )
      gB_acc_list <- lapply( c("BAL", "NS", "TS"), function(mm){                 # mm <- "NS"
        
        nn_df <- dplyr::filter( mgdat, substrate == mm)
        oo_feat <- mgfeat_ra[ , rownames(nn_df)]
        
        ## need to get BAL also
        bal_df <- dplyr::filter( mgdat, substrate == "BAL" )
        bal_feat <- mgfeat_ra[ , rownames(bal_df)]
        
        
        ## here as above, but we match mm culture work with bal_feat abundance 
        ## so find the corresponding SUBJECT's sample in BAL
        pp_dat <- t( sapply( micro_gu[ -c(1,2,6)], function(aa){         # aa <- micro_gu[4]
          
          if( !any( grepl(aa, mgtax[ rownames(oo_feat), "g"] )) ){
            print("missing from feature abundance")
            c( "micro" = aa, "sens" = NA, "spec" = NA, "ppv" = NA, "npv" = NA, "mean_ra" = NA )
          }else{
            if( sum(grepl( aa, mgtax[ rownames(oo_feat), "g"] )) == 1 ){                             ## changed for kraken2 - manage case of 1-row matrices
              bb_df <- mgfeat_ra[ grepl(aa, rownames(mgfeat_ra)), ] 
            }else{
              bb_df <- colSums( mgfeat_ra[ grepl(aa, rownames(mgfeat_ra)), ] )
            }
            
            # get samp, get samp_BAL, compare abund/culture, flatten if complicated
            # some samples have apparently complicated origin - eg SC136C2 / SC136C3 - no TS, possible duplication?
            true_pos <- sum( sapply( as.character(nn_df$participant), function(bbb){  #  bbb <- "SC136C2" bbb <- as.character(nn_df$participant)[3] 
              ccc <- as.character(dplyr::filter( nn_df, participant == bbb)$sample) 
              ccc_bal <- as.character(dplyr::filter( bal_df, participant == bbb)$sample)
              if( length(ccc)==0 | length(ccc_bal) ==0){ NA }else{
                ddd <- bb_df[ ccc_bal ] > thresh & grepl( aa, nn_df[ ccc, "micro_unif" ])
                if( length(ddd)==0){FALSE}else if( sum(ddd) == 0){FALSE}else if( sum(ddd)>0 ){TRUE}
              }
            }), na.rm = TRUE)
            
            true_neg <- sum( sapply( as.character(nn_df$participant), function(bbb){  #  bbb <- "SC136C2" bbb <- as.character(nn_df$participant)[3] 
              ccc <- as.character(dplyr::filter( nn_df, participant == bbb)$sample) 
              ccc_bal <- as.character(dplyr::filter( bal_df, participant == bbb)$sample)
              if( length(ccc)==0 | length(ccc_bal) ==0){ NA }else{
                ddd <- bb_df[ ccc_bal ] < thresh & !grepl( aa, nn_df[ ccc, "micro_unif" ])
                if( length(ddd)==0){FALSE}else if( sum(ddd) == 0){FALSE}else if( sum(ddd)>0 ){TRUE}
              }
            }), na.rm = TRUE)
            
            false_pos <- sum( sapply( as.character(nn_df$participant), function(bbb){  #  bbb <- "SC136C2" bbb <- as.character(nn_df$participant)[3] 
              ccc <- as.character(dplyr::filter( nn_df, participant == bbb)$sample) 
              ccc_bal <- as.character(dplyr::filter( bal_df, participant == bbb)$sample)
              if( length(ccc)==0 | length(ccc_bal) ==0){ NA }else{
                ddd <- bb_df[ ccc_bal ] < thresh & grepl( aa, nn_df[ ccc, "micro_unif" ])
                if( length(ddd)==0){FALSE}else if( sum(ddd) == 0){FALSE}else if( sum(ddd)>0 ){TRUE}
              }
            }), na.rm = TRUE)
            
            false_neg <- sum( sapply( as.character(nn_df$participant), function(bbb){  #  bbb <- "SC136C2" bbb <- as.character(nn_df$participant)[3] 
              ccc <- as.character(dplyr::filter( nn_df, participant == bbb)$sample) 
              ccc_bal <- as.character(dplyr::filter( bal_df, participant == bbb)$sample)
              if( length(ccc)==0 | length(ccc_bal) ==0){ NA }else{
                ddd <- bb_df[ ccc_bal ] > thresh & !grepl( aa, nn_df[ ccc, "micro_unif" ])
                if( length(ddd)==0){FALSE}else if( sum(ddd) == 0){FALSE}else if( sum(ddd)>0 ){TRUE}
              }
            }), na.rm = TRUE)
            
            mean_ra <- sum(bb_df)
            
            c( "micro" = aa , 
               "sens" = (true_pos/ (true_pos + false_neg)) ,
               "spec" = (true_neg/ (true_neg + false_pos)) ,
               "ppv" = (true_pos/ (true_pos + false_pos)) ,
               "npv" = (true_neg/ (true_neg + false_neg)),
               "mean_ra" = mean_ra
            ) 
            
          }
        }) )
        cbind( "substrate" = mm, pp_dat)
      })
      
      
      gB_acc_datish <-as.data.frame( do.call("rbind", gB_acc_list), stringsAsFactors = FALSE )
      gB_acc_dat <- data.frame( t( apply( gB_acc_datish, 1, function(aa){ gsub("NaN", "0", aa) })), stringsAsFactors = FALSE )
      gB_acc_dat$spec <- as.numeric( as.character(gB_acc_dat$spec))
      gB_acc_dat$sens <- as.numeric( as.character(gB_acc_dat$sens))
      gB_acc_dat$ppv <- as.numeric( as.character(gB_acc_dat$ppv))
      gB_acc_dat$npv <- as.numeric( as.character(gB_acc_dat$npv))
      gB_acc_dat$mean_ra <- as.numeric( as.character(gB_acc_dat$mean_ra))
      head(gB_acc_dat)
      
      
      head( gB_acc_dat_m <- reshape2::melt( gB_acc_dat, id.vars = c("micro", "substrate")) )
      gB_acc_dat_m$value <- as.numeric(gB_acc_dat_m$value)
      
    
##  S P E C I E S   v e r s u s    B A L   ===============================================================
    
    micro_su <- micro_uniq
    sB_acc_list <- lapply( c("BAL", "NS", "TS"), function(mm){                 # mm <- "NS"
      
      nn_df <- dplyr::filter( mgdat, substrate == mm)
      oo_feat <- mgfeat_ra[ , rownames(nn_df)]
      
      ## need to get BAL also
      bal_df <- dplyr::filter( mgdat, substrate == "BAL" )
      bal_feat <- mgfeat_ra[ , rownames(bal_df)]
      
      
      ## here as above, but we match mm culture work with bal_feat abundance 
      ## so find the corresponding SUBJECT's sample in BAL
      pp_dat <- t( sapply( micro_su[ -c(1,2,6)], function(aa){         # aa <- micro_su[4]
        
        # if( !any(grepl(aa, rownames(oo_feat))) ){
        if( !any( grepl( gsub(" ", "_", aa), mgtax[ rownames(oo_feat), "s"] )) ){
          print("missing from feature abundance")
          c( "micro" = aa, "sens" = NA, "spec" = NA, "ppv" = NA, "npv" = NA, "mean_ra" = NA)
        }else{
          bb_feat <- mgfeat_ra[ grepl(gsub(" ","_",aa),mgtax[ rownames(oo_feat),"s"]) ,  ]
          if( is.null(dim(bb_feat)) ){
            bb_df <- bb_feat 
          }else{ 
            bb_df <- colSums( mgfeat_ra[ grepl(aa, rownames(mgfeat_ra)),  ] )
          }
        
        # get samp, get samp_BAL, compare abund/culture, flatten if complicated
        # some samples have apparently complicated origin - eg SC136C2 / SC136C3 - no TS, possible duplication?
        true_pos <- sum( sapply( as.character(nn_df$participant), function(bbb){  #  bbb <- "SC136C2" bbb <- as.character(nn_df$participant)[3] 
          ccc <- as.character(dplyr::filter( nn_df, participant == bbb)$sample) 
          ccc_bal <- as.character(dplyr::filter( bal_df, participant == bbb)$sample)
          if( length(ccc)==0 | length(ccc_bal) ==0){ NA }else{
            ddd <- bb_df[ ccc_bal ] > thresh & grepl( aa, nn_df[ ccc, "micro_unif" ])
            if( length(ddd)==0){FALSE}else if( sum(ddd) == 0){FALSE}else if( sum(ddd)>0 ){TRUE}
          }
        }), na.rm = TRUE)
        
        true_neg <- sum( sapply( as.character(nn_df$participant), function(bbb){  #  bbb <- "SC136C2" bbb <- as.character(nn_df$participant)[3] 
          ccc <- as.character(dplyr::filter( nn_df, participant == bbb)$sample) 
          ccc_bal <- as.character(dplyr::filter( bal_df, participant == bbb)$sample)
          if( length(ccc)==0 | length(ccc_bal) ==0){ NA }else{
            ddd <- bb_df[ ccc_bal ] < thresh & !grepl( aa, nn_df[ ccc, "micro_unif" ])
            if( length(ddd)==0){FALSE}else if( sum(ddd) == 0){FALSE}else if( sum(ddd)>0 ){TRUE}
          }
        }), na.rm = TRUE)
        
        false_pos <- sum( sapply( as.character(nn_df$participant), function(bbb){  #  bbb <- "SC136C2" bbb <- as.character(nn_df$participant)[3] 
          ccc <- as.character(dplyr::filter( nn_df, participant == bbb)$sample) 
          ccc_bal <- as.character(dplyr::filter( bal_df, participant == bbb)$sample)
          if( length(ccc)==0 | length(ccc_bal) ==0){ NA }else{
            ddd <- bb_df[ ccc_bal ] < thresh & grepl( aa, nn_df[ ccc, "micro_unif" ])
            if( length(ddd)==0){FALSE}else if( sum(ddd) == 0){FALSE}else if( sum(ddd)>0 ){TRUE}
          }
        }), na.rm = TRUE)
        
        false_neg <- sum( sapply( as.character(nn_df$participant), function(bbb){  #  bbb <- "SC136C2" bbb <- as.character(nn_df$participant)[3] 
          ccc <- as.character(dplyr::filter( nn_df, participant == bbb)$sample) 
          ccc_bal <- as.character(dplyr::filter( bal_df, participant == bbb)$sample)
          if( length(ccc)==0 | length(ccc_bal) ==0){ NA }else{
            ddd <- bb_df[ ccc_bal ] > thresh & !grepl( aa, nn_df[ ccc, "micro_unif" ])
            if( length(ddd)==0){FALSE}else if( sum(ddd) == 0){FALSE}else if( sum(ddd)>0 ){TRUE}
          }
        }), na.rm = TRUE)
        
        mean_ra <- sum(bb_df)
        
        c( "micro" = aa , 
           "sens" = (true_pos/ (true_pos + false_neg)) ,
           "spec" = (true_neg/ (true_neg + false_pos)) ,
           "ppv" = (true_pos/ (true_pos + false_pos)) ,
           "npv" = (true_neg/ (true_neg + false_neg)),
           "mean_ra" = mean_ra
        ) 
        
        }
      }) )
      cbind( "substrate" = mm, pp_dat)
    })
    
    
    sB_acc_datish <-as.data.frame( do.call("rbind", sB_acc_list), stringsAsFactors = FALSE )
    sB_acc_dat <- data.frame( t( apply( sB_acc_datish, 1, function(aa){ gsub("NaN", "0", aa) })), stringsAsFactors = FALSE )
    sB_acc_dat$spec <- as.numeric( as.character(sB_acc_dat$spec))
    sB_acc_dat$sens <- as.numeric( as.character(sB_acc_dat$sens))
    sB_acc_dat$ppv <- as.numeric( as.character(sB_acc_dat$ppv))
    sB_acc_dat$npv <- as.numeric( as.character(sB_acc_dat$npv))
    sB_acc_dat$mean_ra <- as.numeric( as.character(sB_acc_dat$mean_ra))
    head(sB_acc_dat)
    
    
    head( sB_acc_dat_m <- reshape2::melt( sB_acc_dat, id.vars = c("micro", "substrate")) )
    sB_acc_dat_m$value <- as.numeric(sB_acc_dat_m$value)
    
    
## NMDS of Jaccardian culture  ----------------------------------------------- 
    
    # use bug dat and micro_uniq to create a matrix
    head(bug_dat)
    micro_uniq
    
    micro_uniq %in% bug_dat$micro
    unique(bug_dat$micro) %in% micro_uniq
    
    set.seed(2205)
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
    p_s$nom <- gsub("^(\\w)\\w* ", "\\1. ", rownames(p_s))

    ggplot( p_m, aes(fill = substrate, colour = substrate, x = NMDS1, y = NMDS2)) + 
      geom_text_repel(aes(label = sample_proc), size = 3, max.overlaps = 30, min.segment.length = 0, segment.linetype = 3) + #, position = position_jitterdodge(jitter.width = 0.4, jitter.height = 0.4)) +
      geom_path(aes(group = parti), colour = "grey30", size = 0.1, lty = 1) +
      geom_label(data = p_s, aes(label = nom, fill = NULL), size = 5, colour = "grey40", alpha  =0.5) #, position = position_jitterdodge(jitter.width = 1, jitter.height = 1))
    
    
      
      
##   T A B L E   =======================================================================    
      
      knitr::kable(cbind(
        # g_acc_dat[ , 1:2],
        "sens" = aggregate( sens ~ micro, data=g_acc_dat, FUN = mean),
        "spec" = aggregate( spec ~ micro, data=g_acc_dat, FUN = mean)[,2],
        "ppv" = aggregate( ppv ~ micro, data=g_acc_dat, FUN = mean)[,2],
        "npv" = aggregate( npv ~ micro, data=g_acc_dat, FUN = mean)[,2],
        "wrong RA" = aggregate( mean_ra ~ micro, data=g_acc_dat, FUN = mean)[,2] )
        )
      
      
##   P L O T S    =======================================================================    
      
      # dev.off()
      
      
      ## ABUNDANCE in the metagenomics
      # dev.new()
      # dev.new()
      plot4 <- ggplot( mp_m, aes( Var2, value)) +
        coord_flip() +
        facet_grid(  decis ~ substrate ) +
        geom_col(aes(fill = Var2)) +
        labs(
          x = "cultured microbes",
          y = "summed %",
          title = "abundance (WGS)"
        ) +
        # margin(0,0,0,2, unit =  "pt") +
        # scale_fill_manual( values = nr_col) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 0.9,  vjust = 0.9),
          legend.position = "none",
          strip.text.y = element_text(angle = 0),
        )
      
      
      ## OCCURRENCE in the metagenomics
      
      # dev.new()
      plot5 <- ggplot( mp_bin_m, aes( Var2 , value)) +
        coord_flip() +
        facet_grid(  decis ~ substrate ) +
        geom_col(aes(fill = Var2)) +
        # scale_fill_manual( values = nr_col) +
        labs(
          x = "",
          y = "incidence",
          title = "incidence (WGS)"
        ) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 0.9,  vjust = 0.9),
          axis.text.y = element_blank(),
          legend.position = "none",
          strip.text.y = element_text(angle = 0),
        )
      
      
      ## occurrence in the culture-work
      
      # dev.new()
      plot6 <- ggplot( metad_melt, aes( micro) ) +
        coord_flip() +
        facet_grid(  decis ~ substrate ) +
        geom_bar(aes(fill = micro)) +
        # scale_fill_manual( values = nr_col) +
        labs(
          x = "",
          y = "incidence",
          title = "incidence (culture-work)"
        ) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 0.9,  vjust = 0.9, size = 12),
          # axis.text.y = element_blank(),
          legend.position = "none",
          strip.text.y = element_text(angle = 0),
        )
      
      
      plot7 <- ggplot( bug_dat, aes( micro, incidence)) +
        coord_flip() +
        facet_grid(  . ~ substrate ) +
        geom_col(aes(fill = micro)) +
        labs(
          x = "cultured microbes",
          y = "summed %",
          title = "abundance (WGS)"
        ) +
        # margin(0,0,0,2, unit =  "pt") +
        # scale_fill_manual( values = nr_col) +
        theme(
          axis.text.y = element_text(size = 12, face = "italic"),
          axis.text.x = element_text(angle = 45, hjust = 0.9,  vjust = 0.9, size = 12),
          legend.position = "none",
          strip.text.y = element_text(angle = 0),
          # axis.text.x = element_text(angle = 45, hjust = 0.9,  vjust = 0.9, size = 12),
          # legend.position = "none",
          # strip.text.y = element_text(angle = 0),
        )
      # ggsave( plot7, filename="vis/mtu__paedcf__micro_culture_abundance.png", device = "png", width = 8, height = 10)
      
      
      plot8 <- ggplot( acc_dat_m, aes( micro, value)) +
        coord_flip() +
        facet_grid( variable ~ substrate ) +
        geom_col(aes(fill = micro)) +
        labs(
          x = "cultured species",
          y = "% correspondence with metagenomic assignment",
          title = "sensitivity of culture-work (w.r.t. metagenomic sequencing: 1% threshold)"
        ) +
        # margin(0,0,0,2, unit =  "pt") +
        # scale_fill_manual( values = nr_col) +
        theme(
          axis.text.y = element_text(size = 12, face = "italic"),
          axis.text.x = element_text(angle = 45, hjust = 0.9,  vjust = 0.9, size = 12),
          legend.position = "none",
          strip.text.y = element_text(angle = 0),
        )
      # ggsave( plot8, filename="vis/mtu__paedcf__micro_culture_accuracy_0.00.png", device = "png", width = 8, height = 10)
      # ggsave( plot8, filename="vis/mtu__paedcf__micro_culture_accuracy_0.01.png", device = "png", width = 8, height = 10)

            
      plot9 <- ggplot( acc_dat, aes( sens, spec, fill = micro )) +
        coord_fixed(xlim = c(0,1), ylim = c(0,1)) +
        geom_abline(intercept = c(0,0), slope = 1, colour = "grey50", lty = 2) +
        facet_grid( . ~substrate ) +
        geom_point(shape = 21, aes(size = mean_ra), alpha = 0.5 ) +
        labs(
          x = "sensitivity\n(true culture positive / total MGX positive)",
          y = "specificity\n(true culture negative / total MGX negative)",
          title = paste0("species-level sensitivity/specificity (MGX detection >", thresh, " threshold)")
        ) +
        # margin(0,0,0,2, unit =  "pt") +
        # scale_fill_manual( values = nr_col) +
        theme(
          axis.text.y = element_text(size = 12, face = "italic"),
          axis.text.x = element_text(angle = 45, hjust = 0.9,  vjust = 0.9, size = 12),
          legend.position = "bottom",
          strip.text.y = element_text(angle = 0),
        )
      # ggsave( plot9, filename="vis/mtu__paedcf__micro_culture_accuracy_0.00.png", device = "png", width = 8, height = 10)
      # ggsave( plot9, filename="vis/mtu__paedcf__micro_culture_accuracy_0.01.png", device = "png", width = 8, height = 10)

            
      plot10 <- ggplot( acc_dat, aes( ppv, npv, fill = micro )) +
        coord_fixed(xlim = c(0,1), ylim = c(0,1)) +
        geom_abline(intercept = c(0,0), slope = 1, colour = "grey50", lty = 2) +
        facet_grid( . ~substrate ) +
        geom_point(shape = 21, aes(size = mean_ra), alpha = 0.5 ) +
        geom_text_repel(aes( label = gsub("(.).* (.*)", "\\1. \\2", micro, perl = TRUE) ),
                        direction = "both", max.time = 2, max.iter = Inf,
                        size = 4, segment.color = "grey80") +
        labs(
          x = "positive predictive value\n(true culture positive / all culture positive)",
          y = "negative predictive value\n(true culture negative / all culture negative)",
          title = paste0("species-level diagnostic value (MGX detection >", thresh, " threshold)")
        ) +
        # margin(0,0,0,2, unit =  "pt") +
        # scale_fill_manual( values = nr_col) +
        theme(
          axis.text.y = element_text(size = 12, face = "italic"),
          axis.text.x = element_text(angle = 45, hjust = 0.9,  vjust = 0.9, size = 12),
          legend.position = "bottom",
          strip.text.y = element_text(angle = 0),
        )
      # ggsave( plot10, filename="vis/mtu__paedcf__micro_culture_accuracy_0.00.png", device = "png", width = 8, height = 10)
      # ggsave( plot10, filename="vis/mtu__paedcf__micro_culture_accuracy_0.01.png", device = "png", width = 8, height = 10)
      # ggsave( plot10, filename="vis/mtu__paedcf__micro_culture_accuracy_0.025.png", device = "png", width = 8, height = 10)

      
      
      plot11 <- ggplot( g_acc_dat, aes( sens, spec, fill = micro )) +
        coord_fixed(xlim = c(0,1.3), ylim = c(0,1.3)) +
        geom_abline(intercept = c(0,0), slope = 1, colour = "grey65", lty = 2) +
        facet_grid( . ~substrate ) +
        geom_point(shape = 21, aes(size = mean_ra), alpha = 0.5 ) +
        geom_text_repel(aes( label = gsub("(.).* (.*)", "\\1. \\2", micro, perl = TRUE) ),
                        direction = "both", max.time = 2, max.iter = Inf,
                        size = 3, fontface = "italic", segment.color = "grey80", max.overlaps = 100) +
        labs(
          x = "sensitivity\n(true culture positive / total MGX positive)",
          y = "specificity\n(true culture negative / total MGX negative)",
          title = paste0("genus sensitivity/specificity at ", thresh)
        ) +
        # margin(0,0,0,2, unit =  "pt") +
        # scale_fill_manual( values = nr_col) +
        theme(
          axis.text.y = element_text(size = 12, face = "italic"),
          axis.text.x = element_text(angle = 45, hjust = 0.9,  vjust = 0.9, size = 12),
          legend.position = "bottom",
          strip.text.y = element_text(angle = 0),
        )

            
      plot12 <- ggplot( g_acc_dat, aes( ppv, npv, fill = micro )) +
        coord_fixed(xlim = c(0,1.3), ylim = c(0,1.3)) +
        geom_abline(intercept = c(0,0), slope = 1, colour = "grey65", lty = 2) +
        facet_grid( . ~substrate ) +
        geom_point(shape = 21, aes(size = mean_ra), alpha = 0.5 ) +
        geom_text_repel(aes( label = gsub("(.).* (.*)", "\\1. \\2", micro, perl = TRUE) ),
                        direction = "both", max.time = 2, max.iter = Inf,
                        size = 3, fontface = "italic", segment.color = "grey80") +
        labs(
          x = "positive predictive value\n(true culture positive / all culture positive)",
          y = "negative predictive value\n(true culture negative / all culture negative)",
          title = paste0("genus diagnostic ability at ", thresh)
        ) +
        # margin(0,0,0,2, unit =  "pt") +
        # scale_fill_manual( values = nr_col) +
        theme(
          axis.text.y = element_text(size = 12, face = "italic"),
          axis.text.x = element_text(angle = 45, hjust = 0.9,  vjust = 0.9, size = 12),
          legend.position = "bottom",
          strip.text.y = element_text(angle = 0),
        )
      
      plot_11_12 <- ggpubr::ggarrange( plot11, plot12, ncol = 1, nrow = 2, common.legend = TRUE, legend = "right")
      plot_11_12
      # ggsave( plot_11_12, filename="vis/mtu__paedcf__micro_culture_SS_PNPV.png", device = "png", width = 15, height = 10)

      
      
    # diagnose w.r.t. BAL
      
      
      plot13 <- ggplot( gB_acc_dat, aes( sens, spec, fill = micro )) +
        coord_fixed(xlim = c(0,1.05), ylim = c(0,1.05)) +
        geom_abline(intercept = c(0,0), slope = 1, colour = "grey65", lty = 2) +
        facet_grid( . ~substrate ) +
        geom_point(shape = 21, aes(size = mean_ra), alpha = 0.5 ) +
        geom_text_repel(aes( label = gsub("(.).* (.*)", "\\1. \\2", micro, perl = TRUE) ),
                        direction = "both", max.time = 2, max.iter = Inf,
                        size = 3, fontface = "italic", segment.color = "grey80", max.overlaps = 100) +
        labs(
          x = "sensitivity\n(true culture positive / total MGX positive in BAL)",
          y = "specificity\n(true culture negative / total MGX negative in BAL)",
          title = paste0("genus sensitivity/specificity to BAL abundance at ", thresh)
        ) +
        # margin(0,0,0,2, unit =  "pt") +
        # scale_fill_manual( values = nr_col) +
        theme(
          axis.text.y = element_text(size = 12, face = "italic"),
          axis.text.x = element_text(angle = 45, hjust = 0.9,  vjust = 0.9, size = 12),
          legend.position = "bottom",
          strip.text.y = element_text(angle = 0),
        )
      
      
      plot14 <- ggplot( gB_acc_dat, aes( ppv, npv, fill = micro )) +
        coord_fixed(xlim = c(0,1.3), ylim = c(0,1.3)) +
        geom_abline(intercept = c(0,0), slope = 1, colour = "grey65", lty = 2) +
        facet_grid( . ~substrate ) +
        geom_point(shape = 21, aes(size = mean_ra), alpha = 0.5 ) +
        geom_text_repel(aes( label = gsub("(.).* (.*)", "\\1. \\2", micro, perl = TRUE) ),
                        direction = "both", max.time = 2, max.iter = Inf,
                        size = 3, fontface = "italic", segment.color = "grey80") +
        labs(
          x = "positive predictive value\n(true culture positive / all culture positive in BAL)",
          y = "negative predictive value\n(true culture negative / all culture negative in BAL)",
          title = paste0("genus ability to diagnose BAL at ", thresh)
        ) +
        # margin(0,0,0,2, unit =  "pt") +
        # scale_fill_manual( values = nr_col) +
        theme(
          axis.text.y = element_text(size = 12, face = "italic"),
          axis.text.x = element_text(angle = 45, hjust = 0.9,  vjust = 0.9, size = 12),
          legend.position = "bottom",
          strip.text.y = element_text(angle = 0),
        )
      
      plot_13_14 <- ggpubr::ggarrange( plot13, plot14, ncol = 1, nrow = 2, common.legend = TRUE, legend = "right")
      plot_13_14
      # ggsave( plot_13_14, filename="vis/mtu__paedcf__micro_culture_SS_PNPV_wrt-BAL.png", device = "png", width = 15, height = 10, bg = "white")
      

    ## spec BAL      
      
      (plot15 <- ggplot( sB_acc_dat, aes( ppv, npv, fill = micro )) +
        coord_fixed(xlim = c(0,1.3), ylim = c(0,1.3)) +
        geom_abline(intercept = c(0,0), slope = 1, colour = "grey65", lty = 2) +
        facet_grid( . ~substrate ) +
        geom_point(shape = 21, aes(size = mean_ra), alpha = 0.5 ) +
        geom_text_repel(aes( label = gsub("(.).* (.*)", "\\1. \\2", micro, perl = TRUE) ),
                        direction = "both", max.time = 2, max.iter = Inf,
                        size = 3, fontface = "italic", segment.color = "grey80", max.overlaps = 100) +
        labs(
          x = "positive predictive value\n(true culture positive / all culture positive in BAL)",
          y = "negative predictive value\n(true culture negative / all culture negative in BAL)",
          title = paste0("genus ability to diagnose BAL at ", thresh)
        ) +
        # margin(0,0,0,2, unit =  "pt") +
        # scale_fill_manual( values = nr_col) +
        theme(
          axis.text.y = element_text(size = 12, face = "italic"),
          axis.text.x = element_text(angle = 45, hjust = 0.9,  vjust = 0.9, size = 12),
          legend.position = "bottom",
          strip.text.y = element_text(angle = 0),
        ))

      
      
      
      
      # plot3   # NMDS of pathogens
      # plot4   # abundance of pathogen species
      # plot5   # incidence of pathogen (genera?)
      # plot6   # incidence in culture work
      # plot7   # abundance of pathogen species  -  difference to 4?
      # plot8   # sens, spec, ppv, npv, meanRA of pathogens at threshold% with MGX-
      # plot9   # sens-spec plot
      # plot10  # ppv-npv plot
      # # dev.new()
      # plot_11_12  #  accuracy & diagnostic ability of MGX v. Culture
      # # dev.new()
      # plot_13_14  #  accuracy & diagnostic ability of MGX v. BAL culture
      
      
      # ## mgx abundance, MGX incidence, culture incidence
      # ggpubr::ggarrange( plot4, plot5, plot6, ncol = 3, nrow = 1, widths = c(0.5, 0.25, 0.25))
      # 
      # 
      # ## mgx abundance, culture incidence
      # ggpubr::ggarrange( plot4 +
      #        # scale_x_discrete(position = "top") +
      #        # scale_y_reverse() +
      #        labs(x = "") +
      #        theme( axis.text.y = element_text(colour = "grey10", size =12, face = "italic")), plot6, ncol = 2, nrow = 1, widths = c(0.65, 0.35))
      # 
      #   
      # ## mgx incidence, culture incidence
      # (plot_micro <- ggpubr::ggarrange( plot5 +
      #                      # scale_x_discrete(position = "top") +
      #                      # scale_y_reverse() +
      #                      labs(x = "") +
      #                      theme( axis.text.y = element_text(colour = "grey10", size =12, face = "italic")), plot6, ncol = 2, nrow = 1, widths = c(0.65, 0.35)) )
      
      # ggsave( plot_micro, filename="vis/mtu__paedcf__micro_culture_incidence.png", device = "png", width = 12, height = 10)
      
      
      
      # save.image("output/mtu__paedcf__microbiology_data.RData")
      
      