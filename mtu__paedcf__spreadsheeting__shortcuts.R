

## mtu__paedcf__spreadsheeting :: copies of copies of copies of copies of 

      # rm(list=ls())
      set.seed(2205)


## libs  =========

      library(vegan)
      library(ggplot2)
      # library(ggpubr)  #  what now
      library(patchwork)

      source("analysis/background_code/R__fns_jfg/fn_definitions.R")     # this is just using a filtering fn

      print("  + + +   mixed Metagenomics, copied from PROJECT run  ::  kraken2:conf=0.10, hits=5, bracken:--t = 50    + + +")
      print("  + + +    | |                                                                     + + +")
      print("  + + +    '---->  negative  controls  empty  /  aspecific decontam only,          + + +")
      print("  + + +      '->   unpaired reads re-incorporated, even if mainly human            + + +")

      
# ## recycled ------------------------------------------------------------------
      clr <- function(aa){    print(paste("samples must be ROWS ; dim =", dim(aa)[1], "x", dim(aa)[2]))
        t(apply(aa, 1, function(bb) log(bb) - log(mean(bb)) ))
      }


## Kraken2 + unpaired ====================================================================================================================

  # Kraken2
      krapre_k2_brack_file="input/mtu__paedcf_krakenStnd_abundances.tsv"
      krapre_k2_mpa_file="input/mtu__paedcf__krakenStnd_taxonomy.tsv"
    # secondary, almost all just human sequences
      krapre_k2_brack_file__unp="input/mtu__paedcf_krakenStnd_abundances_unpaired.tsv"


    ## Re-arrange data   -----------------------------------

      pre_k2dat <- read.table(file = krapre_k2_brack_file, header=TRUE, quote = "", sep = '\t')
      pre_k2dat[1:10, 1:10]

      prenames_pre_k2dat <- paste0("k2_", stringr::str_pad(pre_k2dat[,2], width=7, pad=0))
      ## non-unique names -  assume only one dupe of each
      prenames_pre_k2dat[ duplicated(pre_k2dat[,2])] <- gsub('k2_', 'k2-A_', prenames_pre_k2dat[ duplicated(pre_k2dat[,2])])
      rownames(pre_k2dat) <- prenames_pre_k2dat

      # counts
      pre_k2_counts <- pre_k2dat[ , grep('_num', colnames(pre_k2dat))]
      colnames(pre_k2_counts) <- gsub( "_S\\d*.bracken_num", "", colnames(pre_k2_counts), perl = TRUE )

      # relab%
      pre_k2_relab <- pre_k2dat[ , grep('_frac', colnames(pre_k2dat))]
      colnames(pre_k2_relab) <- gsub( "_S\\d*.bracken_frac", "", colnames(pre_k2_relab), perl = TRUE )


    ## secondary arrangements with unpaired abundances
      pre_k2unpdat <- read.table(file = krapre_k2_brack_file__unp, header=TRUE, quote = "", sep = '\t')
      pre_k2unpdat[1:10, 1:10]

      prenames_pre_k2unpdat <- paste0("k2_", stringr::str_pad(pre_k2unpdat[,2], width=7, pad=0))
      ## non-unique names -  assume only one dupe of each
      prenames_pre_k2unpdat[ duplicated(pre_k2unpdat[,2])] <- gsub('k2_', 'k2-A_', prenames_pre_k2unpdat[ duplicated(pre_k2unpdat[,2])])
      rownames(pre_k2unpdat) <- prenames_pre_k2unpdat

      # counts
      pre_k2unp_counts <- pre_k2unpdat[ , grep('_num', colnames(pre_k2unpdat))]
      colnames(pre_k2unp_counts) <- gsub( "_S\\d*.bracken_num", "", colnames(pre_k2unp_counts), perl = TRUE )

      # relab%
      pre_k2unp_relab <- pre_k2unpdat[ , grep('_frac', colnames(pre_k2unpdat))]
      colnames(pre_k2unp_relab) <- gsub( "_S\\d*.bracken_frac", "", colnames(pre_k2unp_relab), perl = TRUE )


    ## congeal again. legacy material has samples are COLS
      dim(pre_k2_counts)
      dim(pre_k2unp_counts)
      length(dem_bugs <- unique (rownames(pre_k2_counts), rownames( pre_k2unp_counts)))
      length(dem_samps <- unique (colnames(pre_k2_counts), colnames( pre_k2unp_counts)))


      dim( k2_counts <-  t( sapply( dem_bugs, function(aa){   # aa <- dem_bugs[40]
        # sapply( dem_samps, function(aaa){
        colSums(rbind( pre_k2_counts[ aa, dem_samps ], pre_k2unp_counts[ aa, dem_samps ]), na.rm = TRUE)
      })) )
      sum( pre_k2_counts, pre_k2unp_counts ) == sum( k2_counts )


## tax -start with data from the abundance table names   ====================================

      tax_a <- data.frame(taxon=pre_k2dat[ ,1], k2_id=pre_k2dat[ ,2], row.names = rownames(pre_k2dat) , stringsAsFactors = FALSE)
      tax_a[,1] <- gsub('\\/','', tax_a[,1])
      tax_a[,1] <- gsub('ALOs_','ALOs-', tax_a[,1])
      # tax_a[,1] <- gsub("'","", tax_a[,1])      # none?
      head(tax_a)

      ## NOTE use quote="" to avoid nightmare of special characters like '"`/,| etc in names
      tax_b <- read.table( krapre_k2_mpa_file, sep='\t', header=FALSE, fill = TRUE, stringsAsFactors = FALSE, quote="")
      head(tax_b)


    ## regularise taxonomic ranks (not all have K:P:C:O:F:G:S etc.)   -----------------

      ranks <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__")
      tax_b <- t(apply(tax_b, 1, function(aa){  # aa <- tax_b[1,]
        sapply(ranks, function(aaa){ # aaa <- "c__"
          bbb <- grep(aaa, unlist(aa), value=TRUE)
          ifelse( length(bbb) == 0, "unkn.", bbb)
        })
      }))
      dim(tax_b) ; str(tax_b) ; head(tax_b)


    ## more tax issues due to names: doubled ranks:

      prior_index <- (unlist( lapply(1:nrow(tax_b), function(aa){ if( sum( grepl("unkn.", tax_b[aa,])) == 6){aa} }) ) - 1)
      ## compare:
      # tax_b[ sort( c(prior_index, prior_index+1)) , ]
      # tax_b[ prior_index , ]
      # tax_b[ prior_index+1 , ]
      # amend, by doubling up (will remove duplicates later)
      tax_b[ prior_index , "s__"] <- tax_b[ (prior_index+1) , "s__"]
      # row is only useful if has a species: cannot have a species, and hit 6, And be worth keeping
      tax_b <- tax_b[ !apply(tax_b, 1, function(aa) sum(grepl("unkn.", unlist(aa)) ) == 6) , ]

      colnames(tax_b) <- c("k", "p", "c", "o", "f", "g", "s")

      tax_a[ , "taxon"] <- gsub(" ", "_", paste0("s__", tax_a[ , "taxon"]))
      head(tax_a) ; head(tax_a[ , "taxon"]) ;  dim(tax_a)
      head(tax_b) ;   head(tax_b[ , "s"]) ;  dim(tax_b)


    ## stitch  &  name   ----------------------------------

      tax_c <- merge(tax_b, tax_a, by.x="s", by.y="taxon")[, c(2:7,1,8)]
      head(tax_c) ; dim(tax_c)
      prenames_tax <- paste0("k2_", stringr::str_pad(tax_c[,"k2_id"], width=7, pad=0))
      rownames(tax_c) <- prenames_tax
      head(tax_c)


    ## check all on the same page   -----------------------

      all(rownames(pre_k2dat) == rownames(pre_k2_counts) & rownames(pre_k2dat) == rownames(pre_k2_relab))
      all(rownames(pre_k2dat) %in% rownames(tax_c))

      kraken_ids <- sort(rownames(pre_k2dat))
      pre_k2dat <- pre_k2dat [kraken_ids , ]
      pre_k2_counts <- pre_k2_counts [kraken_ids , ]
      pre_k2_relab <- pre_k2_relab [kraken_ids , ]
      tax_c <- tax_c[kraken_ids , ]


    ## remove if no values (kraken2 set to give all taxa)

      pos_abund <- rownames(pre_k2_relab)[rowSums(pre_k2_relab) > 0]
      pos_abund_c <- colnames(pre_k2_relab)[colSums(pre_k2_relab) > 0]
      pre_k2_relab <- pre_k2_relab[ pos_abund , pos_abund_c ]
      pre_k2_counts <- pre_k2_counts[ pos_abund , pos_abund_c ]
      pre_k2_tax <- tax_c[ pos_abund , ]


##   T I E   T O   M E T A D A T A     =============================================================

      m30 <- as.data.frame(readxl::read_xlsx("input/SAINTS_metadata_updated_infl_jfg.xlsx", sheet = 1) )
      # View(m30)
      colnames(m30)

      # missing data
      colnames(pre_k2_counts)[ !( gsub("-.*", "", colnames(pre_k2_counts), perl = TRUE) %in% gsub( "_", ".", m30$`Subject Code`) ) ]
      # no missing samples!
      m30$`Subject Code`[ !( m30$`Subject Code` %in% gsub("-.*", "", colnames(pre_k2_counts), perl = TRUE)) ]


    ## some points miss representation for specific samples

      metad <- data.frame(
        "sample" = colnames(pre_k2_counts),
        "participant" = sapply(colnames(pre_k2_counts), function(aa){ gsub( "\\..*", "", aa)  }) ,
        "substrate" = sapply(colnames(pre_k2_counts), function(aa){ gsub(".*\\.(\\w*)", "\\1", aa) }),
        m30[ match( sapply(colnames(pre_k2_counts), function(aa){ gsub( "\\..*", "", aa)  }) , m30$`Subject Code`) , ]# aa <- colnames(pre_k2_counts)[4]
      )
      metad$countTotal <- colSums(pre_k2_counts[ , rownames(metad)])
      # remove the tag at the end..
      metad$substrate <- gsub("[a,b]$", "", metad$substrate, perl = TRUE)

      # View(metad)   #  note 5 empty rows at base

      metad$BMI.z.Score <- as.numeric(metad$BMI.z.Score)

      metad$proph_bin <- ifelse( metad$Prohphylaxis.ABs == "No", "none", "prophylaxis")
      metad$modul_bin <- ifelse( metad$Modulator...Start.Date == "No", "none", "modulator")
      head(metad)


      ## straighten the microflora
      metad$micro.BAL <-  unlist( lapply( metad$micro.BAL, function(aa){  #  aa <-  metad$micro.BAL[14]
        bb_vec <- gsub("^\\s", "", unlist(strsplit(aa, split = ",")), perl = TRUE)
        paste0( sort(bb_vec), collapse = ",") }))
      metad$micro.NS <-  unlist( lapply( metad$micro.NS, function(aa){  #  aa <-  metad$micro.BAL[14]
        bb_vec <- gsub("^\\s", "", unlist(strsplit(aa, split = ",")), perl = TRUE)
        paste0( sort(bb_vec), collapse = ",")  }))
      metad$micro.TS <-  unlist( lapply( metad$micro.TS, function(aa){  #  aa <-  metad$micro.BAL[14]
        bb_vec <- gsub("^\\s", "", unlist(strsplit(aa, split = ",")), perl = TRUE)
        paste0( sort(bb_vec), collapse = ",")  }))


      ## incorrect. rethink - only one micro profile per sample.

      metad$micro <- gsub(", ", ",",
                          unlist( lapply( rownames(metad), function(aa){  # aa <- rownames(metad)[86]
                            bb_df <- m30[ grep( gsub("\\..*", "", aa), m30$`Subject Code`) , ]
                            if( nrow(bb_df) == 0){
                              "NA"
                            }else if( grepl("BAL", aa)){
                              bb_df$micro.BAL
                            }else if( grepl("NS", aa)){
                              bb_df$micro.NS
                            }else if( grepl("TS", aa)){
                              bb_df$micro.TS
                            }else{
                              "NA"
                            }
                          }) )
      )

      # View(metad[ , c(1, 3, 17:19,25)])


      ## straighen out all the numeric inflamation data
      metad$infl.IL.8 <- as.numeric( as.character( metad$infl.IL.8))
      metad$infl.NE <- as.numeric( as.character( metad$infl.NE))
      metad$infl.TCC <- as.numeric( as.character( metad$infl.TCC))

      # View(metad)
      mgdat <- metad


    ## dissims - here for rownames and dissim   ------------------

      pre_k2_ra_bcdist <- as.dist(vegdist( t(pre_k2_relab[ , rownames(mgdat)]), method = "bray"))
      pre_k2_ra_jacdist <- as.dist(vegdist( t(pre_k2_relab[ , rownames(mgdat)]), method = "jacc"))


    ## back to mgdat ----

      mgdat$bc_wardD2 <- hclust( pre_k2_ra_bcdist, method = "ward.D2")$order
      mgdat$jac_wardD2 <- hclust( pre_k2_ra_bcdist, method = "ward.D2")$order

      mgdat$hclust <- mgdat$jac_wardD2


    ## check - do all match?   -----------------------

      dim(pre_k2_counts)
      dim(pre_k2_relab)
      dim(pre_k2_tax)
      dim(mgdat)

      # look
      pre_k2_counts[1:10, 1:10]
      pre_k2_relab[1:10, 1:10]

      # View(pre_k2_counts)
      # hist( log10(pre_k2_counts ))


## transformations etc   ---------------------------------------------------

    # CLR    -------------
      ## zCompositions assumes samples are COLUMNS
      # zPatterns(pre_k2_counts, 0)
      # heatmap( apply(pre_k2_counts[ , mgdat$bc_wardD2 ] > 0, 1, as.numeric), Rowv = NULL )

      dim(pre_k2_counts)
      print("  + + +   NOTE :: CZM used, and z.warn set to 0.9999 to avoid culling data    + + +")
      pre_k2_cmult <- zCompositions::cmultRepl( pre_k2_counts, method = "CZM", z.warning = 0.9999, )
      dim(pre_k2_cmult)
      pre_k2_clr <- apply( pre_k2_cmult, 2, function(aa){ log(aa) - mean(log(aa)) })
      dim(pre_k2_clr)


## write out   ---------------------------------------------------

    ## general metadata
      dim(mgdat)
      # dir.create("output/tsv")
      # write.table(mgdat, "output/tsv/mtu__paedcf__dat-88-31.tsv", sep = "\t")
      saveRDS( mgdat, "output/mtu__paedcf__dat-88-31.RDS")

    ## tsv output
      # write.table(pre_k2_counts, "output/tsv/mtu__paedcf__kraken.10.5.20.50-88__num.tsv", sep = "\t")
      # write.table(pre_k2_relab, "output/tsv/mtu__paedcf__kraken.10.5.20.50-88__frac.tsv", sep = "\t")
      # write.table(pre_k2_tax, "output/tsv/mtu__paedcf__kraken.10.5.20.50-88__tax.tsv", sep = "\t")

    ## RDS for R output - note samples are ROWS
      saveRDS( t(pre_k2_counts), "output/mtu__paedcf__kraken.10.5.20.50-88__num.RDS")
      saveRDS( t(pre_k2_relab), "output/mtu__paedcf__kraken.10.5.20.50-88__frac.RDS")
      saveRDS( t(pre_k2_ra_bcdist), "output/mtu__paedcf__kraken.10.5.20.50-88__BCdist.RDS")
      saveRDS( t(pre_k2_ra_jacdist), "output/mtu__paedcf__kraken.10.5.20.50-88__Jaccdist.RDS")
      saveRDS( t(pre_k2_clr), "output/mtu__paedcf__kraken.10.5.20.50-88__clr.RDS")
      saveRDS( pre_k2_tax, "output/mtu__paedcf__kraken.10.5.20.50-88__tax.RDS")


## read in   ---------------------------------------------------

      mgdat <- readRDS("output/mtu__paedcf__dat-88-31.RDS")
      mgfeat <- t(readRDS("output/mtu__paedcf__kraken.10.5.20.50-88__num.RDS"))
      mgfeat_ra <- t(readRDS("output/mtu__paedcf__kraken.10.5.20.50-88__frac.RDS"))
      mgtax <- readRDS("output/mtu__paedcf__kraken.10.5.20.50-88__tax.RDS")
      mgtax[ , "kid"] <- rownames(mgtax)
      
      mgfeat_clr <- t(readRDS("output/mtu__paedcf__kraken.10.5.20.50-88__clr.RDS"))
      # mgfeat_dist_bray <- readRDS("output/mtu__paedcf__kraken.10.5.20.50-88__BCdist.RDS")
      # mgfeat_dist_jac <- readRDS("output/mtu__paedcf__kraken.10.5.20.50-88__Jaccdist.RDS")
      mgfeat_bc <- readRDS("output/mtu__paedcf__kraken.10.5.20.50-88__BCdist.RDS")
      mgfeat_ja <- readRDS("output/mtu__paedcf__kraken.10.5.20.50-88__Jaccdist.RDS")

      
    ## emulation of Kaiju
      print("  + + +    be wide :   k2id rownames replaced with mgtax[ k2id , \"s\"]      + + +")
      rownames(mgfeat) <- mgtax[ rownames(mgfeat) , "s"]      
      rownames(mgfeat_ra) <- mgtax[ rownames(mgfeat_ra) , "s"]      
      rownames(mgfeat_clr) <- mgtax[ rownames(mgfeat_clr) , "s"]      
      rownames(mgtax) <- mgtax[ rownames(mgtax) , "s"]      
      
    
    ## v.0.3.5   -    remove the humans reads   ---------------------------------
      print("  + + +                                                                            + + +")
      print("  + + +    < ! >    v 0.3.5  ---    death to all humans' reads :   nhsfeat / _ra   + + +")
      
      ## standby problems.com      
      nhsfeat <- mgfeat[ -grep( "Homo", rownames(mgfeat)) , ]
      nhsfeat_ra <- mgfeat_ra[ -grep( "Homo", rownames(mgfeat_ra)) , ]
      nhsfeat_clr <- mgfeat_clr[ -grep( "Homo", rownames(mgfeat_clr)) , ]     ## not actually necessary
      
      # bc_no_Hsap <- as.dist(vegdist( t(nhsfeat_ra[ , rownames(mgdat)]), method = "bray"))
      # jac_no_Hsap <- as.dist(vegdist( t(nhsfeat_ra[ , rownames(mgdat)]), method = "jacc"))
      # euc_no_Hsap <- as.dist(vegdist( t(nhsfeat_clr[ , rownames(mgdat)]), method = "euc"))
      # # 
      # saveRDS(bc_no_Hsap,  "output/mtu__paedcf__kraken.10.5.20.50-88__noHsap__BCdist.RDS")
      # saveRDS(jac_no_Hsap, "output/mtu__paedcf__kraken.10.5.20.50-88__noHsap__Jaccdist.RDS")
      # saveRDS(euc_no_Hsap, "output/mtu__paedcf__kraken.10.5.20.50-88__noHsap__CLR-Eucdist.RDS")
      print("  + + +    < ! >     --- --- ---    though note manus_3.5 redoes this step...      + + +")
      print("  + + +                                                                            + + +")
       bc_no_Hsap <- readRDS("output/mtu__paedcf__kraken.10.5.20.50-88__noHsap__BCdist.RDS")
      jac_no_Hsap <- readRDS("output/mtu__paedcf__kraken.10.5.20.50-88__noHsap__Jaccdist.RDS")
      euc_no_Hsap <- readRDS("output/mtu__paedcf__kraken.10.5.20.50-88__noHsap__CLR-Eucdist.RDS")
      
      
      
## plotting info  =========
      
    ## colours for variables
      
      those_effing_cols <- c("#FBA90A","#9A0794","#4D854C")
      names(those_effing_cols) <- c("BAL", "MMS", "OPS")
      
      those_relic_cols <- c("#FBA90A","#9A0794","#4D854C")
      names(those_effing_cols) <- c("BAL", "NS", "TS")
      
      
      ## recycled ----------------------------------------------------------------------------
      ab_col <- ( c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a", "#FF1493",
                    "#b15928", "#737f33", "#8B008B", "#32fbd8", "#fdbf6f",
                    RColorBrewer::brewer.pal(5, "Spectral"),
                    "#b2df8a", "#fb9a99", "#d9e627", "#EE82EE", "#DEB887",
                    "#a6cee3",
                    RColorBrewer::brewer.pal(4,'Accent')
      ) )
      
      # colours, factors, whathavveyou  
      length(nr_col <- c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a", "#FF1493",
                         "#b15928", "#737f33", "#8B008B", "#32fbd8", "#fdbf6f",
                         RColorBrewer::brewer.pal(9, "Spectral")[],
                         "#b2df8a", "#fb9a99", "#d9e627", "#EE82EE", "#DEB887",
                         "#a6cee3" ))

      length(va_cols <- c("#20DE8B", "#CCDE8B", "#FFDE8B", "#FFA88B", "#FF6A8B", "#FF6AD5", "#C874AA", "#C774E7", "#AD8CFF", "#966BFF", "#90CFFF",
                          "#296656", "#569874", "#7EC488", "#A997AB", "#532E57", "#F9897B", "#D7509F", "#F9247E", "#AE1357", "#661246", "#9239F6",
                          "#903495", "#6F3460", "#4A354F", "#D20076", "#FF0076", "#FF4373", "#FF6B58", "#F8B660", "#7FD4C1", "#30BFDD", "#8690FF",
                          "#ACD0F4", "#F7C0BB", "#FBCFF3", "#65323E", "#FE7F9D", "#FFC0CB", "#75D8D5", "#09979B", "#063B41", "#7B556C", "#86486F",
                          "#F1956E", "#EB9B60", "#864A42", "#FAD36A", "#F5FFBF", "#EFDB9C" ))

            
## theme variables    -------------------------------------------------------------------
      
      print("  + + +   applying theme_update in spreadsheeting_shortcuts   + + +    ")
      theme_set( theme_minimal())
      theme_update(
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.line = element_line(colour = "grey80", size = 0.2),
        #
        legend.box = NULL,  # controls the direction of multiple legend boxes
        legend.position = "right",
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.75, "cm"),
        #
        panel.border = element_rect(colour = "grey20", fill=NA, linewidth=0.8),
        panel.spacing.y = unit(2, "lines"),
        panel.grid.major.x = element_line(colour = "lightskyblue3", linewidth = 0.1),
        panel.grid.major.y = element_line(colour = "lightskyblue3", linewidth = 0.1),
        panel.grid.minor.x = element_line(colour = "lightskyblue3", linewidth = 0.1),
        panel.grid.minor.y = element_line(colour = "lightskyblue3", linewidth = 0.1),
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
        
        text = element_text(colour = "grey20"),     # "mono", 
        NULL
      )

