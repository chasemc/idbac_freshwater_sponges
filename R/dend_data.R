
micromonospora_boundaries <- c(624, 684)
ids <- labels(protein_peak_dend$dendrogram)[micromonospora_boundaries[1]:micromonospora_boundaries[2]]
not_micro <- labels(protein_peak_dend$dendrogram)[!labels(protein_peak_dend$dendrogram) %in% ids]

dend_data <- function(dend,
                      ids,
                      pool){

  if (!all(ids %in% labels(dend))) {
    stop()
  }

  if (!all(labels(dend) %in% ids)){
    # get all labels not in ids
    to_prune <- labels(dend)[!labels(dend) %in% ids]
    # remove leaves not specified in ids
    dend <- dendextend::prune(dend,  to_prune)
}

  pb <- txtProgressBar(min = 0, max = length(ids), style = 3)
  message("retriveing small mol spectra")

    small_data <- lapply(seq_along(ids), function(x){
      setTxtProgressBar(pb, x)
      small_data <- IDBacApp::idbac_get_spectra(pool = pool,
                                                sampleIDs = ids[[x]],
                                                type = "small",
                                                MALDIquant = F)


      return(lapply(small_data, function(x){
        x[ , mass := as.integer(mass)]
        x <- x[ , list(intensity = mean(intensity))   , mass]
        x[mass > 200 ,  , ][mass < 2000 , , ]
      }))

    })
    close(pb)

    small_data <- unlist(small_data, recursive = FALSE)
    small_data <- split(small_data, names(small_data))
    small_data <- small_data[ match(labels(dend), names(small_data)) ]

    small_data <- lapply(small_data, function(x){
      #combine all replicates into one data.table
      x <- do.call(rbind, x)
      # mean intensity at unit mass resolution across all replicates
      x[ , list(intensity = mean(intensity))   , mass]
    })

    small_data  <- lapply(seq_along(small_data), function(x) {
      cbind(small_data[[x]],
            isolate = x)
    })

    small_datas <- do.call(rbind, small_data)




    pp <- ggplot(small_datas,
                 aes(mass,
                     as.factor(isolate))) +
      geom_raster(aes(fill = intensity),
                  interpolate = TRUE)+
      theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()) +
      ylab("Intensity") +
      xlab("m/z") +
      theme(axis.text.x = element_text(face = "italic")) +
      theme(legend.position = "none") +
      scale_fill_distiller(palette = "Greys",direction = 1,  trans="sqrt", limits=c(0,1))

    # prot --------------------------------------------------------------------
    pb <- txtProgressBar(min = 0, max = length(ids), style = 3)
    message("retriveing protein spectra")

    small_data <- lapply(seq_along(ids), function(x){
      setTxtProgressBar(pb, x)
      small_data <- IDBacApp::idbac_get_spectra(pool = pool,
                                                sampleIDs = ids[[x]],
                                                type = "protein",
                                                MALDIquant = F)
      return(lapply(small_data, function(x){
        x[ , mass := as.integer(mass)]
        x <- x[ , list(intensity = mean(intensity))   , mass]
        x[mass > 4000 ,  , ][mass < 15000 , , ]
      }))

    })
    close(pb)

    small_data <- unlist(small_data, recursive = FALSE)
    small_data <- split(small_data, names(small_data))
    small_data <- small_data[ match(labels(dend), names(small_data)) ]

    small_data <- lapply(small_data, function(x){
      #combine all replicates into one data.table
      x <- do.call(rbind, x)
      # mean intensity at unit mass resolution across all replicates
      x[ , list(intensity = mean(intensity))   , mass]
    })

    small_data  <- lapply(seq_along(small_data), function(x) {
      cbind(small_data[[x]],
            isolate = x)
    })

    small_data <- do.call(rbind, small_data)

    pp2 <- ggplot(small_data,
                  aes(mass,
                      as.factor(isolate))) +
      geom_raster(aes(fill = intensity),
                  interpolate = TRUE)+
      theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()) +
      ylab("Intensity") +
      xlab("m/z") +
      theme(axis.text.x = element_text(face = "italic")) +
      theme(legend.position = "none") +
      scale_fill_distiller(palette = "Greys",direction = 1,  trans="sqrt", limits=c(0,100))


    # Dendrogram plotting -----------------------------------------------------
    mets <- IDBacApp::idbac_get_metadata(strainID = labels(dend),
                                         metadataColumn = "sponge_identifier",
                                         pool = pool)
    colnames(mets) <- c("taxa", "sponge")
    # ggplot dendrogram portion
    p <- ggtree(ape::as.phylo(dend),
                ladderize = FALSE)

    p %<+% mets + geom_tippoint(aes(color=sponge), size=3)




    # get 16s genus information
    genera_16s <- IDBacApp::idbac_get_metadata(strainID = labels(dend),
                                               metadataColumn = "genus",
                                               pool = pool)[, 2]

    # combine 16s genera data with MALDI-TOF reference db matches
    genera_maldi <- search_result$genus[match(labels(dend),
                                              search_result$query)]
    genera_maldi[!is.na(genera_maldi)] <- paste0(genera_maldi[!is.na(genera_maldi)], "*")
    genera_maldi[is.na(genera_maldi)] <- ""
    genera <- genera_16s
    genera[is.na(genera)] <- ""
    genera <- paste0(genera, genera_maldi)

    # Add genera labels to plot
    p <- p +
      lapply(unique(genera),
             function(x){
               annotate("text",
                        x = .85,
                        y = which(genera == x),
                        label = x,
                        size = 3,
                        angle = 0,
                        fontface = 'italic')
             }) +
      coord_cartesian(xlim = c(0, 1))


    cowplot::plot_grid(p,  pp, pp2, ncol = 3, align = "h", axis = "b")
}
