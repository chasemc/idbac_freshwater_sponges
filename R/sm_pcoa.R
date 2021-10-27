
sm_pcoa <- function(ids,
                    small_peak_data,
                    pool) {

  sm <- IDBacApp::createFuzzyVector(massStart = 200,
                                    massEnd = 3000,
                                    massList = lapply(small_peak_data, function(x) x@mass),
                                    intensityList = lapply(small_peak_data, function(x) x@intensity),
                                    ppm = 300)

  sm <- sm[, labs %in% ids]

  labs <- colnames(sm)

  pcoa_result <- ape::pcoa(as.dist(1 - coop::cosine(sm)))

  pco <- cbind.data.frame(labs, pcoa_result$vectors[ , 1:2],
                          stringsAsFactors = FALSE,
                          IDBacApp::idbac_get_metadata(strainID =labs,
                                                       metadataColumn = "sponge_identifier",
                                                       pool = pool)[, 2])



  colnames(pco) <- c("isolate", "PCoA1", "PCoA2", "sponge")


  ggplot(data = pco,
         aes(x =  PCoA1,
             y = PCoA2,
             shape = sponge,
             color = sponge)) +
    geom_point(alpha = 0.7, size = 2.5) +
    # stat_ellipse(type = 't',
    #              size = .5,
    #              show.legend = F,
    #              alpha = .7) +
    xlab(paste0("PCoA1 (",
                round(pcoa_result$values$Relative_eig[[1]] * 100, 2),
                "%)")) +
    ylab(paste0("PCoA2 (",
                round(pcoa_result$values$Relative_eig[[2]] * 100, 2),
                "%)")) +
    scale_color_manual(values = c("#1b9e77",
                                  "#7570b3"))



}
