
bob_the_two_plot_builder <- function(pool = pool, input_a, input_b, protein_dist, dend){

  spectra <- idbac_get_spectra(pool = pool,
                               sampleIDs = c(input_a, input_b),
                               type = "protein",
                               MALDIquant = T)

  spectra <- MALDIquant::averageMassSpectra(spectra,
                                            labels = names(spectra))

  #Ensure expected order
  spectra <- spectra[order(c(input_a, input_b))]


  biggy <- rbind(data.table(mass = spectra[[1]]@mass,
                            intensity = spectra[[1]]@intensity,
                            Isolate = input_a),
                 data.table(mass = spectra[[2]]@mass,
                            intensity = -spectra[[2]]@intensity,
                            Isolate = input_b))

  biggy <- biggy[order(mass), ]


  spectra <- idbac_get_spectra(pool = pool,
                               sampleIDs = c(input_a, input_b),
                               type = "small",
                               MALDIquant = T)

  spectra <- MALDIquant::averageMassSpectra(spectra,
                                            labels = names(spectra))

  #Ensure expected order
  spectra <- spectra[order(c(input_a, input_b))]

  smallz <- rbind(data.table(mass = spectra[[1]]@mass,
                             intensity = spectra[[1]]@intensity,
                             Isolate = input_a),
                  data.table(mass = spectra[[2]]@mass,
                             intensity = spectra[[2]]@intensity,
                             Isolate = input_b))

  smallz <- smallz[order(mass), ]

  # Square root transform spectrum to reduce differences in intensities
  smallz <- smallz[ , intensity := sqrt(intensity)]

  # Flip input_b spectrum into negative intensity
  smallz <- smallz[Isolate == input_b, intensity := -intensity]


  lca <- lca_height(dend = dend,
                    input_a = input_a,
                    input_b = input_b)


  lca <- paste0("Join height on dendrogram: ",
                round(lca, digits = 2))

  a <- ggplot(biggy) +
    geom_line(aes(x = mass,
                  y = intensity,
                  color = Isolate)) +
    scale_color_manual(values=c("#000000", "#afafaf")) +
    theme_minimal() +
    labs(main="",
         x = "Daltons",
         y = "Intensity") +
    #annotate("text", x = 15000, y = 75, label = lca) +
    theme(axis.title.x = element_blank()) +
    scale_x_continuous(breaks=seq(2000, 20000, 3000))


  b <- ggplot(smallz) +
    geom_line(aes(x = mass,
                  y = intensity,
                  color = Isolate)) +
    scale_color_manual(values=c("#000000", "#afafaf")) +
    theme_minimal() +
    labs(main="",
         x = "Daltons",
         y = "Intensity (square root)") +
    scale_x_continuous(breaks=seq(0, 3000, 500))



  p <- ggpubr::ggarrange(create_dend_plot(input_a = input_a,
                                          input_b = input_b,
                                          dend = dend),
                         ggpubr::ggarrange(
                           a,
                           b,
                           ncol = 1,
                           nrow = 2, labels = c("b", "c")), ncol = 2, nrow = 1, widths = c(1.5,2), labels="a")
  annotate_figure(p,
                  top = ggpubr::text_grob(paste0("Comparing Isolates ",
                                                 input_a,
                                                 " and ",
                                                 input_b),
                                          color = "black",
                                          face = "bold",
                                          size = 14)
  )

}




create_dend_plot <- function(input_a, input_b, dend){

  subtrees <- dendextend::partition_leaves(dend)

  find_ancestors <- function(input_labels) {
    which(sapply(subtrees, function(x) input_labels %in% x))
  }

  ancestors <- lapply(c(input_a, input_b), find_ancestors)
  # find largest common


  dend <- ape::as.phylo(dend)

  internal_node <- ape::getMRCA(dend, c(input_a, input_b))


  gg_data <- ggtree::groupClade(dend, internal_node)
  ###
  ids <- c(input_a, input_b)
  new_labs <- sapply(ids,
                 function(x){
                   grep(x, labels(gg_data))
                 },
                 USE.NAMES = F
  )

  labels(gg_data) = rep("", length(labels(gg_data)))

  labels(gg_data)[new_labs] <- ids

  ###


  p <- ggtree(gg_data,
              aes(color = group), ladderize = F) +
    theme(legend.position = 'none') +
    scale_color_manual(values=c("black", "red"))  + geom_tiplab(size=2.5)


  #p2 <- viewClade(p, internal_node) + geom_tiplab(size=2.5) + coord_cartesian(clip="off")

  cowplot::plot_grid(p)

}








lca_height <- function(dend,
                       input_a,
                       input_b) {

  subtrees <- dendextend::partition_leaves(dend)

  find_ancestors <- function(input_labels) {
    which(sapply(subtrees, function(x) input_labels %in% x))
  }

  ancestors <- lapply(c(input_a, input_b), find_ancestors)
  # find largest common
  lca <- ancestors[[1]][max(which(ancestors[[1]] %in% ancestors[[2]]))]

  # which column represents "evolutionary" height
  index <-  which(sapply(as.data.frame(dendextend::get_nodes_xy(dend)), max) != attributes(dend)$members)[[1]]
  dendextend::get_nodes_xy(dend)[lca, index]
}


































bob_the_mini_plot_builder <- function(input_a, input_b){

    spectra <- idbac_get_spectra(pool = pool,
                                 sampleIDs = c(input_a, input_b),
                                 type = "protein",
                                 MALDIquant = T)

  spectra <- MALDIquant::averageMassSpectra(spectra,
                                            labels = names(spectra))

  #Ensure expected order
  spectra <- spectra[order(c(input_a, input_b))]

  one <- data.table(mass = spectra[[1]]@mass,
                    intensity = spectra[[1]]@intensity,
                    Isolate = input_a)
  one[, mass := round(mass, 0)]
  one <- one[ , list(intensity = sum(intensity), Isolate = unique(Isolate)), mass]
  two <- data.table(mass = spectra[[2]]@mass,
                    intensity = -spectra[[2]]@intensity,
                    Isolate = input_b)
  two[, mass := round(mass, 0)]
  two <- two[ , list(intensity = sum(intensity), Isolate = unique(Isolate)), mass]


  biggy <- rbind(one,
                 two)

  biggy <- biggy[order(mass), ]


  spectra <- idbac_get_spectra(pool = pool,
                               sampleIDs = c(input_a, input_b),
                               type = "small",
                               MALDIquant = T)

  spectra <- MALDIquant::averageMassSpectra(spectra,
                                            labels = names(spectra))

  #Ensure expected order
  spectra <- spectra[order(c(input_a, input_b))]
  one <- data.table(mass = spectra[[1]]@mass,
                    intensity = spectra[[1]]@intensity,
                    Isolate = input_a)
  one[, mass := round(mass, 0)]
  one <- one[ , list(intensity = sum(intensity), Isolate = unique(Isolate)), mass]
  two <- data.table(mass = spectra[[2]]@mass,
                    intensity = spectra[[2]]@intensity,
                    Isolate = input_b)
  two[, mass := round(mass, 0)]
  two <- two[ , list(intensity = sum(intensity), Isolate = unique(Isolate)), mass]


  smallz <- rbind(one, two)

  smallz <- smallz[order(mass), ]

  # Square root transform spectrum to reduce differences in intensities
  smallz <- smallz[ , intensity := sqrt(intensity)]

  # Flip input_b spectrum into negative intensity
  smallz <- smallz[Isolate == input_b, intensity := -intensity]



  a <- ggplot(biggy) +
    geom_line(aes(x = mass,
                  y = intensity,
                  color = Isolate), size = .25) +
    scale_color_manual(values=c("#000000", "#afafaf")) +
    theme_minimal() +
    labs(main="",
         x = "Daltons",
         y = "Intensity") +
    #annotate("text", x = 15000, y = 75, label = lca) +
    theme(axis.title.x = element_blank()) +
    scale_x_continuous(breaks=seq(2000, 20000, 3000))


  b <- ggplot(smallz) +
    geom_line(aes(x = mass,
                  y = intensity,
                  color = Isolate), size = .25) +
    scale_color_manual(values=c("#000000", "#afafaf")) +
    theme_minimal() +
    labs(main="",
         x = "Daltons",
         y = "Intensity (square root)") +
    scale_x_continuous(breaks=seq(0, 3000, 500))



  p <- ggpubr::ggarrange(a,
                         b,
                         ncol = 1,
                         nrow = 2)

  annotate_figure(p,
                  top = ggpubr::text_grob(paste0("Comparing ",
                                                 input_a,
                                                 " and ",
                                                 input_b),
                                          color = "black",
                                          face = "bold",
                                          size = 14)
  )

}


