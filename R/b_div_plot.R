
b_div_plot <- function(p){

  ggpubr::ggarrange(p[[1]],
                    p[[2]],
                    p[[3]],
                    p[[4]], align = "hv")

}




b_div <- function(meta,
                  meta_names,
                  p_dist_temp){

  lapply(meta_names,
         function(x){
           vegan::betadisper(p_dist_temp, as.factor(meta[, get(x)]))
         })

}


b_div_plot_calc <- function(meta,
                            meta_names,
                            p_dist_temp){


  beta_div <- b_div(meta,
                    meta_names,
                    p_dist_temp)

  lapply(seq_along(beta_div),
         function(x){

           b <- cbind.data.frame(x = beta_div[[x]]$vectors[, 1],
                                 y = beta_div[[x]]$vectors[, 2],
                                 color = beta_div[[x]]$group)

           perc1 <- round((beta_div[[x]]$eig[[1]] / sum(beta_div[[x]]$eig)) * 100, 1)
           perc2 <- round((beta_div[[x]]$eig[[2]] / sum(beta_div[[x]]$eig)) * 100, 1)


           legend_name <- meta_names[[x]]
           if (legend_name == "sponge_identifier") legend_name <- "sponge"

           ggplot(data = b,
                  aes(x = x,
                      y = y,
                      shape = color,
                      color = color)) +
             geom_point(alpha = 0.7, size = 2.5) +
             stat_ellipse(type = 't',
                          size = .5,
                          show.legend = F,
                          alpha = .7) +
             xlab(paste0("PCoA1 (", perc1, "%)")) +
             ylab(paste0("PCoA2 (", perc2, "%)")) +
             scale_color_manual(values = c("#1b9e77",
                                           "#7570b3"))

         })


}
