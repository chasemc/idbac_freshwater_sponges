
cut_analysis <- function(input_dist, input_dend, start_seq, stop_seq, by_seq){

  # Loop through all cut heights
  ww <- lapply(seq(start_seq, stop_seq, by_seq), function(y){
    # Convert dist to matrix
    input_dist <- as.matrix(input_dist)
    # cut tree
    h <- cutree(input_dend, h = y)
    # Loop through all groups
    w <- lapply(1:max(h),
                function(x){

                  a <- which(h == x)
                  a <- input_dist[colnames(input_dist) %in% names(a), colnames(input_dist) %in% names(a)]
                  a[lower.tri(a, diag = F)] <- NA

                  if (is.null(nrow(a)) || nrow(a) == 0) {
                    # if only one sample in a group...
                    a <- cbind.data.frame(group = x,
                                          Var1 = names(h)[[x]],
                                          Var2 = names(h)[[x]],
                                          value = 1)
                  } else{
                    # if > 1 sample in a group
                    a <- suppressWarnings(cbind(group = x,
                                                melt(a,
                                                     fill = T)))
                  }
                  as.data.table(a)
                })

    cbind(iter = y,
          do.call(rbind,
                  w))
  })

  ww <- do.call(rbind, ww)

  # remove upper triangle of distance matrix
  ww <- ww[!is.na(value), ]


  dplyr::tbl(pool, "metadata") %>%
    dplyr::select(c("strain_id", "sponge_identifier")) %>%
    collect() %>%
    return(.) -> sample_meta


  ww <- merge(ww, sample_meta, by.x = "Var1", by.y = "strain_id")
  colnames(ww)[colnames(ww) == "sponge_identifier"] <- "var1_sponge"
  ww <- merge(ww, sample_meta, by.x = "Var2", by.y = "strain_id")
  colnames(ww)[colnames(ww) == "sponge_identifier"] <- "var2_sponge"
  ww[, mixed := var1_sponge != var2_sponge ]

  ww

}
