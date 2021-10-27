

recolor_mans2 <- function(pool,
                          tsv_path,
                          json_path,
                          save_path,
                          save_name,
                          png,
                          dend = NULL,
                          k = NULL,
                          meta_column,
                          display_labels = FALSE){

  IDBacApp:::.checkPool(pool = pool)

  # MAN edge table (source, target, weight)
  # with added column for edge length in python
  z <- data.table::fread(tsv_path)

  if (is.null(k)) {

    # get metadata
    b <- IDBacApp::idbac_get_metadata(strainID = z$source,
                                      metadataColumn = meta_column,
                                      pool = pool)
    meta <- unique(b[, 2])
    # convert selected meta column to integer representation
    b[, 2] <- as.integer(factor(b[, 2], levels = sort(unique(meta))))

    # bind integer rep to MAN edge table
    z <- cbind(z,
               group = b[, 2])
  } else {
    meta <- 1:k

    selected_names <- unique(z$source)
    # remove all non-selected names from dendrogram
    dend <- dendextend::prune(dend, labels(dend)[!labels(dend) %in% selected_names])
    # cut protein dendrogram
    # creates named integer vector, integers represent groups
    b <- cutree(dend, k = k)

    # bind integer rep to MAN edge table
    z <- cbind(z,
               group = as.numeric(b)[match(z$source,
                                           names(b))])
  }
  # Define color palette, IDBac's starts colorblind friendly
  cols <- RColorBrewer::brewer.pal(7, "Set1")


  con= file("recolor_mans2", open = "a")
  writeLines(paste0(save_name, ": ", max(dend$height)), con = con)
  close(con)

  # setup plotly plot + legend
  fig <- plotly::plotly_empty(data = data.frame(x = rep(2, length(unique(meta))),
                                                y = seq(1,
                                                        length(unique(meta)),
                                                        0.25)[1:length(unique(meta))],
                                                cc = as.factor(sort(unique(meta)))),
                              x = ~x,
                              y = ~y,
                              type = 'scatter',
                              mode = 'text',
                              text = as.factor(sort(unique(meta))),
                              color = ~cc,
                              colors = cols[1:length(unique(meta))]) %>%
    layout(xaxis = list(title = '',
                        range = c(0, 4)),
           yaxis = list(title = '',
                        range = c(0,4))) %>%
    layout(showlegend = FALSE) %>%
  layout(plot_bgcolor='transparent') %>%
    layout(paper_bgcolor='transparent') %>% config(displayModeBar = F)



  # read plotly json
  a <- readr::read_file(json_path)
  a <- jsonlite::parse_json(a,
                            simplifyVector = TRUE,
                            simplifyDataFrame = FALSE,
                            simplifyMatrix = FALSE)
  # make all nodes in the network grey
  a$data[[2]]$marker$color[a$data[[2]]$marker$color == "red"] <- "grey"

  #get which nodes are sampels (larger size)
  # not text match in case someone has sample-id same as m/z
  s_mark <- which(a$data[[2]]$marker$size == max(a$data[[2]]$marker$size))


  grps <- z[match(a$data[[2]]$text[s_mark], z$source), group]

  a$data[[2]]$marker$color[s_mark] <- cols[grps]
  a$data[[2]]$marker$size <- a$data[[2]]$marker$size /2
  #removes all labels
  if (!display_labels) {
    a$data[[2]]$text <- ""
  }

  p <- plotly::plotly_build(a) %>%
    layout(annotations = list(x = .5,
                              y = 1.05,
                              text =  sprintf("<i>%s</i>", save_name),
                              showarrow = F,
                              xref = 'paper',
                              yref = 'paper'))

  p <- subplot(p,
               fig,
               widths = c(0.85,
                          0.15)) %>% layout(plot_bgcolor='transparent') %>%
    layout(paper_bgcolor='transparent') %>% config(displayModeBar = F)
  return(p)
}






recolor_mans_dend <- function(pool,
                          tsv_path,
                          json_path,
                          save_path,
                          save_name,
                          png,
                          dend = NULL,
                          k = NULL,
                          meta_column,
                          display_labels = FALSE){



  IDBacApp:::.checkPool(pool = pool)

  # MAN edge table (source, target, weight)
  # with added column for edge length in python
  z <- data.table::fread(tsv_path)

  meta <- 1:k

  selected_names <- unique(z$source)
  # remove all non-selected names from dendrogram
  dend <- dendextend::prune(dend, labels(dend)[!labels(dend) %in% selected_names])
  # cut protein dendrogram
  # creates named integer vector, integers represent groups
  b <- cutree(dend, k = k)

  # bind integer rep to MAN edge table
  z <- cbind(z,
             group = as.numeric(b)[match(z$source,
                                         names(b))])



  # Define color palette, IDBac's starts colorblind friendly
  cols <- RColorBrewer::brewer.pal(7, "Set1")[1:length(unique(meta))]

  # Change color of dendrogram
  dendextend::color_branches(dend, clusters = b[match(labels(dend), names(b))],
                             col = cols[unique(b[match(labels(dend), names(b))])]) %>%
    dendextend::set_labels(ifelse(display_labels,labels(dend), "")) %>%
    ggplot() %>%
    ggplotly() %>%
    hide_legend() %>%
    return() -> dend_plot

  con= file("recolor_mans_dend", open = "a")
  writeLines(paste0(save_name, ": ", max(dend$height)), con = con)
  close(con)
  saveRDS(dend, paste0("recolor_mans_dend_", save_name, ".rds"))

# setup plotly plot + legend
fig <- plotly::plotly_empty(data = data.frame(x = rep(2, length(unique(meta))),
                                              y = seq(1,
                                                      length(unique(meta)),
                                                      0.25)[1:length(unique(meta))],
                                              cc = as.factor(sort(unique(meta)))),
                            x = ~x,
                            y = ~y,
                            type = 'scatter',
                            hoverinfo = 'text',
                            mode = 'text',
                            text = as.factor(sort(unique(meta))),
                            color = ~cc,
                            colors = cols) %>%
  layout(xaxis = list(title = '',
                      range = c(0, 4)),
         yaxis = list(title = '',
                      range = c(0,4))) %>%
  layout(showlegend = FALSE) %>%
  layout(plot_bgcolor = 'rgb(250, 250, 250)') %>%
  layout(paper_bgcolor = 'rgb(250, 250, 250)')



# read plotly json
a <- readr::read_file(json_path)
a <- jsonlite::parse_json(a,
                          simplifyVector = TRUE,
                          simplifyDataFrame = FALSE,
                          simplifyMatrix = FALSE)
# make all nodes in the network grey
a$data[[2]]$marker$color[a$data[[2]]$marker$color == "red"] <- "grey"

#get which nodes are samples (larger size)
# not text match in case someone has sample-id same as m/z
s_mark <- which(a$data[[2]]$marker$size == max(a$data[[2]]$marker$size))


grps <- z[match(a$data[[2]]$text[s_mark], z$source), group]

a$data[[2]]$marker$color[s_mark] <- cols[grps]

#change marker size
a$data[[2]]$marker$size <- a$data[[2]]$marker$size

# plot large nodes on top layer
temp <- a$data[[2]]
a$data[[3]] <- temp

large_node_index <- which(temp$marker$size == max(temp$marker$size))
a$data[[3]]$text <- a$data[[3]]$text[large_node_index]
a$data[[3]]$x <- a$data[[3]]$x[large_node_index]
a$data[[3]]$y <- a$data[[3]]$y[large_node_index]
a$data[[3]]$marker$color <- a$data[[3]]$marker$color[large_node_index]
a$data[[3]]$marker$size <- a$data[[3]]$marker$size[large_node_index]


lil_node_index <- which(temp$marker$size == min(temp$marker$size))
a$data[[2]]$text <- a$data[[2]]$text[lil_node_index]
a$data[[2]]$x <- a$data[[2]]$x[lil_node_index]
a$data[[2]]$y <- a$data[[2]]$y[lil_node_index]
a$data[[2]]$marker$color <- a$data[[2]]$marker$color[lil_node_index]
a$data[[2]]$marker$size <- a$data[[2]]$marker$size[lil_node_index]



p <- plotly::plotly_build(a) %>%
  layout(annotations = list(x = .5,
                            y = 1.05,
                            text =  sprintf("<i>%s</i>", save_name),
                            showarrow = F,
                            xref = 'paper',
                            yref = 'paper'))


#removes small node labels
  a$data[[2]]$text <- ""

#removes large node labels
if (!display_labels) {
  a$data[[3]]$text <- ""
}



p <- plotly::plotly_build(a) %>%
  layout(annotations = list(x = .5,
                            y = 1.05,
                            text =  sprintf("<i>%s</i>", save_name),
                            showarrow = F,
                            xref = 'paper',
                            yref = 'paper'))

subplot(p, dend_plot) %>% layout(plot_bgcolor='transparent') %>%
  layout(paper_bgcolor='transparent') %>% config(displayModeBar = F)


}






recolor_mans_sponge_dend <- function(pool,
                          tsv_path,
                          json_path,
                          save_path,
                          save_name,
                          png,
                          dend = NULL,
                          k = NULL,
                          meta_column,
                          display_labels = FALSE){



  IDBacApp:::.checkPool(pool = pool)

  # MAN edge table (source, target, weight)
  # with added column for edge length in python
  z <- data.table::fread(tsv_path)

  # get metadata
  b <- IDBacApp::idbac_get_metadata(strainID = z$source,
                                    metadataColumn = "sponge_identifier",
                                    pool = pool)
  meta <- unique(b[, 2])
  # convert selected meta column to integer representation
  b[, 2] <- as.integer(as.factor(b[, 2]))

  # bind integer rep to MAN edge table
  z <- cbind(z,
             group = b[, 2])





  meta <- 1:k

  selected_names <- unique(z$source)
  # remove all non-selected names from dendrogram
  dend <- dendextend::prune(dend, labels(dend)[!labels(dend) %in% selected_names])
  # cut protein dendrogram
  # creates named integer vector, integers represent groups





  # Define color palette, IDBac's starts colorblind friendly
  cols <- RColorBrewer::brewer.pal(7, "Set1")[1:length(unique(meta))]

  # Change color of dendrogram
  dendextend::color_branches(dend, clusters = b$sponge_identifier[match(labels(dend), z$source)],
                             col = cols[unique(b$sponge_identifier[match(labels(dend), z$source)])]) %>%
    dendextend::set_labels("") %>%
    ggplot() %>%
    ggplotly() %>%
    hide_legend() %>%
    return() -> dend_plot


  con= file("recolor_mans_sponge_dend", open = "a")
  writeLines(paste0(save_name, ": ", max(dend$height)), con = con)
  close(con)



  # setup plotly plot + legend
  fig <- plotly::plotly_empty(data = data.frame(x = rep(2, length(unique(meta))),
                                                y = seq(1,
                                                        length(unique(meta)),
                                                        0.25)[1:length(unique(meta))],
                                                cc = as.factor(sort(unique(meta)))),
                              x = ~x,
                              y = ~y,
                              type = 'scatter',
                              mode = 'text',
                              text = as.factor(sort(unique(meta))),
                              color = ~cc,
                              colors = cols) %>%
    layout(xaxis = list(title = '',
                        range = c(0, 4)),
           yaxis = list(title = '',
                        range = c(0,4))) %>%
    layout(showlegend = FALSE) %>%
    layout(plot_bgcolor='transparent') %>%
    layout(paper_bgcolor='transparent') %>% config(displayModeBar = F)



  # read plotly json
  a <- readr::read_file(json_path)
  a <- jsonlite::parse_json(a,
                            simplifyVector = TRUE,
                            simplifyDataFrame = FALSE,
                            simplifyMatrix = FALSE)
  # make all nodes in the network grey
  a$data[[2]]$marker$color[a$data[[2]]$marker$color == "red"] <- "grey"

  #get which nodes are samples (larger size)
  # not text match in case someone has sample-id same as m/z
  s_mark <- which(a$data[[2]]$marker$size == max(a$data[[2]]$marker$size))


  grps <- z[match(a$data[[2]]$text[s_mark], z$source), group]

  a$data[[2]]$marker$color[s_mark] <- cols[grps]

  #change marker size
  a$data[[2]]$marker$size <- a$data[[2]]$marker$size

  # plot large nodes on top layer
  temp <- a$data[[2]]
  a$data[[3]] <- temp

  large_node_index <- which(temp$marker$size == max(temp$marker$size))
  a$data[[3]]$text <- a$data[[3]]$text[large_node_index]
  a$data[[3]]$x <- a$data[[3]]$x[large_node_index]
  a$data[[3]]$y <- a$data[[3]]$y[large_node_index]
  a$data[[3]]$marker$color <- a$data[[3]]$marker$color[large_node_index]
  a$data[[3]]$marker$size <- a$data[[3]]$marker$size[large_node_index]


  lil_node_index <- which(temp$marker$size == min(temp$marker$size))
  a$data[[2]]$text <- a$data[[2]]$text[lil_node_index]
  a$data[[2]]$x <- a$data[[2]]$x[lil_node_index]
  a$data[[2]]$y <- a$data[[2]]$y[lil_node_index]
  a$data[[2]]$marker$color <- a$data[[2]]$marker$color[lil_node_index]
  a$data[[2]]$marker$size <- a$data[[2]]$marker$size[lil_node_index]




  #removes small node labels
  if (!display_labels) {
    a$data[[2]]$text <- ""
  }
  #removes large node labels
  if (!display_labels) {
    a$data[[3]]$text <- ""
  }



  p <- plotly::plotly_build(a)

  subplot(p, dend_plot) %>% layout(plot_bgcolor='transparent') %>%
    layout(paper_bgcolor='transparent') %>% config(displayModeBar = F)


}


