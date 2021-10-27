library(IDBacApp)

data_directory <- here::here("data",
                             "raw_data",
                             "two_sponge")
a <- list.dirs(data_directory,
               recursive = F,
               full.names = T)

file_path <- here::here("data",
                        "sqlite")
temp_dir <- tempdir()
temp_dir <- file.path(temp_dir,
                      "mzml")
dir.create(temp_dir)

IDBacApp::idbac_create(fileName = "two_sponge",
                       filePath = file_path)

for (i in seq_along(a)){
  print(i)
  sampleMap <- file.path(a[[i]], "layout.tsv")
  sampleMap <- read.delim(sampleMap, header = F, stringsAsFactors = F)

  temp_dim <- dim(sampleMap)

  #fill empty columns
  if (temp_dim[[2]] < 24) {
    sampleMap <- cbind(sampleMap,
                       do.call(cbind,
                               lapply(1:(24 - temp_dim[[2]]), function(x) rep(NA, nrow(sampleMap)))))
  }
  #fill empty rows
  if (temp_dim[[1]] < 16) {
    sampleMap <- cbind(sampleMap, do.call(rbind, lapply(1:(16 - temp_dim[[1]]), function(x) rep(NA, ncol(sampleMap)))))
  }


  IDBacApp::db_from_bruker(dataDirectory = a[[i]],
                           fileName = "two_sponge",
                           filePath = file_path,
                           sampleMap = sampleMap,
                           tempDir = temp_dir)
  file.remove(list.files(temp_dir))


}





# Associate strain id with location on plate
b <- read.delim(here::here("data/raw_data/two_sponge/2018_05_23/names_association.tsv"), sep = "\t", header = F)
d <- IDBacApp:::map384Well()
d <- as.data.frame(d[1:6, 1:8])
b <- unlist(b)
names(b) <- as.character(unlist(d))

for(i in seq_along(b)) {

  pool::poolWithTransaction(two_sponge_pool, function(conn){
    statement <- DBI::dbSendStatement(conn, 'UPDATE spectra
                                      SET strain_id = $new_id
                                      WHERE strain_id = $old_id;')

    statement <- DBI::dbBind(statement,
                             list(new_id = b[[i]],
                                  old_id = names(b[i])))
    res = DBI::dbGetRowsAffected(statement)
    DBI::dbClearResult(statement)
    return(res)
  })

}



for(i in seq_along(b)) {

  pool::poolWithTransaction(two_sponge_pool, function(conn){
    statement <- DBI::dbSendStatement(conn, 'UPDATE metadata
                                      SET strain_id = $new_id
                                      WHERE strain_id = $old_id;')

    statement <- DBI::dbBind(statement,
                             list(new_id = b[[i]],
                                  old_id = names(b[i])))
    res = DBI::dbGetRowsAffected(statement)
    DBI::dbClearResult(statement)
    return(res)
  })

}



# Delete samples that were run on the same day as sponge samples, but aren't part of this study
# Also rename part of a plate from 2018-05-23


# There's no function in IDBAc, yet, to delete samples from the database
# This was intentional, as the GUI allows people to create new DB with subsets of data
# I am considering adding the functionality, mostly for use on the commandline.


to_delete <- c(
  'ESP1-2; EDP1', # removed because part of a different experiment
  'ESP1-EDP2', # removed because part of a different experiment
  'ESP10-2; EDP1', # removed because part of a different experiment
  'ESP10-EDP2', # removed because part of a different experiment
  'ESP11-2; EDP1', # removed because part of a different experiment
  'ESP11-EDP2', # removed because part of a different experiment
  'ESP12-2; EDP1', # removed because part of a different experiment
  'ESP12-EDP2', # removed because part of a different experiment
  'ESP13-2; EDP1', # removed because part of a different experiment
  'ESP13-EDP2', # removed because part of a different experiment
  'ESP14-EDP2', # removed because part of a different experiment
  'ESP15-2; EDP1', # removed because part of a different experiment
  'ESP15-EDP2', # removed because part of a different experiment
  'ESP16-2; EDP1', # removed because part of a different experiment
  'ESP17-2; EDP1', # removed because part of a different experiment
  'ESP18-2; EDP1', # removed because part of a different experiment
  'ESP19-2; EDP1', # removed because part of a different experiment
  'ESP2-2; EDP1', # removed because part of a different experiment
  'ESP2-EDP2', # removed because part of a different experiment
  'ESP20-EDP1', # removed because part of a different experiment
  'ESP21-EDP1', # removed because part of a different experiment
  'ESP22-EDP1', # removed because part of a different experiment
  'ESP24-EDP1', # removed because part of a different experiment
  'ESP25-EDP1', # removed because part of a different experiment
  'ESP26-EDP1', # removed because part of a different experiment
  'ESP27-EDP1', # removed because part of a different experiment
  'ESP28-EDP1', # removed because part of a different experiment
  'ESP29-EDP1', # removed because part of a different experiment
  'ESP3-2; EDP1', # removed because part of a different experiment
  'ESP3-EDP2', # removed because part of a different experiment
  'ESP30-EDP1', # removed because part of a different experiment
  'ESP33-EDP1', # removed because part of a different experiment
  'ESP34-EDP1', # removed because part of a different experiment
  'ESP35-EDP1', # removed because part of a different experiment
  'ESP36-EDP1', # removed because part of a different experiment
  'ESP37-EDP1', # removed because part of a different experiment
  'ESP38-EDP1', # removed because part of a different experiment
  'ESP39-EDP1', # removed because part of a different experiment
  'ESP4-2; EDP1', # removed because part of a different experiment
  'ESP4-EDP2', # removed because part of a different experiment
  'ESP40-EDP1', # removed because part of a different experiment
  'ESP41-EDP1', # removed because part of a different experiment
  'ESP43-EDP1', # removed because part of a different experiment
  'ESP44-EDP1', # removed because part of a different experiment
  'ESP5-2; EDP1', # removed because part of a different experiment
  'ESP5-EDP2', # removed because part of a different experiment
  'ESP6-EDP2', # removed because part of a different experiment
  'ESP7-2; EDP1', # removed because part of a different experiment
  'ESP7-EDP2', # removed because part of a different experiment
  'ESP8-2; EDP1', # removed because part of a different experiment
  'ESP8-EDP2', # removed because part of a different experiment
  'ESP9-2; EDP1', # removed because part of a different experiment
  'ESP9-EDP2', # removed because part of a different experiment
  'ESP23', # removed because part of a different experiment
  'M', # removed blanks and controls
  'P', # removed blanks and controls
  'PEP', # removed blanks and controls
  'BTS', # removed blanks and controls
  'EXTRACTS', # removed because part of a different experiment
  '', # removed names missing in sample spreadsheets (other experiments)
  '996', # missing metadata
  '998', # missing metadata
  '994', # missing metadata
  '993', # missing metadata
  '995', # only two replicates
  "891", # only two replicates
  "958", # only two replicates
  "993") # only two replicates


two_sponge_pool <- idbac_connect(fileName = "two_sponge",
                                 filePath = file_path)[[1]]

pool::poolWithTransaction(two_sponge_pool, function(conn){
  statement <- DBI::dbSendStatement(conn, 'DELETE FROM spectra WHERE strain_id=$sample_id')
  statement <- DBI::dbBind(statement, list(sample_id = to_delete))
  bb=DBI::dbGetRowsAffected(statement)
  DBI::dbClearResult(statement)
  return(bb)
})


pool::poolWithTransaction(two_sponge_pool, function(conn){
  statement <- DBI::dbSendStatement(conn, 'DELETE FROM metadata WHERE strain_id=$sample_id')
  statement <- DBI::dbBind(statement, list(sample_id = to_delete))
  bb=DBI::dbGetRowsAffected(statement)
  DBI::dbClearResult(statement)
  return(bb)
})

