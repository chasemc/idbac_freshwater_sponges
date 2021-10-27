
Note: This github repository contains none of the data needed to run these scripts, for that see: https://doi.org/10.5281/zenodo.5123348

 

# Reproducibility

TO reproduce the figures/analyses from this work, download the entire "sponge_manuscript" directory (warning- this is a few GBlarge)

A docker image is available from dockerhub, but also available in this repo:


To build the image from scratch open a terminal at `sponge_meta/docker` and run the following command (warning: this is not fast): 
```
docker build -f dockerfile --tag chasemc2/sponge_manuscript .
```


To build notebook open a terminal at `sponge_meta`:

``` 
docker run --user `id -u` --rm -ti -v${PWD}:/top -w/top/analyses chasemc2/sponge_manuscript Rscript -e 'rmarkdown::render("analysis.Rmd",  params = list(conditions = "all") ,  output_file = "all_isolates.html")' 
docker run --user `id -u` --rm -ti -v${PWD}:/top -w/top/analyses chasemc2/sponge_manuscript Rscript -e 'rmarkdown::render("analysis.Rmd", params = list(conditions = "matched"),  output_file = "matched_isolates.html")' 
docker run --user `id -u` --rm -ti -v${PWD}:/top -w/top/analyses chasemc2/sponge_manuscript Rscript -e 'rmarkdown::render("16s_analysis.Rmd")' 
```
