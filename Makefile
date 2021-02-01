scripts/frechet_dist.html: scripts/frechet_dist.Rmd
	cd scripts; Rscript -e 'rmarkdown::render("frechet_dist.Rmd")'
