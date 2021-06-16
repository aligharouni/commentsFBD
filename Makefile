scripts/frechet_dist.html: scripts/frechet_dist.Rmd
	cd scripts; Rscript -e 'rmarkdown::render("frechet_dist.Rmd")'

curveBP.pdf: curveBP.tex scripts/cent_plot.pdf
	texi2dvi -p curveBP.tex

scripts/cent_plot.pdf: scripts/CurveBoxplot.R
	cd scripts; R CMD BATCH --vanilla CurveBoxplot.R
