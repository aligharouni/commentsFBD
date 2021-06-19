scripts/frechet_dist.html: scripts/frechet_dist.Rmd
	cd scripts; Rscript -e 'rmarkdown::render("frechet_dist.Rmd")'

curveBP.pdf: curveBP.tex cent_plot.tex
	texi2dvi -p curveBP.tex

cent_plot.tex: scripts/CurveBoxplot.R
	cd scripts; R CMD BATCH --vanilla CurveBoxplot.R
	mv scripts/cent_plot.tex scripts/cent_plot_ras1.png .

