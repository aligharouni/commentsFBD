
scripts/frechet_dist.html: scripts/frechet_dist.Rmd
	cd scripts; Rscript -e 'rmarkdown::render("frechet_dist.Rmd")'

curveBP.pdf: curveBP.tex scripts/cent_plot.pdf
	texi2dvi -p curveBP.tex

scripts/cent_plot.pdf: scripts/cent_plot.tex
	cd scripts; pdflatex cent_plot.tex

scripts/cent_plot.tex: scripts/CurveBoxplot.R
	cd scripts; R CMD BATCH --vanilla CurveBoxplot.R
	
## There must a better way to do this, just for now;
## adding the L-2 norm to the cent_plot
scripts/cent_plot2.pdf: scripts/cent_plot2.tex
	cd scripts; pdflatex cent_plot2.tex

scripts/cent_plot2.tex: scripts/CurveBoxplot2.R
	cd scripts; R CMD BATCH --vanilla CurveBoxplot2.R

## New submission to PeerJ Life & Environment
curveBP_peerJ.pdf: curveBP_peerJ.tex scripts/cent_plot.pdf scripts/cent_plot2.pdf
	texi2dvi -p curveBP_peerJ.tex
