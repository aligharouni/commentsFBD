## New submission to PeerJ Life & Environment
curveBP_peerJ.pdf: curveBP_peerJ.tex figs/cent_plot.pdf
	texi2dvi -p curveBP_peerJ.tex

figs/cent_plot.pdf: scripts/cent_plot.tex
## must compile in scripts/ dir because of aux files
	mkdir -p figs
	cd scripts; pdflatex --output-directory=../figs cent_plot.tex 

## Appendix plot, with L_2 norm result	(not sure where this is used?)
figs/cent_plot2.pdf: scripts/cent_plot2.tex
	cd scripts; pdflatex  --output-directory=../figs cent_plot2.tex

scripts/cent_plot.tex: scripts/CurveBoxplot.R
	cd scripts; R CMD BATCH --vanilla CurveBoxplot.R

scripts/cent_plot2.tex: scripts/CurveBoxplot.R
	cd scripts; R CMD BATCH --vanilla CurveBoxplot.R

### all this stuff is older; will need GH archaeology to build it

scripts/frechet_dist.html: scripts/frechet_dist.Rmd
	cd scripts; Rscript -e 'rmarkdown::render("frechet_dist.Rmd")'

curveBP.pdf: curveBP.tex scripts/cent_plot.pdf
	texi2dvi -p curveBP.tex
