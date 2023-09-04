pkgs <- c("tidyverse","latex2exp","fda","roahd","directlabels",
          "colorspace","ggrepel","Cairo","ggrastr","tikzDevice")
i1 <- installed.packages()
pkgs <- setdiff(pkgs, rownames(i1))
install.packages(pkgs, repos = "https://cloud.r-project.org")
