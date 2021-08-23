library(tidyverse); theme_set(theme_bw())

ensemble_J <- read_csv("./data/juul1.csv", col_names = FALSE) ## the columns are the trajs
long_ensemble <- ensemble_J %>% as.matrix() %>% reshape2::melt() %>%
  rename(tvec="Var1",grp="Var2") %>% mutate(across(tvec, ~.-1))

####
envdat <- readRDS("envdat.rds")

envdat <- (envdat
    %>% mutate(across(method, ~ifelse(. == "fda_roahd", "fda_J2_all", .)))
    %>% separate(method, into = c("pkg", "J", "nsamp"))
)

with(envdat, table(pkg, nsamp, J))

## from Juul_work.rmd
## want J to vary *slowest* but come first in the label
nmvec <- expand.grid(e=c("lwr","upr"),
                     q=c(50,90),
                     J=c(2, 5, 10, 20, 30, 40, 50)) %>%
    apply(1, function(x) paste(rev(x), collapse="_")) %>%
    trimws()
juul1 <- read_csv("data/juul_boundary1.csv", col_names = FALSE) %>%
    setNames(nmvec) %>%
    mutate(data.tvec = 0:(n()-1)) %>%
    pivot_longer(-data.tvec) %>%
    ## make it look like envdat
    separate(name, into = c("J","quantile","bound")) %>%
    filter(bound=="upr", quantile=="90") %>%
    transmute(pkg = "Juul", J = paste0("J",J), data.tvec = data.tvec,
              data.upr = value, nsample = 1000)
              

## replicate fda data in the J50 facet, for comparison
envdat <- bind_rows(envdat,
(filter(envdat, pkg == "fda") %>%
 mutate(pkg = "fda_fake", J = "J50")))
 


ggplot(envdat, aes(data.tvec, data.upr)) +
    geom_line(aes(colour=pkg, linetype = nsamp)) +
    facet_wrap(~J) +
    scale_colour_brewer(palette="Dark2")

juul1B <- bind_rows(juul1 ,
                    filter(envdat, pkg == "fda")) %>%
  mutate(across(J, ~as.numeric(gsub("J","",.))))

ggplot(juul1B, aes(data.tvec, data.upr)) +
    geom_line(aes(colour=J, group=interaction(J,pkg), linetype=pkg))


## conclusions: 

## * in general Juul results may be shifted by 1 index point relative to AG results
## (unimportant)
## * for J=50, AG & Juul agree very well
## * for J=2, AG & Juul agree pretty well (I'm not too worried about discrepancies here)
