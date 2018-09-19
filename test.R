
source("~/programs/fragxpress/pkg/R/read_efg.R")

base.path <- file.path("/home/raim/programs/fragxpress/data",
                       "Elektroferogramme")

## DNA libraries
fg.path <- file.path(base.path,"Coli Seq Samples")

fgs <- read_efg(path=fg.path, type="aati", verb=1)

samples <- c("5B_160818_1Step",
             "5F_170818_PhCh_Son")


plot_efg(fgs, sample=samples[1], main=samples[1])
plot_efg(fgs, sample=samples[2], main=samples[2])

## Supercoiling
fg.path <- file.path(base.path,"CQ Salima Serie")

fgs <- read_efg(path=fg.path, type="aati", verb=1)

samples <- c("pUC19_topA_0.1ng.uL.1000CQ",
             "high.sensitivity.ladder")

plot_efg(fgs, sample=samples[1], main=samples[1], xlim=c(1200,3000), log="x")
plot_efg(fgs, sample=samples[2], main=samples[2], xlim=c(100,6000),log="x")

