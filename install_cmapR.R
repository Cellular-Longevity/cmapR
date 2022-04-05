## Installing known dependencies: data.table & rhdf5
if (!require('data.table')) install.packages("data.table",repos = "http://cran.us.r-project.org")
if (!require('BiocManager')) install.packages("BiocManager",repos = "http://cran.us.r-project.org")

BiocManager::install(version = '3.14', ask=F)
if (!require('rhdf5')) BiocManager::install(c("rhdf5"),version = "3.14")
if (!require('prada')) BiocManager::install(c("prada"),version = "3.14")
if (!require('matrixStats')) BiocManager::install('matrixStats',version = "3.14")

if (!require('devtools')) install.packages("devtools", repos = "http://cran.us.r-project.org")
devtools::install(pkg = '.')

