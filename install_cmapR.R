## Installing known dependencies: data.table & rhdf5
if (!require('data.table')) install.packages("data.table",repos = "http://cran.us.r-project.org")

if (!require('BiocManager')) install.packages("BiocManager",repos = "http://cran.us.r-project.org")
if (!require('rhdf5')) BiocManager::install(c("rhdf5"))
if (!require('prada')) BiocManager::install(c("prada"))

if (!require('devtools')) install.packages("devtools", repos = "http://cran.us.r-project.org")
devtools::install(pkg = '.')
