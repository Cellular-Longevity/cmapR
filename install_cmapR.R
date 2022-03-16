## Installing known dependencies: data.table & rhdf5
if (!require('data.table')) install.packages("data.table")

if (!require('BiocManager')) install.packages("BiocManager")
if (!require('rhdf5')) BiocManager::install(c("rhdf5"))
if (!require('prada')) BiocManager::install(c("prada"))

if (!require('devtools')) install.packages("devtools")
devtools::install(pkg = '.')
