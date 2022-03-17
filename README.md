# cmapR: modified by Loyal

Tools and utilities for loading and analyzing methylation data. The GCT object and .gctx file format are a way to keep data and associated row and column metadata together at all times, including operations for reading, writing and subsetting these data. 

For a tutorial of the available functions in this package, see `tutorial/loyal_basics.R`

### Install instructions
Create a conda environment for the package by installing R and the dependencies. Here I use mamba for increased speed in creating the environment. You can also do this in an environment that already exists or has R installed. 

```
mamba create -n cmapr -c conda-forge -c r -c bioconda r-base=4.1.2
conda activate cmapr
```

Clone the github repo to a convenient location, and install the package locally by running the install R script.
```
git clone https://github.com/Cellular-Longevity/cmapR.git
cd cmapR
Rscript install_cmapR.R
```

In order to use the functions pulling data from AWS, you need to have the `AWS_ACCESS_KEY_ID`, `AWS_DEFAULT_REGION` and `AWS_SECRET_ACCESS_KEY` environment variables configured. You can set these by running `aws configure` or modifying the `~/.aws/credentials` file. 

You should then be all set up! Please file an issue with any problems installing or using this package.

### Citation information
This repo was originally developed as cmap/cmapR. "Parsing and utility functions for analyzing CMap data. To learn more about the CMap project at the Broad Institute, please visit [clue.io](https://clue.io)." If you use GCTx and/or cmapR in your work, please cite [Enache et al.](https://www.biorxiv.org/content/early/2017/11/30/227041)
