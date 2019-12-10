# curated_adipo_descriptor

Code for the Adipocyte paper (Modelling the gene expression and the DNA-binding
in the 3T3-L1 differentiating adipocytes)

## Setting up the docker environment

The analysis was run on a [docker](https://hub.docker.com/r/bcmslab/adiporeg/)
image based on the the latest **rocker/verse**.
Other R packages were added to the image and were made available as an image 
that can be obtained and launched on any local machine running
[docker](https://hub.docker.com/r/bcmslab/adiporeg/).

```bash
$ docker pull bcmslab/adiporeg:latest
$ docker run -it bcmslab/adiporeg:latest bash
```

## Obtaining the source code

The source code is hosted publicly on this repository in the form of a research
compendium. This includes the scripts to reproduce the figures and tables in 
this manuscript. Another public repository 
[curated_adipo_descriptor](https://github.com/BCMSLab/curated_adipo_descriptor)
contains the scripts to download, prepare and analyze the data.

```bash
$ git clone http://github.com/BCMSLab/auto_adipo_diff
```

## Runing the analysis

In the directory `auto_adipo_diff`, run `make`

```bash
$ cd auto_adipo_diff
$ make all
```

## Details of the R environment
The version of **R** that was used to perform this analysis is the 3.5.2
(2018-12-20) on `x86\_64-pc-linux-gnu`.

## More

This manuscript was published under the title [Modelling the gene expression and the DNA-binding in the 3T3-L1 differentiating adipocytes](https://www.tandfonline.com/doi/full/10.1080/21623945.2019.1697563)
