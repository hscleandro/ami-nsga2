# Installation and usage documentation of AMI-NSGA2 (2020)
![](https://img.shields.io/badge/last%20edited-10--09--202020-yellow.svg)
![](https://img.shields.io/badge/author-Leandro%20Corrêa-blue.svg)

## CONTENTS OF THIS FILE
* [Introduction](#introduction)
* [Software requirements](#software-requirements)
* [Installation](#installation)
* [Input files](#input-files)
* [Execution](#Exectuion)


## Introduction
This document is about the evolutionary algorithm (EA) for the identification of active modules in biological networks called AMI-NSGA2. Using the concept of the popular algorithm NSGA-II  was developed a new method that applies multiple connected sub-graphs as candidate solution and news adapted genetics opperators (crossover and mutation) using a singular approach based in the connection among the genes in the biological network.

## Software Requirements
* Linux, (32 bit or 64 bit) or MAC OS X Mavericks (64 bit recommended)
* R version >= 3.2.1.

## Installation
Install the R dependences exectuing the script packages.R

## Input files
AMI-NSGA2 supports two types of files: (1) network structure represented by a list of edges or adjacency matrix; (2) A weight-list of differentially expressed genes and their corresponding statistics. In the benchmark/output there are two files containing examples of an interaction network (sim1_connection_list.edgelist) and a weight-list of genes (sim1_weight.csv) in addition to a file containing the simulated artificial active module in this network (sim1_hittrue.csv). 

*Warning*:  The interaction network must not present disconnection between the vertices. That is, all vertices must be parts of a single connected component. Another important point, the weight-list must contain the specific name of the gene in the first column, and its signal variation value (range of 0 to 1) in the second. Both files (the interaction network and weight-list) must contain the same genes name.

## Execution

1.) In the terminal, go to the aminsga2 folder:
```
cd ~/aminsga2/
```
2.) Configure the config.json file to set the hyperparameters or use  "--help" to aid you do this manually.


3.) Running this command:
```
Rscript aminsga2.R --net ~/network_interatcion_file.edgelist --profile weight_list_file.csv --output /outpu_directory/
```

## Help

For information about the AMI-NSGA2 options run: 
```
Rscript aminsga2.R --help
```


## Built With

* [R](https://www.r-project.org/) - Main computational language used
* [python](https://www.python.org/) - Support computational language used for benchmark


## Versioning

AMINSGA-2 version 1.0. Last edited: 10/09/2020

## Authors

* **Leandro Corrêa** - *developer, researcher, software architect* - [personal web page](https://hscleandro.wixsite.com/professional)

## License

This project is licensed under the GNU general public licence - see the [LICENSE](https://github.com/hscleandro/AMINSGA2/LICENSE) file for details

## Acknowledgments

* CNRS [CNRS](http://http://www.cnrs.fr/)
* I3S [I3S](https://www.i3s.unice.fr/en)
