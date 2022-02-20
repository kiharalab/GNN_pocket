# GNN_pocket

<a href="https://github.com/marktext/marktext/releases/latest">
   <img src="https://img.shields.io/badge/GNN_pocket-v1.0.0-green">
   <img src="https://img.shields.io/badge/platform-Linux%20-green">
   <img src="https://img.shields.io/badge/Language-python3-green">
   <img src="https://img.shields.io/badge/Language-C-green">
   <img src="https://img.shields.io/badge/dependencies-tested-green">
   <img src="https://img.shields.io/badge/licence-GNU-green">
</a>      <br>
GNN_pocket  is a tool to predict the pocket regions in a protein structure, where the residues that are at the edges of a pocket get tagged as active.

Copyright (C) 2022 Yuanyuan Zhang, Xiao Wang, Charles Christoffer, & Daisuke Kihara, and Purdue University.

License: GPL v3 for academic use. (For commercial use, please contact us for different licensing.)

Contact: Daisuke Kihara (dkihara@purdue.edu)

## Installation

### 1. [`Install git`](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)

### 2. Clone the repository in your machine

```
git clone https://github.com/kiharalab/GNN_pocket
cd GNN_pocket
```
### 3. Build dependencies.

You have two options to install dependency on your computer, install requirements.txt or install dependencies one by one:

#### 3.1 Install with pip and python(Ver 3.6.9).

##### 3.1.1[`install pip`](https://pip.pypa.io/en/stable/installing/).

##### 3.1.2  Install dependency in command line.

```
pip install -r requirements.txt 
```

If you encounter any errors, you can install each library one by one:

```
click==8.0.3
matplotlib==3.4.3
numpy==1.20.3
pandas==1.3.4
scikit_learn==1.0.2
scipy==1.7.1
torch==1.10.1
tqdm==4.62.3
```


## Usage

```
python3 main.py -h:
  --mode                Running Mode
  --gpu                 Specify the gpu to use
  --test_dir            Specify the directory of original dataset, the structure under this dataset 
                        should be "./id/structure.pqr"
  --test_odir           Specify the directory that you want to save the processed data
```


### Running Example

```
python main.py --mode=1 --gpu=1 --test_dir="./dataset/test" --test_odir="./dataset/test_processed"
```

### Running Time

We strongly recommend you to run on a machine with a SSD, the overall running time for an example would be around 2 minutes with a SSD.

But the running time would be much longer without SSD, it is possible to over 10 minutes in some cases with large structures.


## Output file

The output files are saved as final_pred which is at the same directory.
For example:
```
cd ./final_pred
cd 23
```
structure.pqr: the prediction is at column 9, which 1.000 means pocket 1, 0.000 means not a pocket atom.
## Visualization
Two examples about pocktes, the first one is example 23 from test set, the second one is example 142 from validation set.
<figure class="half">
    <img src="https://github.com/kiharalab/GNN_pocket/blob/master/examples/pred_23.png" width="500"/>
    <img src="https://github.com/kiharalab/GNN_pocket/blob/master/examples/pred_142.png" width="500"/>
</figure>
## Referece

- [1] [tools:pqr2pdb](https://github.com/hleonov/bin/blob/master/pqr2pdb.py)
- [2] [tools:VisGrid](https://kiharalab.org/VisGrid/)
- [3] [tools:ghecom](https://pdbj.org/ghecom/)
- [4] [Li, Bin, et al. "Characterization of local geometry of protein surfaces with the visibility criterion." *Proteins: Structure, Function, and Bioinformatics* 71.2 (2008): 670-683.](https://onlinelibrary.wiley.com/doi/abs/10.1002/prot.21732?casa_token=vfhZxgyYvAUAAAAA:FwaQlnRvl1z-jSwfHW_DLb7yjwRn5FXZklyxIyY18mWaAfEjVBDXnv86aLr32z6Jtj8jya4VILWVA-Y)
- [5] [Kawabata, Takeshi. "Detection of multiscale pockets on protein surfaces using mathematical morphology." *Proteins: Structure, Function, and Bioinformatics* 78.5 (2010): 1195-1211.](https://onlinelibrary.wiley.com/doi/abs/10.1002/prot.22639?casa_token=gZElgIqLTX4AAAAA:ioK8V7Ajbzk1Vl9qsvpz8nsg6_6fKZanm6zdoKh9_QM6TJ-hQmVdVpWPvThhwJle9FShzSLkOz00VQc)
- [6] [Kipf, Thomas N., and Max Welling. "Semi-supervised classification with graph convolutional networks." *arXiv preprint arXiv:1609.02907* (2016).](https://arxiv.org/pdf/1609.02907.pdf)
- [7][Wang, Xiao, Sean T. Flannery, and Daisuke Kihara. "Protein docking model evaluation by graph neural networks." *Frontiers in Molecular Biosciences* 8 (2021): 402.](https://www.frontiersin.org/articles/10.3389/fmolb.2021.647915/full)
