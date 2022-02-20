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

## Usage

```
python3 main.py -h:
  --mode                Running Mode
  --gpu                   specify the gpu to use
  --test_dir            specify the directory of dataset
  --test_odir          specify the directory that you want to save the processed data
```


### Running Example

```
python main.py --mode=1 --gpu=1 --test_dir="../dataset/test" --test_odir="../dataset/test_processed"
```

## Output file

id_structure.pqr: the prediction is at column 9, which 1.000 means pocket 1, 0.000 means not a pocket atom.
## Visualization
Two examples about pocktes, the first one is example 23 from test set, the second one is example 103 from validation set.
<figure class="half">
    <img src="https://github.com/Zhang038/GNN_pocket/blob/master/example/23_pred.png" width="500"/>
    <img src="https://github.com/Zhang038/GNN_pocket/blob/master/example/103_pred.png" width="500"/>
</figure>
