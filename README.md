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

Copyright (C) 2021 Yuanyuan Zhang, Xiao Wang, Charles Christoffer, & Daisuke Kihara, and Purdue University.

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
  --gpu                   specify the gpu to use
  --test_dir            specify the directory of dataset
  --test_odir          specify the directory that you want to save the processed data
```

### Running Example

```
python main.py --mode=1 --gpu=1 --test_dir="../dataset/test" --test_odir="../dataset/test_processed"```T
```

## Output file

Structure.pqr: the prediction is at column 9, which 1.000 means pocket 1, 0.000 means not a pocket atom.

## Visualization

1 example in test set for visualization: 23_structure.pqr
![image](https://github.com/Zhang038/GNN_pocket/blob/master/example/23_pred.png)
