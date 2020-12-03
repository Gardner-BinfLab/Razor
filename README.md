## Razor
Razor is a tool to detect signal peptides for *eukaryotic* protein sequences. In addition to signal peptide detection, we also detect:

 - If the signal peptide carries toxic protein.
 - Whether the signal peptide is from fungi.

## Installation
#### Prerequisite
 - Python 3.6+

Download/Clone the source code to your device. In the source code directory, execute these commands:

```
pip3 install -r requirements.txt
```
It is highly recommended to use a virtual environment ```venv``` and install the dependencies to that environment. If you are interested in a [webserver version of this tool](https://tisigner.com/razor), please check [TISIGNER_ReactJS repository](https://github.com/Gardner-BinfLab/TISIGNER-ReactJS).

## Usage

Description of available options:
```
usage: Razor [-h] [-v] -f FASTAFILE [-o OUTPUT] [-m MAXSCAN] [-n NCORES]

A tool to detect signal peptide

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         Show program's version number and exit.
  -f FASTAFILE, --fastafile FASTAFILE
                        Input fasta file
  -o OUTPUT, --output OUTPUT
                        Output file name.
  -m MAXSCAN, --maxscan MAXSCAN
                        Check for cleavage site upto this residue. Default: 80
  -n NCORES, --ncores NCORES
                        Number of cores to use. Default: 1/4 of total cores.
  -q QUIET, --quiet QUIET
                        Do not show warnings. (yes/no). Default: yes

(c) Authors
```

 - ```m``` is the maximum length upto which we scan for the possible cleavage site. For example: ```-m 50``` means we will scan upto 50th residue for the presence of a cleavage site. By default, upto 80th residue is scanned.

 - ```n``` is the number of CPU cores we will use for the computation. This will be turned off if number of sequences is less than 100. Above that, we will use one fourth of your available CPU cores by default.

Sample usage:
```
python3 razor.py -f example_fasta.fa
```

## Description of results
Razor detects signal peptide in the given sequence. If signal peptide is found, it also checks if the signal peptide carries toxic proteins or is from fungi. This is done using 5 random forest model at each detection step. Consequently, we have 5 scores for each step. These scores are described below:

|                                       | Signal peptide | Toxin               | Fungi               |
|---------------------------------------|----------------|---------------------|---------------------|
| Scores from 5 models                  | Y_score        | Toxin_Scores        | Fungi_Scores        |
| Prediction from 5 models              | True/False     | True/False          | True/False          |
| Final scores (Median of scores above) | SP_score       | Toxin_scores_Median | Fungi_scores_Median |

#### Cleavage site identification
Possible cleavage site is the residues where the C-score is maximum. There will be 5 probable cleavage sites form 5 models. The location of the median of these max C-scores is regarded as the final cleavage site. If all of the signal peptide predictions are False, the final cleavage site will be 0 regardless of the values in possible cleavage sites.

Final cleavage site is labelled as ```Cleavage after residue``` in the results.

### Example results file
For an example signal peptide: [Q07310](https://www.uniprot.org/uniprot/Q07310) with cleavage after 27th residue, the result file looks like this. This has a SP_score of 0.75, with 4 out of 5 models returning True in ```SP_prediction```. Looking at predictions for fungi and toxin, we are certain that it does not have any toxic proteins and is not of fungi origin.


| Accession | Sequence                  | Y_score                        | SP_Prediction                   | Max_C                         | Probable Cleavage after | Cleavage after residue | SP_score | Fungi_Scores                   | Fungi_Prediction                    | Fungi_scores_Median | Toxin_Scores                   | Toxin_Prediction                    | Toxin_scores_Median |
|-----------|---------------------------|--------------------------------|---------------------------------|-------------------------------|-------------------------|------------------------|----------|--------------------------------|-------------------------------------|---------------------|--------------------------------|-------------------------------------|---------------------|
| Q07310    | MSFTLHSVFFTLKVSSFLGSLV... | [0.81, 0.75, 0.34, 0.75, 0.76] | [True, True, False, True, True] | [0.87, 0.81, 0.54, 0.8, 0.81] | [27, 27, 27, 27, 27]    | 27                     | 0.75     | [0.06, 0.07, 0.14, 0.09, 0.06] | [False, False, False, False, False] | 0.07                | [0.07, 0.03, 0.04, 0.16, 0.02] | [False, False, False, False, False] | 0.04                |


## Cite
If you find Razor useful, please cite the following paper:
 - Bikash K Bhandari, Paul P Gardner, Chun Shen Lim. (2020). Annotating eukaryotic and toxin-specific signal peptides using Razor. bioRxiv. DOI:[10.1101/2020.11.30.405613](https://doi.org/10.1101/2020.11.30.405613)
