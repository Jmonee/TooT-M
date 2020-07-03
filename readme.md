# TooT-M
This tool predicts transporter proteins.
 
Input: proteins sequences in Fasta format
Output: the predicted location,1=membrane, 0=non-membrane



The training dataset  can be found [here](https://tootsuite.encs.concordia.ca/datasets/membrane/)


## FOLDERS
There are a number of folders that support the running of TooT-M and its outputs.

### intermediate_files
Contains the homology details needed to extract the features. Details of the  `psiBlast hits` for each sequence is found here.

intermediate_files/Compostions: Contains the extracted `MSA_TPAAC` features of test set

### db
Contains the database to be used when performing BLAST.


### src
The scripts needed to use the tool.

## HOW TO USE
 - This tool requires that `BLAST` be pre-installed
 - Usage: `Rscript src/TooT_M_V1_11.R -query=<input> [-TooTM=<TooTMdir>] [-out=<outdir>] [-db=<path to db>][-topcons2=<TOPCONS2 results path>]`
  - `<input>` is your sequence input file in fasta format
  - `<out>` is the output directory where you want the predicted 	results, formatted as csv
  - `<TooTMdir>` is the directory where the base TooT-M files 	are located
  - `<path to db> is the directory where the database is stored`
   - `<TOPCONS2 results path> is the directory where  TOPCONS2 results (.txt) is stored; results can be obtained by there available website [here](https://tootsuite.encs.concordia.ca/datasets/membrane/)  or by downloading the tool as described [here](https://github.com/ElofssonLab/TOPCONS2)`
 - Pse-PSSM features (lambda in (0,1,...49) of each sequence in the test set is  found under [intermediate_files/Compositions/Testing/](intermediate_files/Compositions/Testing/)
