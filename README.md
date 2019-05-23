# Fun CANDIY

Used to host data and methods related to the psychoactives paper published [here](https://chemrxiv.org/articles/Computational_Chemoproteomics_to_Understand_the_Role_of_Selected_Psychoactives_in_Treating_Mental_Health_Indications/6148940).

### Manifest

* README.md

This file.

* all\_disorders\_10.csv.gz
* all\_disorders\_25.csv.gz
* all\_disorders\_40.csv.gz
* all\_disorders\_100.csv.gz

These files contain all CANDO predictions for all indications. They are obtained from transcribing the respective canpredict files obtained from the Samudrala group to the CSV format and are used for the original analysis presented in the paper above. Any further analysis of psychoactive and/or other known drugs should be done using these files.

* all\_new.protein\_compound\_interact\_real.fpt.xz

The 46784x3733 matrix file as obtained from the Samudrala group. 

* compound\_column\_cando\_id\_mapping.tab

A file which maps a compound to its CANDO ID. Listed as obtained from the Samudrala group.

* compound\_indication\_mapping.tab.filtered.v0

A file which maps the **known** treatments of a compound to the indications that it treats.

* process\_fpt.R
* Canpredict.jl

A set of scripts which read the above fpt file and generate the files `all_disorders_10.csv` `all_disorders_25.csv` `all_disorders_40.csv` and `all_disorders_100.csv`. You must run *process\_fpt.R* before running *Canpredict.jl*. Note that the *process\_fpt.R* script produces a modified version of the FPT matrix which is used by *Canpredict. jl* to produce the aforementioned CSV files since the original binary is not available.

* all\_compounds.smi

All compounds present in the CANDO platform.

* cando\_proteins.lst

All proteins present in the CANDO platform.

* psychoactives.tsv

The category, CANDO column ID, CANDO compound ID, and name of all compounds considered to be *psychoactive* in this work.

* cando\_mental\_disorder.lst

The MESH IDs for all mental health indications studied in this work.

* extract.R
* ranker.R
* utils.R

The code required to calculate the compound rank and indication rank.

* random\_analysis.R

The code required to repeat the above scripts after randomizing the compound predictions.

* graph.R

The code required to generate figures 1-5 and figure S1.

* chorddiagram.R

The code required to generate figure 6.

* table.R

The code required to generate the tables in the paper's supporting information.

* main.R

Runs the above R code.


