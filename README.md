# ILM-MLE

This repo contains code to find the Maximum likelihood estimates (MLE) of parameters for an epidemic individual-level model (ILM) using Iterated Filtering in the package [POMP](https://kingaa.github.io/pomp/). These parameters could be the coefficients of susceptibility or infectivity scores defined at the individual level. The population level data (observed case counts) is used to try to fit the parameters. Right now, the individual level data are exposure interactions between susceptible and infected individuals, together with the relevant covariates for each individual, simulated with an agent-based model in [GAMA](https://gama-platform.github.io/covid19).

## Usage

The boarding.Rmd notebook is provided as a simple example of MLE for an epidemic model using POMP. To evaluate parameter fitting interactively, use the gama.Rmd notebook. The properties including IF parameters are stored in the .properties files located inside the gama folder.

To run a long MLE calculation, the gama.R script can be used. The accompaning shell script helps launching the R process:

```run.sh gama normal.properties```

where the first argument is the folder containing the data, and the second the name of the properties file to use

## Preparing data

The python script parse.py will read in contact and case data from GAMA, and output the corresponding files needed in POMP. It requires two arguments: the folder where the data files generated with GAMA are located, and the maximum time for the resulting case data. 

```python parse.py -i gama2 -t 120```

The flag ```-d``` can also be added to print out contact debug information. This script generates three output files in the specified folder: case_counts.csv containing the new cases per each time unit, and the contact information needed to calculate beta in the SIER model in POMP in two separate files, indices and contacts.

In order to test that the resulting files were properly parsed, the test.sh script can be used to compare the contacts printed by the parse script with the debug options with the contacts as read from the indices and contacts files by the test script:

```./test.sh <dir> <a0> <a1> <b0> <b1>```

The argument ```<dir>``` is the directory where the parsed data files are saved, the arguments a0 through b1 are the coefficients used to generate the GAMA data, and are optional. If provided, the value of beta for each time will calculated as printed as well.

## Cleaning-up

The IF function in POMP creates bake and stew files that are stored in output/bake. The advantage of keeping those files is that the expensive MLE calculations do not need to be run again when those files are already generated, but if the parameters change, they should be removed to generate new ones. Do that with the ```clean.sh``` shell script. Adding the ```all``` argument it will delete all the output files, including plots and parameters.