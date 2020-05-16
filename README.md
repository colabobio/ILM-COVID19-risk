# ILM-MLE

This repo contains code to find the Maximum likelihood estimates (MLE) of parameters for an epidemic individual-level model (ILM) using Iterated Filtering in the package [POMP](https://kingaa.github.io/pomp/). These parameters could be the coefficients of susceptibility or infectivity scores defined at the individual level. The population level data (observed case counts) is used to try to fit the parameters. Right now, the individual level data are exposure interactions between susceptible and infected individuals, together with the relevant covariates for each individual, simulated with an agent-based model in [GAMA](https://gama-platform.github.io/covid19).

## Usage

The boarding.Rmd notebook is provided as a simple example of MLE for an epidemic model using POMP. To evaluate parameter fitting interactively, use the gama.Rmd notebook. The properties including IF parameters are stored in the .properties files located inside the gama folder.

To run a long MLE calculation, the gama.R script can be used. The accompaning shell script helps launching the R process:

```run.sh gama normal.properties```

where the first argument is the folder containing the data, and the second the name of the properties file to use

## Preparing data

The python script parse.py will read in contact and case data from GAMA, and output the corresponding files needed in POMP. It can be run without any arguments. In order to test that the resulting files were properly parsed, the test.sh script can be run.