# intravar

Manuscript, Supporting Information, and code used in "The effect of intraspecific variation and heritability on community pattern and robustness". In this repository, you will find:
* **intravar_final.pdf**: final manuscript file, with the supporting information attached
* **sim\_QG\_LV.nb**: code for simulating the quantitative genetic Lotka-Volterra model
* **sim\_QG\_LV\_2D.nb**: code for simulating the quantitative genetic Lotka-Volterra model with two trait dimensions
* **sim\_QG\_LV\_all.m**: code to generate all our data using the quantitative genetic Lotka-Volterra model
* **get_statistics.R**: code for obtaining species richness, trait pattern, and community robustness information
* **sim\_SK\_LV.R**: code to simulate model with hypergeometric model of inheritance

Documentation for each code file in the repository can be found below.


## sim\_QG\_LV.nb

This Mathematica 10.1 notebook simulates a single instance of the quantitative genetic Lotka-Volterra model, and plots the results. One can make a movie of the dynamics unfolding.

In solving the differential equations, we have introduced a weak Allee effect into the Lotka-Volterra dynamics: species whose density falls below a predetermined threshold experience a smooth decline of their per capita growth rates to zero. This is to prevent densities from accidentally dropping to negative values due to numerical error. As long as the densities of extant species are substantially (i.e., several orders of magnitude) greater than the threshold density, this will not make a difference to the outcome of the dynamics.

#### Input

There is no special input: the notebook can be executed as is. To change parameters, simply modify the variables under the "Define parameters" header.

#### Output

All output is generated within the notebook itself (time series, plots, and movie).

#### To run

Import the notebook in Mathematica and execute.


## sim\_QG\_LV\_2D.nb

This Mathematica 10.1 notebook simulates a single instance of the quantitative genetic Lotka-Volterra model with two trait dimensions, and plots the results.

In solving the differential equations, we have introduced a weak Allee effect into the Lotka-Volterra dynamics: species whose density falls below a predetermined threshold experience a smooth decline of their per capita growth rates to zero. This is to prevent densities from accidentally dropping to negative values due to numerical error. As long as the densities of extant species are substantially (i.e., several orders of magnitude) greater than the threshold density, this will not make a difference to the outcome of the dynamics.

#### Input

There is no special input: the notebook can be executed as is. To change parameters, simply modify the variables under the "Define parameters" header.

#### Output

All output is generated within the notebook itself (time series and plots).

#### To run

Import the notebook in Mathematica and execute.


## sim\_QG\_LV\_all.m

This Mathematica 10.1 script generates the data in all our simulation sets for all parameter values, for the quantitative genetic Lotka-Volterra model.

In solving the differential equations, we have introduced a weak Allee effect into the Lotka-Volterra dynamics: species whose density falls below a predetermined threshold experience a smooth decline of their per capita growth rates to zero. This is to prevent densities from accidentally dropping to negative values due to numerical error. As long as the densities of extant species are substantially (i.e., several orders of magnitude) greater than the threshold density, this will not make a difference to the outcome of the dynamics.

#### Input

There is no special input: the script can be executed as is. However, near the top of the script, be sure to change the variable `[outfile]`, which determines where and under what name the output will be saved.

#### Output

An R-compatible data frame with the rows corresponding to individual simulations, and the following columns:
* **svar**: level of intraspecific variability (0: no variation; 1: low levels; 2: mixed levels; 3: high levels of intraspecific variation)
* **hsquare**: heritability (either 0, 0.1, or 0.5)
* **omega**: competition width (either 0.1, 0.15, or 0.2)
* **bshape**: shape of intrinsic growth function *b*(*z*) (1: rectangular growth;  function; 2: quadratic growth function; 3: triangular growth function)
* **replicate**: which simulation set does the given simulation belong to (from 0 to 99)
* **sigma1**, **sigma2**, ...: *S* columns with the *i*th one containing the intraspecific standard deviation of the *i*th species
* **n1**, **n2**, ...: *S* columns with the final densities of the species
* **mu1**, **mu2**, ...: *S* columns with the final mean trait values of the species

The data frame gets saved into `[outfile]`, defined near the top of the script.

#### To run

Import the script in Mathematica, then execute "Run Package".

Note that running the code as is may take several weeks. One way to speed things up is to simulate just a handful of simulation sets at a time, by adjusting the range of sets (second line from the bottom in the script). One can then run several processes in parallel, each simulating a different subrange of the simulation sets, and then concatenate the results into a single file. In that case, be sure not to copy the output file headers except at the very top of the concatenated file.


## get_statistics.R

This R script receives the data frame generated by the previous Mathematica script (sim\_QG\_LV\_all.nb), generates the statistics used in the box plots of the manuscript for each simulation, attaches them as new columns to the existing data frame, and saves this updated table into another file.

#### Input

A file containing a data frame with whose structure is that of the output of sim\_QG\_LV\_all.nb.

#### Output

An R-compatible data frame written to a file, with the following columns:
* **svar**: level of intraspecific variability (0: no variation; 1: low levels; 2: mixed levels; 3: high levels of intraspecific variation)
* **hsquare**: heritability (either 0, 0.1, or 0.5)
* **omega**: competition width (either 0.1, 0.15, or 0.2)
* **bshape**: shape of intrinsic growth function *b*(*z*) (1: rectangular growth;  function; 2: quadratic growth function; 3: triangular growth function)
* **replicate**: which simulation set does the given simulation belong to (from 0 to 99)
* **sigma1**, **sigma2**, ...: *S* columns with the *i*th one containing the intraspecific standard deviation of the *i*th species
* **n1**, **n2**, ...: *S* columns with the final densities of the species
* **mu1**, **mu2**, ...: *S* columns with the final mean trait values of the species
* **richness**: the number of coexisting species in the final state
* **CV**: coefficient of variation of adjacent species trait differences
* **pval**: p-values of the coefficients of variation
* **gmean**: average community robustness
* **gsd**: robustness heterogeneity

#### To run

Invoke, from the command line,

    Rscript get_statistics.R [infile] [outfile]

where the command line arguments are:                                                                   
- `[infile]`: the name of the input data file (the one generated by the Mathematica script)
- `[outfile]`: name of file to which modified data should be saved


## sim\_SK\_LV.R

This R script simulates the finite-locus eco-evolutionary dynamics with the Shpak-Kondrashov hypergeometric model of inheritance (Shpak and Kondrashov 1999, Evolution 53: 600-604). This is the version with species-specific environmental variation already implemented. Species' environmental standard deviations are independently and uniformly sampled from the interval `[0.01, 0.05]`. To change this, alter the min and max limits in line 36 of the script.

#### Input

Name of file to save results in, the number of simulation steps, the size (in units of time) of each step, the number of initial species, and the number of loci.

#### Output

An (*S* x *G*) table, where S is the number of species and G is the number of genotypes. The (*i*,  *j*)th entry is the density of species *i*'s genotype *G* at the end of the simulation.


#### To run

Make sure the R packages "tensor" and "pracma" are installed. Then run the script from the command line by invoking

    Rscript SK_LV_model.R [outfile] [steps] [dt] [species] [loci]

where the command line arguments are as follows:
* `[outfile]`: name of file to save results in
* `[steps]`: the number of time steps to run the program for (to reach equilibrium, this may need to be large, e.g., 1e7)
* `[dt]`: the size of each time step (e.g., 0.02)
* `[species]`: the number of initial species (e.g., 51)
* `[loci]`: the number of distinct loci contributing to the quantitative trait (*n* loci result in 2*n* + 1 distinct genotypes)
