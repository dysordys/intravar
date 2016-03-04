# intravar

Manuscript, code, and full supporting material for "The effect of intraspecific variation and heritability on community pattern and robustness". Documentation for each file in the repository can be found below.


### sim\_QG\_LV.nb

This Mathematica 10.1 script generates the data in all our simulation sets for the quantitative genetic Lotka-Volterra model. It simulates the model for all parameters and simulation sets. The output is an R-compatible data frame with the rows corresponding to individual simulations, and the following columns:
* **svar**: level of intraspecific variability
  - **svar** = 0: no intraspecific variation
  - **svar** = 1: low levels of intraspecific variation
  - **svar** = 2: mixed levels of intraspecific variation
  - **svar** = 3: high levels of intraspecific variation
* **hsquare**: heritability (either 0, 0.1, or 0.5)
* **omega**: competition width (either 0.1, 0.15, or 0.2)
* **bshape**: shape of intrinsic growth function *b*(*z*)
  - **bshape** = 1: rectangular growth function
  - **bshape** = 2: quadratic growth function
  - **bshape** = 3: triangular growth function
* **replicate**: which simulation set does the given simulation belong to (from 0 to 99)
* **sigma1**, **sigma2**, ...: *S* columns with the *i*th one containing the intraspecific standard deviation of the *i*th species
* **n1**, **n2**, ...: *S* columns with the final densities of the species
* **mu1**, **mu2**, ...: *S* columns with the final mean trait values of the species

In solving the differential equations, we have introduced a weak Allee effect into the Lotka-Volterra dynamics: species whose density falls below a predetermined threshold experience a smooth decline of their per capita growth rates to zero. This is to prevent densities from accidentally dropping to negative values due to numerical error. As long as the densities of extant species are substantially (i.e., several orders of magnitude) greater than the threshold density, this will not make a difference to the outcome of the dynamics.

Note that running the code as is may take several weeks. One way to speed things up is to simulate just a handful of simulation sets at a time, by adjusting the range of sets in line 98. One can then run several processes in parallel, each simulating a different subrange of the simulation sets, and then concatenate the results into a single file. In that case, be sure not to copy the output file headers except at the very top of the concatenated file.


### get_statistics.R

This R script receives the data frame generated by the previous Mathematica script (sim\_QG\_LV.nb), generates the statistics used in the box plots of the manuscript for each simulation, attaches them as new columns to the existing data frame, and saves this updated table into another file.

#### Input

A file containing a data frame with whose structure is that of the output of sim\_QG\_LV.nb.

#### Output

An R-compatible data frame written to a file, with the following columns:
* **svar**: level of intraspecific variability
  - **svar** = 0: no intraspecific variation
  - **svar** = 1: low levels of intraspecific variation
  - **svar** = 2: mixed levels of intraspecific variation
  - **svar** = 3: high levels of intraspecific variation
* **hsquare**: heritability (either 0, 0.1, or 0.5)
* **omega**: competition width (either 0.1, 0.15, or 0.2)
* **bshape**: shape of intrinsic growth function *b*(*z*)
  - **bshape** = 1: rectangular growth function
  - **bshape** = 2: quadratic growth function
  - **bshape** = 3: triangular growth function
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
- [infile]: the name of the input data file (the one generated by the Mathematica script)
- [outfile]: name of file to which modified data should be saved


### sim\_SK\_LV.R

This script simulates the finite-locus eco-evolutionary dynamics with the Shpak-Kondrashov hypergeometric model of inheritance (Shpak and Kondrashov 1999, Evolution 53: 600-604). This is the version with species-specific environmental variation already implemented. Species' environmental standard deviations are independently and uniformly sampled from the interval [0.01, 0.05]. To change this, alter the min and max limits in line 36 of the script.

#### Input

Name of file to save results in, the number of simulation steps, the size (in units of time) of each step, the number of initial species, and the number of loci.

#### Output

An (*S* x *G*) table, where S is the number of species and G is the number of genotypes. The (*i*,  *j*)th entry is the density of species *i*'s genotype *G* at the end of the simulation.


#### To run

Make sure the R packages "tensor" and "pracma" are installed. Then run the script from the command line by invoking

    Rscript SK_LV_model.R [outfile] [steps] [dt] [species] [loci]

where the command line arguments are as follows:
* [outfile]: name of file to save results in
* [steps]: the number of time steps to run the program for (to reach equilibrium, this may need to be large, e.g., 1e7)
* [dt]: the size of each time step (e.g., 0.02)
* [species]: the number of initial species (e.g., 51)
* [loci]: the number of distinct loci contributing to the quantitative trait (e.g., 25 means there are 25 loci and therefore 2 * 25 + 1 = 51 distinct genotypes)
