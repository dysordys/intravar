# intravar

Manuscript, code, and full supporting material for "The effect of intraspecific variation and heritability on community pattern and robustness"

## sim\_QG\_LV.nb

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
* **sigma1**, **sigma2**, ...: *S* columns with the *i*th one containing the intraspecific standard deviation of the *i*th species (in our case, *S* = 51)
* **n1**, **n2**, ...: *S* columns with the final densities of the species
* **mu1**, **mu2**, ...: *S* columns with the final mean trait values of the species

In solving the differential equations, we have introduced a weak Allee effect into the Lotka-Volterra dynamics: species whose density falls below a predetermined threshold experience a smooth decline of their per capita growth rates to zero. This is to prevent densities from accidentally dropping to negative values due to numerical error. As long as the densities of extant species are substantially (i.e., several orders of magnitude) greater than the threshold density, this will not make a difference to the outcome of the dynamics.

Note that running the code as is may take several weeks. One way to speed things up is to simulate just a handful of simulation sets at a time, by adjusting the range of sets in line 98. One can then run several processes in parallel, each simulating a different subrange of the simulation sets, and then concatenate the results into a single file. In that case, be sure not to copy the output file headers except at the very top of the concatenated file.
