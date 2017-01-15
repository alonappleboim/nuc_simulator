Simulators Summary:
===================
The simulator is split into several modules, spread across the different folders:
The Working Simulator - contains the main scripts that are called when a simulation runs. Those are the functions for extracting the different sites given a DNA sequence, for translating the sites into eviction, assembly and sliding rates and the functions for running the Monte Carlo simulation itself, given the rates.
Features And Data Comparison - contains the scripts that are used to extract the features from a nucleosome distribution that came from a simulation. In the end, the only feature that was really used in the project was the likelihood, which is given by the get_NFR_features function. The function for comparing the simulation to the experimental data (which returns the likelihood) is Compare_Sum_To_Data.
Sequence Preference Model - contains the scripts and tests used to create the null model of the simulator (which creates the rates assuming that there is no RSC). The useful function in this folder is run_null_model_simulation_from_genome, which replaces the usual function run_simulation_from_genome (which takes RSC into account).
Dynamics - contains the functions for running a dynamic simulation and reviewing the results. Also, the functions for getting the correct time expansion constant are in this folder (getTimeExpansion, dynamicsCostFunction...)


=====|||*|||=====


The Possible Parameters For A Static Simulation:
================================================
parameters that I actually change have an "***" before them...
you can see the defaults of all the parameters in the function run_simulation_from_genome.

***n_steps - the number of steps in the simulation
gen_len - the length of the genome in the simulation (at least 2501)
linker_len - the linker_len parameter (the width of the sigmoid added to the ends of the nucleosome footprint)
nuc_width - the width of the nucleosome (before adding the sigmoids)
REB1_width - the width of the REB1 trans factor
ABF1_width - the width of the ABF1 trans factor
RAP1_width - the width of the RAP1 trans factor
slide_len - the amount of bps in a single slide
poly_rate - the rate of the left sliding that is due to the polymerase
poly_pos - a vector indicating the positions of the polymerase (for example: 1200:2800)
REB1_a_rate - the REB1 assembly rate
REB1_e_rate - the REB1 eviction rate
ABF1_a_rate - the ABF1 assembly rate
ABF1_e_rate - the ABF1 eviction rate
RAP1_a_rate - the RAP1 assembly rate
RAP1_e_rate - the RAP1 eviction rate
nuc_base_a_rate - the global starting assembly rate of the nucs, before changes
                  from the binding sites
nuc_base_e_rate - the global starting eviction rate of the nucs, before changes
                  from the binding sites
nuc_base_r_rate - the global starting right rate of the nucs, before changes
                  from the binding sites
nuc_base_l_rate - the global starting left rate of the nucs, before changes
                  from the binding sites
***TF_evic_intensity - the factor that multiplies the convolution of the
                       TF sites on the nuc eviction rate
***RSC_evic_intensity - the factor that multiplies the convolution of the
                        PolyAT sites on the nuc eviction rate
***RSC_evic_length - the length of the convolution of the
                     PolyAT sites on the nuc eviction rate
***RSC_slide_intensity - the factor that multiplies the convolution of the
                         PolyAT sites on the nuc sliding rates
***RSC_slide_length - the length of the convolution of the
                      PolyAT sites on the nuc sliding rates

					  
The Simulation Stages:
======================
There are a few stages to a simulation (where its only input is a 2501-bp-long sequence of A/C/G/T):
1. Extracting the Poly(dA:dT) and binding sites from the sequence (using Extract_Sites_From_Gene).
2. Generate the left, right, eviction and assembly rates, given the different sites (using generate_rates_from_sites).
3. Running the simulation, given the sliding, eviction and sliding parameters (using Simulate).

Extracting The Sites:
=====================
For now, the extraction is not too sophisticated. For each TF, it searches for its motif in the sequence using PWMs, where the stronger the correlation - the stronger the TF will effect the rates. At this moment the TF part of the simulation has yet to be tested thoroughly and so the simulations that were run did not take into account the TFs...
Poly(dA:dT) patterns are also extracted, where any sequence of 5 As (or more) or 5 Ts (or more) are taken into account.

Generating The Rates:
=====================
What I do, is take the base rates and then add to them convolutions of the site vectors with different windows, that are defined by the simulation parameters:
1. TF sites effect the nucleosome eviction rate. The convolution is between the sites and a step vector which is the length of the TF plus the length of the nucleosome. The step vector is multiplied by the "TF_evic_intensity" parameter.
2. Poly(dA:dT) sites effect eviction and sliding. Here we also do a convolution with a step vector, but this time the vector can have different lengths (added to the base length) and different intensities, all dictated by the "RSC_evic_intensity", "rsc_evic_slide_intensity_ratio", "RSC_evic_slide_length_ratio" and "RSC_slide_length" parameters. The sliding parameters change so that PolyA pushes right and PolyT pushes left.
3. The Polymerase also effects the different rates in its location (downstream of the TSS). it adds a left sliding rate (the "Poly_Rate" parameter) and also weakens the effects of RSC and the TFs. It also increases the rate of TF eviction. All of these haven't been explored too much yet because I'm focusing on the NFR in my simulations so far. Note that for the simulator at the moment, the Polymerase doesn't effect the rates (although it can).

The Simulation Itself:
======================
The simulation keeps a state that changes over time - it is a binary vector, where 1s signify nucleosome centers. Each place that has no nucleosome center in it (and there is room for another nucleosome to enter) has an assembly rate. Every nucleosome also has an eviction rate, left and right sliding rates. The full procedure of the simulation is documented in the "Simulate" function. In the end it returns the sum of all the states that were during the simulation (representing the time that each bp had a nuc center on it), along with other values which are less important.

Comparing The Simulation To The Data:
=====================================
Once we have the vector of centers from the simulation, we can compare it to the vector of reads from the experiment. Right now there is one major way to compare - using a likelihood feature.
The simulation is normalized and then multiplied by the amount of total reads. Then, we assume that that represents the average of nucs that we expect to see in each bp, and that the distribution is a poisson distribution. We check to see what is the chance that, given that distribution, we got the experiment data (looking only at the NFR region and in log scale), and that is the likelihood feature. All this is done in the "Compare_Sum_To_Data" function.
There were also attempts to make other features like +1 position, NFR width and so on, but they were put aside for now.



=====|||*|||=====

The Dynamic Simulator:
======================
The dynamic simulator basically runs two static simulations - one with parameters that are optimal for the t=0 minutes part of the experimental data, and another for our null model (no RSC model). The simulator runs a normal simulation, then passes the final state of that simulation to a new simulation which has the parameters of the null model.
The resulting state history is like a state history of a gene that had RSC kicked out of it in the middle. We then take the simulation time and find the best time expansion coefficient, that when multiplied by the simulation time reproduces a timeline that correlated to the different timestamps of the experimental data.
