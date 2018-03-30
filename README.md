
# SBI-project

A python package for protein complex modeling from protein pairwise interactions.

## *Table of Contents*:
<!-- TOC depthFrom:1 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->

- [SBI-project](#sbi-project)
	- [*Table of Contents*:](#table-of-contents)
	- [Objective](#objective)
	- [Theoretical background](#theoretical-background)
	- [Input](#input)
		- [Steps](#steps)
			- [Classifying the chains](#classifying-the-chains)
			- [Macro-complex assembly](#macro-complex-assembly)
			- [Evaluation](#evaluation)
			- [Output](#output)
			- [Complementary software](#complementary-software)
	- [Tutorial](#tutorial)
		- [Command line arguments:](#command-line-arguments)
	- [Requirements](#requirements)
	- [Analysis](#analysis)
		- [[1gzx - Hemoglobin](https://www.rcsb.org/structure/1gzx)](#1gzx-hemoglobinhttpswwwrcsborgstructure1gzx)
		- [[3kuy - Nucleosome](https://www.rcsb.org/structure/3kuy)](#3kuy-nucleosomehttpswwwrcsborgstructure3kuy)
		- [[1g65 - proteasome](https://www.rcsb.org/structure/1g65)](#1g65-proteasomehttpswwwrcsborgstructure1g65)
		- [[5vox - Yeast V-ATPase](https://www.rcsb.org/structure/5vox)](#5vox-yeast-v-atpasehttpswwwrcsborgstructure5vox)
		- [[5oom - native assembly intermediate of the human mitochondrial ribosome with unfolded interfacial rRNA](https://www.rcsb.org/structure/5oom)](#5oom-native-assembly-intermediate-of-the-human-mitochondrial-ribosome-with-unfolded-interfacial-rrnahttpswwwrcsborgstructure5oom)
		- [Limitations](#limitations)
	- [Acknowledgement](#acknowledgement)

<!-- /TOC -->

## Objective

This python package takes a set of PDB files of pairwise protein interactions
and returns a PDB file with the multi protein complex that is formed with these
interactions.

## Theoretical background
TODOS

## Input

The input is a set of PDB files holding pairs of proteins interacting. This
input does not have to be perfect:

* There may be many protein complexes that can be formed using the interactions provided.
* The sequences of the same protein appearing in several PDB files may not be the same although they have to be very similar.
* The same pair of proteins can interact in different ways.
* A pair of proteins provided by the user does not have to end up in the protein macro-complex.

The user can also introduce some information about the output that wants to obtain.
For instance, the stoichiometry of the final complex. They can also provide a
subunit limit to define the size of polymeric complexes which otherwise would be
virtually endless.

### Steps

#### Classifying the chains

Because chains may not be labeled in a consistent way one of the firsts steps in
the process is labeling all the chains that we can find in the input.
Registering how they are called in their PDB files if they are the first of their type
or changing the labels in a consistent way. At this point, we are able to identify
the interactions each protein establishes.

####  Macro-complex assembly

Once we have processed and classified all the pairwise interactions from the input set
we start constructing the macro-complex using superimposition.

The approach adopted in this package is based on a recursive function which will
attempt to add as many interactions as each chain forming the macro-complex at
that stage can.

By doing so we start with all the different chains the user has found in the
input and recursively add chains on top of them. Each node of the recursive tree has an
identifier indicating the interactions occurring at that stage, this identifier is saved
for further usage. This identifier enables us to assess if a macro-complex at a specific
node has already been processed in a previous node and, therefore,
stop that branch.

Before we add a new chain to the macro-complex we check that it is not clashing with
the other chains already in the structure. Furthermore, once we add it we also register
the interactions this new chain has with the surrounding ones so as to not attempt to
superimpose that interaction in nodes to come thus, reducing processing demand.

Once the recursive function has finished we are able to build the macro-complex/es
obtained from the identifiers at the final nodes of the recursive tree. This enables
us to work only with one structure to which we add and remove chains every time we move
up/down the tree, therefore, minimizing the memory usage of the computer.

Lastly, we optimize the model using **modeler conjugate gradient**. This function
tweaks the sidechains so as to minimize the overall energy of the macro-complex.

#### Evaluation
*Do we have to evaluate the energy of the model??*
In order to validate the model, we should analyze using a prose like software to
assess the energy of the interacting surfaces that were not in the starting set.


#### Output

The output of this function are PDB files of the macro-complex/es built from the
pairwise interactions. If the user specified passing an already build structure
it returns a new directory named XXXX_all_interactions containing the PDB files
of the pairwise interactions.

#### Complementary software

This package includes complementary functions which can be set by the user.
* Allows the user to pass PDB interaction files of more than 2 proteins.
* Allows the user to pass a full protein complex to be broken into pairwise interaction PDB files.
* Allows the user to pass the stoichiometry so as to specify the components of the macro-complex.
* Allows the user to pass a number of subunits he wants the model to have in case the macro-complex is theoretically infinite.

## Tutorial

In the following section there is a little tutorial in order to make the program
of easy use and for the user to get a grasp of all the functionalities available .

### Command line arguments:
* -i --input: Input a directory containing only the pdb files with the interactions the user wants to process or a macro-complex pdb file if the break option is set.
* -v --verbose: Increase output verbosity, by default is is set to False.
* -k --subunit_limit: Subunit threshold if the protein can theoretically be limitless, by default it is set to False.
* -opt -- optimization: Indicate if you want to optimize the model, by default it is True.
* -int --intensive: Indicate if you want to find all possible structures or just the first one found, by default it will return the first one found, on the other hand, if it is intensive the programm will attempt to find all the possible structures.
* -br --break: Indicate if you want to return all the pairwise interactions or just one of each type. If you pass 'all' the program will output all the pairwise interactions found and if you pass 'unique' it will only return one interaction of each type.
* -st --stoichiometry: Parameter defining the stoichiometry of the complex, the program will assume the different interaction pdb files passed are the interactions forming the macro-complex and will attempt to build the complex using each interaction passed once.

1. Default settings

 *python3 macrocomplex_builder.py -i interactions_3kuy*
 * interactions_3kuy is the path to the directory containing all the interaction pdb files that are to be used in the construction of the macro-complex

2.  Assembly specifications

 *python3 macrocomplex_builder.py -i interactions_3kuy -opt -int*

 This command is the same as before but this time specifying you want to perform an intensive search and to optimize models obtained.
 * -opt: specification for the final model to be optimized.
 * -intensive: specification for the program to return all the possible structures found.

 *python3 macrocomplex_builder.py -i interactions_1TUB -k 100*

 This command in turn specifies the number of subunits the macro-complex has to have. This feature is set for proteins such as tubuline that would otherwise be endless.
 * -k 100: indicates the macro-complex has to have a maximum of 100 subunits.

 *python3 macrocomplex_builder-py -i interaction_3kuy -st*

 This command allows the user to pass the stoichiometry he desires to the model. To do so the program will first detect all the non redundant chains and then return to the user one by one these chains. The user will then have 2 options: i) set a number, which will correspond to the number of chains of that type he wishes the final model to have and ii) press enter, which will be as if no stoichiometry is set for that chain, the program will try to fit as many as it can in the final structure.

3. Get interactions

 *python3 macrocomplex_builder.py -i interactions_3kuy -br all*

 This command is different form the previous ones. In this case you pass a fully built macro-complex and the output will be a directory with all the pairwise interactions found in the structure.

 *python3 macrocomplex_builder.py -i interactions_3kuy -br unique*

 This command, as the previous one, returns the interactions forming the complex but in this case without any redundancies. It will return a directory with non-redundant interactions.

## Requirements

In order to run this package with all its functionalities the user must have several programs:
* Python 3
* Python modules
  *  [Modeller version 9.19](https://salilab.org/modeller/download_installation.html)
  *  Biopython
  *  numpy
  *  gzip
  *  re
  *  copy
  *  argparse

## Analysis

In the following section we are going to discuss some examples of inputs-outputs
and how it worked for each one.
### [1gzx - Hemoglobin](https://www.rcsb.org/structure/1gzx)

This is an example of complex formed by 2 different chains that form 3 types of
interactions.This program is able to properly fully reconstruct this complex with
the interactions given as we can see in the images below. Our program, thus, has
no problem dealing with this type of interactions. We can see there are no how there
are minimal to no differences between the superimposed original 1gzx and the ones
the program has built.


| <img src="https://github.com/Amoralpa11/SBI-project/blob/complex_breaker/img/1gzx_original.png" width="200" height="200"> | <img src="https://github.com/Amoralpa11/SBI-project/blob/complex_breaker/img/1gzx_built.png" width="200" height="200"> | <img src="https://github.com/Amoralpa11/SBI-project/blob/complex_breaker/img/1gzx_built_optimized.png" width="200" height="200"> |
| :---: | :---: | :---: |
| *Original* | *Built-No Optimization* | *Optimization*|

### [3kuy - Nucleosome](https://www.rcsb.org/structure/3kuy)

This is an example of complex formed by 10 different chains, 8 protein and 2 DNA.
In the following images we can see the original 3kuy, the 3kuy we built superimposed
the real 3kuy and the optimized 3kuy with, again, the original 3kuy superimposed.
No flaws are seen in the built complexes, therefore, we can say that the
program can work with DNA as well as with protein chains with relatively simple
structures.

| <img src="https://github.com/Amoralpa11/SBI-project/blob/complex_breaker/img/3kuy_original.png" width="200" height="200"> | <img src="https://github.com/Amoralpa11/SBI-project/blob/complex_breaker/img/3kuy_built.png" width="200" height="200"> | <img src="https://github.com/Amoralpa11/SBI-project/blob/complex_breaker/img/3kuy_built_optimized.png" width="200" height="200"> |
| :---: | :---: | :---: |
| *Original* | *Built-No Optimization* | *Optimization*|

### [1g65 - proteasome](https://www.rcsb.org/structure/1g65)

When we run a more complex model, 30 protein chains, the program can still rebuild
it. Due to eliminating redundancy we can see that some chains are slightly shifted
with respect to the original molecule but the complex is till built without missing
any chain

| <img src="https://github.com/Amoralpa11/SBI-project/blob/complex_breaker/img/1g65_original.png" width="200" height="200"> | <img src="https://github.com/Amoralpa11/SBI-project/blob/complex_breaker/img/1g65_built.png" width="200" height="200"> | <img src="https://github.com/Amoralpa11/SBI-project/blob/complex_breaker/img/1g65_built_optimized.png" width="200" height="200"> |
| :---: | :---: | :---: |
| *Original* | *Built-No Optimization* | *Optimization*|

### [5vox - Yeast V-ATPase](https://www.rcsb.org/structure/5vox)

When we run a more complex model with a particular chain composition, knowing the
chain stoichiometry is necessary. If the stoichiometry is not set it will fit all
the chains it can in the model even if the original model doesn't have them. If it
is set then it will respect the number of chains but it won't necessarily respect
the position some of these chains are in the original model. If its a lateral
chain/chain compound, for example then it will place it in the first place it can.

| <img src="https://github.com/Amoralpa11/SBI-project/blob/complex_breaker/img/5vox_original.png" width="150" height="150"> | <img src="https://github.com/Amoralpa11/SBI-project/blob/complex_breaker/img/5vox_built.png" width="150" height="150"> | <img src="https://github.com/Amoralpa11/SBI-project/blob/complex_breaker/img/5vox_built_optimized.png" width="150" height="150"> | <img src="https://github.com/Amoralpa11/SBI-project/blob/complex_breaker/img/5vox_built_st.png" width="150" height="150"> | <img src="https://github.com/Amoralpa11/SBI-project/blob/complex_breaker/img/5vox_built_optimized_st.png" width="150" height="150"> |
| :---: | :---: | :---: | :---: | :---: |
| *Original* | *Built-No Optimization* | *Optimization*| *Built-Not Opt-St* | *Built-Opt-St* |

### [5oom - native assembly intermediate of the human mitochondrial ribosome with unfolded interfacial rRNA](https://www.rcsb.org/structure/5oom)

5omm is another example of this programs robustness. As we can see in the images
below the program was able to reconstruct this complex compounded of 53 chains. It
did not place 2 small lateral alpha-helix probably because it didn't have the
interaction file since these helices in the pdb files are made out of "unkown"
residues and atoms. Therefore, our complex breaker didn't identify the interaction
and was unable to incorporate it in the complex building process. Still, the
program was able to rebuild it taking into account both protein and RNA chains.

| <img src="https://github.com/Amoralpa11/SBI-project/blob/complex_breaker/img/5oom_original.png" width="200" height="200"> | <img src="https://github.com/Amoralpa11/SBI-project/blob/complex_breaker/img/5oom_built.png" width="200" height="200"> | <img src="https://github.com/Amoralpa11/SBI-project/blob/complex_breaker/img/5oom_built_optimized.png" width="200" height="200"> |
| :---: | :---: | :---: |
| *Original* | *Built-No Optimization* | *Optimization*|

### Limitations

The main limitations of this program are the following:

* When performing the reduction of repeated interactions in the input by similarity we set certain thresholds of similarity. These thresholds can lead to consider subunits formed by the same chains and that have the same interaction but that are slightly rotated, so as to fit in a barrel for example, to be deprecated. In these cases we only gather one of these interactions and ultimately when building the complex, barrel in this case, the last subunit may not be added due to clashes derived from an accumulation of little rotation mistakes in the rest of subunits forming the barrel.

* With the recursive approach the amount of processing time when many chains with many interactions are passed and an intensive search is specified can be very large. In these cases we recommend to use the option -simp of the program which returns the first structure found, therefore, reducing processing time but with the inconvenient that the structure returned may not be the desired by the user.

* We would have liked to give the user a little more control over the parameters used to define clashes and interactions allowing him to pass a file with the specifications of the parameters. For example, if he wanted to use CA or CB to measure distances, minimum distance to consider a clash/interaction, etc...

## Acknowledgement
All the protein visualization was done using [Chimera.](https://www.cgl.ucsf.edu/chimera/)
All the code of this project was programmed in [Python 3.6.3.](https://www.python.org/downloads/release/python-363/)
