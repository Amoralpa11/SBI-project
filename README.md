# SBI-project
A python package for protein complex modelling

## Objective

This python package takes a set of pdb files of pairwise protein interactions
and returns a pdb file with the multiplotein complex that is formed with these
interactions.

## Input

The input, is a set of pdb files holding pairs of proteins interacting. This
input does not have to be perfect:

* There may be many protein complexes that can be formed using the interactions provided.
* The sequences of the same protein appearing in several pdb files may not be the same although they have to be very
similar.
* The same pair of proteins can interact in different ways
* A pair of proteins provided by the user does not have to end up in the protein macro-complex 
* etc.

The user can also introduce some information about the output that wants to
obtain, For instance, the the stoichiometry of the final
complex. They can also provide a subunit limit to define the size of
polymeric complexes which otherwise would be virtually endless.

## Steps

### Classifying the chains

Because chains may not be labeled in a consistent way one of the firsts steps in
the process is labelling all the chains that we can find in input.
Registering how they are called in their pdb files if they are the first of their type
or changing the labels in a consistent way. At this point we are able to identify
the interactions each protein stablishes.
**DO WE???**
**We also check which interactions are possible at the same time.**

###  Macro-complex assembly

Once we have processed and classified all the pairwise interactions from the input set
we start constructing the macro-complex using superimposition.

The approach adopted in this package is based on a recursive function which will
attempt to add as many interactions as each chain forming the macro-complex at 
that stage can.

By doing so we will start with all the different chains the user has found in the
input and recursively build on top of them. Each node of the recursive tree has an
identifier indicating the interactions occuring at that stage, this identifier is 
saved for further usage. This identifier enables us to assess if a macro-complex at 
a specific node has already been processed in a previous node and, therefore, 
stop that branch. 

Before we add a new chain to the macro-complex we check that it is not clashing with 
the other chains already in the structure. Furthermore, once we add it we also register
the interactions this new chain has with the surrounding ones so as to not attempt to
superimpose that interaction in nodes to come thus, reducing processing demand.

Once the recursive function has finished we are able to build the macro-complex/es 
obtained from the identifiers at the final nodes of the recursive tree. This enables
us to work only with one structure throughout the entire thus minimizing the memory
usage of the computer.


If we find a pair of exclusive interactions we should keep an individual thread
for the two possible interactions and either return all the produced complexes
or just the ones that meets the user specifications if set.

### Evaluation

In order to validate the model we should analyse using a prosa like software to
asses the energy of the interacting surfaces that where not in the starting set.


## Complementary software

It would be interesting to make a program that takes a molecular complex and
generates the pdb files with the pairwise interactions.



