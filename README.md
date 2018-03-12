# SBI-project
A python package for protein complex modelling

## Objective

This python package takes a set of pdb files of pairwise protein interactions
and returns a pdb file with the multiplotein complex that is formed with these
interactions.

## Input

The input, is a set of pdb files holding pairs of proteins interacting. This
input does not have to be perfect:

* there may be many protein complexes that can be formed using the interactions provided.
* The sequences of the same protein appearing in several pdb files may not be the same although they have to be very
similar.
* The same pair of proteins can interact in different ways
* etc.

The user can also introduce some information about the output that wants to
obtain, For instance, the user may indicate the stoichiometry of the final
complex. They can also provide a limit of subunits limiting the size of
polymeric complexes that can be virtually endless.

## Steps

### Classifying the chains

Because they may not be labeled in a consistent way one of the firsts steps in
the process would be labelling all the chains that we can find in input.
Registering how are they called in their pdb files or even changing those labels
in a consistent way. At this point we would be able of identify the interactions
that stablish every protein.

We would also need to check which interactions are possible at the same time.

###  assembling

Once we have processed and classified all the interactions inside the input set
we can start constructing the complex.

One possible approach would be choosing a monomer, randomly or based in a
criteria, like the number of different interactions that can have at the same
time. This monomer would be a core, around which we would start adding proteins
in the same way they interact in the pdb file.

Every time we add a protein we would have to check if there is a crash with
another protein already added.

If we find a pair of exclusive interactions we should keep an individual thread
for the two possible interactions and either return all the produced complexes
or just the ones that meets the user specifications if set.

### Evaluation

In order to validate the model we should analyse using a prosa like software to
asses the energy of the interacting surfaces that where not in the starting set.


## Complementary software

It would be interesting to make a program that takes a molecular complex and
generates the pdb files with the pairwise interactions.



