%*
===================================================================
Author: Julián María Galindo Álvarez

Given a RNA sequence in "input/rna_ss_input.lp" encoded in predicates 
of the form seq(sequence_index,sequence_base), this program maximizes
the number of contacts calculated in rna_ss_prediction_base_edited.lp
by minimizing the negative value of the contacts (Nussinov energy).

Based on 'Exploring Life through Logic Programming: Logic Programming in 
Bioinformatics -  RNA secondary structure prediction' available in 
https://computerscience.nmsu.edu/_files/documents/TR-CS-NMSU-2014-10-24.pdf
===================================================================
*%

#include "rna_ss_prediction_base.lp".

%% minimize energy function
#minimize{E: energy1(E)}.		%R13_1

%% Nussinov Energy E1
energy1(-C) :- contacts(C).	    %R13_1_1

#show energy1/1.
