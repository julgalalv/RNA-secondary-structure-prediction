%*
===================================================================
Author: Julián María Galindo Álvarez

Given a RNA sequence in "input/rna_ss_input.lp" encoded in predicates 
of the form seq(sequence_index,sequence_base), this program minimizes
the following energy function:

	E2 = (c1*n-|P|) + |AU-0.35|P|/2| + |CG-0.53|P|/2|

where |P| is the number of contacts (number of bases with pairings
, the same as 2 * number of pairings) calculated in 
rna_ss_prediction_base_edited.lp, n is the length of the sequence and
AU, CG are the number of AU, CG pairings. In this model, c1 is a 
learnable parameter.

Based on 'Exploring Life through Logic Programming: Logic Programming in 
Bioinformatics -  RNA secondary structure prediction' available in 
https://computerscience.nmsu.edu/_files/documents/TR-CS-NMSU-2014-10-24.pdf
===================================================================
*%

#include "rna_ss_prediction_base.lp".

%% minimize energy function
#minimize{E : energy2(E)}.										%R13_2

%% length of the sequence
total(N) :- N=#count{X,Y : seq(X,Y)}.							%R14

%% number of AU and CG pairings
au(N) :- N=#count{A,B:pairing(A,B),seq(A,a),seq(B,u)}.			%R15
cg(N) :- N=#count{A,B:pairing(A,B),seq(A,c),seq(B,g)}.			%R16

%% search grid for parameter c1 
c(40;50;60;70;80;100).											%R17
1{c1(C) : c(C)}1.												%R18

%% cost function 
energy2(E) :- total(N), contacts(C), P = C/2,					%R19
			  au(AU), cg(CG), c1(C1),
			  E = |C1*N-100*C| 
			      + |100*AU - 35*P| 
				  + |100*CG - 53*P|.
			  	   
			      
#show c1/1.
#show energy2/1.