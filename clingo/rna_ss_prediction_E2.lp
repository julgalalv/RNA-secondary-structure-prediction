%*
===================================================================
Author: Julián María Galindo Álvarez

Given a RNA sequence in "input/rna_ss_input.lp" encoded in predicates 
of the form seq(sequence_index,sequence_base), this program maximizes
the followinw energy function:

	E2 = c1(n-|P|) + c2|AU-0.35|P|| + c3|CG-0.53|P||

where |P| is the number of contacts calculated in rna_ss_prediction_base_edited.lp,
n is the length of the sequence, AU and CG are the numer of AU and CG 
pairings and c1, c2, c3 are some suitable weigths.

Based on 'Exploring Life through Logic Programming: Logic Programming in 
Bioinformatics -  RNA secondary structure prediction' available in 
https://computerscience.nmsu.edu/_files/documents/TR-CS-NMSU-2014-10-24.pdf
===================================================================
*%

#include "rna_ss_prediction_base.lp".

%% maximize energy function
#maximize{E : energy2(E)}.										%R13_2

%% length of the sequence
total(N) :- N=#count{X,Y : seq(X,Y)}.							%R14

%% number of AU and CG pairings

au(N) :- N=#count{A,B:pairing(A,B),seq(A,a),seq(B,u)}.			%R15
cg(N) :- N=#count{A,B:pairing(A,B),seq(A,c),seq(B,g)}.			%R16

%% energy function wit unitary weigths
energy2(E) :- C1=1, C2=1, C3=1,									%R17
			 total(N), contacts(C), au(AU), cg(CG),
			 P = C/2,
			 E = C1 * 100 * (N-P) + C2 * |100*AU - 35*P| +
			     C3 * |100*CG - 53*P|.
				 
#show energy2/1.