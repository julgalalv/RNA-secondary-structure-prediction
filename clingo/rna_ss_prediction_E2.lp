%*
===================================================================
Author: Julián María Galindo Álvarez

Given a RNA sequence in "input/rna_ss_input.lp" encoded in predicates 
of the form seq(sequence_index,sequence_base), this program maximizes
the followinw energy function:

	E2 = c1(n-|P|) + c2|AU-0.35|P|| + c3|CG-0.53|P|| + c4|GU-0.12|P||

where |P| is the number of contacts calculated in rna_ss_prediction_base_edited.lp,
n is the length of the sequence, AU, CG and GU are the numer of AU, CG and GU 
pairings and c1, c2, c3, c4 are some suitable weigths.

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

%% energy function wit unitary weigths
energy2(E) :- total(N), contacts(C), au(AU), cg(CG), 
			  P = C/2, GU = P - AU - CG,
			  E = 100*(N-C) 
			  	  + |100*AU - 35*P| 
			      + |100*CG - 53*P| 
				  + |100*GU - 12*P|.
				 
#show energy2/1.