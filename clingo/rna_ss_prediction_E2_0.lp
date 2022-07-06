%*
===================================================================
Author: Julián María Galindo Álvarez

Implementation of cost function E2 found in 'Exploring Life through 
Logic Programming: Logic Programming in Bioinformatics -  RNA 
secondary structure prediction' available in 
https://computerscience.nmsu.edu/_files/documents/TR-CS-NMSU-2014-10-24.pdf
===================================================================
*%

#include "rna_ss_prediction_base.lp".

%% minimize energy function
#maximize{E : energy2(E)}.										%R13_2

%% length of the sequence
total(N) :- N=#count{X,Y : seq(X,Y)}.							%R14

%% number of AU and CG pairings
au(N) :- N=#count{A,B:pairing(A,B),seq(A,a),seq(B,u)}.			%R15
cg(N) :- N=#count{A,B:pairing(A,B),seq(A,c),seq(B,g)}.			%R16


%% cost function 
energy2(E) :- C1=1, C2=1, C3=1,                                  
			  total(N), contacts(C), au(AU), cg(CG),
			  E = C1 * (N-C/2) + C2 * |100*AU - 35*C| +
				 C3 * |100*CG - 53*C|.	
			  	   
			      
#show energy2/1.