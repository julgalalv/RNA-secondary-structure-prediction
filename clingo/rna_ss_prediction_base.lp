%*
===================================================================
Author: Julián María Galindo Álvarez

Given a RNA sequence in "input/rna_ss_input.lp" encoded in predicates 
of the form seq(sequence_index,sequence_base), this program computes 
the models of possible base pairings of the bases ACGU (A-U,C-G,U-G).
This file is ment to be used in rna_ss_prediction_E1.lp and 
rna_ss_prediction_E2.lp to maximize the corresponding energy function 
in order to predict the secondary structure of the RNA sequence.

Based on 'Exploring Life through Logic Programming: Logic Programming in 
Bioinformatics -  RNA secondary structure prediction' available in 
https://computerscience.nmsu.edu/_files/documents/TR-CS-NMSU-2014-10-24.pdf
===================================================================
*%

#const allow_pseudo_knots=false.
pseudo_knots(allow_pseudo_knots).

#include "input/rna_ss_input.lp".

%% domain predicates
sequence_index(X) :- seq(X,_).                                 %R1
sequence_base(B) :- seq(_,B).                                  %R2
%%% Definition of the pairing function
0 {pairing(X,Y):sequence_index(Y)} 1 :- sequence_index(X).     %R3

%% the pairing is injective and symmetric
:- sequence_index(X1),sequence_index(X2),sequence_index(Y),    %R4
   X1<X2,pairing(X1,Y),pairing(X2,Y).	
pairing(B,A):- pairing(A,B),                                   %R5
			      sequence_index(A), 
			      sequence_index(B).
		   
%% wrong associations
wrong(X,X):- sequence_base(X).                                 %R6
wrong(a,c). wrong(a,g). wrong(c,u).                            %R7
:-wrong(B1,B2),seq(X1,B1),seq(X2,B2),pairing(X1,X2).           %R8

%% each position can have at most one pairing	
%% a base cannot match with itself and the successive				
:- sequence_index(X1), pairing(X1,X1).                         %R9
:- sequence_index(X1), sequence_index(X2),                     %R10
   X1=X2+1, pairing(X1,X2).


%% Optional constraint: no pseudo-knots
:- sequence_index(X1), sequence_index(X2),                     %R11
   sequence_index(X3), sequence_index(X4),
   X1<X3,X3<X2,X2<X4,
   pairing(X1,X2),pairing(X3,X4),
   pseudo_knots(false).

%% Number of contacts (Nussinov Energy E1)
contacts(C):- C = #count{A,B : pairing(A,B)}.                  %R12
			    
#show contacts/1.
#show pairing/2.
#show seq/2.


