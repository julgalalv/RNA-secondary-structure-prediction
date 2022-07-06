
seq(1,a). seq(2,g). seq(3,u). seq(4,c). seq(5,c). seq(6,a).		
																
%% domain predicates                                            
sequence_index(X) :- seq(X,_).									
sequence_base(B) :- seq(_,B).									
%%% Definition of the pairing function                          
0 {pairing(X,Y):sequence_index(Y)} 1 :- sequence_index(X).		
																
%% the pairing is injective and symmetric                       
:-sequence_index(X1;X2;Y),X1<X2,pairing(X1,Y),pairing(X2,Y).	
pairing(B,A):- sequence_index(A;B),pairing(A,B).				
																
%% wrong associations                                           
wrong(X,X):- sequence_base(X).									
wrong(a,c). wrong(a,g). wrong(c,u).								
:-wrong(B1,B2),seq(X1,B1),seq(X2,B2),pairing(X1,X2).			
																
%% each position can have at most one pairing	                
%% a base cannot match with itself and the successive			
:- sequence_index(X1), pairing(X1,X1).							
:- sequence_index(X1), sequence_index(X2), 						
   X1=X2+1, pairing(X1,X2).                                     
																
																
%% Optional constraint: no pseudo-knots                         
:- sequence_index(X1;X2;X3;X4), X1<X3,X3<X2,X2<X4,			
   pairing(X1,X2),pairing(X3,X4).                               
																
%% Nussinov Energy E1                                           
contacts(C/2):- C = #count{pairing(A,B)}.		 			
#maximize[contacts(C)=C].											
																
#maximize[energy(E)=E].											
total(N) :- N=#count{ seq(X,Y)}.                                
au(N) :- N=#count{pairing(A,B):seq(A,a):seq(B,u)}.            
cg(N) :- N=#count{pairing(A,B):seq(A,c):seq(B,g)}.              
energy(E) :- C1=1, C2=1, C3=1,                                  
			 total(N), contacts(C), au(AU), cg(CG),
			 E = C1 * (N-C/2) + C2 * #abs(100*AU - 35*C) +
				 C3 * #abs(100*CG - 53*C).			    

				   