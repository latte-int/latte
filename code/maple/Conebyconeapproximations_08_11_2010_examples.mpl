### Examples split out from main file

$include "Conebyconeapproximations_08_11_2010.mpl";

#########################################################################
## SOME EXAMPLES OF WEIGHTED EHRHART.


Delta:=[[0,0],[5/28,0],[5/28,5/14]];

#######################################################################################


##########



CompleteEhrhartweighted(n,Delta,[1,1],1);

######################################################################################
#I verified that it takes integral values; I could verify the complete ehrhart against the one
# computed with S^L; L:=[]; but the two program looks that I cannot run them together... (I #suppose because of linalg)
#One example:
knapsack:=proc(n) local ze, S,j,zej; ze:=[seq(0,i=1..n)];
    S:=[ze];
    for j from 1 to n do zej:=subsop(j=1,ze);
        S:=[op(S),zej/j];
    od;
end:
Sourat:=proc(n) local ze, S,j,zej; ze:=[seq(0,i=1..n)];
    S:=[ze];
    for j from 1 to n do zej:=subsop(j=j,ze);
        S:=[op(S),zej/(j+1)];
    od;
end:
#############  DIMENSION 4 ###################################################################
R4 := [[0, 1, 1/3, 2], [1, 0, 1/3, 1], [1, 0, 0, 0], [2/3, 1, 0, 2/3], [0, 1/2, 2, 0]];
#simplify(TopEhrhartweighted(n,R4,[1,1,1,1],0,4));
#simplify(TopEhrhartweighted(n,R4,[1,1,1,1],0,3));
#simplify(TopEhrhartweighted(n,R4,[1,1,1,1],0,2));
#TopEhrhartweighted(n,R4,[1,1,1,1],1,5);
#TopEhrhartweighted(n,R4,[1,1,1,1],1,4);
#TopEhrhartweighted(n,R4,[1,1,1,1],1,3);

#############  DIMENSION 5 ###################################################################
#R5:=[[0, 4/5, 0, 0, 3/4], [1, 1, 3, 5/6, 4], [0, 5, 3/5, 3/5, 0], [1/5, 4/5, 1/3, 5/6, 0], [1/2, 1/2, 1/4, 0, 1/5], [5/2, 4/5, 0, 1/2, 1/5]];#random_rational_simplex(6,5);
#TopEhrhartweighted(n,R5,[1,1,1,1,1],0,5);
#TopEhrhartweighted(n,R5,[1,1,1,1,1],0,4);
#TopEhrhartweighted(n,R5,[1,1,1,1,1],0,3);
#TopEhrhartweighted(n,R5,[1,1,1,1,1],1,6);
#TopEhrhartweighted(n,R5,[1,1,1,1,1],1,5);
#TopEhrhartweighted(n,R4,[1,1,1,1,1],1,4);

######################################################################################

#EXAMPLES FOR THE REAL EHRHART POLYNOMIAL versus integers:

CompleteEhrhartweighted(n,Delta,[1,1],0);
CompleteEhrhartweighted_real(n,Delta,[1,1],0);

lattice3:=[[0,0,0],[1,0,0],[1,2,0],[0,0,1]];

divlattice:=proc(div); [[0,0,0],[div*1,0,0],[div*1,div*2,0],[0,0,div*1]];end:

CompleteEhrhartweighted(n,divlattice(1),[1,1,1],0);
CompleteEhrhartweighted_real(n,divlattice(1),[1,1,1],0);
CompleteEhrhartweighted(n,divlattice(1/2),[1,1,1],0);
CompleteEhrhartweighted_real(n,divlattice(1/2),[1,1,1],0);

# I
