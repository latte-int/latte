 #---------------------------------------------------------
#This is MAPLE code to generate random $k$-simplices
# inside the Birkhoff polytope B_n
# THe user can prescribe in advance forbidden positions in the
# permutation matrices with the idea of forcing some simplices
# to share some facet with B_n.
#
# INPUT:
# n:=Birkhoff size of matrices; k:=dimension of simplex;
# F:=list of forbidden positions given as a list of 2-tuples
#--------------------------------------------------------

# GLOBAL INITIALIZATION
with(combinat):
with(LinearAlgebra):
with(linalg):
interface(quiet=true):
randomize():

#-------------------------------------------
# This subroutine transfers a list of lists OR a matrix into a long n^2 vector.
#------------------------------------------------------------------------------
straighten:=proc(M)
local a,JJ,h;
JJ:=convert(M,listlist);
h:=[];
for a in JJ do
  h:=[op(h),op(a)];
od;
return(h);
end proc:

#------------------------------
# Takes the Birkhoff number $n$ and a list of coordinates in the
# $n$-by-$n$ matrix which MUST be set to 0.  Outputs a straightened
# vector of a permuated identity matrix which satisfies the required 0 spaces.
#-----------------------------------------
RestrictedPositionRandPerm:=proc(n,F)
local roll,i, P, In, frow, fcol, prow, pcol, forbidden, temp, count;

count:=1;

# a list of the number of available spots at an index.
# each index corresponds to a row of P
# initialize to n, as nothing is forbidden currently.
In:=[n $j=1..n];

# a list for each row containing what the numbers are.  -1 means
# it's not yet set, 0 and 1 mean it's set to 0 and 1, respectively.
# initialize to unset
P:=[]:
for i from 1 to n do
	P:=[op(P),[-1 $j=1..n]];
od;

# begin reading in the forbidden spaces
for i from 1 to nops(F) do
	forbidden:=op(i, F);

	# the row coordinate
	frow:=op(1, forbidden);

	# the col coordinate
	fcol:=op(2, forbidden);

	# decrement the correct position in In to reflect a space being set
	In:=subsop(frow = (op(frow, In) - 1), In);

	#remove the possible column from the row's possibilities
	P:=subsop(frow=subsop(fcol = 0, op(frow, P)),P);
od;


# begin greedily choosing random positions in the row with the
# least number of possibilities, removing it from the other rows
# if we get stuck, the prescribed subspace is impossible to achieve
# (i.e. 3 rows with only the same 2 spaces available for 1's)

while count<n+1 do
	prow:=1;
    while op(prow,In) < 1 do
        prow:=prow+1;
    od;
	# find the row with the least number of possibilities and save it to row
	for i from prow+1 to nops(In) do
		if op(i, In) > 0 then
			if op(i, In) < op(prow, In) then
				prow:=i;
			fi; #if
		fi; #if
	od; #for


	# make sure there is an available position in the row, otherwise notify and quit
	if op(prow, In)=0 then
		print("Why would you do this?  It's impossible!, choose different forbidden positions");
		quit;
	fi; #if

    In:=subsop(prow=0,In);
     roll:=rand(1..n):

	# begin checking random numbers until an index is one
	# we can set and set it to 1, turn everything else into 0 and save the column so we can
	# make that column 0 for the rest of the rows
	temp:= roll();
	pcol:=-1;
	while pcol < 0 do
		if op(temp, op(prow, P))=-1 then
			pcol:=temp;
		else
            temp:=roll();
		fi;
	od; # while
	for i from 1 to n do
		if i=pcol then
          P:=subsop(prow=subsop(i = 1, op(prow, P)),P);
		else
          P:=subsop(prow=subsop(i = 0, op(prow, P)),P);
		fi;
	od; # for
	for i from 1 to n do
		if op(i, In)>0 then
          if op(pcol,op(i,P))=-1 then
             In:=subsop(i = (op(i, In) - 1), In);
          fi;
            P:=subsop(i=subsop(pcol = 0, op(i, P)),P);
		fi; # if
	od; # for
	count:=count+1;
od;

# returns the list of lists P
return(P);
end proc:


#----------------------------------------------------------
# Roughly speaking, the permutations are generated randomly, first
# from among those that respect the forbidden positions F, while
# this is possible. Then general random permutations are allowed to
# fill the matrix until we reach the desired rank r.
#
# Finally to known when it is time to quit looking for points in that face
# We use the fact that in a k-dimensional face of B_n the following
# inequality (n-1)^2-k<=nops(F) is satisfied. That help us decide when
# we need to add more points from outside this face.
#
#-----MAIN ROUTINE--------------------------

gensimplexB_n:=proc(n,k,F)
local V,II,simp,r,U,J,simp2;

	if k> (n-1)^2+1 then
		printf("This is impossible: %d >= (%d -1)^2 + 1 = %d\n", k, n, (n-1)^2 +1);
		quit;
	fi;
	II:=IdentityMatrix(n);
	simp2:=[straighten(II)];
	J:=SubMatrix(IdentityMatrix(n),[$1..(n-1)],[$1..(n-1)]):
	simp:=[straighten(J)]:

	Rank(Matrix(simp));
	r:=1:
	while r<k do
		if (r <=(n-1)^2-nops(F)) and nops(F)>0 then
			U:=Matrix(RestrictedPositionRandPerm(n,F));
		else
			 U:=SubMatrix(IdentityMatrix(n),[$1..n],randperm(n));
		fi;
		V:=SubMatrix(U,[$1..(n-1)],[$1..(n-1)]);
		simp2:=[op(simp2),straighten(U)];
		if r<Rank(Matrix(simp2)) then
			simp:=[op(simp),straighten(V)]:
			r:=r+1;
		fi;
	od;

	return(matrix(simp));
end proc:

gensimplexB_n_nonProjected:=proc(n,k,F, dataBaseFileName)
local V,II,simp,r,U,J,simp2;
local nonProjectedSimplex, numberIterations;

	numberIteration:=0; # the number of times we tried to find a NEW simplex that was not already in the database.

	if k> (n-1)^2+1 then
		printf("This is impossible: %d >= (%d -1)^2 + 1 = %d\n", k, n, (n-1)^2 +1);
		quit;
	fi;




	II:=IdentityMatrix(n);
	simp2:=[straighten(II)];
	J:=SubMatrix(IdentityMatrix(n),[$1..(n-1)],[$1..(n-1)]):
	simp:=[straighten(J)]:

	Rank(Matrix(simp));
	nonProjectedSimplex:=[];
	r:=1:
	while r<k do
		if (r <=(n-1)^2-nops(F)) and nops(F)>0 then
			U:=Matrix(RestrictedPositionRandPerm(n,F));
		else
			 U:=SubMatrix(IdentityMatrix(n),[$1..n],randperm(n));
		fi;
		V:=SubMatrix(U,[$1..(n-1)],[$1..(n-1)]);
		simp2:=[op(simp2),straighten(U)];
		if r<Rank(Matrix(simp2)) then
			simp:=[op(simp),straighten(V)]:
			nonProjectedSimplex:=[op(nonProjectedSimplex), straighten(U)];
			r:=r+1;
		fi;
	od;

	return([simp, nonProjectedSimplex]);
end proc:

#--------------EXAMPLE-----------------------
# This is a Chan-Robbins-Yuen face of B_5 and we are sampling
# from it.

#Zeros:=[[1,3],[1,4],[1,5],[2,4],[2,5],[3,5]];
#gensimplexB_n(8,50,[]);
