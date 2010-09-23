with(linalg):
with(LinearAlgebra):
with(numapprox,laurent):
read("testingSL/forGregtests.mpl"); #load the ehrhart functions

#Input: simplexDim: the amb. dim of the simplex.
#Output: returns a list of simplexDim+1 vectors in R^(simplexDim).
create_random_simplex:=proc(simplexDim)
	local i, j, M, checkRankMatrix:

	do
		M:=RandomMatrix(simplexDim+1, simplexDim);
		
		checkRankMatrix:=RandomMatrix(simplexDim, simplexDim);
		for i from 1 to simplexDim do
			for j from 1 to simplexDim do
				checkRankMatrix[i, j] := M[i, j] - M[simplexDim+1, j];
			od;
		od;
		#print("got here");
		#print(checkRankMatrix, "=checkRandMatrix");
		
		if Rank(checkRankMatrix) = simplexDim then:
			break;
		fi;
	end do;
	#print(M, "=M");
	
	return(convert(M, listlist));
end:

#Input: n+1 points in R^n
#Output: the facet equations that define the simplex of the n+1 points.
simplex_to_hyperplanes:=proc(simplex)
local i, b_ax, equations, n_points;

	equations:= [];
	for i from 1 to nops(simplex) do:
		n_points:= subsop(i = NULL, simplex);
		b_ax:= facets_equation(n_points, simplex[i]); # 0 < b - a*x
		equations:=[ op(b_ax), op(equations)];
	od;
	return(equations);
end:

#Input:facet equations in maple syntax.
#Output: Writes to file fileName the latte-style facet equations.
write_facets_to_file:=proc(equations, fileName, simplexDim)
	local i, j, filePtr, lcmDenom, M:
	
	filePtr:=fopen(fileName, WRITE, TEXT);
	
	#print(equations, "=equations");
	#print(fileName, "=fileName");
	
	
	
	#sorry, I had to do a bunch of converting because I couldn't pull the elements in the origional structure.
	M:=convert(equations, Matrix);
	#print(M, "=M1");
	M:=convert(M, list);
	#print(M, "=M3");
	#print(simplexDim, "=simplexDim");
	M:=matrix(simplexDim+1, simplexDim+1, M);
	#print(M, "=M2");
	#print(M[1], "=M2[1,1]");
	
	#write the latte file.
	fprintf(filePtr,"%d %d\n", simplexDim+1, simplexDim+1);
	for i from 1 to simplexDim+1 do:
		lcmDenom:=1;
		
		for j from 1 to simplexDim+1 do:
			lcmDenom:=lcm(lcmDenom, denom(M[i, j]));
		od:

		for j from 1 to simplexDim+1 do:
			#print(M[i, j], "M[i, j]");
			fprintf(filePtr, "%d ", M[i, j]*lcmDenom);
		od:
		fprintf(filePtr, "\n");
	od:
	close(filePtr);
end:

#input:filename and list of simplices
#output: writes each simplex to a new line of the file in maple syntax.
write_simplex_to_file:=proc(simplexList, fileName)
	local filePtr, i:
	filePtr:=fopen(fileName,WRITE,TEXT):
	
	for i from 1 to nops(simplexList) do:	
		writeline(filePtr, convert(simplexList[i], string));
	od;	
	close(filePtr);
end:

#Input:List of n points of a simplex and the last point.
#Output: the Facet of the equation containing the n points with the sign such that the last point is in the halfspace. 
#Description: Feeding the n points of a n-simplex and the 1 point that
# is not in their facet it generates ONE equation of the simplex in LattE format
#for those points.
facets_equation:=proc(L,notinL)
	local x,aux1,aux2,M,i;
	M:=Matrix([op(1,L)]);
	for i from 2 to nops(L) do
		M:=stackmatrix(M,op(i,L));
	od:
	M:=stackmatrix(M,[x[j] $j=1..nops(L)]);
	M:=transpose(M);
	M:=stackmatrix(M,[1 $j=1..nops(L)+1]);
	aux1:=det(M):
	for i from 1 to nops(notinL) do
		aux1:=subs(x[i]=op(i,notinL),aux1):
	od;
	#print(sign(aux1),"hi I am here");
	aux2:=sign(aux1)*det(M);
	for i from 1 to nops(notinL) do
		aux2:=subs(x[i]=0,aux2):
	od;
	return(genmatrix({aux2*x[0]+sign(aux1)*det(M)},[x[j] $j=0..nops(L)]));
end:

#input:filename and list of simplices
#output: writes each simplex to a new line of the file in maple syntax.
write_simplex_to_file:=proc(simplexList, fileName)
	local filePtr, i:
	filePtr:=fopen(fileName,WRITE,TEXT):
	
	for i from 1 to nops(simplexList) do:	
		writeline(filePtr, convert(simplexList[i], string));
	od;	
	close(filePtr);
end:

test_sl_ehrhart:=proc(mydim, myfilename):
	local myi, myCC:=[], mysimplex:
	randomize():
	mysimplex:=create_random_simplex(mydim):
	myfileName:="testingSL/testingEhrhartDim9Test1Top3":
	write_simplex_to_file(mysimplex,cat(myfileName,".simp")):
	write_facets_to_file(simplex_to_hyperplanes(mysimplex),cat(myfileName,".latte"),mydim):

	for myi from 0 to 2 do:
		myCC:=[op(myCC),[coeff_dminusk_Eh(mysimplex,myi)]]:
	od:
	myCC;
end: