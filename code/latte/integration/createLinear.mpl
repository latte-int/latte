with(LinearAlgebra):
random_linearform_given_degree_dimension_maxcoef_maxterm:=proc(m,d,maxcoef,maxterm)
	random_linearform_given_degree_dimension_maxcoef_componentmax_maxterm(m,d,maxcoef,maxcoef,maxterm);
end:
#This function takes in degree, dimension and some bounds, creates a random linear form
random_linearform_given_degree_dimension_maxcoef_componentmax_maxterm:=proc(m,d,maxcoef,componentmax,maxterm)
local temp,l,termcount,i,j,R;
R:=rand(maxterm);
termcount:=R() + 1;
l:=[];
for i from 1 to termcount do
temp:=[];
for j from 1 to d do
R:=rand(componentmax);
temp:=[R(),op(temp)];
od;
temp:=[m,temp];
R:=rand(maxcoef);
temp:=[R()-trunc(maxcoef/2),temp];
l:=[op(l),temp];
od;
end:

#This function can take in a specified multiplicity and produce a desired simplex plus linear form. But it only produces one term.
#For example, if we want a linear form of power 4, dimension 3, maxcoef 10, vertexmax 10, k=2 (which means we want two different inner products), the multiplicity is a list [2,2] (which means one inner product with multiplicity 2 and the other 2 as well).
random_linform_given_multiplicity:=proc(m,d,maxcoef,vertexmax,k,multiplicity)
local l,temp,i,j,str,t,simplex,mat,R,singular,inputFile;
l:=[];str:=[];R:=rand(vertexmax);
if k=1 then for i from 1 to d do l:=[op(l),1];od;
else
for j from 1 to multiplicity[1]-1 do l:=[op(l),0];od;
for i from 2 to k do
  t:=R();
  for j from 1 to multiplicity[i] do l:=[op(l),t];od;
od;
end if;
singular:=1;
while (singular=1) do
  simplex:=[];
  for i from 1 to d do
    temp:=[];
    for j from 1 to d do
    if k=1 then 
      temp:=[op(temp),2*R()];
    else 
      temp:=[op(temp),R()];
    end if; 
    od;
    simplex:=[op(simplex),temp];
  od;
  mat:=convert(simplex,Matrix);
  #print(mat);
  if (Determinant(mat)<>0) then singular:=0; end if;
od;
temp:=[];
for i from 1 to d do temp:=[op(temp),0];od;
if k=1 then temp[1]:=1;temp[2]:=1;fi;
simplex:=[convert(Transpose(convert(temp,Vector)).mat/2,'list'),op(simplex)];
mat:=Determinant(mat).MatrixInverse(mat);
l:=convert(mat.convert(l,Vector),'list');
str:=[m,l];
R:=rand(maxcoef);
str:=[R()-trunc(maxcoef/2),str];
str:=[str,simplex];
#writeline(inputFile,convert(str[1],string));
#writeline(inputFile,convert([op(2..nops(str),str)],string));
end:

# This test function tests each dimension, from dim1 to dim2 and for each one, tests 3 kinds of multiplicities of inner products. First one is uniform, for example d=7, k=2, take 2 different innerproducts, with multiplicity 4 and 4. Second one is identical, i.e. there is only one inner product with multiplicity d+1. Third one is all distinct, i.e. all distinct innerproducts. We assume that results are correct.
generate_random_linear_with_random_multiplicity:=proc(dim1,dim2,degreemax,coefmax,vertexmax,input)
local k,p,i,j,t,l,d,R,m,maxcoef,maxver;
#number of linear forms in category 1
local uniformFormCount,num,str;
local inputFile;
uniformFormCount:=10;
inputFile:=fopen(input,WRITE,TEXT);
for d from dim1 to dim2 do
  R:=rand(degreemax);
  m:=R();
  R:=rand(coefmax);
  maxcoef:=R();
  R:=rand(vertexmax);
  maxver:=R();
  k:=1;num:=d-1;
  if num>uniformFormCount then num:=uniformFormCount;fi;  
  for i from 1 to num do
    R:=rand(d-k-num+i);
    k:=k+1+R();
    #print(k," different inner products with multiplicity"):
    p:=[];l:=d+1;
    for j from 1 to k-1 do
      R:=rand(l-1-k+j);
      t:=1+R();
      p:=[op(p),t];
      l:=l-t;
    od;
    p:=[op(p),l];
    #print(p);
    str:=random_linform_given_multiplicity(m,d,maxcoef,maxver,k,p);    
    #print(str);
    writeline(inputFile,convert(str[1],string));
    writeline(inputFile,convert(str[2],string));
  od;
  #print("dimension ",d, "identical products testing");
  p:=[d+1];
  #print(p);
  str:=random_linform_given_multiplicity(m,d,maxcoef,maxver,1,p);
  #print(str);
  #writeline(inputFile, str);
  writeline(inputFile, convert(str[1],string));
  writeline(inputFile, convert(str[2],string));
  #print("dimension ",d, "distinct products testing"):
  p:=[];
  for j from 1 to d+1 do p:=[op(p),1];od;  
  #print(p);
  str:=random_linform_given_multiplicity(m,d,maxcoef,maxver,d+1,p);
  #print(str);
  #writeline(inputFile, str);
  writeline(inputFile, convert(str[1],string));
  writeline(inputFile, convert(str[2],string));
od;
close(inputFile);
end:
#random_linform_given_multiplicity(2,2,1,10,2,[2,1]);
generate_random_linear_with_random_multiplicity(2,10,8,20,20,"integration/randomForms.txt");
