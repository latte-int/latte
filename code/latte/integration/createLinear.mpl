with(LinearAlgebra):
#This function takes in degree, dimension and some bounds, creates a random linear form
random_linearform_given_degree_dimension_maxcoef_componentmax_maxterm:=proc(m,d,maxcoef,componentmax,maxterm)
local temp,l,termcount,i,j,R;
R:=rand(maxterm);
termcount:=R();
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
random_linform_given_multiplicity:=proc(m,d,maxcoef,vertexmax,k,multiplicity)
local l,temp,i,j,str,t,simplex,mat,R,singular;
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
  print(mat);
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
str:=[op(str),simplex];
end:

#random_linearform_given_degree_dimension_maxcoef_componentmax_maxterm(4,3,50,6,10);
#random_linform_given_multiplicity(4,4,50,30,1,[5]);


