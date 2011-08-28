with(LinearAlgebra):
#This function takes in degree, dimension and some bounds, creates a random linear form in list form.
random_linearform_given_degree_dimension_maxcoef_componentmax_maxterm:=proc(m,d,maxcoef,componentmax,maxterm, rationalCoeff)
	local temp,l,termcount,i,j,R, R_denom;
	R:=rand(maxterm); #max number of powers of linear forms.
	R_denom:=rand(100);
	termcount:=R() + 1;
	R:=rand(componentmax); #now R will give the coeff. inside the linear form: (coeff.)^power.
	l:=[];
	for i from 1 to termcount do
		temp:=[];
			for j from 1 to d do
				
				temp:=[R(),op(temp)];
				#temp:=[1, op(temp)];
			od;
		temp:=[m,temp];
		R:=rand(maxcoef);
		if ( rationalCoeff = 0) then 
			temp:=[R()-trunc(maxcoef/2),temp];
		else
			temp:=[(R()-trunc(maxcoef/2))/ (R_denom() + 1),temp];
		fi;
		l:=[op(l),temp];
	od;
	#printf("Random Linear Form:\n\n");
	#print(l);
	l;
end:

random_linearform_given_degree_dimension_maxcoef_componentmax_maxterm_nonhomogen:=proc(m,d,maxcoef,componentmax,maxterm, rationalCoeff)
	local l;
	 
	l:=random_linearform_given_degree_dimension_maxcoef_componentmax_maxterm(m,d,maxcoef,componentmax,maxterm, rationalCoeff);
	l:=[op(l), [rand(maxcoef)()*(-1)^(rand(2)())/(rand(maxcoef)()+1), [0, [seq(1, i=1..d)] ]
	           ] 
	   ]; #add c(1+...+1)^0
	return l;
end;

# @parm pairwiseDegree: max power of 1 linear form
# @parm dimension: variable count
# @parm maxNumberTerms: in the sum (not including the constant.
# @return random non-homogenous sum of products of powers of linear forms in list notation.
random_linearFormProduct_nonhomogen:=proc(pairwiseDegree, dimension, maxNumberTerms)
	local answer;
	local numberOfTerms, randNum, oneProduct;

	numberOfTerms:= rand(1..maxNumberTerms)();
	randNum := rand(100);
	answer:=[];
		
	for i from 1 to numberOfTerms do
	
		 numberProducts:=ceil(stats[random, uniform[0.5, 1.3]]()*dimension + 1); #in one terms. this number is D : l_1^m_1 ... l_D^m_D
		 
		 oneProduct:=[];
		 for j from 1 to numberProducts do
		 	 oneProduct:=[op(oneProduct), [rand(pairwiseDegree)()+1, [seq(randNum(), k=1..dimension)]]];
		 end; #j
		 
		 oneProduct:=[ randNum()/(randNum()+1), oneProduct];
	
		answer:=[op(answer), oneProduct]; 
	end; #for i

	#now add the constent term.
	
	answer:=[op(answer), [randNum()+1,[[0,[seq(1, k=1..dimension)]]]]]; 
	
	return answer;
end;

# @parm plfList (product of linear form) list
# @parm dimension: var count
# @return: maxple expression for the plf.
convert_linearFormProductList_to_product_maple_expression:=proc(plfList, dimension)
	local answer;
	local term, onePLF, k;

	answer:=0;

	for term in plfList do  #term is [c, [[p, [l]], [p, [l]], ...]]
		mapleTerm:=term[1]; #coefficient.
		
		for onePLF in term[2] do
			mapleTerm:=mapleTerm* (sum(x[k]*onePLF[2][k], k=1..dimension))^onePLF[1]
		end; #onePLF
		
		answer:=answer + mapleTerm;
	
	end;

	return answer;
end;



#
# Given a linear form in list form [ [coeff, [power, [linear form]]] ], we convert it to an
# expression in the variables x[1],..., x[d].
#
#@parm l : the linear form in list form
#@parm d : number of variables
#@return : maple expression
convert_linearFormList_to_maple_expression:=proc(l, d)
	local lformList;
	local lformMaple;
	local coeff;
	local power;
	local ans, i;
	
	ans:=0;
	for lformList in l do
	
		lformMaple:=0;
		
		coeff:=lformList[1];
		power:=lformList[2][1];
		for i from 1 to d do
			lformMaple:= lformMaple + lformList[2][2][i]*x[i];
		end;
		lformMaple:= coeff*(lformMaple)^power;
		
		ans:=ans + lformMaple
	end; #for each linear form in the list

	return ans;
end:


#
# Given a linear form in list form [ [coeff, [power, [linear form]]] ], we convert it 
# to a product of powers of linear fomrs and then convert this to an 
# expression in the variables x[1],..., x[d] 
#
#@parm l : the linear form in list form
#@parm d : number of variables
#@return : maple expression
convert_linearFormList_to_product_maple_expression:=proc(l, d)
	local lformList;
	local lformMaple;
	local coeff;
	local power;
	local ans, i;
	
	ans:=1;
	for lformList in l do
	
		lformMaple:=0;
		
		coeff:=lformList[1];
		power:=lformList[2][1];
		for i from 1 to d do
			lformMaple:= lformMaple + lformList[2][2][i]*x[i];
		end;
		lformMaple:= coeff*(lformMaple)^power;
		
		ans:=ans * lformMaple
	end; #for each linear form in the list

	return ans;
end:

del_space:=proc(s)
local i,str;
str:="";
for i from 1 to length(s) do
if not s[i]=" " then str:=cat(str,s[i]);fi;
od;
RETURN(str);
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
str:=[[R()-trunc(maxcoef/2),str]];
str:=[str,simplex];
end:

# This test function tests each dimension dim, tests 3 kinds of multiplicities of inner products. First one is uniform, for example d=7, k=2, take 2 different innerproducts, with multiplicity 4 and 4. Second one is identical, i.e. there is only one inner product with multiplicity d+1. Third one is all distinct, i.e. all distinct innerproducts. We assume that results are correct.
generate_random_linear_with_random_multiplicity:=proc(dim,coefmax,vertexmax)
local k,p,i,j,t,l,d,R,m,s,myDegree,maxcoef,maxver;
local randomFormCount,num,str,timeout;
local inputFile,benchmarks2,markFile;
randomFormCount:=10;
  R:=rand(coefmax);
  maxcoef:=1+R();
  R:=rand(vertexmax);
  maxver:=1+R();
  k:=1;num:=dim-2;
  if num>randomFormCount then num:=randomFormCount;fi;  
  for i from 1 to num do
    R:=rand(dim-k-num+i);
    k:=k+1+R();
    p:=[];l:=dim+1;
    for j from 1 to k-1 do
      R:=rand(l-k+j);
      t:=1+R();
      p:=[op(p),t];
      l:=l-t;
    od;
    p:=[op(p),l];
    benchmarks2:=fopen("integration/benchmarks.txt",APPEND,TEXT);  
    writeline(benchmarks2,"");
    writeline(benchmarks2,del_space(convert(p,string)));
    fprintf(benchmarks2,"                                   ");
    close(benchmarks2);
    for myDegree from 1 to 300 do
      inputFile:=fopen("integration/multipoly.txt",WRITE,TEXT);
      for s from 1 to 50 do
        str:=random_linform_given_multiplicity(myDegree,dim,maxcoef,maxver,k,p);        
        writeline(inputFile,convert(str[1],string));
        writeline(inputFile,convert(str[2],string));
      od;
      close(inputFile);
      system("./integrate integration/multipoly.txt integration/results.txt");
      markFile:=fopen("integration/mark.txt",READ,TEXT);
      timeout:=fscanf(markFile,"%d");
      close(markFile);
      if timeout=[0] then
      if myDegree=2 then myDegree:=4;fi;
      if myDegree=5 then myDegree:=9;fi;
      if myDegree=10 then myDegree:=19;fi;
      if myDegree=20 then myDegree:=29;fi;
      if myDegree=30 then myDegree:=39;fi;
      if myDegree=40 then myDegree:=49;fi;
      if myDegree=50 then myDegree:=99;fi;
      if myDegree=100 then myDegree:=199;fi;
      if myDegree=200 then myDegree:=299;fi;
      else myDegree:=301;
      end if;
    od;
  od;
  p:=[dim+1];
  benchmarks2:=fopen("integration/benchmarks.txt",APPEND,TEXT):  
  writeline(benchmarks2,"");
  writeline(benchmarks2,del_space(convert(p,string)));
  fprintf(benchmarks2,"                                   ");
  close(benchmarks2);
  for myDegree from 1 to 300 do
  inputFile:=fopen("integration/multipoly.txt",WRITE,TEXT);
  for s from 1 to 50 do
    str:=random_linform_given_multiplicity(myDegree,dim,maxcoef,maxver,1,p);
    writeline(inputFile, convert(str[1],string));
    writeline(inputFile, convert(str[2],string));
  od;
  close(inputFile);
  system("./integrate integration/multipoly.txt integration/results.txt");
  markFile:=fopen("integration/mark.txt",READ,TEXT);
  timeout:=fscanf(markFile,"%d");
  close(markFile);
  if timeout=[0] then
  if myDegree=2 then myDegree:=4;fi;
  if myDegree=5 then myDegree:=9;fi;
  if myDegree=10 then myDegree:=19;fi;
  if myDegree=20 then myDegree:=29;fi;
  if myDegree=30 then myDegree:=39;fi;
  if myDegree=40 then myDegree:=49;fi;
  if myDegree=50 then myDegree:=99;fi;
  if myDegree=100 then myDegree:=199;fi;
  if myDegree=200 then myDegree:=299;fi;
  else myDegree:=301;
  end if;
  od;
  p:=[];
  for j from 1 to dim+1 do p:=[op(p),1];od;
  benchmarks2:=fopen("integration/benchmarks.txt",APPEND,TEXT):
  writeline(benchmarks2,"");
  writeline(benchmarks2,del_space(convert(p,string)));
  fprintf(benchmarks2,"                                   ");
  close(benchmarks2);
  for myDegree from 1 to 300 do
  inputFile:=fopen("integration/multipoly.txt",WRITE,TEXT);
  for s from 1 to 50 do
    str:=random_linform_given_multiplicity(myDegree,dim,maxcoef,maxver,dim+1,p);
    writeline(inputFile, convert(str[1],string));
    writeline(inputFile, convert(str[2],string));
  od;
  close(inputFile);
  system("./integrate integration/multipoly.txt integration/results.txt");
  markFile:=fopen("integration/mark.txt",READ,TEXT);
  timeout:=fscanf(markFile,"%d");
  close(markFile);
  if timeout=[0] then
  if myDegree=2 then myDegree:=4;fi;
  if myDegree=5 then myDegree:=9;fi;
  if myDegree=10 then myDegree:=19;fi;
  if myDegree=20 then myDegree:=29;fi;
  if myDegree=30 then myDegree:=39;fi;
  if myDegree=40 then myDegree:=49;fi;
  if myDegree=50 then myDegree:=99;fi;
  if myDegree=100 then myDegree:=199;fi;
  if myDegree=200 then myDegree:=299;fi;
  else myDegree:=301;
  end if;
  od;
end:
#generate_random_linear_with_random_multiplicity(2,10,8,20,20,"integration/randomForms.txt");

#local benchmarks, result:
#local myDim, myDegree, timedOut:
#benchmarks:=fopen("integration/benchmarks.txt",WRITE,TEXT):
#writeline(benchmarks, "Average integration time of a linear form with specified multiplicity with a sample size of 50"):
#writeline(benchmarks, ""):
#writeline(benchmarks, ""):
#writeline(benchmarks, "                                                                Degree                                                            "):
#close(benchmarks):
#for myDim from 70 to 1000 do
#benchmarks:=fopen("integration/benchmarks.txt",APPEND,TEXT):
#writeline(benchmarks,""):
#writeline(benchmarks,""):
#fprintf(benchmarks, "      %3.3s                        |       1         2         5         10        20        30        40        50        100       200       300",convert(myDim,string)):
#writeline(benchmarks,""):
#  close(benchmarks):
#  generate_random_linear_with_random_multiplicity(myDim,20,20):
#  if myDim=70 then myDim:=79:fi:
#  if myDim=80 then myDim:=89:fi:
#  if myDim=90 then myDim:=99:fi:
#  if myDim=100 then myDim:=199:fi:
#  if myDim=200 then myDim:=299:fi:
#  if myDim=300 then myDim:=399:fi:
#  if myDim=400 then myDim:=499:fi:
#  if myDim=500 then myDim:=599:fi:
#  if myDim=600 then myDim:=699:fi:
#  if myDim=700 then myDim:=799:fi:
#  if myDim=800 then myDim:=899:fi:
#  if myDim=900 then myDim:=999:fi:
#od:

