polynomial_to_sparsepoly := proc(p, dimension)
  local variables,coefficients, monomials, exponents; 
  variables := [ seq(x[i], i=1..dimension) ];     
  #print (variables);     
  coefficients := coeffs(p, variables, 'monomials');
  #print (coefficients);
  #print (monomials);
  exponents := map( monomial -> map( variable -> degree(monomial, variable),
                            variables),
                     [monomials]);
  #print(exponents);
  zip((c,m) -> [c,m],
      [coefficients], 
      exponents);
end:

lattice_random_simplex:=proc(d,N) local R,U;
  R := rand(N):
  U:=proc()[seq(R(),i=1..d)] end proc:
  [ seq(U(), i=1..d+1) ];
end:

random_linearform_given_degree_dimension_maxcoef_componentmax_maxterm:=proc(m,d,maxcoef,componentmax,maxterm)
local temp,l,termcount,i,j,R;
R:=rand(maxterm);
termcount:=R()+1;
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

random_sparse_homogeneous_polynomial_with_degree:=proc(N,d,M,r) 
  local p, R;
  ## Give up if too large polynomials requested
  if (r > 500000) then
    error "Too large a polynomial requested"
  fi;
  R := rand(N);
  p := randpoly([ seq(x[i], i=1..d) ], 
                homogeneous, degree = M, terms = r, coeffs = proc() R() + 1; end);
  polynomial_to_sparsepoly(p, d);
end:

test_integration:=proc(polyCount, bigConstant, numTerms, dimension, myDegree, decomposing)
  global filename:
  local errors, wrong:
  local myMonomials, mySimplices, myLinForms, mapleLinForms, myResults, mapleResults:
  local curForms, curTerm, curSet:
  local myIndex, formIndex, i, j:
  local myTime, temp, intTime, L:
  local inputFile, outputFile, errorFile:
  
  print(randomGen(bigConstant, dimension, myDegree, numTerms)):
  #get polynomials
  myTime:=0:
  intTime:=0:
  for myIndex from 1 to polyCount do
    mySimplices[myIndex]:=lattice_random_simplex(dimension, bigConstant):
    if decomposing = 1 then
      myMonomials[myIndex]:=random_sparse_homogeneous_polynomial_with_degree(bigConstant, dimension, myDegree, numTerms);
    else
      #print(myDegree, dimension, bigConstant, bigConstant, numTerms);
      mapleLinForms[myIndex]:=random_linearform_given_degree_dimension_maxcoef_componentmax_maxterm(myDegree, dimension, bigConstant, bigConstant, numTerms);
      #print(random_linearform_given_degree_dimension_maxcoef_componentmax_maxterm(myDegree, dimension, bigConstant, bigConstant, numTerms));
    end if:
  od:

  #write to file
  if decomposing = 1 then
    inputFile:=fopen("integration/check_in.tmp",WRITE,TEXT):
    for i from 1 to polyCount do
      writeline(inputFile, convert(myMonomials[i], string));
      writeline(inputFile, StringTools[DeleteSpace](convert(mySimplices[i], string)));
    od:
    close(inputFile):

    system("./BurstTest integration/check_in.tmp integration/check_out.tmp polynomial"):

    outputFile:=fopen("integration/mark.txt",READ,TEXT):
    curTerm:= readline(outputFile):
    if (curTerm = "1") then
      print("Integration timed out.");
      close(outputFile):
      return -1:
    end if:
    close(outputFile):
  else
    inputFile:=fopen("integration/check_in.tmp",WRITE,TEXT):
    for i from 1 to polyCount do
      writeline(inputFile, convert(mapleLinForms[i], string));
      writeline(inputFile, StringTools[DeleteSpace](convert(mySimplices[i], string)));
    od:
    close(inputFile):

    system("./BurstTest integration/check_in.tmp integration/check_out.tmp linear"):

    outputFile:=fopen("integration/mark.txt",READ,TEXT):
    curTerm:= readline(outputFile):
    if (curTerm = "1") then
   #   print("Integration timed out.");
      close(outputFile):
      return -1:
    end if:
    close(outputFile):
  end if:
end:

draw_table5:=proc()
global filename:
local benchmarks, result:
local myDim, myDegree, timedOut:
benchmarks:=fopen("integration/burstbench1.txt",APPEND,TEXT):
writeline(benchmarks, "Average integration time of one power of a linear form over random simplices (average time over 50 random forms)"):
writeline(benchmarks, ""):
writeline(benchmarks, "                                            Degree                                     "):
writeline(benchmarks, "      _________________________________________________________________________________"):
writeline(benchmarks, "   n  |       2        10        20         50       100       300      1000"):
writeline(benchmarks, "_______________________________________________________________________________________"):
fprintf(benchmarks, "%4.4s  |", "1000"):
close(benchmarks):
for myDim from 1000 to 1000 do
  timedOut:= 0:
  for myDegree from 2 to 1000 do
    #samplesize, bigConstant, numTerms, dimension, myDegree, decomposing
    if timedOut = 0 then
      result:= test_integration(50, 100, 1, myDim, myDegree, 0):
      if result = -1 then
       timedOut:= 1:
      end if:
    else
      benchmarks:=fopen("integration/burstbench1.txt",APPEND,TEXT):
      fprintf(benchmarks, "%10s", "--"):
      close(benchmarks):
    end if:
    if result = -1 then
      myDegree:=1000:
    end if:
    if myDegree = 2 then
      myDegree:= 9:
    elif myDegree = 10 then
      myDegree:=19:
    elif myDegree = 20 then
      myDegree:=49:
    elif myDegree = 50 then
      myDegree:=99:
    elif myDegree = 100 then
      myDegree:=299:
    elif myDegree = 300 then
      myDegree:=999:
    end if:
  od:
  if myDim = 10 then
      myDim:= 19:
  elif myDim = 20 then
    myDim:=49:
  elif myDim = 50 then
    myDim:=99:
  elif myDim = 100 then
    myDim:=299:
  elif myDim = 300 then
    myDim:=399:
  elif myDim = 400 then
    myDim:=999:
  end if:
  benchmarks:=fopen("integration/burstbench1.txt",APPEND,TEXT):
  writeline(benchmarks, ""):
  fprintf(benchmarks, "%4.4s  |", convert(myDim + 1, string)):
  close(benchmarks):
od:
end:

draw_table10:=proc()
global filename:
local benchmarks, result:
local myDim, myDegree, timedOut:
benchmarks:=fopen("integration/burstbench1.txt",APPEND,TEXT):
writeline(benchmarks, "Average integration time of a random monomial of prescribed degree by decomposition into linear forms (average time over 50 random forms)"):
writeline(benchmarks, ""):
writeline(benchmarks, "     ______________________________________________________________________________________________________________"):
writeline(benchmarks, "  n  |       1         2         5         10        20        30        40        50         100       200      300"):
writeline(benchmarks, "___________________________________________________________________________________________________________________"):
fprintf(benchmarks, "%3.3s  |", "6"):
close(benchmarks):
for myDim from 6 to 50 do
  timedOut:= 0:
  for myDegree from 1 to 300 do
    #samplesize, bigConstant, numTerms, dimension, myDegree, decomposing
    if timedOut = 0 then
      result:= test_integration(50, 1000, 1, myDim, myDegree, 1):
      if result = -1 then
       timedOut:= 1:
      end if:
    else
      benchmarks:=fopen("integration/burstbench1.txt",APPEND,TEXT):
      fprintf(benchmarks, "%10s", "--"):
      close(benchmarks):
    end if:
    if myDegree = 2 then
      myDegree:= 4:
    elif myDegree = 5 then
      myDegree:=9:
    elif myDegree = 10 then
      myDegree:=19:
    elif myDegree = 20 then
      myDegree:=29:
    elif myDegree = 30 then
      myDegree:=39:
    elif myDegree = 40 then
      myDegree:=49:
    elif myDegree = 50 then
      myDegree:=99:
    elif myDegree = 100 then
      myDegree:=199:
    elif myDegree = 200 then
      myDegree:=299:
    end if:
  od:
  if myDim = 8 then
      myDim:= 9:
  elif myDim = 10 then
    myDim:=14:
  elif myDim = 15 then
    myDim:=19:
  elif myDim = 20 then
    myDim:=29:
  elif myDim = 30 then
    myDim:=39:
  elif myDim = 40 then
    myDim:=49:
  end if:
  benchmarks:=fopen("integration/burstbench1.txt",APPEND,TEXT):
  writeline(benchmarks, ""):
  fprintf(benchmarks, "%3.3s  |", convert(myDim + 1, string)):
  close(benchmarks):
od:
end:

draw_table5():
draw_table10():


