$include "m-knapsack.mpl";
 
randomize(12345);
#randomize();

  
testSame:=proc(a, b, numTValues)
 	local dif;
 	
 	for i from 0 to numTValues do:
 		dif:=eval(subs({MOD=latteMod, T=i, t=i}, a - b));
 		
 		if ( dif <> 0) then
            print("Error: Results are not the same");
	 		print("Maple: ", a);
 			print("latte: ", b);
 	
 			print("ops");
 			print("T=",i,dif);
 			quit;
 		fi;
 	od;
 end:
  
 
 
nextWeeklyIncreasingSequence:=proc(startingList, upperLimit)
	local L, n, N, i;
	L:=startingList;
	if( L[1] > upperLimit) then
		return [];
	fi;
	n:=nops(L);
	N:=n;
	while ( n>0 and L[n] >= upperLimit ) do
		n:= n-1;
	end;
	
	if (n = 0) then
		return [];
	fi;
	L[n] := L[n]+1;
	for i from n+1 to N do
		L[i] := L[n];
	od;
	
	
	if( L[n] > upperLimit and n > 1) then
		L[n-1]:=L[n-1]+1;
		for i from n to N do
			L[i] = L[n-1];
		od;
	fi;
	return L;
end:	

nextRandomSequence:=proc(L, k)
	N:=L;
	
	r:=rand(k);
	for i from 1 to nops(L) do
		N[i] := 1 + r();
	end;
	

	return N;
end:


#@param termType: if 1, then only coptue the kth coeff, if 2, then comptue the first k+1 terms.
cmpLatteNminusK:=proc(L, k, termType)
	local typeStr;
	
	if (termType = 1) then 
		#call latte	
		system(cat("../latte/top-ehrhart-knapsack -f knapsackEquation.latte -o knapsackEquation.latte.topKehrhart.mpl --gcd-polynomial 0 -k ",k+1, " > /dev/null 2>&1"));
		read("knapsackEquation.latte.topKehrhart.mpl");
	
		#call maple
		mapleNminusK:=coeff_Nminusk_knapsack(L, t, T, k);

		testSame(coeff(mapleNminusK, T, nops(L)-1-k), parse(cat("coeff",nops(L)-1,"minus",k), statement), ilcm(op(L)) );
		parse(cat("coeff",nops(L)-1,"minus",k,":= '","coeff",nops(L)-1,"minus",k,"'"), statement); #clear the variable we just read in	
	end;
		
	if (termType = 2) then
		#call latte	
		system(cat("../latte/top-ehrhart-knapsack -f knapsackEquation.latte -o knapsackEquation.latte.topKehrhart.mpl --gcd-polynomial 0 --all-k ",k+1, " > /dev/null  2>&1"));
		read("knapsackEquation.latte.topKehrhart.mpl");
	
		#call maple
		mapleNminusK:=knapsackKTerms(L, t, T, k);

		testSame(mapleNminusK, eval(parse("topKPolynomial", statement)), ilcm(op(L)) );
		topKPolynomial:='topKPolynomial';#clear the variable we just read in
	end;
end:

testLatteNminusK:=proc(termType)

	for dim from 5 to 5 do
		L:= [ seq(1, i = 1..dim) ];
		counter = 1;
		
		counter:=1;
		while (counter < 20) do
			if ( L = []) then
				break;
			end;
			if ( igcd(op(L)) = 1) then
			
				print("start testing", L);
				#write knapsack to file.
				knapFP:=fopen("knapsackEquation.latte", WRITE, TEXT);
				fprintf(knapFP,"%d ", nops(L));
				for i from 1 to nops(L) do
					fprintf(knapFP,"%d ", L[i]);
				od;
				fprintf(knapFP, "\n");
				fclose(knapFP);
				
				for k from 0 to nops(L)-1 do					
					cmpLatteNminusK(L, k, termType);				
				od;			
				
				print("Finished testing", L);
				counter:=counter +1;
			end;

			#L:=nextWeeklyIncreasingSequence(L, 10);
			L:=nextRandomSequence(L, 10);			
			f:=fopen("LkFile.txt", APPEND, TEXT);
			fprintf(f,"%a, %d\n", L, k);
			fclose(f);
		end; #while	
	od;


end:
 
 
testLatteNminusK(1): #just the N-k term
#testLatteNminusK(2): #the N, N-1, ..., N-k term

