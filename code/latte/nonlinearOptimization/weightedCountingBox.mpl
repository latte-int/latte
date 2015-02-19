with(LinearAlgebra):

myLessThanEq:=proc(lb, ub)
	local n, i;
	n := nops(lb);
	for i from 1 to n do
		if ( lb[i] > ub[i]) then
			return false;
		end;
	end;
	return true;
end:

#@param filename: fileName for a file that has a latte-style polynomial on the first line
#@return [n, poly], where n is the number of variables x[1]..x[n], and poly is a maple-style polynomail in x[i]
getPolynomial:=proc(fileName)
	local line, n, poly;
	line:=readline(fileName);
	line:=parse(line);
	n:=nops(line[1][2]);
	
	poly:=0;
	for i from 1 to nops(line) do
		poly:=poly+ line[i][1]*product(x[j]^(line[i][2][j]),j=1..n);
	end;
	
	return [n,poly];
	
end:



#@param lb, up: bounds
#@param c: coefficient
#@param l: linear form list, ex[1,3,4,2]
#@param p: power for the linear form, ex (1*x[1]+3*x[2]+4*x[3]+2*x[4])^p
#@return \sum_{x \in box} <x, l>^p
weightedBoxCount:=proc(lb, ub, c, l, p)
	local current, ans, i;
	ASSERT(nops(lb) = nops(ub));

	current := lb;
	
	ans := 0;
	while ( myLessThanEq(current, ub)) do
		#print(current);
		
		ans := ans + DotProduct(Vector(current), Vector(l))^p;
		
		current[1] := current[1] + 1;
	
		i := 1;
		while (current[i] > ub[i] and i < nops(lb)) do
			current[i] := lb[i];
			current[i+1] := current[i+1] + 1;
			i := i + 1;
		od;
	
	od;
	ans := ans*c;
end:

#@param fx: maple expression in x[i], and maybe s.
#@returns \sum_{x \in box} fx(x). fx(x) can be anything, power of LF or (f+x)^k
weightedBoxCountExpression:=proc(lb, ub, fx)
	local current, ans, i;
	ASSERT(nops(lb) = nops(ub));

	current := lb;
	
	ans := 0;
	while ( myLessThanEq(current, ub)) do
		#print(current);
		
		ans := ans + subs({seq(x[i]=current[i],i=1..nops(lb))}, fx);
		
		current[1] := current[1] + 1;
	
		i := 1;
		while (current[i] > ub[i] and i < nops(lb)) do
			current[i] := lb[i];
			current[i+1] := current[i+1] + 1;
			i := i + 1;
		od;
	
	od;
	ans := simplify(expand(ans));
end:


#@param fx: maple expression in x[i], and maybe s.
#@returns \sum_{x \in box} fx(x). fx(x) can be anything, power of LF or (f+x)^k
maxValueExpression:=proc(lb, ub, fx)
	local current, ans, i, maxVal, maxFval;
	ASSERT(nops(lb) = nops(ub));

	current := lb;
	
	maxVal := lb;
	maxFval := subs({seq(x[i]=lb[i],i=1..nops(lb))}, fx);
	while ( myLessThanEq(current, ub)) do
		#print(current);
		
		
		ans := subs({seq(x[i]=current[i],i=1..nops(lb))}, fx);
		if ( ans > maxFval) then
			maxFval := ans;
			maxVal := current;
		end;
		
		current[1] := current[1] + 1;
	
		i := 1;
		while (current[i] > ub[i] and i < nops(lb)) do
			current[i] := lb[i];
			current[i+1] := current[i+1] + 1;
			i := i + 1;
		od;
	
	od;
	return [maxFval, maxVal];
end:



#@param fx: maple expression in x[i] only
#@param lb, ub: bounds arrays
#@returns [min fx, max fx, avg fx, std dev fx] 
polynomialStats:=proc(lb, ub, fx)
	local current, tmax, tmin, i;
	ASSERT(nops(lb) = nops(ub));

	current := lb;
	fxMax:= subs({seq(x[i]=current[i],i=1..nops(lb))}, fx);
	fxMin:=fxMax;
	
	N:=product(ub[i] - lb[i]+1, i=1..nops(lb));
	
	
	sumk1:=0;
	sumk2:=0;
	while ( myLessThanEq(current, ub)) do
		#print(current);
		
		sumk1_part := subs({seq(x[i]=current[i],i=1..nops(lb))}, fx);
		sumk2 := sumk2 + subs({seq(x[i]=current[i],i=1..nops(lb))}, fx^2);
		sumk1 := sumk1 + sumk1_part;
		fxMax:=max(fxMax,sumk1_part);
		fxMin:=min(fxMin,sumk1_part);
		
		current[1] := current[1] + 1;
	
		i := 1;
		while (current[i] > ub[i] and i < nops(lb)) do
			current[i] := lb[i];
			current[i+1] := current[i+1] + 1;
			i := i + 1;
		od;
	
	od;
	
	fxAvg:= evalf(sumk1/N);
	fxSD:= evalf(sqrt(sumk2/N - (fxAvg^2)));
	
	
	print("min f", fxMin);
	print("max f", fxMax);
	print("avg f", fxAvg);
	print("SD  f", fxSD);
	
	return [fxMin, fxMax, fxAvg, fxSD];
end:


