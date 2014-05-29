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
end;


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
	
end;

sumPolynomial:=proc(lb, ub, poly, point)
	local current, ans, i;
	
	current := lb;
	ans:=0;
	while ( myLessThanEq(current, ub)) do
		#print(current);
		
		ans := ans + eval(subs({seq(x[i]=current[i],i=1..n)}, poly));
		
		current[1] := current[1] + 1;
	
		i := 1;
		while (current[i] > ub[i] and i < nops(lb)) do
			current[i] := lb[i];
			current[i+1] := current[i+1] + 1;
			i := i + 1;
		od;
	
	od;
	return ans;
end;




#lb, up: lower and upperbound
#np[1] = number of vars,
#np[2] = polynomial in x[i], i =1..n
#k to take the k-norm
weightedPolynomialLP:=proc(lb, ub, np, k)
	local current, ans, i;
	ASSERT(nops(lb) = nops(ub));

	current := lb;
	n:=np[1]:
	p:=np[2]:
	pk:=(p)^k:
	
	ans := 0;
	minVal:=eval(subs({seq(x[i]=lb[i],i=1..n)}, p));
	maxVal:=minVal;
	while ( myLessThanEq(current, ub)) do
		#print(current);
		
		ans := ans + eval(subs({seq(x[i]=current[i],i=1..n)}, pk));
		t:=eval(subs({seq(x[i]=current[i],i=1..n)}, p));
		minVal:=min(minVal, t);
		maxVal:=max(maxVal, t);
		
		current[1] := current[1] + 1;
	
		i := 1;
		while (current[i] > ub[i] and i < nops(lb)) do
			current[i] := lb[i];
			current[i+1] := current[i+1] + 1;
			i := i + 1;
		od;
	
	od;
	printf("minValue of f %E\n",minVal);
	printf("maxValue of f %E\n", maxVal);
	printf("sum f^k       %E\n", ans);
	printf("||f||_k       %E\n", evalf(ans^(1/k)));
	printf(" k            %E\n", k);
	
end;

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
end;

sPolynomial:=proc(lb, ub, fx, k)
	local current, ans, i;
	ASSERT(nops(lb) = nops(ub));

	current := lb;
	
	ans := 0;
	while ( myLessThanEq(current, ub)) do
		#print(current);
		
		ans := ans + subs({seq(x[i]=current[i],i=1..nops(lb))}, (fx)^k);
		
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

polynomialRange:=proc(lb, ub, fx)
	local current, tmax, tmin, i;
	ASSERT(nops(lb) = nops(ub));

	current := lb;
	tmax:= subs({seq(x[i]=current[i],i=1..nops(lb))}, fx);
	tmin:=tmax;
	
	ans := 0;
	while ( myLessThanEq(current, ub)) do
		#print(current);
		
		t := subs({seq(x[i]=current[i],i=1..nops(lb))}, fx);
		tmax:=max(tmax,t);
		tmin:=min(tmin,t);
		
		current[1] := current[1] + 1;
	
		i := 1;
		while (current[i] > ub[i] and i < nops(lb)) do
			current[i] := lb[i];
			current[i+1] := current[i+1] + 1;
			i := i + 1;
		od;
	
	od;
	print("min=", tmin);
	print("max=", tmax);
end:
#weightedBoxCount([4, 10], [5, 13], 1, [1,2], 10);

