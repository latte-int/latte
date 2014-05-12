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



#weightedBoxCount([4, 10], [5, 13], 1, [1,2], 10);

