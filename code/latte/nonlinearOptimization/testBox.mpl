read("weightedCountingBox.mpl");


toMapleFromLatte:=proc(pxList)

	px:=0;
	for m in pxList do:
		oneTerm:=m[1];
		for i from 1 to nops(m[2]) do:
			oneTerm:= oneTerm* x[i]^(m[2][i]);
		od;
		
		px:= px + oneTerm;
	od;
	
	return px;
end:


testRandomForm:=proc(dim)

	randomize();
	r:=rand(-10..10):
	
	u:=[]:
	l:=[]:
	linForm:=[]:
	for i from 1 to dim do:
		x:=r():
		l:=[x, op(l)]:
		u:=[x+rand(10)()+1, op(u)]:
		linForm:=[(rand(10)()+1)*(-1)^(rand(0..1)()), op(linForm)]:
	end;

	linForm[rand(1..dim)()] := 0;
	linForm[rand(1..dim)()] := 0;

	fptr:=fopen("boxTesting.box", WRITE, TEXT);
	fprintf(fptr, "%d\n", dim);
	for i from 1 to dim do:
		fprintf(fptr, "%d %d\n", l[i], u[i]);
	end;
	fclose(fptr);
		
	fptr:=fopen("boxTesting.linForms", WRITE, TEXT);

	c:=rand(3)()+1;
	M:=rand(10)();
	fprintf(fptr, "[[%d, [%d, [", c, M);
	for i from 1 to dim do
		fprintf(fptr, "%d", linForm[i]);
		if i <> dim then
			fprintf(fptr, ",");
		end;
	end;
	fprintf(fptr,"]]] ]\n");
	fclose(fptr);


	
	system("../boxOpt --count --boxFile=boxTesting.box --linear-forms=boxTesting.linForms");
	#system("./boxOpt boxTesting.hrep --linFile=boxTesting.linForms");
	print("Just called latte");
	cans:=weightedBoxCount(l, u, c, linForm, M):
	printf("Final count: %d\n", cans);
	

end:



testRandPolynomial:=proc(dim)


	randomize();
	r:=rand(-10..10):
	
	u:=[]:
	l:=[]:
	pxList:=[]:
	for i from 1 to dim do:
		x:=r():
		l:=[x, op(l)]:
		u:=[x+rand(10)()+1, op(u)]:
	end;
	x:='x';


	fptr:=fopen("boxTesting.box", WRITE, TEXT);
	fprintf(fptr, "%d\n", dim);
	for i from 1 to dim do:
		fprintf(fptr, "%d %d\n", l[i], u[i]);
	end;
	fclose(fptr);

	for j from 1 to 1 do:
		oneTerm:=[];
		for i from 1 to dim do:
			oneTerm:=[rand(0..2)(), op(oneTerm)]:
		end;
		pxList:= [[rand(1..5)(), oneTerm], op(pxList)];
	od;

	

	fptr:=fopen("boxTesting.poly", WRITE, TEXT);
	px:=toMapleFromLatte(pxList);
	fprintf(fptr,"%a\n", pxList);
	fclose(fptr);
	
	print(pxList);
	print(px);
	
	system("../boxOpt --count --boxFile=boxTesting.box --monomials=boxTesting.poly --k 1");
	print("Just called latte");
	cans:=weightedBoxCountExpression(l, u, px^1):
	printf("Final count: %d  (maple)\n", cans);	
	

end:


#weightedBoxCount([4, 10], [5, 13], 1, [1,2], 10);
#testRandomForm(5);
#testRandPolynomial(5);


fff:= x[1]^2 * x[6]^1 * x[9]^1 * x[14]^1 - 4 * x[1]^2 *  x[7]^1 * x[15]^1 * x[16]^1 + 3 * x[2]^1 * x[4]^1 * x[13]^1 * x[15]^2+3 *x[2]^1 * x[5]^2 * x[6]^1 * x[16]^1 - 8 * x[2]^1 * x[5]^2 * x[9]^1* x[19]^1;
a:=weightedBoxCountExpression([121, 168, 137, 83, 40, 101, 43, 65, 193, 83, 121, 121, 139, 124, 125, 128, 106, 176, 193, 78],  [121, 168, 137, 83, 40, 101, 43, 65, 193, 83, 121, 125, 140, 125, 125, 128, 106, 176, 193, 78], (fff - 16275962210)^3);
evalf(16275962210  + a^(1/3));
maxValueExpression([121, 168, 137, 83, 40, 101, 43, 65, 193, 83, 121, 121, 139, 124, 125, 128, 106, 176, 193, 78],  [121, 168, 137, 83, 40, 101, 43, 65, 193, 83, 121, 125, 140, 125, 125, 128, 106, 176, 193, 78], (fff - 0)^1);


