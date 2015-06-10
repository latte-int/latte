read("weightedCountingBox.mpl"):


toMapleFromLatte:=proc(pxList)
	local px,m,oneTerm,i;
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
	local r,u,l,linForm,i,x,fptr,c,M,cans;
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



testRandPolynomialSum:=proc(dim)
	local r,u,l,pxList,i,t,fptr,j,oneTerm,px,cans,x; 
	randomize();
	r:=rand(-10..10):
	
	u:=[]:
	l:=[]:
	pxList:=[]:
	for i from 1 to dim do:
		t:=r():
		l:=[t, op(l)]:
		u:=[t+rand(10)()+1, op(u)]:
	end;



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
		pxList:= [[rand(1..5)()/rand(1..10)(), oneTerm], op(pxList)];
	od;

	

	fptr:=fopen("boxTesting.poly", WRITE, TEXT);
	px:=toMapleFromLatte(pxList);
	fprintf(fptr,"%a\n", pxList);
	fclose(fptr);
	
	print(pxList);
	print(px);
	
	#system("../boxOpt --count --boxFile=boxTesting.box --monomials=boxTesting.poly --k 1");
	system("../boxOpt --opt-lf --boxFile=boxTesting.box --monomials=boxTesting.poly --k 3");
	print("Just called latte");
	cans:=weightedBoxCountExpression(l, u, (px+s)^3):
	printf("Final count: %a  (maple)\n", sort(cans));	
	

end:

testRandPolynomialIntegration:=proc(dim)
	local r,u,l,psList,i,t,fptr,j,oneTerm,px,cans;
	randomize();
	r:=rand(-10..10):
	
	u:=[]:
	l:=[]:
	pxList:=[]:
	for i from 1 to dim do:
		t:=r():
		#l:=[t, op(l)]:
		#u:=[t+rand(10)()+1, op(u)]:
		l:=[-1, op(l)];
		u:=[1, op(u)];
	end;


	fptr:=fopen("boxTesting.box", WRITE, TEXT);
	fprintf(fptr, "%d\n", dim);
	for i from 1 to dim do:
		fprintf(fptr, "%d %d\n", l[i], u[i]);
	end;
	fclose(fptr);

	for j from 1 to 2 do:
		oneTerm:=[];
		for i from 1 to dim do:
			oneTerm:=[rand(0..2)(), op(oneTerm)]:
		end;
		pxList:= [[rand(1..5)()/rand(1..10)(), oneTerm], op(pxList)];
	od;

	

	fptr:=fopen("boxTesting.poly", WRITE, TEXT);
	px:=toMapleFromLatte(pxList);
	fprintf(fptr,"%a\n", pxList);
	fclose(fptr);
	
	print(pxList);
	print(px);
	
	
	system("../boxOpt --opt-cont-ns --boxFile=boxTesting.box --monomials=boxTesting.poly --k 3");
	print("Just called latte");
	print((px+3)^3, [seq(x[i] = l[i]..u[i],i=1..dim)]);
	cans:=int((px+s)^3, [seq(x[i] = l[i]..u[i],i=1..dim)]);
	printf("Final count: %a  (maple)\n", evalf(sort(cans)));	
	
end:

#weightedBoxCount([4, 10], [5, 13], 1, [1,2], 10);
#testRandomForm(5);
#testRandPolynomialSum(5);
testRandPolynomialIntegration(5);






