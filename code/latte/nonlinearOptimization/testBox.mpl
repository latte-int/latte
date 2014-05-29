read("weightedCountingBox.mpl");






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
		
	#fptr:=fopen("boxTesting.hrep", WRITE, TEXT);
	#fprintf(fptr, "%d %d\n", 2*dim, dim+1);
	#for i from 1 to dim do:
	#	fprintf(fptr, "%d ", u[i]);
	#	for j  from 1 to dim do:
	#		if j = i then
	#			fprintf(fptr, "-1 ");
	#		else
	#			fprintf(fptr, "0 ");
	#		end;
	#	end;
	#	fprintf(fptr, "\n");
	#	
	#	fprintf(fptr, "%d ", -1*l[i]);
	#	for j  from 1 to dim do:
	#		if j = i then
	#			fprintf(fptr, "1 ");
	#		else
	#			fprintf(fptr, "0 ");
	#		end;
	#	end;
	#	fprintf(fptr, "\n");		
	#end;
	#fclose(fptr);
	
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
	

end;



testPolynomial:=proc()

	np:=getPolynomial("poly2.poly");
	print(np);
	ans:=weightedPolynomialLP([4, 10], [5, 13], np, 5);
	
	evalf(ans);

end;


#weightedBoxCount([4, 10], [5, 13], 1, [1,2], 10);
#testRandomForm(5);
testPolynomial();
