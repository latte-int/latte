read("/home/mathman/ProgrammingSoftware/sourceFiles/latte-4ti2-merged/dest/share/latte-int/Conebyconeapproximations_08_11_2010.mpl"):
Delta := [[0, 0, 0, 0]/1,
[0, 0, 0, 1]/2,
[0, 0, 1, 0]/2,
[0, 5, 5, 0]/2,
[5, 7, 1, 0]/2]:
Form := [0, 0, 0, 0]:
Exponent := 0:
k := 3:
#printTopEhrhartweightedPoly(n,Delta,Form,Exponent,2):

seed:=randomize();

cor:=printIncrementalEhrhartweightedPoly(n, N, Delta, Form, Exponent):


quit;


#my:=printIncrementalEhrhartweightedPoly2(n, Delta, Form, Exponent, 4):









quit;
p:=3;
for m from Exponent+4 to 0 by -1 do	
	theError:=expand(cor[m+1]-my[m+1]):
	for i from 0 to p-1 do
		printf("%a:%a=%a\t", m,i, eval(subs(MOD=modp, subs(n=i, theError))));	
	od;
end:

print("ok, here we go:");

printf("cor=%a\n", cor);
printf("my =%a\n", my);


quit;

Simplex:=tangentConesToSimplex(simpleCones, 4);


#cones:=SimplexToTangentCones(Simplex);


correct:=0;
#printTopEhrhartweightedPoly_real(n,N, Simplex, linearForms[1][2][2], linearForms[1][2][1],2);
#correct:=printIncrementalEhrhartweightedPoly_real(n, N, Simplex, linearForms[1][2][2], linearForms[1][2][1]):
 
my:=printIncrementalEhrhartPolynomial(n,N,simpleCones,linearForms,4, true, 7):

print("n=",n,"N=", N);
p:=3;
for d from 4+linearForms[1][2][1] to 0 by -1 do

	theError:=coeff(correct - my, N, d):
	
	for i from 0 to p do
		printf("\t%a:%a:%a",d,i, eval(subs({n=1,MOD=latteMod},theError))); 
	end:

end:

