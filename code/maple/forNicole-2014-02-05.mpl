
with(combinat):

#  4 LINES ; VERTICES AND CONES;
# 
mu:=[-x1-b1,-x2-b2,x1+x2-b3,-x1+x2-b4];
Vertices:=proc() local VV,d,ss,sv; VV:={};
for d from 1 to 6 do
ss:=solve({mu[choose(4,2)[d][1]],mu[choose(4,2)[d][2]]},{x1,x2}); print(ss);
sv:=subs(ss,[x1,x2]);
VV:={op(VV),s[choose(4,2)[d]]=sv};
od;
VV;
end:
Vertices();
# The cone $C_{[1,2]}$  is the cone with edges $(1,0),(0,1)$,
# The cone $C_{[1,3]}$  is the cone with edges $(0,-1),(1,-1)$,
# The cone $C_{[1,4]}$ is the cone with edges $(0,-1),(1,1)$,
# The cone $C_{[2,3]}$ is the cone with edges $(-1,0),(-1,1)$,
#The cone $C_{[2,4]}$ is the cone with edges $(1,0),(1,1)$,
#The cone  $C_{[3,4]}$ is the cone with edges $(-1,-1),(1,-1)$:



L:=[[0,1]];
C12:=[[0,1],[1,0]]:s12:=[-b[1],-b[2]];
C13:=[[0,-1],[1,-1]]:s13:=[-b[1],b[1]+b[3]];
C14:=[[0,-1],[1,1]]:s14:=[-b[1], -b[1]+b[4]];
C23:=[[-1,0],[-1,1]]: s23:=[b[2]+b[3],-b[2]];
C24:=[[1,0],[1,1]]: s24:=[-b[2]-b[4],-b[2]];
C34:=[[-1,-1],[1,-1]]: s34:=[-(1/2)*b[4]+(1/2)*b[3], (1/2)*b[4]+(1/2)*b[3]];
# 
E12:=eval(subs({EXP=exp,TODD=Todd},subs(T=1,tfunction_SL(1,s12,C12,[],[t*x1,t*x2]))));
E13:=eval(subs({EXP=exp,TODD=Todd},subs(T=1,tfunction_SL(1,s13,C13,[],[t*x1,t*x2]))));
E14:=eval(subs({EXP=exp,TODD=Todd},subs(T=1,tfunction_SL(1,s14,C14,[],[t*x1,t*x2]))));
E23:=eval(subs({EXP=exp,TODD=Todd},subs(T=1,tfunction_SL(1,s23,C23,[],[t*x1,t*x2]))));
E24:=eval(subs({EXP=exp,TODD=Todd},subs(T=1,tfunction_SL(1,s24,C24,[],[t*x1,t*x2]))));
E34:=eval(subs({EXP=exp,TODD=Todd},subs(T=1,tfunction_SL(1,s34,C34,[],[t*x1,t*x2]))));




V12:=eval(subs({EXP=exp,TODD=Todd},subs(T=1,tfunction_SL(1,s12,C12,C12,[t*x1,t*x2]))));
V13:=eval(subs({EXP=exp,TODD=Todd},subs(T=1,tfunction_SL(1,s13,C13,C12,[t*x1,t*x2]))));
V14:=eval(subs({EXP=exp,TODD=Todd},subs(T=1,tfunction_SL(1,s14,C14,C12,[t*x1,t*x2]))));
V23:=eval(subs({EXP=exp,TODD=Todd},subs(T=1,tfunction_SL(1,s23,C23,C12,[t*x1,t*x2]))));
V24:=eval(subs({EXP=exp,TODD=Todd},subs(T=1,tfunction_SL(1,s24,C24,C12,[t*x1,t*x2]))));
V34:=eval(subs({EXP=exp,TODD=Todd},subs(T=1,tfunction_SL(1,s34,C34,C12,[t*x1,t*x2]))));



M12:=eval(subs({EXP=exp,TODD=Todd},subs(T=1,tfunction_SL(1,s12,C12,L,[t*x1,t*x2]))));
M13:=eval(subs({EXP=exp,TODD=Todd},subs(T=1,tfunction_SL(1,s13,C13,L,[t*x1,t*x2]))));
M14:=eval(subs({EXP=exp,TODD=Todd},subs(T=1,tfunction_SL(1,s14,C14,L,[t*x1,t*x2]))));
M23:=eval(subs({EXP=exp,TODD=Todd},subs(T=1,tfunction_SL(1,s23,C23,L,[t*x1,t*x2]))));
M24:=eval(subs({EXP=exp,TODD=Todd},subs(T=1,tfunction_SL(1,s24,C24,L,[t*x1,t*x2]))));
M34:=eval(subs({EXP=exp,TODD=Todd},subs(T=1,tfunction_SL(1,s34,C34,L,[t*x1,t*x2]))));





# CHAMBER TAU1;  VERTICES ARE [2,3],[2,4],[3,4];

V1:=simplify(coeff(series(V23+V24+V34,t=0),t,0)); factor(V1);;
S1:=simplify(coeff(series(E24+E23+E34,t=0),t,0));
simplify(S1-V1);
M1:=simplify(coeff(series(M23+M24+M34,t=0),t,0));
simplify(M1-V1);
## DOES NOT HAVE LINEAR TERMS ?? IS THIS CORRECT;

# CHAMBER TAU2; VERTICES ARE [1,2],[1,4],[2,3],[3,4];


Vol2:=simplify(coeff(series(V12+V14+V23+V34,t=0),t,0));
Sum2:=subs(frac=fpart,simplify(coeff(series(E12+E14+E23+E34,t=0),t,0)));
simplify(Sum2-Vol2);
Mix2:=simplify(coeff(series(M12+M14+M23+M34,t=0),t,0));
simplify(Mix2-Vol2);

# CHAMBER TAU3; VERTICES [1,2],[1,3],[1,4]




