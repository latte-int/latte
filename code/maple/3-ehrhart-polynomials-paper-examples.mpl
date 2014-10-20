# Examples from paper "Three Ehrhart Quasi-Polynomials"

read("RealBarvinok-mars-exemples-2014-03-10.mpl"):

############################################################

"##### Example 1.3 - Numbers of lattice points in a semi-rational rectangle dilated by a real number t. Figure 3. #####";

vertices := [[0, 0], [sqrt(2), 0], [sqrt(2), 1], [0, 1]];
simple_vertex_cones := [[[1, 0], [0, 1]], 
                        [[-1, 0], [0, 1]], 
                        [[-1, 0], [0, -1]], 
                        [[1, 0], [0, -1]]]:
L := []:       # we are counting lattice points.
ell := [1, 0]: # an arbitrary linear form.
M := 0:        # we compute the sum of ell^M.
reg := [1, 1]: # a regular linear form.
S := expand(add(ttruncatedSL('t', vertices[i], simple_vertex_cones[i], L, ell, reg, M), i=1..nops(vertices))):
E2 := coeff(S, t, 2);
E1 := coeff(S, t, 1);
E0 := factor(coeff(S, t, 0));

############################################################

"##### Examples 2.8/2.10/2.15/2.29/2.34/2.36 - Parametric polytope in dimension 2 from 4 hyperplanes (lines) #####";

# - chamber tau2 (rectangle) "was computed with our Maple program". Figure 7, 8.

## from file forNicole-2014-02-05.mpl:

with(combinat):

# Computation of the parametric vertices from 4 lines:
# 
"Parametric hyperplane arrangement (shown in Example 2.8) given by:";
mu:=[-x1-b[1],-x2-b[2],x1+x2-b[3],-x1+x2-b[4]];

Vertices:=proc(mu, Bases) local VV,d,ss,sv; 
    VV:={};
    for d from 1 to 6 do
        ss:=solve({mu[Bases[d][1]],mu[Bases[d][2]]},{x1,x2}); print(ss);
        sv:=subs(ss,[x1,x2]);
        VV:={op(VV),s[Bases[d]]=sv};
    od;
    VV;
end:

Bases := choose(4,2):

"Parametric vertices (shown in Example 2.10):";

Vertices(mu, Bases);

# Corresponding vertex cones:
#
# (We do not include code for computing these cones; let us just write
# them down by hand.  The vertices below have been computed above.)
#
# The cone $C_{[1,2]}$  is the cone with edges $(1,0),(0,1)$,
# The cone $C_{[1,3]}$  is the cone with edges $(0,-1),(1,-1)$,
# The cone $C_{[1,4]}$ is the cone with edges $(0,-1),(1,1)$,
# The cone $C_{[2,3]}$ is the cone with edges $(-1,0),(-1,1)$,
# The cone $C_{[2,4]}$ is the cone with edges $(1,0),(1,1)$,
# The cone  $C_{[3,4]}$ is the cone with edges $(-1,-1),(1,-1)$:

C[[1,2]]:=[[0,1],[1,0]]: s[[1,2]]:=[-b[1],-b[2]]:
C[[1,3]]:=[[0,-1],[1,-1]]: s[[1,3]]:=[-b[1],b[1]+b[3]]:
C[[1,4]]:=[[0,-1],[1,1]]: s[[1,4]]:=[-b[1], -b[1]+b[4]]:
C[[2,3]]:=[[-1,0],[-1,1]]: s[[2,3]]:=[b[2]+b[3],-b[2]]:
C[[2,4]]:=[[1,0],[1,1]]: s[[2,4]]:=[-b[2]-b[4],-b[2]]:
C[[3,4]]:=[[-1,-1],[1,-1]]: s[[3,4]]:=[-(1/2)*b[4]+(1/2)*b[3], (1/2)*b[4]+(1/2)*b[3]]:

# Admissible chambers

"Bases of admissible chambers (shown in Example 2.15):";

Taus:=[[],[],[]]:
Taus[1] := [[2,3],[2,4],[3,4]];
Taus[2] := [[1,2],[1,4],[2,3],[3,4]];
Taus[3] := [[1, 2], [1, 3], [2, 3]]; # according to paper
#Taus[3]:= [[1,2],[1,3],[1,4]];  # according to file forNicole-2014-02-05.mpl

# Subspaces

"Subspaces L to use (in Example 2.29):";

V_space := [[0,1],[1,0]];                     #full space 
L_space := [[0,1]];
O_space := [];                                #trivial space;

# Computations

for L in [V_space, L_space, O_space] do
    for bas in Bases do
        E[L][bas] := coeff(series(eval(subs({EXP=exp,TODD=Todd},subs(T=1,tfunction_SL(1,s[bas],C[bas],L_space,[t*x1,t*x2])))), t=0), t, 0);
    od;
    for tau from 1 to nops(Taus) do
        E[L][tau] := simplify(add(E[L][bas], bas in Taus[tau]));
    od;
od:

# Output

for tau from 1 to nops(Taus) do

    "Chamber tau", tau;

    "Volume";

    factor(E[V_space][tau]);

    "Sum";

    E[O_space][tau];

    "Intermediate sum";

    E[L_Space][tau];

od:

# Example 5.4 - 4-dimensional simplex with vertices
# [4,6,4,3],[5,7,9,1],[5,7,3,7],[6,8,3,9],[2,1,8,0]. Table 1.

# Example 5.5 - triangle with vertices [1, 1], [1, 2], [2, 2]. Figure 11.

mars := [[1,1], [1, 2], [2, 2]];
#mars0:=FBbonfrac(t,mars,0,[1,5,7,2],0);
#mars1:=FBbonfrac(t,mars,1,[1,5,7,2],0);
#mars2:=FBbonfrac(t,mars,2,[1,5,7,2],0):
#mars3:=FBbonfrac(t,mars,3,[1,5,7,2],0):
#mars4:=FBbonfrac(t,mars,4,[1,5,7,2],0):
#### plot([mars4,mars3,mars2,mars1,mars0], t=0.6...0.7,color=[black,red,blue,violet,green], axis[1]=[gridlines=[30, thickness=1, subticks=false]],axis[2]=[gridlines=[6, thickness=1, subticks=false]],thickness=1);


# Example 5.6 - 3-dimensional simplex with vertices
# [0,1,1],[4,2,1],[1,1,2],[1,2,4]. Figures 12, 13.

#====================
# unrelated tests, moved here from RealBarvinok-mars-exemples-2014-03-10.mpl

#### #testfrac := [[2, 1, 4], [2, 8, 7], [4, 3, 2], [2, 8, 2]]:
#### #mars:=[[6, 9, 5, 1], [3, 5, 4, 0], [7, 4, 9, 1], [1, 3, 7, 2], [8, 9, 1, 7]];

#### marsdim3:=[[1,1,1],[4,2,1],[1,1,2],[2,2,2] ];
#### M:=[[0,1,1],[4,2,1],[1,1,2],[1,2,4] ];
#### mars0:=FBbonfrac(t,marsdim3,0,[1,5,7,2],0);
#### M0:=FBbonfrac(t,M,0,[1,5,7,2],0);
#### M1:=FBbonfrac(t,M,1,[1,5,7,2],0):
#### M2:=FBbonfrac(t,M,2,[1,5,7,2],0):
#### M3:=FBbonfrac(t,M,3,[1,5,7,2],0):
#### mars1:=FBbonfrac(t,marsdim3,1,[1,5,7,2],0);M1:=FBbonfrac(t,M,1,[1,5,7,2],0);
#### mars2:=FBbonfrac(t,marsdim3,2,[1,5,7,2],0):
#### mars3:=FBbonfrac(t,marsdim3,3,[1,5,7,2],0):eval(subs(t=9/2,mars3));
#### L:=[seq(seq(k+n/10,n=0..9),k=1..3)];

#### plot([mars3,mars2,mars1,mars0], t=1...4.5,color=[black,red,blue,green], axis[1]=[gridlines=[10, thickness=1, subticks=false]],axis[2]=[gridlines=[60, thickness=1, subticks=false]],thickness=2, adaptive=20,sample=[0,1/2,1,3/2,2,5/2,3,7/2,4,9/2]);
#### plot([M3,M2,M1,M0], t=1.3...3.9,color=[black,red,blue,green], axis[1]=[gridlines=[20, thickness=1, subticks=false]],axis[2]=[gridlines=[60, thickness=1, subticks=false]],thickness=2, adaptive=20,sample=L);
#### L:=[seq(seq(k+n/12,n=0..11),k=0..1)];
#### plot([M2,M1], t=1.3...2.9,color=[red,blue,green], axis[1]=[gridlines=[20, thickness=1, subticks=false]],axis[2]=[gridlines=[60, thickness=1, subticks=false]],thickness=2, adaptive=20,sample=L);
#### eval(subs(t=5/2-0.1,M3));
#### plot([M3], t=1.6...2.1,color=[black,red,blue,green], axis[1]=[gridlines=[10, thickness=1, subticks=true]],axis[2]=[gridlines=[60, thickness=1, subticks=false]],thickness=2, adaptive=10,sample=L);

