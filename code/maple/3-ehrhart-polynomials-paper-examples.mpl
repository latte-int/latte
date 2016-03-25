# Examples from paper "Three Ehrhart Quasi-Polynomials"

with(plots):
read("RealBarvinok-mars-exemples-2014-03-10.mpl"): # this also loads Conebycone...

plotfont := [Times, roman, 12]:

############################################################

printf("##### Example 1.3 - Numbers of lattice points in a semi-rational rectangle dilated by a real number t. Figure 3. #####\n\n");

vertices := [[0, 0], [sqrt(2), 0], [sqrt(2), 1], [0, 1]];
simple_vertex_cones := [[[1, 0], [0, 1]], 
                        [[-1, 0], [0, 1]], 
                        [[-1, 0], [0, -1]], 
                        [[1, 0], [0, -1]]]:
L := []:       # we are counting lattice points.
ell := [1, 0]: # an arbitrary linear form.
M := 0:        # we compute the sum of ell^M.
reg := [1, 1]: # a regular linear form.
S := expand(add(ttruncatedSL('t', 'T', vertices[i], simple_vertex_cones[i], L, ell, reg, M), i=1..nops(vertices))):
E2 := coeff(S, T, 2);
E1 := coeff(S, T, 1);
E0 := factor(coeff(S, T, 0));

# Plot Figure 3.
BasisPlots := Array(1 .. 3, 1 .. 1):
XRange := -1.5 .. 2.5:
BasisPlotsTicks := [-1=-1, 0=0, 1=1, 2=2]:
BasisPlots[1, 1] := plot(t, 't' = XRange, color = "Blue",
                         tickmarks = [BasisPlotsTicks, [-1=-1, 1=1]],
                         axis = [gridlines = [linestyle = dot]],
                         font = plotfont):
YRange := 0..1.1:
BasisPlots[2, 1] := plot([t-floor(t)], 't' = XRange, YRange, color = "Green",
                         discont = true, symbolsize = 15, symbol = solidcircle,
                         tickmarks = [BasisPlotsTicks, [1=1]],
                         axis = [gridlines = [linestyle = dot]],
                         font = plotfont):
BasisPlots[3, 1] := plot([sqrt(2)*t-floor(sqrt(2)*t)], 't' = XRange, YRange, color = "Red",
                         discont = true, symbolsize = 15, symbol = solidcircle,
                         tickmarks = [BasisPlotsTicks, [1=1]],
                         axis = [gridlines = [linestyle = dot]],
                         font = plotfont):
display(BasisPlots);

Ex13Volume := E2*t^2:
Ex13LatticePoints := evaluateEhrhart(S, 't'):
Ex13Plot := plot([Ex13LatticePoints, Ex13Volume], 't' = 0 .. 3.4,
                 color = ["Black", "Blue"],
                 discont = true, symbolsize = 15, symbol = solidcircle,
                 tickmarks = [[0 = 0, 1/sqrt(2) = "", 1 = 1, sqrt(2) = sqrt(2), 2 = 2,
                               (3/2)*sqrt(2) = "", 2*sqrt(2) = "", 3 = 3],
                              [0, 1, 2, 4, 6, 9, 12, 15, 20]],
                 axis = [gridlines = [linestyle = dot]],
                 font = plotfont):
display(Ex13Plot);

############################################################

printf("##### Examples 2.8/2.10/2.15/2.29/2.34/2.36 - Parametric polytope in dimension 2 from 4 hyperplanes (lines) #####\n\n");

# - chamber tau2 (rectangle) "was computed with our Maple program". Figure 7, 8.

## from file forNicole-2014-02-05.mpl:

with(combinat):

# Computation of the parametric vertices from 4 lines:
# 
printf("## Parametric hyperplane arrangement (shown in Example 2.8) given by:\n");
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

printf("## Parametric vertices (shown in Example 2.10):\n");

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

printf("## Bases of admissible chambers (shown in Example 2.15):\n");

Taus:=[[],[],[]]:
Taus[1] := [[2,3],[2,4],[3,4]];
Taus[2] := [[1,2],[1,4],[2,3],[3,4]];
Taus[3] := [[1, 2], [1, 3], [2, 3]]; # according to paper
#Taus[3]:= [[1,2],[1,3],[1,4]];  # according to file forNicole-2014-02-05.mpl

# Subspaces

printf("## Subspaces L to use (in Example 2.29 and following):\n");

V_linspace := [[0,1],[1,0]];                     #full space 
L_linspace := [[0,1]];
O_linspace := [];                                #trivial space;

# Computations

printf("## Computations of parametric volume, sum and intermediate sums in Example 2.29:\n");

for L in [V_linspace, L_linspace, O_linspace] do
    for bas in Bases do
        E[L][bas] := coeff(series(eval(subs({EXP=exp,TODD=Todd},tfunction_SL(1,1,s[bas],C[bas],L,[t*x1,t*x2]))), t=0), t, 0);
    od;
    for tau from 1 to nops(Taus) do
        E[L][tau] := simplify(add(E[L][bas], bas in Taus[tau]));
    od;
od:

# Output

for tau from 1 to nops(Taus) do
    printf("## Chamber tau_%d ##\n", tau):
    "Volume" = factor(E[V_linspace][tau]);
    "Sum" = E[O_linspace][tau];
    "Intermediate sum" = E[L_linspace][tau];
od;

#

printf("## Specialization to multiples of the following vector from tau_2 (Example 2.34):\n");

b_0 := [0, 0, 5, 3];
tau := 2;

Ex234Volume := collect(eval(subs({b = expand(t * b_0)}, E[V_linspace][tau])), t);
Ex234Sum := collect(eval(subs({b = expand(t * b_0)}, E[O_linspace][tau])), t);
Ex234IntermediateSum := collect(eval(subs({b = expand(t * b_0)}, E[L_linspace][tau])), t);
## could as well compute without first calculating the fully
## parametric sum!!

printf("## Figure 7:\n");

Ex234Plots := Array(1..2):
Ex234Plots[1] := plot([value(Ex234Sum), value(Ex234IntermediateSum), Ex234Volume],
                      't' = 0 .. 1,
                      color = ["Black", "Red", "Blue"],
                      discont = true, symbolsize = 15, symbol = solidcircle,
                      axis = [gridlines = [linestyle = dot]],
                      tickmarks = [[0 = 0, 1/5 = 1/5, 1/3 = 1/3, 2/5 = 2/5, 3/5 = 3/5,
                                    2/3 = 2/3, 4/5 = 4/5, 1 = 1],
                                   [1, 2, 3, 5, 8, 9, 13, 19]],
                      font = plotfont):
Ex234Plots[2] := plot([value(Ex234Sum), value(Ex234IntermediateSum), Ex234Volume],
                      't' = 1 .. 2.4,
                      color = ["Black", "Red", "Blue"],
                      discont = true, symbolsize = 15, symbol = solidcircle,
                      axis = [gridlines = [linestyle = dot]],
                      tickmarks = [[1 = 1, 6/5 = 6/5, 4/3 = 4/3, 7/5 = 7/5, 8/5 = 8/5,
                                    5/3 = 5/3, 2 = 2, 11/5 = 11/5, 7/3 = 7/3, 12/5 = 12/5],
                                   [19 = 19, 24 = 24, 26 = 26, 32 = 32, 39 = 39,
                                    41 = 41, 49 = 49, 60 = 60, 69 = 69, 72 = 72]],
                      font = plotfont):
display(Ex234Plots);

printf("## Specialization to multiples of the following irrational vector from tau_2 (Example 2.36):\n");

b_I := [0, 0, 3*sqrt(2), 3];
tau := 2;

Ex236Volume := collect(eval(subs({b = expand(t * b_I)}, E[V_linspace][tau])), t);
Ex236Sum := collect(eval(subs({b = expand(t * b_I)}, E[O_linspace][tau])), t);
Ex236IntermediateSum := collect(eval(subs({b = expand(t * b_I)}, E[L_linspace][tau])), t);

printf("## Figure 8:\n");

Ex236Plots := Array(1..2):
# simplify for plotting purposes; otherwise too many discontinuity
# markers are shown.
Ex236Sum:=simplify(evaluateEhrhart(Ex236Sum, 't')):
Ex236IntermediateSum:=simplify(evaluateEhrhart(Ex236IntermediateSum, 't')):
Ex236Plots[1] := plot([Ex236Sum, Ex236IntermediateSum, Ex236Volume],
                      't' = 0 .. 1,
                      color = ["Black", "Red", "Blue"],
                      discont = true, symbolsize = 15, symbol = solidcircle,
                      axis = [gridlines = [linestyle = dot]],
                      tickmarks = [[sqrt(2)/6 = sqrt(2)/6, 1/3 = 1/3, sqrt(2)/3 = sqrt(2)/3, 2/3 = 2/3, sqrt(2)/2 = sqrt(2)/2, 2*sqrt(2)/3 = 2*sqrt(2)/3, 1 = 1],
                                   [0, 1, 2, 3, 5, 6, 9, 13, 14]],
                      font = plotfont):
# must plot separately, otherwise a discontinuity marker is displayed
# in the intermediate sum graph.
Ex236Ticks := [1, 5*sqrt(2)/6, 4/3, sqrt(2), 7*sqrt(2)/6, 5/3, 8*sqrt(2)/6]:
Ex236Ticks := map(x -> (x=x), Ex236Ticks):
Ex236DiscontinuousPlot := plot([Ex236Sum],
                               't' = 1 .. 2,
                               color = ["Black"],
                               discont = true, symbolsize = 15, symbol = solidcircle,
                               axis = [gridlines = [linestyle = dot]],
                               tickmarks = [Ex236Ticks,
                                            [14, 18, 20, 26, 32, 34, 41]],
                               font = plotfont):
Ex236ContinuousPlot := plot([Ex236IntermediateSum, Ex236Volume],
                            't' = 1 .. 2,
                            color = ["Red", "Blue"],
                            discont = false,
                            axis = [gridlines = [linestyle = dot]],
                            tickmarks = [Ex236Ticks,
                                         [14, 18, 20, 26, 32, 34, 41]],
                            font = plotfont):
Ex236Plots[2] := display({Ex236DiscontinuousPlot, Ex236ContinuousPlot}):
display(Ex236Plots);

############################################################

printf("##### Approximating polynomials for the integer dilation of a 4-dimensional lattice simplex with the following vertices (Example 5.4): #####\n\n");

Simplex4 := [[4,6,4,3],[5,7,9,1],[5,7,3,7],[6,8,3,9],[2,1,8,0]];

for k from 0 to 4 do
    printf("## k = %d:\n", k):
    "Cone-by-cone" = collect((cone_by_cone(t, Simplex4, [0, 0, 0, 0], 0, k)), t);
    "Full-Barvinok" = collect((fullbarvinok(t, Simplex4, k, [0, 0, 0, 0], 0) assuming t::integer), t);
od;

############################################################

printf("##### Example 5.5 - approximating quasipolynomials for real dilations of triangle with the following vertices. Figure 11. #####\n\n");

mars := [[1,1], [1, 2], [2, 2]];

for k from 0 to 2 do
    printf("## k = %d:\n", k):
    cbc := collect((cone_by_cone_real(t, mars, [0, 0], 0, k)), t):
    Ex55ConeByCone[k] := cbc;
    fb := collect(fullbarvinok(t, mars, k, [0, 0], 0), t):
    Ex55FullBarvinok[k] := fb;
od;

printf("## Evaluations for k=2:\n"):

for t in [1/2, Pi/6, 1, 1-10^(-6), 1+10^(-6)] do  
    ## Skipping the example t = (17/10)^(1/3)/2, because 
    ## of a weird Maple bug in our installation of Maple 14, 
    ## same as this one: http://www.mapleprimes.com/questions/130220-Maple-15-X86-64-LINUX-Fedora-16-Log2
    "Cone-by-cone"(t) = simplify(value(cbc));
    "Full-Barvinok"(t) = simplify(value(fb));
od;
t := 't';

# Plot Figure 11.
Figure11 := Array(1..2):
subplots := {}:
# Can't just use "plot" with a list of functions here;
# or discontinuity markers will have the wrong colors.
for k from 0 to 2 do
    Ex55ConeByCone[k] := simplify(evaluateEhrhart(Ex55ConeByCone[k], 't'));
    subplots := subplots union {plot(Ex55ConeByCone[k], 't' = 0 .. 5.4,
                                     color = ["Green", "Red", "Black"][k+1],
                                     discont = true, symbolsize = 15, symbol = solidcircle,
                                     axis = [gridlines = [linestyle = dot]],
                                     font = plotfont)};
od:
Figure11[1] := display(subplots):
subplots := {}:
for k from 0 to 2 do
    Ex55FullBarvinok[k] := simplify(evaluateEhrhart(Ex55FullBarvinok[k], 't'));
    subplots := subplots union {plot(Ex55FullBarvinok[k], 't' = 0 .. 5.4,
                                     color = ["Green", "Red", "Black"][k+1],
                                     discont = true, symbolsize = 15, symbol = solidcircle,
                                      axis = [gridlines = [linestyle = dot]],
                                     font = plotfont)};
od:
Figure11[2] := display(subplots):
display(Figure11);

############################################################


printf("##### Example 5.6 - approximating quasipolynomials for real dilations of 3-dimensional simplex with the following vertices (Figures 12, 13): #####\n");


Simplex56 := [[0,1,1],[4,2,1],[1,1,2],[1,2,4]];

for k from 0 to 3 do
    printf("## k = %d:\n", k):
    cbc := collect((cone_by_cone_real(t, Simplex56, [0, 0, 0], 0, k)), t):
    Ex56ConeByCone[k] := cbc;
    fb := collect(fullbarvinok(t, Simplex56, k, [0, 0, 0], 0), t):
    Ex56FullBarvinok[k] := fb;
od;

# Plot Figures 12 and 13.
subplots := {}:
# Can't just use "plot" with a list of functions here;
# or discontinuity markers will have the wrong colors.
for k from 0 to 3 do
    Ex56ConeByCone[k] := simplify(evaluateEhrhart(Ex56ConeByCone[k], 't'));
    subplots := subplots union {plot(Ex56ConeByCone[k], 't' = 1.3 .. 3.9,
                                     color = ["Green", "Blue", "Red", "Black"][k+1],
                                     discont = true, symbolsize = 15, symbol = solidcircle,
                                     axis = [gridlines = [linestyle = dot]],
                                     font = plotfont)};
od:
Figure12 := display(subplots):
display(Figure12);
subplots := {}:
for k from 0 to 3 do
    Ex56FullBarvinok[k] := simplify(evaluateEhrhart(Ex56FullBarvinok[k], 't'));
    subplots := subplots union {plot(Ex56FullBarvinok[k], 't' = 1.3 .. 3.9,
                                     color = ["Green", "Blue", "Red", "Black"][k+1],
                                     discont = true, symbolsize = 15, symbol = solidcircle,
                                     axis = [gridlines = [linestyle = dot]],
                                     font = plotfont)};
od:
Figure13 := display(subplots):
display(Figure13);


########################### END ############################






## The following commented-out code could be used for plotting above.

############################################################

#mars0:=FBbonfrac(t,mars,0,[1,5,7,2],0);
#mars1:=FBbonfrac(t,mars,1,[1,5,7,2],0);
#mars2:=FBbonfrac(t,mars,2,[1,5,7,2],0):
#mars3:=FBbonfrac(t,mars,3,[1,5,7,2],0):
#mars4:=FBbonfrac(t,mars,4,[1,5,7,2],0):
#### plot([mars4,mars3,mars2,mars1,mars0], t=0.6...0.7,color=[black,red,blue,violet,green], axis[1]=[gridlines=[30, thickness=1, subticks=false]],axis[2]=[gridlines=[6, thickness=1, subticks=false]],thickness=1);

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

