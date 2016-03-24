## Abandoned attempt to have a 
## driver file that creates the figures for the paper "3 Ehrhart Polynomials".
with(plots):
#p := plot(x^2,x=-5..5):
read("3-ehrhart-polynomials-paper-examples.mpl"):
#### Figure 3 - Example 1.3
plotopts:="portrait,color,leftmargin=0,bottommargin=0,width=15cm,height=10cm":
plotsetup(eps, plotoutput="rectangle-irrationnel.eps", plotoptions=plotopts):
print(display(Ex13Plot)); # This pushes it out 
#### Unfortunately this method does not create plots as nice as the
#### worksheet: No gridlines! No typeset axis labels!!
.....
