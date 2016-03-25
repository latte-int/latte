## Driver file that creates the figures for the paper "3 Ehrhart
## Polynomials".
##
## To create quality figures, this file needs to be run in a
## worksheet.  If run from the command line, this method does not
## create nice plots: No gridlines! No typeset axis labels!!

with(plots):
read("3-ehrhart-polynomials-paper-examples.mpl"):
plotopts:="portrait,color,leftmargin=0,bottommargin=0,width=15cm,height=10cm,noborder,font=[Times, roman, 12]":  ### doesn't seem to be used....
#### Figure 3 - Example 1.3
plotsetup(eps, plotoutput="rectangle-irrationnel-basis.eps", plotoptions=plotopts):
print(display(BasisPlots)); # This pushes it out 
plotsetup(eps, plotoutput="rectangle-irrationnel.eps", plotoptions=plotopts):
print(display(Ex13Plot)); # This pushes it out 
#### Figure 7
plotsetup(eps, plotoutput="SlicedDilatedtetra0-1.eps", plotoptions=plotopts):
print(display(Ex234Plots[1]));
plotsetup(eps, plotoutput="SlicedDilatedtetra1-2.eps", plotoptions=plotopts):
print(display(Ex234Plots[2]));
#### Figure 8
plotsetup(eps, plotoutput="SlicedDilatedtetra-irrational0-1.eps", plotoptions=plotopts):
print(display(Ex236Plots[1]));
plotsetup(eps, plotoutput="SlicedDilatedtetra-irrational1-2.eps", plotoptions=plotopts):
print(display(Ex236Plots[2]));
#### Figure 11
plotsetup(eps, plotoutput="pic-en-dim2-conebycone.eps", plotoptions=plotopts):
print(display(Figure11[1]));
plotsetup(eps, plotoutput="pic-en-dim2-barvinok.eps", plotoptions=plotopts):
print(display(Figure11[2]));
#### Figures 12 and 13
plotsetup(eps, plotoutput="conebycone-simplexDim3v3.eps", plotoptions=plotopts):
print(display(Figure12));
plotsetup(eps, plotoutput="fullbarvinok-simplexDim3v3.eps", plotoptions=plotopts):
print(display(Figure13));
####
plotsetup(default);
