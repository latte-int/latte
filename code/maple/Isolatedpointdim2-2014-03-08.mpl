with(LinearAlgebra):with(combinat):L:=[[1,1], [1, 2], [2, 2]];
 F :=[x-t, y-x, 2*t-y];
subs({x=L[1][1],y=L[1][2],t=1},F[3]);
subs({x=L[2][1],y=L[2][2],t=1},F[3]);
subs({x=L[3][1],y=L[3][2],t=1},F[3]);

numberS:=proc(t) local nn,x,y,z;nn:={}; 
for x from 0 to 20 do 
for y from 0 to 20 do 
  if x-t>=0 and  y-x>=0 and 
2*t-y>=0  
then nn:={op(nn),[x,y]}; 
fi;;od;od;nn;end:
numberS(1-1/100);numberS(1);numberS(1+1/100); 

subs({x=2,y=2,t=1+eps},F); #thus [1,1] remains left; [2,2] remains right;

