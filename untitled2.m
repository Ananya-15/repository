clear 
clc 

fid=fopen('StreamFunctionVal.txt','r'); 
formatSpec="%f"; 
psi=fscanf(fid,formatSpec); 

Lx=1; 
Ly=1; 
Nx=120; 
Ny=120; 
dx=Lx/(Nx-1); 
dy=Ly/(Ny-1); 

psi=reshape(psi,[120,120]); 
contour(0:dx:Lx,0:dy:Ly,psi)