%
% Finite element solution for steady state heat conduction equation % 
% with structured rectangular element
%
clc;
clear all;
close all;
T0=100;
a=1;
conductivity=25;
nelem=16;
nnode=25;
globalK =zeros(nnode,nnode) ;

B=[1 2 7 6;
   2 3 8 7;
   3 4 9 8;
   4 5 10 9;
   6 7 12 11;
   7 8 13 12;
   8 9 14 13;
   9 10 15 14;
   11 12 17 16;
   12 13 18 17;
   13 14 19 18;
   14 15 20 19;
   16 17 22 21;
   17 18 23 22;
   18 19 24 23;
   19 20 25 24];


k=conductivity/6*[4 -1 -2 -1;
        -1 4 -1 -2;
        -2 -1 4 -1;
        -1 -2 -1 4];




%assemble global matrix
for i=1:nelem
    for irow=1:4
        for icol=1:4
            globalK(B(i,irow),B(i,icol))=  globalK(B(i,irow),B(i,icol))+ ...
                k(irow,icol);
        end;
    end;
end;

zeroNodeIndex=[1 2 3 4 5 10 15 20 25];
topBoundaryNodeIndex=[21 22 23 24];
topBoundaryNodeValues=[T0;T0*cos(pi/8);T0*cos(pi/4);T0*cos(3*pi/8)];
rhsVector=-1*sum(globalK(:,topBoundaryNodeIndex(1,:) )*topBoundaryNodeValues,2);

%remove rows and columns

deleteColumn=cat(2,zeroNodeIndex,topBoundaryNodeIndex);
globalK(:,deleteColumn(1,:) )=[]; % delete columns
globalK(deleteColumn(1, :), :)=[]; % delete rows
rhsVector(deleteColumn(1,:))=[];

feSol=globalK\rhsVector;

count=1;
for i=1:nnode
    zeroIndex=find(zeroNodeIndex==i);
    nonZeroIndex=find(topBoundaryNodeIndex==i);
    if zeroIndex>0
        solutionVector(i)=0.0;
    elseif nonZeroIndex>0
        solutionVector(i)=topBoundaryNodeValues(nonZeroIndex);
    else
        solutionVector(i)=feSol(count);
        count=count+1;
    end;
end

 x=[0 1 2 3 4 0 1 2 3 4 0 1 2 3 4 0 1 2 3 4 0 1 2 3 4];
 y=[0 0 0 0 0 1 1 1 1 1 2 2 2 2 2 3 3 3 3 3 4 4 4 4 4];
 [xp,yp]=meshgrid(0:a/10:4*a,0:a/10:4*a);
 T=griddata(x,y,solutionVector,xp,yp);
 contourf(xp,yp,T);
 title('Temperature distribution using rectangular elements');
 xlabel('x(m)');
 ylabel('y(m)');
 colorbar;



