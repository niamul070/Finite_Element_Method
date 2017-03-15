%
% Finite element solution for steady state heat conduction equation % 
% with structured right triangular element
%
clc;
clear all;
close all;
T0=100;
a=1;
conductivity=25;
nelem=32;
nnode=25;
globalK =zeros(nnode,nnode) ;

B = [7, 6, 1;
1, 2, 7;
8, 7, 2;
2, 3, 8;
9, 8, 3;
3, 4, 9;
10, 9, 4;
4, 5, 10;
12, 11, 6;
6, 7, 12;
13,12, 7;
7, 8, 13;
14,13, 8;
8, 9, 14;
15, 14, 9;
9, 10, 15;
17, 16, 11;
11, 12, 17;
18, 17, 12;
12, 13, 18;
19, 18, 13;
13, 14, 19;
20, 19, 14;
14, 15, 20;
22, 21, 16;
16, 17, 22;
23, 22, 17;
17, 18, 23;
24, 23, 18;
18, 19, 24;
25, 24, 19;
19, 20, 25];


k=conductivity/2*[1 -1 0;
                  -1 2 -1;
                  0 -1 1];




%assemble global matrix
for i=1:nelem
    for irow=1:3
        for icol=1:3
            globalK(B(i,irow),B(i,icol))=  globalK(B(i,irow),B(i,icol))+ ...
                k(irow,icol);
        end;
    end;
end;

zeroNodeIndex=[1 2 3 4 5 10 15 20 25];
topBoundaryNodeIndex=[21 22 23 24];
topBoundaryNodeValues=[T0;T0*cos(pi/8);T0*cos(pi/4);T0*cos(3*pi/8)];
rhsVector=-sum(globalK(:,topBoundaryNodeIndex(1,:) )*topBoundaryNodeValues,2);

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
end;

 x=[0 1 2 3 4 0 1 2 3 4 0 1 2 3 4 0 1 2 3 4 0 1 2 3 4];
 y=[0 0 0 0 0 1 1 1 1 1 2 2 2 2 2 3 3 3 3 3 4 4 4 4 4];
 [xp,yp]=meshgrid(0:a/10:4*a,0:a/10:4*a);
 T=griddata(x,y,solutionVector,xp,yp);
 contourf(xp,yp,T);
 title('Temperature distribution using triangular elements');
 xlabel('x(m)');
 ylabel('y(m)');
 colorbar;



