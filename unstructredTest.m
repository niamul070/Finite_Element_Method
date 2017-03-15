%
% this is a test solution with sturctured triangles for unstructured grid
% code.
%
clc;
clear all;
close all;
T0=100;
a=1;
conductivity=25;
fid=fopen('tridomain.data');
line=fgetl(fid);
temp=fscanf(fid,'%d');
nelem=temp(1);
nnode=temp(2);
line=fgetl(fid);  %bface
connectivity=zeros(nelem,4);
nodePoints=zeros(nnode,3);
beta=zeros(nelem,3);
gamma=zeros(nelem,3);
k=zeros(3,3,nelem);
area=zeros(nelem,1);
globalK =zeros(nnode,nnode) ;


for i=1:nelem
    line=fgetl(fid);
    temp=sscanf(line,'%d')';
    connectivity(i,1)=temp(1);
    connectivity(i,2)=temp(2);
    connectivity(i,3)=temp(3);
    connectivity(i,4)=temp(4);
end

line=fgetl(fid); %coordinates of the nodePoints

for i = 1: nnode
    tline = fgetl(fid);
    temp = sscanf(tline, '%f')';
%     nodePoints(i,1) = temp(2);
%     nodePoints(i,2)=temp(3);
%     nodePoints(i,3)=temp(4);
    nodePoints(i,1) = temp(1);
    nodePoints(i,2)=temp(2);
    nodePoints(i,3)=temp(3);
end


%calculate element matrix and store in k 3d arrayt
for n=1:nelem
    xi=nodePoints(connectivity(n,2),1);
    xj=nodePoints(connectivity(n,3),1);
    xk=nodePoints(connectivity(n,4),1);
    yi=nodePoints(connectivity(n,2),2);
    yj=nodePoints(connectivity(n,3),2);
    yk=nodePoints(connectivity(n,4),2);
    beta(n,2)=yk-yi;
    beta(n,3)=yi-yj;
    beta(n,1)=yj-yk;
    gamma(n,1)=-(xj-xk);
    gamma(n,2)=-(xk-xi);
    gamma(n,3)=-(xi-xj);
    area(n)=polyarea([xi,xj,xk],[yi,yj,yk]);
    
end

for e=1:nelem
    for i=1:3
        for j=1:3
            k(i,j,e)=conductivity/(4*area(e))*(beta(e,j)*beta(e,i)+gamma(e,j)*gamma(e,i));
%            fprintf('i=%d j=%d out=%d\n',i,j,e);
        end
    end
end

B=connectivity(:,2:4);

%assemble global matrix
for i=1:nelem
    for irow=1:3
        for icol=1:3
            globalK(B(i,irow),B(i,icol))=  globalK(B(i,irow),B(i,icol))+ ...
                k(irow,icol,i);
        end;
    end;
end;
%initialize variables for ease of calculation
i=1;
count1=1;
count2=1;
left=0;
right=4.0;
top=4.0;
bottom=0;

%define boundary nodes 
for n=1:nnode
    if nodePoints(n,1)==right ||nodePoints(n,2)==bottom
        zeroNodeIndex(count1)=n;
        count1=count1+1;
    end
    if nodePoints(n,2)==top  %top boundary,y=4
      topBoundaryNodeIndex(count2)=n;
      topBoundaryNodeValues(count2,1)=T0*cos(pi*nodePoints(n,1)/(8*a));
      count2=count2+1;
    end
    
end

%right hand vector
rhsVector=-sum(globalK(:,topBoundaryNodeIndex)*topBoundaryNodeValues,2);

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

figure(1);
[xp,yp]=meshgrid(0:a/10:4*a,0:a/10:4*a);
T=griddata(nodePoints(:,1),nodePoints(:,2),solutionVector,xp,yp);
contourf(xp,yp,T);
title('Temperature contour for structured triangle with unstructured grid');
xlabel('x(m)');
ylabel('y(m)');
colorbar;

