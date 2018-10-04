
function writeNodeData2D(inN)

Globals2D;

N = inN;
Nfp = N+1;
Nfaces = 3;

%% Nodal data
[r,s] = Nodes2D(N);
[r,s] = xytors(r,s);
Np = length(r);

% find all the nodes that lie on each edge
NODETOL = 1e-8;
faceNodes1   = find( abs(s+1) < NODETOL)';
faceNodes2   = find( abs(r+s) < NODETOL)';
faceNodes3   = find( abs(r+1) < NODETOL)';
FaceNodes  = [faceNodes1;faceNodes2;faceNodes3]';

V = Vandermonde2D(N, r, s);
MM = inv(transpose(V))/V;
[Dr,Ds] = Dmatrices2D(N, r, s, V);
LIFT = Lift2D(N, FaceNodes, r, s);

fname = sprintf('triangleN%02d.dat', N);
fid = fopen(fname, 'w');

writeFloatMatrix(fid, r, 'Nodal r-coordinates');
writeFloatMatrix(fid, s, 'Nodal s-coordinates');
writeFloatMatrix(fid, Dr, 'Nodal Dr differentiation matrix');
writeFloatMatrix(fid, Ds, 'Nodal Ds differentiation matrix');
writeFloatMatrix(fid, MM, 'Nodal Mass Matrix');

writeIntMatrix(fid, FaceNodes'-1, 'Nodal Face nodes');
writeFloatMatrix(fid, LIFT, 'Nodal Lift Matrix');

%% Plotting data
%compute equispaced nodes on equilateral triangle
[plotR,plotS] = EquiNodes2D(N+4);
plotNp = length(plotR);
plotEToV = delaunay(plotR,plotS)-1;
plotNelements = size(plotEToV,1);
[plotR,plotS] = xytors(plotR,plotS);

%check triangulation
before = plotNelements;
sk = 0;
for e=1:plotNelements
  v1 = plotEToV(e,1)+1;
  v2 = plotEToV(e,2)+1;
  v3 = plotEToV(e,3)+1;

  x1 = plotR(v1);
  x2 = plotR(v2);
  x3 = plotR(v3);

  y1 = plotS(v1);
  y2 = plotS(v2);
  y3 = plotS(v3);

  plotA = (x2-x1)*(y3-y1) - (y2-y1)*(x3-x1);
  if(abs(plotA)>1e-5) sk = sk+1; plotEToV(sk,:) = [v1-1,v2-1,v3-1]; end
end
plotNelements = sk;
plotEToV = plotEToV(1:sk,:);
after = plotNelements;

plotInterp = Vandermonde2D(N, plotR,plotS)/V;

writeFloatMatrix(fid, plotR, 'Plotting r-coordinates');
writeFloatMatrix(fid, plotS, 'Plotting s-coordinates');
writeFloatMatrix(fid, plotInterp, 'Plotting Interpolation Matrix');
writeIntMatrix(fid, plotEToV, 'Plotting triangulation');

%% Cubature data
%volume cubature
[cubr,cubs,cubw] = Cubature2D(3*N);
cInterp = Vandermonde2D(N, cubr, cubs)/V;
Ncub = length(cubr);

cV = Vandermonde2D(N, cubr, cubs);
cV'*diag(cubw)*cV;

[cVr,cVs] = GradVandermonde2D(N, cubr, cubs);
cubDrT = V*transpose(cVr)*diag(cubw);
cubDsT = V*transpose(cVs)*diag(cubw);
cubProject = V*cV'*diag(cubw); %% relies on (transpose(cV)*diag(cubw)*cV being the identity)

writeFloatMatrix(fid, cubr, 'Cubature r-coordinates');
writeFloatMatrix(fid, cubs, 'Cubature s-coordinates');
writeFloatMatrix(fid, cubw, 'Cubature weights');

writeFloatMatrix(fid, cInterp, 'Cubature Interpolation Matrix');
writeFloatMatrix(fid, cubDrT, 'Cubature Weak Dr Differentiation Matrix');
writeFloatMatrix(fid, cubDsT, 'Cubature Weak Ds Differentiation Matrix');
writeFloatMatrix(fid, cubProject, 'Cubature Projection Matrix');

%surface cubature
[z,w] = JacobiGQ(0,0,ceil(3*N/2));
Nfi = length(z);

ir = [z,-z,-ones(Nfi,1)];
is = [-ones(Nfi,1), z, -z];
iw = [w,w,w];

sV = Vandermonde2D(N, ir(:), is(:));
sInterp = sV/V;
interp = [sInterp(1:Nfi,FaceNodes(:,1));sInterp(Nfi+1:2*Nfi,FaceNodes(:,2));sInterp(2*Nfi+1:3*Nfi,FaceNodes(:,3))];

%integration node lift matrix
iLIFT = V*V'*sInterp'*diag(iw(:));

writeFloatMatrix(fid, interp, 'Cubature Surface Interpolation Matrix');
writeFloatMatrix(fid, iLIFT, 'Cubature Surface Lift Matrix');

fclose(fid)

end
