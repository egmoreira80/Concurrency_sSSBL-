function [emd] = estimating_EMD(sim,sol,vertices,faces)
nEigs          = 100;
FEM            = firstOrderFEM(vertices,faces);
[evecs,evals]  = eigs(FEM.laplacian,FEM.vtxInnerProds,nEigs,'sm'); 
evalsd         = diag(evals);
index          = find(abs(evalsd)<10^(-12)*abs(evalsd(end)));
evecs(:,index) = [];
structure      = precomputeEarthMoversADMM(vertices,faces,evecs);
[emd,J]        = earthMoversADMM(vertices, faces, sim, sol, structure);