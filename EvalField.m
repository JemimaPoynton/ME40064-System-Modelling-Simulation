function [val] = EvalField(mesh, eID, xipt)
% evaluates the diffusion coefficient at eID. Assumes linear basis
% functions
psi0 = EvalBasis(0,xipt);
psi1 = EvalBasis(1,xipt);

val = mesh.DVec(eID)*psi0 + mesh.DVec(eID + 1)*psi1;
