function globalMassMat = massMat_global(mesh)
% function massMat_global calculates the global mass matrix for linear
% basis functions.

%% Initialise Global Matrix
globalMassMat = zeros(mesh.ngn, mesh.ngn);

%% Calculate local matrices and place in correct position
for eID = 1:mesh.ne
    localMassMat = massMat_local(mesh, eID);
    globalMassMat(eID:eID+1, eID:eID+1) = localMassMat + globalMassMat(eID:eID+1, eID:eID+1);
end