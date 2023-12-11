function globalMassMat = massMat_global_qb(mesh)
% function massMat_global_qb calculates the global mass matrix for quadratic
% basis functions.

%% Initialise global matrix
globalMassMat = zeros(2*mesh.ngn-1, 2*mesh.ngn-1);

%% Calculate local matrices and place in correct position
for eID = 1:mesh.ne
    localMassMat = massMat_local_qb(mesh, eID);
    globalMassMat(2*eID -1 : 2*eID +1, 2*eID -1 : 2*eID +1) = ... 
        localMassMat + globalMassMat(2*eID -1 : 2*eID +1, 2*eID -1 :2*eID +1);
end