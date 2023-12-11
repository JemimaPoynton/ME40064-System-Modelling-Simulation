function [globalMatrix] = stiffMat_global_qb(mesh, lambda)
% function stiffMat_global_qb creates the global stiffness matrix by placing
% the diffusion and reaction local element matrices in the correct
% positions. Uses quadratic basis functions.
%
% lambda: reaction coefficient
% 
% Jemima Poynton 12/23

%% Create Global Matrix
globalMatrix = zeros(2*mesh.ngn-1, 2*mesh.ngn-1);

%% Place diffusion and reaction terms
for eID = 1:mesh.ne
    diffusion = diffusionLEM_qb(eID, mesh);
    globalMatrix(2*eID-1:2*eID+1, 2*eID-1:2*eID+1) = diffusion + globalMatrix(2*eID-1:2*eID+1, 2*eID-1:2*eID+1);
    
    reaction = reactionLEM_qb(lambda, eID, mesh);
    globalMatrix(2*eID-1:2*eID+1, 2*eID-1:2*eID+1) = -reaction + globalMatrix(2*eID-1:2*eID+1, 2*eID-1:2*eID+1);
end


