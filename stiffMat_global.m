function [globalMatrix] = stiffMat_global(mesh, lambda)
% function stiffMat_global creates the global stiffness matrix by placing
% the diffusion and reaction local element matrices in the correct
% positions. Uses linear basis functions.
%
% lambda: reaction coefficient
% 
% Jemima Poynton 11/23

%% Create Global Matrix
globalMatrix = zeros(mesh.ngn, mesh.ngn);

%% Place diffusion and reaction terms
for eID = 1:mesh.ne
    diffusion = diffusionLEM(eID, mesh);
    globalMatrix(eID:eID+1, eID:eID+1) = diffusion + globalMatrix(eID:eID+1, eID:eID+1);

    reaction = reactionLEM(lambda, eID, mesh);
    globalMatrix(eID:eID+1, eID:eID+1) = -reaction + globalMatrix(eID:eID+1, eID:eID+1);
end



