function [mesh, E_idx, De_idx, B_idx] = skinProperties(mesh, layers, DVals, betaVals, gammaVals)
% function skinProperties creates a mesh property for each coefficient (D, beta, gamma)
% with spatially varying properties at each layer. 
%
% layers: layers of the skin in the form [E De B]
% DVals: diffusion coefficient in each layer of the skin [D1 D2 D3]
% betaVals: 1-by-3 matrix defining beta at each layer
% gammaVals: 1-by-3 matrix defining gamma at each layer

x = mesh.nvec;
E = layers(1); De = layers(2); B = layers(3);

[~, E_idx] = min(abs(x-E)); % Find index of boundary
[~, De_idx] = min(abs(x-De));
[~, B_idx] = min(abs(x-B));

mesh.DVec(:,1:E_idx) = DVals(1); % Set up varying D field
mesh.DVec(:,E_idx+1:De_idx) = DVals(2);
mesh.DVec(:, De_idx+1:mesh.ngn) = DVals(3);

mesh.betaVec(:,1:E_idx) = betaVals(1); % Set up varying beta field
mesh.betaVec(:,E_idx+1:De_idx) = betaVals(2);
mesh.betaVec(:, De_idx+1:mesh.ngn) = betaVals(3);

mesh.gammaVec(:,1:E_idx) = gammaVals(1); % Set a varying gamma field
mesh.gammaVec(:,E_idx+1:De_idx) = gammaVals(2);
mesh.gammaVec(:, De_idx+1:mesh.ngn) = gammaVals(3);