function [globalMatrix, globalVector] = applyBCs(globalMatrix, globalVector, BCs, dt, theta, t, mesh)
% function applyBCs applies the boundary conditions (BCs) defined in the
% object BCs to the global matrix and global vectoy, 
% first applying Neumann conditions. 
% The BCs are applied at either side of the mesh.
% This function is capable of handling both time-varying, and time
% invariant terms
% 
% BCs.[type].pos: boolean representing whether a BC is present or not at the [start end]
%                 of the mesh, where [type] is 'neumann'/'dirichlet'
% BCs.[type].val: value of the condition at the [start end] of the mesh,
%                 defined as a gradient for type Neumann.
% t: current time
%
% Jemima Poynton 12/23

neumannVector = zeros(size(globalVector)); % initialise

%% Extract current time BCs 
[~, idx] = min(abs(BCs.neumann.t-t)); % finds index of closest value to t for handling time-varying terms

%% Extract D at boundaries
D = [mesh.DVec(1) mesh.DVec(end)];

%% Calculate Neumann Boundary Condition Vector
if BCs.neumann.pos(1) == 1 % check if at position 1
    neumannVector(1) = -D(1)*BCs.neumann.val(idx,1); % -ve for x=start condition (see derivation)
end

if BCs.neumann.pos(2) == 1 % check if at position 2
    neumannVector(end) = D(2)*BCs.neumann.val(idx,2); % +ve for x=end condition (see derivation)
end

%% Apply BCs
globalVector = globalVector + dt*(theta*neumannVector + (1- theta)*neumannVector);

%% Apply Dirichlet Boundary Conditions  
[~, idx] = min(abs(BCs.dirichlet.t-t));

if BCs.dirichlet.pos(1) == 1  % check if at position 1
    globalMatrix(1,:) = zeros(size(globalMatrix(1,:)));
    globalMatrix(1,1) = 1;

    globalVector(1) = BCs.dirichlet.val(idx,1);
end

if BCs.dirichlet.pos(2) == 1 % check if at position 2
    globalMatrix(end,:) = zeros(size(globalMatrix(end,:)));
    globalMatrix(end,end) = 1;

    globalVector(end) = BCs.dirichlet.val(idx,2);
end


