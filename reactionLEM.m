function LE_matrix = reactionLEM(lambda, eID, mesh)
% functions reactionLEM outputs the local element matrix of the reaction
% coefficient in the form 
%
% [int00 int01;
%  int10 int11]
%
% This function assumes linear basis functions.
%
% eID: local element number
% mesh: mesh data structure
% lambda: reaction coefficient
% LE_matrix: local element matrix
%
% Jemima Poynton 11/23

%% Extract Jacobian and boundary conditions
J = mesh.elem(eID).J;

%% Calculate Integrals
% The integrals can be defined directly

Int00 = (2/3)*lambda*J;
Int01 = (1/3)*lambda*J; 

Int10 = Int01; % Noting that the matrix is always symmetric (at this stage)
Int11 = Int00;

%% Form Matrix
LE_matrix = [Int00  Int01; 
             Int10  Int11];