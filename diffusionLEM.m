function LE_matrix = diffusionLEM(eID, mesh)
% functions diffusionLEM outputs the local element matrix of the diffusion
% coefficient in the form 
%
% [int00 int01;
%  int10 int11]
%
% This function assumes linear basis functions.
%
% eID: local element number
% mesh: mesh data structure
% LE_matrix: local element matrix
%
% Jemima Poynton 12/23

%% Extract Jacobian and boundary conditions
J = mesh.elem(eID).J;

%% Define derivatives
dpsi_0_dxi = EvalBasisGrad(0,[]);
dpsi_1_dxi = EvalBasisGrad(1,[]);

%% Initiate Gaussian Quadrature and Extract D
N = 4;
gq = CreateGQScheme(N);

%% Calculate Integrals
Int00 = 0; Int01 = 0; Int10 = 0; Int11 = 0;

for i = 1:N % integrals defined directly for speed
    D = EvalField(mesh, eID, gq.xipts(i)); % find what D is at this point on the mesh

    Int00 = Int00 + gq.gsw(i)*D*dpsi_0_dxi*dpsi_0_dxi*(1/J);
    Int01 = Int01 + gq.gsw(i)*D*dpsi_0_dxi*dpsi_1_dxi*(1/J);
    Int10 = Int10 + gq.gsw(i)*D*dpsi_1_dxi*dpsi_0_dxi*(1/J);
    Int11 = Int11 + gq.gsw(i)*D*dpsi_1_dxi*dpsi_1_dxi*(1/J);
end

%% Form Matrix
LE_matrix = [Int00  Int01; 
             Int10  Int11];

