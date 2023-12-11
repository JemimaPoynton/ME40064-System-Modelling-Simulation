function LE_matrix = diffusionLEM_qb(eID, mesh)
% functions diffusionLEM outputs the local element matrix of the diffusion
% coefficient in the form 
%
%  [Int00 Int01 Int02; 
%   Int10 Int11 Int12
%   Int20 Int21 Int22];
%
% This function assumes quadratic basis functions.
%
% eID: local element number
% mesh: mesh data structure
% LE_matrix: local element matrix
%
% Jemima Poynton 12/23

%% Extract Jacobian and boundary conditions
J = mesh.elem(eID).J;

%% Initiate Gaussian Quadrature Scheme
N = 3; % !MAYBE CREATE A FUNCTION THAT DETERMINES WHAT N SHOULD BE! 
gq = CreateGQScheme(N);

%% Define derivatives
dxi_dx = 1/J;
dpsi_0_dxi = EvalBasisGrad_QuadraticBasis(0,gq.xipts); % in form [dpsi0_dxi dpsi1_dxi] for [!BASIS FUNCTION!]
dpsi_1_dxi = EvalBasisGrad_QuadraticBasis(1,gq.xipts);
dpsi_2_dxi = EvalBasisGrad_QuadraticBasis(2,gq.xipts);

%% Calculate Integrals
% Noting that the integral contains no xi terms, the integrals can be defined directly
Int00 = 0; Int01 = 0; Int02 = 0;
Int10 = 0; Int11 = 0; Int12 = 0;
Int20 = 0; Int21 = 0; Int22 = 0;

for i = 1:N
    D = EvalField(mesh, eID, gq.xipts(i));
    
    Int00 = Int00 + gq.gsw(i)*D*dpsi_0_dxi(i)*dpsi_0_dxi(i)*(1/J); % Hard coded for speed
    Int01 = Int01 + gq.gsw(i)*D*dpsi_0_dxi(i)*dpsi_1_dxi(i)*(1/J);
    Int02 = Int02 + gq.gsw(i)*D*dpsi_0_dxi(i)*dpsi_2_dxi(i)*(1/J);

    Int10 = Int10 + gq.gsw(i)*D*dpsi_1_dxi(i)*dpsi_0_dxi(i)*(1/J);
    Int11 = Int11 + gq.gsw(i)*D*dpsi_1_dxi(i)*dpsi_1_dxi(i)*(1/J);
    Int12 = Int12 + gq.gsw(i)*D*dpsi_1_dxi(i)*dpsi_2_dxi(i)*(1/J);

    Int20 = Int20 + gq.gsw(i)*D*dpsi_2_dxi(i)*dpsi_0_dxi(i)*(1/J);
    Int21 = Int21 + gq.gsw(i)*D*dpsi_2_dxi(i)*dpsi_1_dxi(i)*(1/J);
    Int22 = Int22 + gq.gsw(i)*D*dpsi_2_dxi(i)*dpsi_2_dxi(i)*(1/J);
end

%% Form Matrix
LE_matrix = [Int00 Int01 Int02; 
             Int10 Int11 Int12
             Int20 Int21 Int22];

