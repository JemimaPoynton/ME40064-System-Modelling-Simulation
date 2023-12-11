function LE_matrix = gammaLEM_qb(eID, mesh)
% functions gammaLEM_qb outputs the local element matrix of the drug degredation
% rate in the form 
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

%% Initialise Gaussian Quadrature
N = 3;
gq = CreateGQScheme(N);

%% Extract Jacobian and boundary conditions
J = mesh.elem(eID).J;
psi0 = EvalBasis_QuadraticBasis(0, gq.xipts);
psi1 = EvalBasis_QuadraticBasis(1, gq.xipts);
psi2 = EvalBasis_QuadraticBasis(2, gq.xipts);

%% Calculate Integrals
Int00 = 0; Int01 = 0; Int02 = 0;
Int10 = 0; Int11 = 0; Int12 = 0;
Int20 = 0; Int21 = 0; Int22 = 0;

for i = 1:N % defined directly for speed
    gammaVal = EvalField_disc(mesh.gammaVec, eID); % evaluate value of beta at eID as a discontinuous basis function

    Int00 = Int00 + gq.gsw(i)*gammaVal*psi0(i)*psi0(i)*J;
    Int01 = Int01 + gq.gsw(i)*gammaVal*psi0(i)*psi1(i)*J;
    Int02 = Int02 + gq.gsw(i)*gammaVal*psi0(i)*psi2(i)*J;

    Int10 = Int10 + gq.gsw(i)*gammaVal*psi1(i)*psi0(i)*J;
    Int11 = Int11 + gq.gsw(i)*gammaVal*psi1(i)*psi1(i)*J;
    Int12 = Int12 + gq.gsw(i)*gammaVal*psi1(i)*psi2(i)*J;

    Int20 = Int20 + gq.gsw(i)*gammaVal*psi2(i)*psi0(i)*J;
    Int21 = Int21 + gq.gsw(i)*gammaVal*psi2(i)*psi1(i)*J;
    Int22 = Int22 + gq.gsw(i)*gammaVal*psi2(i)*psi2(i)*J;
end
%% Form Matrix
LE_matrix = [Int00 Int01 Int02; 
             Int10 Int11 Int12
             Int20 Int21 Int22];