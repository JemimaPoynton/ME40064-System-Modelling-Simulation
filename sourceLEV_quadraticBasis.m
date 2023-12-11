function LEV = sourceLEV_quadraticBasis(f, eID, mesh)
% sourceLEV_quadraticBasis outputs the local element vector for the source term in the form 
%
% [int0;
%  int1]
%
% This function uses quadratic basis functions
%
% f: source term
% eID: local element number
% mesh: mesh data structure
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
Int0 = 0; Int1 = 0; Int2 = 0;

for i = 1:N
    Int0 = Int0 + gq.gsw(i)*f*psi0(i)*J;
    Int1 = Int1 + gq.gsw(i)*f*psi1(i)*J;
    Int2 = Int2 + gq.gsw(i)*f*psi2(i)*J;
end
%% Form Vector
LEV = [Int0;
       Int1;
       Int2];