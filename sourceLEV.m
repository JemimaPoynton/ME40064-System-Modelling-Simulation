function LEV = sourceLEV(f, eID, mesh)
% sourceLEV outputs the local element vector for the source term in the form 
%
% [int0;
%  int1]
%
% f: source term
% eID: local element number
% mesh: mesh data structure
%
% Jemima Poynton 11/23

%% Initialise GQ
N = 3;
gq = CreateGQScheme(N);

%% Extract Jacobian and boundary conditions
J = mesh.elem(eID).J;
psi0 = EvalBasis(0, gq.xipts);
psi1 = EvalBasis(1, gq.xipts);

%% Calculate Integrals
% The integrals can be defined directly
Int0 = 0; Int1 = 0; 

for i = 1:N
    Int0 = Int0 + gq.gsw(i)*f*psi0(i)*J;
    Int1 = Int1 + gq.gsw(i)*f*psi1(i)*J;
end
%% Form Vector
LEV = [Int0;
       Int1];