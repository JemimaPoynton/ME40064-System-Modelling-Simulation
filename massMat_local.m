function massMat_local = massMat_local(mesh, eID)
% function massMat_local calculates the local mass matix for the element in
% the mesh defined by eID. Uses linear basis functions.
%
% Jemima Poynton 12/23

%% Create Quadrature Scheme
N = 4;
gq = CreateGQScheme(N);

%% Calculate Parameters
J = mesh.elem(eID).J;

psi0 = EvalBasis(0, gq.xipts);
psi1 = EvalBasis(1, gq.xipts);

%% Calculate Local Element Mass Matrices
Int00 = 0; Int01 = 0; Int10 = 0; Int11 = 0;

for i = 1:N % Defined directly for speed (to avoid slow for loops)
    Int00 = Int00 + gq.gsw(i)*psi0(i)*psi0(i)*J;
    Int01 = Int01 + gq.gsw(i)*psi0(i)*psi1(i)*J;
    Int10 = Int10 + gq.gsw(i)*psi1(i)*psi0(i)*J;
    Int11 = Int11 + gq.gsw(i)*psi1(i)*psi1(i)*J;
end

%% Form Matrix
massMat_local = [Int00  Int01; 
                 Int10  Int11];


