function [ psi ] = EvalBasis_QuadraticBasis(lnid,xipts)
% function EvalBasis_quadraticBasis evaluates the value of the basis
% function psi_lnid at the different xi points (xipts)

S1 = (-1)^(lnid);
S2 = [-1 1 1];
I1 = [1 0 1];
K = [0.5 1 0.5];

psi = K(lnid+1)*(S1*xipts.^2 + S2(lnid+1)*xipts.^I1(lnid+1));

end

