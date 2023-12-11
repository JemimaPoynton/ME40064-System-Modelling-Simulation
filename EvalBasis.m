function [ psi ] = EvalBasis(lnid,xipts)
% evaluates the basis function psi_lnid at the xi points defined by xipts

sign = (-1)^(lnid+1);
psi = (1 + sign.*xipts)./2;

end

