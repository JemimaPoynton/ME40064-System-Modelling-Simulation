function [ dpsidxi ] = EvalBasisGrad_QuadraticBasis(lnid,xipts)
% function EvalBasisGrad_QuadraticBasis evaluates the gradient (dpsi/dxi) 
% of the basis function psi_lnid at the different xi points (xipts)

    S1 = (-1)^(lnid);
    S2 = lnid - 1;
    K = [0.5 1 0.5];

    dpsidxi = K(lnid+1)*(S1*2*xipts + S2);

end

