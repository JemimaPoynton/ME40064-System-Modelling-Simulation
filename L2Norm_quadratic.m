function L2Norm = L2Norm_quadratic(N, mesh, t, c)
% function L2 Norm calulated the L2 Norm of the data in c by comparing
% against the analytical solution of the transient diffusion equation
% Assumes quadratic basis functions
%
% N: order of gaussian quadrature scheme


gq = CreateGQScheme(N);
esquared = 0;

for i = 1:mesh.ne
    x = mesh.elem(i).x;
    x0 = x(1);
    x2 = x(2);
    x1 = (x2 + x0)/2;
    
    c0 = c(:,2*i -1);
    c1 = c(:,2*i);
    c2 = c(:,2*i +1);
    
    for n = 1:N
        xi = gq.xipts(n);

        psi0 = EvalBasis_QuadraticBasis(0, xi);
        psi1 = EvalBasis_QuadraticBasis(1, xi);
        psi2 = EvalBasis_QuadraticBasis(2, xi);
        
        xxi = x0*psi0 + x1*psi1 + x2*psi2;
        cxi = c0'*psi0 + c1'*psi1 + c2'*psi2;
        
        cexi = TransientAnalyticSoln(xxi,t);

        esquared = esquared + gq.gsw(n)*mesh.elem(i).J*(cexi - cxi).^2;
    end
end

L2Norm = sqrt(abs(esquared));