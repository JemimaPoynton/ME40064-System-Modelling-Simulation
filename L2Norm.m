function L2Norm = L2Norm(N, mesh, t, c)
% function L2 Norm calulates the L2 Norm of the data in c by comparing
% against the analytical solution of the transient diffusion equation
% Assumes linear basis functions
%
% N: order of gaussian quadrature scheme

gq = CreateGQScheme(N);
esquared = 0;

for i = 1:mesh.ne
    x = mesh.elem(i).x;
    x0 = x(1);
    x1 = x(2);
    c0 = c(:,i);
    c1 = c(:,i+1);
    
    for n = 1:N
        xi = gq.xipts(n);

        psi0 = EvalBasis(0, xi);
        psi1 = EvalBasis(1, xi);
        
        xxi = x0*psi0 + x1*psi1;
        cxi = c0'*psi0 + c1'*psi1;
        
        cexi = TransientAnalyticSoln(xxi,t);

        esquared = esquared + gq.gsw(n)*mesh.elem(i).J*(cexi - cxi).^2;
    end
end

L2Norm = sqrt(abs(esquared));