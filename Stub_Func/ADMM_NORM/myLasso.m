function [C, CKSym] =  myLasso (Y, Q,preZ, lambda0, lambda1, rho, alpha, type)
if nargin <8
    type = 1;
end

if type==0
    C = CVX_Lasso(Y, Q, 0, lambda0, lambda1);
elseif type == 2
    C = lasso_Q_norm2(Y, Q, lambda0, lambda1,rho, alpha, preZ);
else
    C = lasso_Q(Y, Q, lambda0, lambda1,rho, alpha, preZ);
end
[CKSym,CAbs] =  BuildAdjacency(thrC(C,.99));
end

