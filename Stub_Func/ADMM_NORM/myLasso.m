function [C, CKSym] =  myLasso (Y, Q,preZ, lambda0, lambda1, rho, alpha, type)
if nargin <8
    type = 1;
end

if type==0
    C = CVX_Lasso(Y, Q, 0, lambda0, lambda1);
elseif type == 2
        C = lasso_Q_norm2(Y, Q, lambda0, lambda1,rho, alpha, preZ);
%     C = lasso_Q_norm2_V2(Y,lambda0, lambda1, Q, preZ );
%     alpha = lambda0/lambda1; gamma1 = .3;
%     C = lasso_Q_norm2_V2(Y,lambda0, lambda1, Q, preZ, alpha,  gamma1 );

    
else
    C = lasso_Q(Y, Q, lambda0, lambda1,rho, alpha, preZ);
end
[CKSym,CAbs] =  BuildAdjacency(thrC(C,.99));
end

