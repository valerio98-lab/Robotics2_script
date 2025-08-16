function Z = skew_null(q)
%% 
% The function is thougth to generalize the skew-symmetric construction of
% a matrix, tipically usefull when you want to create an alternative
% factorization of S(q,q_dot) and you are not dealing with a 3DOF (you
% cannot use the 3x3 skew-symm standard matrix with a >3 DOF)

%Example: 
% q = [1 2 3 4];
% Z = skew_null(q);
% 
% disp(norm(Z*q));

%%
    q=q(:);
    U=null(q.');
    if size(U,2) < 2
        Z = zeros(length(q));
        return 
    end 
    Z = U(:,1)*U(:,2).' - U(:,2)*U(:,1).';
end