function [alpha_min, alpha_max] = bound_of_g(g,q)
    %Function that outputs the upper bound alpha
    % on g 
    %
    %input:
    %- g = the gravity vector
    %- q = a vertical vector q [q1;q2;q3]
    %
    %output: the value of alpha

    dgdq=jacobian(g,q');
    m=simplify(dgdq'*dgdq);
    eigenvalues=eig(m);
    fprintf('\n maximum eigenvalue:'); 
    max_=simplify(max(eigenvalues));
    alpha_min=simplify(sqrt(max_)); 
    disp(alpha_min);

    fprintf('\n minimum eigenvalue:'); 
    max_=simplify(min(eigenvalues));
    alpha_max=simplify(sqrt(max_)); 
    disp(alpha_max); 

end
    

