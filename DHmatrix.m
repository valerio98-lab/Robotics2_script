function T = DHmatrix(alpha_val, a_val, d_val, theta_val, numerical)
% DHmatrix calcola la matrice di Denavit–Hartenberg
% alpha_val, theta_val, a_val, d_val: valori numerici o simbolici
% numerical: true => sostituisce i simboli con i valori; 
%            false => modalità simbolica con indicizzazione da theta (q1→alpha1,a1,d1)

    if numerical
        % --- caso numerico: creo simboli “placeholder” via sym()
        alpha_s = sym('alpha');   % non “syms alpha”
        theta_s = sym('theta');
        a_s     = sym('a');
        d_s     = sym('d');

        T_sym = [ cos(theta_s), -cos(alpha_s)*sin(theta_s),  sin(alpha_s)*sin(theta_s), a_s*cos(theta_s);
                  sin(theta_s),  cos(alpha_s)*cos(theta_s), -sin(alpha_s)*cos(theta_s), a_s*sin(theta_s);
                         0   ,            sin(alpha_s)       ,       cos(alpha_s)         ,   d_s          ;
                         0   ,              0               ,         0                 ,    1           ];

        T = subs(T_sym, [alpha_s, theta_s, a_s, d_s], ...
                        [alpha_val, theta_val, a_val, d_val]);

    else
        % --- caso simbolico dinamico
        if ~isa(theta_val,'sym')
            error('In modalità simbolica θ deve essere un sym, es. q1, q2, …');
        end
        nameTheta = char(theta_val);
        idx = regexp(nameTheta, '\d+$', 'match', 'once');
        if isempty(idx)
            error('Il nome di theta deve terminare con un numero, es. q1, q2, …');
        end

        % creo alphaN, aN, dN
        alphaN = sym(['alpha' idx]);
        aN     = sym(['a'     idx]);
        dN     = sym(['d'     idx]);

        T = [ cos(theta_val), -cos(alphaN)*sin(theta_val),  sin(alphaN)*sin(theta_val), aN*cos(theta_val);
              sin(theta_val),  cos(alphaN)*cos(theta_val), -sin(alphaN)*cos(theta_val), aN*sin(theta_val);
                     0      ,            sin(alphaN)         ,       cos(alphaN)         ,   dN           ;
                     0      ,              0                 ,         0                 ,    1           ];
    end
end
