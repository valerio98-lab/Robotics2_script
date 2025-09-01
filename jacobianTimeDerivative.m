function Jdot = jacobianTimeDerivative(q, dq, Jdata, varargin)
%JACOBIANTIMEDERIVATIVE  Derivata temporale di un Jacobiano (˙J)
%
% La routine riconosce automaticamente se i dati sono **numerici** o
% **simbolici**; in quest’ultimo caso sfrutta il Symbolic Math Toolbox.
%
% -------------------------------------------------------------------------
% Sintassi
%   Jdot = jacobianTimeDerivative(q, dq, Jsym)
%   Jdot = jacobianTimeDerivative(q, dq, Jhandle)
%   Jdot = jacobianTimeDerivative(q, dq, Jhandle, h)
%
% Input
%   q       [n×1]  : configurazione istantanea (double o sym)
%   dq      [n×1]  : velocità articolari ˙q (double o sym)
%   Jsym             matrice simbolica J(q)       (m×n)          oppure
%   Jhandle          function handle @(q)→J(q)    (m×n)
%   h       (opt.) : passo di differenziazione per il caso numerico
%                    default: 1e-6
%
% Output
%   Jdot  [m×n] : derivata temporale ˙J (double o sym a seconda dei dati)
%
% -------------------------------------------------------------------------
% Esempi d’uso
%
% --- 1) Dati NUMERICI con Jacobiano come function handle ---------------
% l1 = 0.4;  l2 = 0.3;
% Jf = @(q)[ -l1*sin(q(1))-l2*sin(q(1)+q(2)),  -l2*sin(q(1)+q(2));
%             l1*cos(q(1))+l2*cos(q(1)+q(2)),   l2*cos(q(1)+q(2)) ];
% q  = [0.5; 0.2];
% dq = [0.3;-0.1];
% Jdot_num = jacobianTimeDerivative(q,dq,Jf);
%
% --- 2) Dati SIMBOLICI con Jacobiano simbolico -------------------------
% syms q1 q2 dq1 dq2 l1 l2 real
% q  = [q1;q2];   dq = [dq1;dq2];
% J  = [ -l1*sin(q1)-l2*sin(q1+q2),  -l2*sin(q1+q2);
%         l1*cos(q1)+l2*cos(q1+q2),   l2*cos(q1+q2) ];
% Jdot_sym = jacobianTimeDerivative(q,dq,J);
% Jdot_sym = simplify(Jdot_sym)
%
% -------------------------------------------------------------------------

%%% --- 3) Esempio con Matrice Generica ----------------------------  

% syms phi theta psi dphi dtheta dpsi real
% q   = [phi; theta; psi];
% dq  = [dphi; dtheta; dpsi];
% 
% Rz = [cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1];
% Ry = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
% Rx = [1 0 0; 0 cos(phi) -sin(phi); 0 sin(phi) cos(phi)];
% R  = Rz*Ry*Rx;
% 
% Rdot = reshape(jacobian(R(:), q)*dq, size(R));  % tua funzione fa lo stesso
% % (Opz.) velocità angolare:  ω^ = Rᵀ Rdot  (skew)
% omega_hat = R.'*Rdot;  % estrai ω con vee(...)
% --------------------------------------------------------  

    %‒‒‒ Parametri --------------------------------------------------------
    if nargin < 4 || isempty(varargin{1}), h = 1e-6; else, h = varargin{1}; end
    q  = q(:);  dq = dq(:);
    n  = numel(q);

    %‒‒‒ Ottieni il Jacobiano corrente -----------------------------------
    if isa(Jdata,'function_handle')
        J = Jdata(q);
    else
        J = Jdata;                 % matrice simbolica già fornita
    end

    %‒‒‒ Modalità SIMBOLICA ---------------------------------------------
    if isa(J,'sym') || isa(q,'sym') || isa(dq,'sym')
        Jdot = sym(zeros(size(J)));
        for k = 1:n
            Jdot = Jdot + diff(J, q(k)) * dq(k);
        end
        Jdot = simplify(Jdot);
        return
    end

    %‒‒‒ Modalità NUMERICA ----------------------------------------------
    [m, nJ] = size(J);
    if nJ ~= n
        error('Dimensione Jacobiano incoerente con q.');
    end

    Jdot = zeros(m,n);
    e    = eye(n);
    for k = 1:n
        Jplus  = Jdata(q + h*e(:,k));
        Jminus = Jdata(q - h*e(:,k));
        dJdq_k = (Jplus - Jminus) / (2*h);  % ∂J/∂q_k
        Jdot   = Jdot + dJdq_k * dq(k);
    end
end
