function [A,B] = linearize_system(f, x, u, x_eq, u_eq, params, pvals, t_sym, t_eq, verbose)
% LINEARIZE_SYSTEM  Linearizza il sistema non lineare dx/dt = f(x,u[,t]) attorno a (x_eq, u_eq[, t_eq])
%
% USO BASE:
%   [A,B] = linearize_system(f, x, u, x_eq, u_eq)
%
% INPUT:
%   f       : vettore simbolico delle dinamiche  dx/dt = f(x,u[,t])
%             Esempio: f = [x2; -k/m*x1 - b/m*x2 + u/m];
%
%   x       : vettore colonna di stati simbolici (stessa dimensione di f)
%             Esempio: x = [x1; x2];
%
%   u       : vettore (o scalare) di ingressi simbolici. Se il sistema è autonomo,
%             passare [] (vuoto).
%             Esempio: u = u;         % scalare
%                      u = [u1; u2];  % vettore
%
%   x_eq    : punto di equilibrio per gli stati (valori numerici, vettore colonna)
%             Esempio: x_eq = [0; 0];
%
%   u_eq    : punto di equilibrio per gli ingressi (valori numerici). Se u=[], passare [].
%             Esempio: u_eq = 0;  oppure  u_eq = [];
%
%   params  : (opzionale) vettore di parametri simbolici presenti in f/x/u
%             Esempio: params = [m, b, k];
%
%   pvals   : (opzionale) valori numerici dei parametri in 'params' (stessa lunghezza)
%             Esempio: pvals  = [2.0, 0.3, 5.0];
%
%   t_sym   : (opzionale) simbolo del tempo se f dipende esplicitamente da t
%             Esempio: syms t real; t_sym = t;
%
%   t_eq    : (opzionale) valore del tempo in cui linearizzare (numerico)
%             Esempio: t_eq = 0;   % o il tempo di interesse
%
%   verbose : (opzionale, default=false) se true, stampa A, B e il check di equilibrio.
%
% OUTPUT:
%   A : matrice di stato linearizzata  A = (∂f/∂x) |_(x_eq,u_eq[,t_eq])
%   B : matrice d’ingresso linearizzata B = (∂f/∂u) |_(x_eq,u_eq[,t_eq]); se u=[] → B ha 0 colonne
%
% ESEMPI RAPIDI:
%   % 1) Sistema SDOF senza parametri simbolici né dipendenza da t
%   syms x1 x2 u real
%   syms k m b real
%   f = [x2; -k/m*x1 - b/m*x2 + u/m];
%   x = [x1; x2];  u = u;
%   [A,B] = linearize_system(subs(f,{k,m,b},{5,2,0.3}), x, u, [0;0], 0);
%
%   % 2) Stesso sistema ma passando parametri separatamente
%   params = [k m b]; pvals = [5 2 0.3];
%   [A,B] = linearize_system(f, x, u, [0;0], 0, params, pvals);
%
%   % 3) Sistema autonomo (u vuoto)
%   u = []; u_eq = [];
%   [A,B] = linearize_system(f, x, u, [0;0], u_eq, params, pvals);
%
%   % 4) Sistema time-varying (dipendente da t)
%   syms t a0 w real
%   f_tv = [x2;
%           -k/m*x1 - b/m*x2 + a0/m*sin(w*t)];
%   params = [k m b a0 w]; pvals = [5 2 0.3 1.0 2*pi];
%   [A,B] = linearize_system(f_tv, x, [], [0;0], [], params, pvals, t, 0, true);

    if nargin < 6, params = []; pvals = []; end
    if nargin < 8, t_sym = [];  t_eq  = [];  end
    if nargin < 10, verbose = false; end

    % Forza forma colonna per evitare concatenazioni errate
    x    = x(:);
    x_eq = x_eq(:);
    if ~isempty(u)
        u    = u(:);
        u_eq = u_eq(:);
    else
        u_eq = [];
    end

    % Jacobiani simbolici
    A_sym = jacobian(f, x);
    if isempty(u)
        B_sym = sym(zeros(numel(f), 0));
    else
        B_sym = jacobian(f, u);
    end

    % Costruisci lista sostituzioni
    subs_syms = [x; u; params(:)];
    subs_vals = [x_eq; u_eq; pvals(:)];
    if ~isempty(t_sym)
        subs_syms = [subs_syms; t_sym];
        subs_vals = [subs_vals; t_eq];
    end

    % Valutazione A,B al punto di linearizzazione
    A_sub = subs(A_sym, subs_syms, subs_vals);
    B_sub = subs(B_sym, subs_syms, subs_vals);

    if isempty(symvar(A_sub)) && isempty(symvar(B_sub))
        A = double(A_sub);
        B = double(B_sub);
    else
        % Se restano simboli (parametri non assegnati), restituiamo vpa
        A = vpa(A_sub);
        B = vpa(B_sub);
        if verbose
            warning('Linearizzazione con parametri non numerici: ritorno A,B in vpa (simboliche).');
        end
    end

    % ---- CHECK EQUILIBRIO (solo informativo) ----------------------------
    % Per definizione, (x_eq, u_eq[, t_eq]) è punto di equilibrio se f(x_eq, u_eq[, t_eq]) = 0.
    % Calcoliamo la "residua" f_eval e la sua norma: se ||f_eval|| <= tol → equilibrio (entro tolleranza).
    if verbose
        f_eval = subs(f, subs_syms, subs_vals);  % f(x_eq, u_eq [, t_eq])
        tol = 1e-9;

        % Proviamo valutazione numerica della norma (se possibile)
        if isempty(symvar(f_eval))
            res_norm = norm(double(f_eval));
            is_eq = res_norm <= tol;
        else
            % Se ancora simbolico, tentiamo vpa; se rimangono simboli indeterminati, non possiamo confermare
            try
                res_norm = double(norm(vpa(f_eval)));
                is_eq = res_norm <= tol;
            catch
                res_norm = NaN;
                is_eq = false; % non confermabile numericamente
            end
        end

        fprintf('Check equilibrio: ||f(x_eq,u_eq[,t_eq])|| = %g  -> %s\n', ...
                res_norm, ternary(is_eq,'OK (equilibrio)','NON equilibrio (o non verificabile numericamente)'));
        fprintf('A:\n'); disp(A);
        fprintf('B:\n'); disp(B);
    end
end

% Utility inline per messaggistica (evita if-else verbosi in fprintf)
function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
