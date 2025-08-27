function sys = sym_ode_to_state(eqs, inputs, outputs, eq_point)
%SYM_ODE_TO_STATE  ODE simboliche -> forma di stato + (opz.) linearizzazione
%
% sys = sym_ode_to_state(eqs, inputs, outputs, eq_point)
%
% INPUT
%   eqs      : cell array di equazioni simboliche (==) in funzione del tempo, p.es.
%              { diff(x(t),2) + a*diff(x(t),1) + b*x(t) == u1(t),  diff(q(t),1) == -c*q(t) + d*x(t) }
%   inputs   : cell array con i simboli d’ingresso time-dependent, p.es. {u1(t), u2(t)}  (può essere {})
%   outputs  : (opz.) cell array di espressioni simboliche per l’uscita, in termini
%              delle variabili originali (es. {x(t), q(t)} o {x(t)-q(t)}).
%              Se omesso, C = I (uscita = stato completo).
%   eq_point : (opz.) struct con i campi:
%                 .x0  (n×1)  valore di equilibrio degli stati rinominati x1..xn (default zeros)
%                 .u0  (m×1)  valore di equilibrio degli ingressi U1..Um (default zeros)
%                 .subs       struct o containers.Map con assegnazioni parametriche (a=..., b=..., ecc.)
%
% OUTPUT (struct sys)
%   .f      : campo vettoriale simbolico f(x,U)
%   .X      : vettore degli stati rinominati [x1;...;xn]
%   .U      : vettore ingressi rinominati [U1;...;Um]
%   .Yorig  : vettore delle variabili di stato originali restituite da odeToVectorField
%   .subsx  : cell array di sostituzioni { Yorig(i) == X(i) }
%   .subsu  : cell array di sostituzioni { inputs{j} == U(j) }
%   .A,.B,.C,.D : (se eq_point passato) Jacobiane valutate all’equilibrio
%   .A_sym,.B_sym,.C_sym,.D_sym : Jacobiane simboliche (non valutate)
%
% NOTE
%   - Tutto resta simbolico (richiede Symbolic Math Toolbox).
%   - Le equazioni devono essere ODE (no DAE). Usare '==' nelle equazioni.

    arguments
        eqs      (1,:) cell
        inputs   (1,:) cell = {}
        outputs  (1,:) cell = {}
        eq_point struct     = struct()
    end

    % 1) ODE -> vector field
    [F, Y] = odeToVectorField(eqs{:});   % Y: stati originali, F: loro derivate prime

    n = numel(Y);
    m = numel(inputs);

    % 2) Rinomino stati/ingressi
    X = sym('x', [n 1], 'real');      % x1..xn
    if m>0
        U = sym('U', [m 1], 'real');  % U1..Um
    else
        U = sym([]);                  % nessun ingresso
    end

    % 3) Sostituzioni verso le nuove variabili
    subsx = arrayfun(@(i) Y(i)==X(i), 1:n, 'uni', 0);
    subsu = {};
    if m>0
        subsu = arrayfun(@(j) inputs{j}==U(j), 1:m, 'uni', 0);
    end

    % Campo vettoriale f(x,U)
    f = subs(F, [Y(:).', cell2sym(inputs)], [X(:).', U(:).']);

    % 4) Uscite (default = stato completo)
    if isempty(outputs)
        C_expr = X.';          % y = x  (C = I, D = 0)
        D_expr = sym(0);
    else
        outs_sym = cellfun(@(g) subs(g, [Y(:).', cell2sym(inputs)], [X(:).', U(:).']), outputs, 'uni', 0);
        C_expr = [outs_sym{:}].';
        D_expr = sym(0);
    end

    % 5) Jacobiane simboliche
    A_sym = jacobian(f, X);
    B_sym = jacobian(f, U);
    C_sym = jacobian(C_expr, X);
    D_sym = jacobian(C_expr, U);

    % 6) Valutazione opzionale in un equilibrio
    have_eq = ~isempty(fieldnames(eq_point));
    if have_eq
        % param substitutions (se presenti)
        paramSubs = struct2cell(eq_point);
        % costruiamo una lista di sostituzioni parametriche (nome campo ≠ x0,u0)
        subsp = {};
        fn = fieldnames(eq_point);
        for k = 1:numel(fn)
            if ~ismember(fn{k}, {'x0','u0'})
                v = eq_point.(fn{k});
                subsp{end+1} = sym(fn{k}) == v; %#ok<AGROW>
            end
        end

        % default x0,u0 = 0
        x0 = isfield(eq_point,'x0') * eq_point.x0 + ~isfield(eq_point,'x0') * sym(zeros(n,1));
        u0 = isfield(eq_point,'u0') * eq_point.u0 + ~isfield(eq_point,'u0') * sym(zeros(m,1));

        A = simplify( subs(A_sym, [X;U;[subsp{:}].'], [x0;u0;repmat(sym(0),0,1)]) );
        B = simplify( subs(B_sym, [X;U;[subsp{:}].'], [x0;u0;repmat(sym(0),0,1)]) );
        C = simplify( subs(C_sym, [X;U;[subsp{:}].'], [x0;u0;repmat(sym(0),0,1)]) );
        D = simplify( subs(D_sym, [X;U;[subsp{:}].'], [x0;u0;repmat(sym(0),0,1)]) );
    else
        A = sym([]); B = sym([]); C = sym([]); D = sym([]);
    end

    % 7) Pack struct
    sys = struct('f',simplify(f), 'X',X, 'U',U, ...
                 'Yorig',Y, 'subsx',{subsx}, 'subsu',{subsu}, ...
                 'A_sym',simplify(A_sym), 'B_sym',simplify(B_sym), ...
                 'C_sym',simplify(C_sym), 'D_sym',simplify(D_sym), ...
                 'A',A, 'B',B, 'C',C, 'D',D);
end

function S = cell2sym(C)
    if isempty(C), S = sym([]); return; end
    S = [C{:}];
end