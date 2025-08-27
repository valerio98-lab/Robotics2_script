function f = equations_to_state(eqns, states, derivs, inputs)
    if nargin < 4, inputs = []; end

    n = numel(states);
    h = numel(derivs);

    % Controlli di consistenza
    assert(mod(n,2)==0, 'Il numero di stati deve essere pari: [pos; vel].');
    assert(h==n/2, 'Il numero di derivate da risolvere deve essere n/2.');
    
    % Risolvi accelerazioni/derivate richieste
    sol = solve(eqns, derivs);  % struct con campi come i simboli in "derivs"

    f = sym(zeros(n,1));
    % Prima metà: posizioni -> derivate sono le velocità (seconda metà di "states")
    for i = 1:h
        f(i) = states(h+i);
    end
    % Seconda metà: velocità -> derivate sono le accelerazioni risolte
    for i = 1:h
        key = char(derivs(i));
        assert(isfield(sol, key), 'solve non ha prodotto il campo "%s".', key);
        f(h+i) = simplify(sol.(key));
    end

    % (Opzionale) sostituisci input non usati o genera funzione numerica:
    % f_num = matlabFunction(f, 'Vars', {states, inputs});
end
