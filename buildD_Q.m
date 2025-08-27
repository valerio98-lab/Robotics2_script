function D = buildD_Q(A, mode)
% buildD - Costruisce D a partire dal nullspace di A
% mode='rows' (default): D è (n-m) x n con righe ortonormali che span(Null(A))
% mode='cols'          : D è n x (n-m) con colonne ortonormali che span(Null(A))

    if nargin < 2, mode = 'rows'; end

    % Base ortonormale del Nullspace(A) nelle COLONNE
    N = null(A);                 % N è n x (n-m), colonne ortonormali, A*N = 0

    switch mode
        case 'rows'
            D = N.';             % righe di D = colonne di N → span(D^T) = Null(A)
        case 'cols'
            D = N;               % colonne di D = colonne di N → span(D)   = Null(A)
        otherwise
            error("mode deve essere 'rows' o 'cols'");
    end
end