function S = jacobianAnalysis(J, tol, vars)
%JACOBIANANALYSIS  Analyse a Jacobian (numeric or symbolic) + optional singularity report.
%   Returned fields
%   --------------
%     S.size           [m n]
%     S.isSquare       logical
%     S.determinant    det(J)  (numeric / symbolic / NaN)
%     S.rank           rank(J)
%     S.nullity        n‑rank
%     S.colBasis       basis of range(J)           (m × rank)
%     S.rowBasis       basis of range(Jᵀ)          (n × rank)
%     S.nullBasis      basis of null(J)            (n × nullity)
%     S.leftNullBasis  basis of null(Jᵀ)           (m × (m‑rank))
%     S.condNumber     2‑norm condition number (numeric J only)


%   S = jacobianAnalysis(J)                – automatic tolerance, no report.
%   S = jacobianAnalysis(J, tol)           – custom tolerance for numeric J.
%   S = jacobianAnalysis(J, tol, vars)     – plus symbolic singularity report.
%
%   *J*      : m‑by‑n Jacobian, numeric or symbolic.
%   *tol*    : numerical rank tolerance (ignored for symbolic J).
%   *vars*   : vector of symbolic variables; if supplied and J is symbolic,
%              the function returns a field S.singularity containing:
%                 .expr   – symbolic condition(s) for loss of rank
%                 .sol    – solutions for vars (via SOLVE)  ❲may be empty❳
%
%   Core diagnostics (always produced)
%   ----------------------------------
%     S.size, S.isSquare, S.determinant, S.rank, S.nullity
%     S.colBasis, S.rowBasis, S.nullBasis, S.leftNullBasis
%     S.condNumber   (numeric)
%
%   © 2025 – MIT Licence.

% ------------------------------------------------------------
% Handle optional inputs
% ------------------------------------------------------------
if nargin < 2 || isempty(tol)
    if ~isa(J,'sym')
        tol = max(size(J))*eps(norm(J,'fro'));
    else
        tol = sym('0');
    end
end
if nargin < 3
    vars = [];
end

validateattributes(J,{'numeric','sym'},{'2d'},mfilename,'J',1);
[m,n]  = size(J);
isSym  = isa(J,'sym');

S             = struct();
S.size        = [m n];
S.isSquare    = (m==n);
S.condNumber  = NaN;         % filled later for numeric

% ------------------------------------------------------------
% Symbolic processing
% ------------------------------------------------------------
if isSym
    % Determinant / rank
    if S.isSquare
        S.determinant = simplify(det(J));
    else
        S.determinant = sym('NaN');
    end
    S.rank    = double(rank(J));
    S.nullity = n - S.rank;

    % Bases via Symbolic Math Toolbox
    S.colBasis      = colspace(J);
    S.rowBasis      = colspace(J.');
    S.nullBasis     = null(J);
    S.leftNullBasis = null(J.');

    % --------------------------------------------------------
    % Singularity report if variables supplied
    % --------------------------------------------------------
    if ~isempty(vars)
        sing = struct();
        r = S.rank;
        if S.isSquare
            sing.expr = simplify(det(J));  % full det
        else
            % Use Gram determinant (det(J*J^T)) – vanishes ⇔ rank<J
            if m <= n
                sing.expr = simplify(det(J*J.'));
            else
                sing.expr = simplify(det(J.'*J));
            end
        end
        try
            sing.sol = solve(sing.expr == 0, vars, 'ReturnConditions', true);
        catch
            sing.sol = [];  % maybe over‑determined or too many vars
        end
        S.singularity = sing;
    end

    % Optional display
    if nargout == 0
        fprintf('Jacobian analysis (symbolic, %dx%d)\n',m,n);
        fprintf('  Rank        : %d\n',S.rank);
        if S.isSquare
            fprintf('  Determinant : %s\n',char(S.determinant)); end
        if isfield(S,'singularity')
            fprintf('  Singularity condition: %s = 0\n',char(S.singularity.expr));
        end
    end
    return;
end

% ------------------------------------------------------------
% Numeric processing (unchanged)
% ------------------------------------------------------------
if S.isSquare
    S.determinant = det(J);
else
    S.determinant = NaN;
end

S.rank    = rank(J, tol);
S.nullity = n - S.rank;

[U,~,V] = svd(J,'econ');
S.colBasis        = U(:,1:S.rank);
S.rowBasis        = V(:,1:S.rank);
S.nullBasis       = V(:,S.rank+1:end);
S.leftNullBasis   = U(:,S.rank+1:end);

if S.isSquare && S.rank == n
    S.condNumber = cond(J);
else
    S.condNumber = Inf;
end

if nargout == 0
    fprintf('Jacobian analysis (numeric, %dx%d)\n',m,n);
    fprintf('  Rank          : %d\n',S.rank);
    if S.isSquare
        fprintf('  Determinant   : %g\n',S.determinant);
        if isfinite(S.condNumber)
            fprintf('  Condition no. : %g\n',S.condNumber);
        else
            fprintf('  Condition no. : Inf (rank-deficient)\n');
        end
    end
    fprintf('  Nullity       : %d\n',S.nullity);
end
