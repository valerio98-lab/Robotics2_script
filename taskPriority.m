function [qdot, taskErrors, info] = taskPriority(Js, rdotList, varargin)

%% ESEMPI: 

% J1 = double(subs(J1, [q1,q2,q3,q4,q5], [0, pi/3, -pi/4, pi/2, -pi/3]));
% J2 = double(subs(J2, [q1,q2,q3], [0, pi/3, -pi/4]));
% J3 = double(subs(J3, [q1,q2], [0, pi/3]));
% 
% 
% % Esempio con 3 task
% r1dot = [2;3];
% r2dot = [2;-0.5];
% r3dot = [-1;0];
% 
% [qdot, errs, info] = taskPriority({J1,J2,J3}, {r1dot,r2dot,r3dot}, ...
%     'Verbose', true, 'ShowProjector', false, ...
%     'SigmaThresh', 1e-3, 'LambdaMax', 1e-1, 'PrintFinalErrors', true);

% qdot  -> velocità di giunto finale
% errs{k} -> errore residuo per il task k (dopo la soluzione)

%%
% PRIORITIZEDIK  Prioritized resolved-rate IK (multi-task) con Null-Space Projection e DLS.
%
%   [qdot, taskErrors, info] = prioritizedIK(Js, rdotList, Name,Value,...)
%
% INPUT
%   Js         : cell array {J1, J2, ..., JK}, Jk is mk-by-n
%   rdotList   : cell array {r1dot, r2dot, ..., rKdot}, rkdot is mk-by-1
%
% Name-Value options:
%   'qdot0'            : n-by-1 vettore iniziale (default: zeros(n,1))
%   'Verbose'          : true/false stampa breve per iterazione (default: false)
%   'ShowProjector'    : true/false stampa P_A,k ad ogni iterazione (default: false)
%   'SigmaThresh'      : soglia per attivare DLS su sigma_min(Ak) (default: 1e-3)
%   'LambdaMax'        : lambda massimo usato in DLS (default: 1e-1)
%   'ReturnLogs'       : true/false salva storia in info (default: true)
%
% OUTPUT
%   qdot        : n-by-1, velocità finale di giunto
%   taskErrors  : cell array {e1, e2, ..., eK} con ek = rdot_k - Jk*qdot
%   info        : struct con log (singular values, lambda, rank, P history, qdot history)
%
% NOTE
%   - DLS: A#_lambda = A' * (A*A' + lambda^2 I)^(-1). Lambda scelto automaticamente
%          in base a sigma_min(A) con una rampa continua tra 0 e LambdaMax.
%   - Con DLS, A#_lambda * A non è un proiettore idempotente perfetto: è un
%     "quasi-proiettore" che stabilizza vicino alle singolarità.

% --------- parsing input ---------
K = numel(Js);
assert(iscell(Js) && iscell(rdotList) && numel(rdotList)==K, ...
    'Js e rdotList devono essere cell array della stessa lunghezza.');

n = size(Js{1},2);
for k = 1:K
    assert(size(Js{k},2)==n, 'Tutte le J_k devono avere n colonne coerenti.');
    assert(isvector(rdotList{k}) && size(rdotList{k},1)==size(Js{k},1), ...
        'rdot{k} deve essere mk-by-1 coerente con J{k}.');
end

p = inputParser;
p.addParameter('qdot0', zeros(n,1), @(x)isnumeric(x)&&isvector(x)&&numel(x)==n);
p.addParameter('Verbose', false, @(x)islogical(x)||ismember(x,[0,1]));
p.addParameter('ShowProjector', false, @(x)islogical(x)||ismember(x,[0,1]));
p.addParameter('SigmaThresh', 1e-3, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('LambdaMax', 1e-1, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('ReturnLogs', true, @(x)islogical(x)||ismember(x,[0,1]));
p.addParameter('PrintFinalErrors', false, @(x)islogical(x)||ismember(x,[0,1]));

p.parse(varargin{:});

qdot = p.Results.qdot0(:);
Verbose = p.Results.Verbose;
ShowP  = p.Results.ShowProjector;
sigThresh = p.Results.SigmaThresh;
lambdaMax = p.Results.LambdaMax;
retLogs   = p.Results.ReturnLogs;

% --------- init ---------
PA = eye(n);               % P_{A,0}
if retLogs
    info.qdot_hist = zeros(n, K+1); info.qdot_hist(:,1) = qdot;
    info.P_hist = cell(1,K+1); info.P_hist{1} = PA;
    info.lambda = zeros(1,K);
    info.rankA  = zeros(1,K);
    info.sigma_min = zeros(1,K);
    info.residual_norm = zeros(1,K);
    info.residual_proj_norm = zeros(1,K);
end

% --------- loop sui task ---------
for k = 1:K
    Jk = Js{k};
    rk = rdotList{k}(:);

    % Jacobiano ristretto (consentito)
    A = Jk * PA;

    % residuo in task space
    b = rk - Jk*qdot;

    % DLS adattiva in base a sigma_min(A)
    if all(A==0)
        lambda = 0; smin = 0; rA = 0;
    else
        s = svd(A, 'econ');
        smin = s(end);
        rA = sum(s > max(eps, 0));  %#ok<NASGU> % rango numerico
        if smin >= sigThresh
            lambda = 0;
        else
            % rampa continua: lambda cresce da 0 a lambdaMax quando smin scende a 0
            lambda = lambdaMax * (1 - smin/max(sigThresh, eps));
        end
    end

    % pseudoinversa (DLS se lambda>0)
    Asharp = dampedPinv(A, lambda);

    % correzione a norma minima dentro lo spazio ammesso
    dq = Asharp * b;

    % update velocità
    qdot = qdot + dq;

    % aggiorna (quasi-)proiettore del null space aumentato
    PA = PA - Asharp*A;  % NB: con DLS non è idempotente perfetto, ma stabilizza

    % logging / stampa
    if retLogs
        info.qdot_hist(:,k+1) = qdot;
        info.P_hist{k+1} = PA;
        info.lambda(k) = lambda;
        info.sigma_min(k) = smin;
        info.rankA(k) = rank(A);
        info.residual_norm(k) = norm(b);
        info.residual_proj_norm(k) = norm(A*Asharp*b); % parte realizzabile
    end

    if Verbose
        fprintf('--- Task %d ---\n', k);
        fprintf('sigma_min(A)=%.3e, lambda=%.3e, rank(A)=%d\n', smin, lambda, rank(A));
        fprintf('||residuo||=%.3e, ||residuo_realizzabile||=%.3e\n', ...
            norm(b), norm(A*Asharp*b));
        fprintf('qdot (dopo task %d):\n', k); disp(qdot.');
        if ShowP
            fprintf('P_A (dopo task %d):\n', k); disp(PA);
        end
    end
end

% --------- errori finali per ciascun task ---------
taskErrors = cell(1,K);
for k = 1:K
    taskErrors{k} = rdotList{k}(:) - Js{k}*qdot;
end


if p.Results.PrintFinalErrors
    for k = 1:K
        fprintf('Errore finale task %d: ||r_kdot - J_k*qdot|| = %.6e\n', k, norm(taskErrors{k}));
    end
end

% --------- info finale ---------
if retLogs
    info.PA_final = PA;
end
end

% ===== helper: pseudoinversa smorzata =====
function Asharp = dampedPinv(A, lambda)
% Damped least-squares pseudoinverse: A'*(A*A' + lambda^2 I)^(-1)
    [m,~] = size(A);
    if lambda > 0
        Asharp = (A') / (A*A' + (lambda^2)*eye(m));
    else
        % Moore-Penrose standard via SVD (più robusta di pinv di default)
        % (per matrici piccole/medie va bene; in produzione si può usare pinv.)
        [U,S,V] = svd(A, 'econ');
        s = diag(S);
        tol = max(size(A)) * eps(max(s));
        sInv = zeros(size(s));
        sInv(s > tol) = 1./s(s > tol);
        Asharp = V*diag(sInv)*U';
    end
end

