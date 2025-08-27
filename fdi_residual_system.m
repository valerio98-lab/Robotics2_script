function mdl = fdi_residual_system(opts)
% FDI residual generator (momentum-based) — with VERBOSE symbolic prints
%
% Modalità:
%  A) SIMBOLICO (M,g simbolici):
%     mdl = fdi_residual_model(struct('M',Msym,'g',gsym,'K',K, ...
%                                     'Fv',Fv,'Fs',Fs,'verbose',true));
%  B) NUMERICO (funzioni con S):
%     mdl = fdi_residual_model(struct('Mfun',@M,'gfun',@g,'Sfun',@S,'K',K, ...
%                                     'Fv',Fv,'Fs',Fs,'verbose',true));
%
% Campi utili:
%  mdl.alpha(q,qd)     % n×1
%  mdl.p(q,qd)         % n×1
%  mdl.init(q0,qd0,r0)
%  mdl.step(state,q,qd,tau,dt)
%  mdl.print(fmt)      % fmt='pretty' (default) o 'latex' — stampa simbolici
%  mdl.sym.*           % struct con M,g,p,alpha,q,qd simbolici (se disponibili)

% --- parsing ---------------------------------------------------------------
if ~isfield(opts,'K'), error('Provide K (diag vector or diagonal matrix)'); end
K = opts.K; if isvector(K), K = diag(K); end
n = size(K,1);
Fv = zeros(n,1); if isfield(opts,'Fv') && ~isempty(opts.Fv), Fv = opts.Fv(:); end
Fs = zeros(n,1); if isfield(opts,'Fs') && ~isempty(opts.Fs), Fs = opts.Fs(:); end
verbose = isfield(opts,'verbose') && opts.verbose;

use_sym_M = isfield(opts,'M') && isa(opts.M,'sym');
use_fun_S = isfield(opts,'Sfun') && ~isempty(opts.Sfun);

mdl = struct(); mdl.K = K; mdl.n = n; mdl.mode = ''; mdl.sym = struct();

if use_sym_M
    % ---------- modalità simbolica (usa alpha_i = -1/2 qd' dM/dq_i qd + g_i + frizioni) ----------
    Msym = opts.M; gsym = opts.g;
    qsym  = sym('q', [n 1], 'real');    %#ok<SYM>
    qdsym = sym('qd',[n 1], 'real');    %#ok<SYM>

    psym = Msym*qdsym;

    alpha_sym = sym(zeros(n,1));
    for i=1:n
        dMi = diff(Msym, qsym(i)); % n×n
        alpha_sym(i) = -sym(1)/2 * (qdsym.' * dMi * qdsym) ...
                       + gsym(i) + Fv(i)*qdsym(i) + Fs(i)*sign(qdsym(i));
    end

    mdl.alpha = matlabFunction(alpha_sym, 'Vars', {qsym,qdsym});
    mdl.p     = matlabFunction(psym,      'Vars', {qsym,qdsym});

    mdl.mode = 'symM';
    mdl.sym.M = Msym; mdl.sym.g = gsym;
    mdl.sym.q = qsym; mdl.sym.qd = qdsym;
    mdl.sym.p = psym; mdl.sym.alpha = alpha_sym;

elseif use_fun_S
    % ---------- modalità numerica (richiede S: C*qd = S(q,qd)*qd) ----------
    Mfun = opts.Mfun; gfun = opts.gfun; Sfun = opts.Sfun;
    mdl.alpha = @(q,qd) ( -Sfun(q,qd).'*qd + gfun(q) + Fv.*qd + Fs.*sign(qd) );
    mdl.p     = @(q,qd) ( Mfun(q)*qd );
    mdl.mode  = 'funS';

else
    error('Provide either {M (sym), g (sym)} or {Mfun,gfun,Sfun}.');
end

% --- init & step -----------------------------------------------------------
mdl.init = @(q0,qd0,r0) init_state(mdl,q0,qd0,r0);
mdl.step = @(state,q,qd,tau,dt) step_fdi(mdl,state,q,qd,tau,dt);

% --- verbose pretty printer ------------------------------------------------
mdl.print = @(fmt) print_symbolics(mdl, fmt);
if verbose, mdl.print('pretty'); end

end

% ================= helpers =================================================

function state = init_state(mdl,q0,qd0,r0)
p0 = mdl.p(q0,qd0);
if nargin<4 || isempty(r0), r0 = zeros(mdl.n,1); end
z0 = p0 + (mdl.K \ r0);
state.r = r0(:);
state.z = z0(:);
end

function [state,r] = step_fdi(mdl,state,q,qd,tau,dt)
alpha = mdl.alpha(q,qd);
state.z = state.z + dt*( tau(:) - alpha - state.r );
r = mdl.K * ( state.z - mdl.p(q,qd) );
state.r = r;
end

function print_symbolics(mdl, fmt)
if nargin<2 || isempty(fmt), fmt = 'pretty'; end
fprintf('\n==== FDI symbolic summary (mode: %s) ====\n', mdl.mode);
if ~isfield(mdl,'sym') || ~isfield(mdl.sym,'q')
    fprintf('(no symbolic data available in this mode)\n'); return;
end
q  = mdl.sym.q;  qd = mdl.sym.qd;
M  = mdl.sym.M;  g  = mdl.sym.g;
p  = mdl.sym.p;  a  = mdl.sym.alpha;

print_block('q (symbols)', q, fmt);
print_block('qd (symbols)', qd, fmt);
print_block('M(q)', M, fmt);
print_block('g(q)', g, fmt);
print_block('p(q,qd) = M(q) qd', p, fmt);
print_block('alpha(q,qd)', a, fmt);

% componenti alpha: separa i contributi
n = length(q);
T = sym(zeros(n,1));  % -1/2 qd' dM/dqi qd
for i=1:n
    dMi = diff(M, q(i));
    T(i) = -sym(1)/2 * (qd.'*dMi*qd);
end
print_block('alpha components: T_i = -1/2 qd^T (dM/dq_i) qd', T, fmt);
print_block('alpha components: g_i(q)', g, fmt);
print_block('alpha components: Fv_i * qd_i', diag_sym(mdl, qd, 'Fv'), fmt);
print_block('alpha components: Fs_i * sgn(qd_i)', diag_sym(mdl, sign(qd), 'Fs'), fmt);

fprintf('=================================================\n\n');
end

function S = diag_sym(mdl,vec,label)
% restituisce un vettore simbolico dei contributi Fv_i*vec_i o Fs_i*vec_i
if ~isfield(mdl,'K') || isempty(mdl.K); S = []; return; end
n = mdl.n;
Fv = zeros(n,1); Fs = zeros(n,1);
if isfield(mdl,'sym') && isfield(mdl.sym,'q') % solo per coerenza stampa
    % carica da struct opts? Non li abbiamo; ricostruiamo simbolici per stampa
    Fv = sym('Fv',[n,1]); Fs = sym('Fs',[n,1]);
else
    Fv = zeros(n,1); Fs = zeros(n,1);
end
if strcmp(label,'Fv')
    S = Fv .* vec;
else
    S = Fs .* vec;
end
end

function print_block(title_str, expr, fmt)
fprintf('\n-- %s --\n', title_str);
switch fmt
    case 'latex'
        disp(latex(simplify(expr)));
    otherwise
        pretty(simplify(expr));
end
end
