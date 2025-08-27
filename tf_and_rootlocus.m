function out = tf_and_rootlocus(A,B,C,D,varargin)
% TF_AND_ROOTLOCUS  (SISO)
% - Se (A,B,C,D) contengono simboli: ritorna G(s) simbolica + N(s), D(s), p_cl(s,K)
% - Se numeriche e PlotRootLocus=false: stampa G(s) numerica
% - Se numeriche e PlotRootLocus=true: disegna rlocus(G)
%
% Uso tipico (flusso che proponi):
%   % --- FASE 1: simbolico ---
%   out = tf_and_rootlocus(A,B,[1 0 0 0],0);  % no numeri -> solo analisi
%   % ora fai Routh su out.pcl ( = den + K*num ), ricavi condizioni su K,...
%   % --- FASE 2: numerico per il plot ---
%   params = struct('M',1,'B',0.8,'K',10,'D',0.5,'Fd',0.1,'Fq',0.2,'Kp',1,'Kd',2);
%   out = tf_and_rootlocus(A,B,[1 0 0 0],0,'Params',params,'PlotRootLocus',true);

    p = inputParser;
    addParameter(p,'Params',struct(),@(s)isstruct(s));
    addParameter(p,'PlotRootLocus',false,@islogical);
    addParameter(p,'Ksymbol',sym('K'),@(x) true);     % simbolo del guadagno
    parse(p,varargin{:});
    params        = p.Results.Params;
    plotRL        = p.Results.PlotRootLocus;
    Ksym          = p.Results.Ksymbol;

    if isvector(C), C = reshape(C,1,[]); end

    isSym = isa(A,'sym') || isa(B,'sym') || isa(C,'sym') || isa(D,'sym');

    % ---------- Caso con simboli: analisi senza plot ----------
    if isSym
        syms s real
        I = sym(eye(size(A,1)));
        Gsym = simplify(C*((s*I - A)^-1)*B + D);
        [num,den] = numden(Gsym);   % N(s), D(s)
        num = expand(num); den = expand(den);

        out = struct();
        out.Gsym = Gsym;
        out.num  = num;
        out.den  = den;
        out.pcl  = expand(den + Ksym*num);  % polinomio caratteristico in retroazione unitaria

        fprintf('G(s) simbolica:\n'); pretty(Gsym); fprintf('\n');
        fprintf('Polinomio chiuso p_cl(s,K) = den(s) + K*num(s):\n'); pretty(out.pcl); fprintf('\n');

        if plotRL
            error(['PlotRootLocus=true richiede valori numerici. ', ...
                   'Valuta i parametri (Params) dopo aver deciso le regioni stabili via Routh.']);
        end
        return
    end

    % ---------- Caso numerico ----------
    % Controllo SISO
    if size(B,2)~=1 || size(C,1)~=1 || ~isscalar(D)
        error('Gestisco solo SISO per il root locus (m=1, p=1, D scalare).');
    end

    % Se il chiamante ha passato Params ma le matrici sono già numeriche, ignoro Params.
    sys = ss(A,B,C,D);
    G   = tf(sys);

    out = struct('sys',sys,'G',G);

    fprintf('Funzione di trasferimento G(s):\n');
    disp(G);

    if plotRL
        figure; rlocus(G); grid on; title('Root Locus del sistema');
    end
end
