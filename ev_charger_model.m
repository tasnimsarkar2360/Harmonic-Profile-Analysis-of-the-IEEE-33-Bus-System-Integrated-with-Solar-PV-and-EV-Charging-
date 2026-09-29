function ev = ev_charger_model(t, cfg)

standalone = (nargin < 1) || isempty(t);
if nargin < 2 || isempty(cfg), cfg = struct(); end
if ~isstruct(cfg)
    error('ev_charger_model:cfg', 'cfg must be a struct (e.g. struct(''P_kW'',250)).');
end
noUserCfg = isempty(fieldnames(cfg)); 
if standalone
    f_tmp = 50;
    if isfield(cfg,'f0') && ~isempty(cfg.f0), f_tmp = cfg.f0; end
    t = (0:1/10000:3/f_tmp)';           % 3 cycles at 10 kHz
end

d.P_kW        = 300;
d.Q_kvar      = 100;
d.chargerType = 'diode6';
d.f0          = 50;
d.Vll_rated   = 12660;
d.phaseA0     = 0;
f = fieldnames(d);
for i = 1:numel(f)
    if ~isfield(cfg, f{i}) || isempty(cfg.(f{i})), cfg.(f{i}) = d.(f{i}); end
end

%  fundamental 
P = cfg.P_kW*1e3;  Q = cfg.Q_kvar*1e3;
S = hypot(P, Q);
Vph_rms = cfg.Vll_rated/sqrt(3);
I1_rms  = S/(3*Vph_rms);
I1_pk   = sqrt(2)*I1_rms;
phi1    = -atan2(Q, P);               

w0  = 2*pi*cfg.f0;
t   = t(:);
ang = w0*t + cfg.phaseA0;

% characteristic harmonics 

switch lower(cfg.chargerType)
    case 'diode6'
        harmOrders = [5  7  11 13 17 19 23 25];
        harmPct    = [23 11 7  5  2.5 2 1.2 1.0];
        harmSign   = [-1 +1 -1 +1 -1 +1 -1 +1];
    case 'pulse12'
        harmOrders = [11 13 23 25];
        harmPct    = [8  5.5 2  1.5];
        harmSign   = [-1 +1 -1 +1];
    case 'active'
        harmOrders = [5   7   11  13  17  19];
        harmPct    = [2.0 1.5 1.0 0.8 0.5 0.4];
        harmSign   = [-1 +1 -1 +1 -1 +1];
    otherwise
        error('ev_charger_model:type', ...
              'Unknown chargerType "%s" (use diode6 | pulse12 | active).', cfg.chargerType);
end
ev.harmonics = table(harmOrders(:), harmPct(:), harmSign(:), ...
    'VariableNames', {'Order','PercentOfFundamental','Sign'});

Iabc = zeros(numel(t), 3);
for ph = 1:3
    phShift = -(ph-1)*2*pi/3;                          
    th = ang + phi1 + phShift;                         
    wave = I1_pk * sin(th);
    for k = 1:numel(harmOrders)
        h  = harmOrders(k);
        Ih = (harmPct(k)/100) * I1_pk;
        if mod(h,3) == 0
            thh = h*(ang + phi1);                      
        else
            thh = h*th;                                
        end
        wave = wave + harmSign(k) * Ih * cos(thh - h*pi/2);   % cos(h(theta - 90deg))
    end
    Iabc(:,ph) = wave;
end

%  parallel R-L-C filter 
hasFilter = isfield(cfg,'filter') && ~isempty(cfg.filter);
if hasFilter
    if ~isfield(cfg,'Zsrc') || numel(cfg.Zsrc) ~= 2
        error('ev_charger_model:Zsrc', ...
              'cfg.Zsrc = [R X] (ohm at f0) is required when cfg.filter is given.');
    end
    flt = cfg.filter;
    Zs  = @(h) cfg.Zsrc(1) + 1j*h*cfg.Zsrc(2);
    Zf  = @(h) rlc_bank_Z(flt, h, cfg.f0);

    Iev1 = I1_rms*exp(1j*phi1);              
    If1  = Vph_rms / Zf(1);                  % filter fundamental current (capacitive)
    Ig1  = Iev1 + If1;                      

    Igrid = zeros(numel(t), 3);
    Ig_h  = 0;                               % sum of squared grid harmonic amplitudes (peak)
    for ph = 1:3
        phShift = -(ph-1)*2*pi/3;
        wave = sqrt(2)*abs(Ig1) * sin(ang + phShift + angle(Ig1));
        for k = 1:numel(harmOrders)
            h  = harmOrders(k);
            Gh = Zf(h) / (Zf(h) + Zs(h));                     % grid share of harmonic h
            Ih = (harmPct(k)/100) * I1_pk * abs(Gh);
            if mod(h,3) == 0
                thh = h*(ang + phi1);
            else
                thh = h*(ang + phi1 + phShift);
            end
            wave = wave + harmSign(k) * Ih * cos(thh - h*pi/2 + angle(Gh));
            if ph == 1, Ig_h = Ig_h + Ih^2; end
        end
        Igrid(:,ph) = wave;
    end
    Sg   = 3*Vph_rms*conj(Ig1);
    ev.Igrid       = Igrid;
    ev.Ifilter     = Igrid - Iabc;
    ev.Pgrid       = real(Sg);
    ev.Qgrid       = imag(Sg);
    ev.I1grid_rms  = abs(Ig1);
    ev.THD_grid    = 100*sqrt(Ig_h) / (sqrt(2)*abs(Ig1));
    ev.PFgrid      = ev.Pgrid / max(abs(Sg), eps);
    ev.filterLoss_W = 3*real(Vph_rms*conj(If1));
else
    ev.Igrid       = Iabc;
    ev.Ifilter     = zeros(size(Iabc));
    ev.Pgrid       = P;
    ev.Qgrid       = Q;
    ev.I1grid_rms  = I1_rms;
    ev.THD_grid    = 100*sqrt(sum((harmPct/100).^2));
    ev.PFgrid      = P/max(S, eps);
    ev.filterLoss_W = 0;
end
ev.PFtrue_grid = ev.PFgrid / sqrt(1 + (ev.THD_grid/100)^2);

ev.Pev = P;  ev.Qev = Q;  ev.Sev = S;
ev.I1_rms = I1_rms;
ev.Iabc   = Iabc;
ev.THD_I  = 100*sqrt(sum((harmPct/100).^2));
ev.PFdisp = P/max(S, eps);
ev.PFtrue = ev.PFdisp / sqrt(1 + (ev.THD_I/100)^2);
ev.cfgUsed = cfg;

% standalone demo output 
if standalone
    if noUserCfg, src = 'default'; else, src = 'the supplied'; end
    fprintf('[ev_charger_model] Standalone run with %s parameters.\n', src);
    fprintf('  Charger: %s | P = %.1f kW | Q = %.1f kVAr | S = %.1f kVA\n', ...
            cfg.chargerType, ev.Pev/1e3, ev.Qev/1e3, ev.Sev/1e3);
    fprintf('  I1 = %.2f A rms | current THD = %.2f %% | displacement PF = %.3f | true PF = %.3f\n', ...
            ev.I1_rms, ev.THD_I, ev.PFdisp, ev.PFtrue);
    if hasFilter
        fprintf('  With filter: grid-side THD = %.2f %%, grid PF = %.3f, P = %.1f kW, Q = %.1f kVAr\n', ...
                ev.THD_grid, ev.PFgrid, ev.Pgrid/1e3, ev.Qgrid/1e3);
    end
    figure('Name','ev_charger_model demo: drawn current waveform');
    plot(t, ev.Iabc, 'LineWidth', 1.1); grid on;
    xlabel('Time (s)'); ylabel('Current (A)'); legend('Ia','Ib','Ic');
    title(sprintf('EV charger current (%s, %.0f kW, THD = %.1f%%)', ...
                  cfg.chargerType, ev.Pev/1e3, ev.THD_I));
end
end

function Z = rlc_bank_Z(flt, h, f0)
% Impedance of the parallel bank of series R-L-C branches at harmonic order h.
    w  = 2*pi*f0;
    Zb = flt.R(:).' + 1j*(h*w*flt.L(:).' - 1./(h*w*flt.C(:).'));
    Z  = 1 / sum(1 ./ Zb);
end