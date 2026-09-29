function pv = pv_cell_model(t, cfg)

if nargin < 1 || isempty(t)
    t = (0:1/10000:0.06)';          % 3 cycles 50 Hz, 10 kHz sample rate
    if nargin < 2, cfg = struct(); end

    fprintf('[pv_cell_model] Interactive mode -- press Enter to accept each default.\n');
    Gin = input('  Irradiance G (W/m^2) [Enter = 1000]: ');
    if ~isempty(Gin), cfg.G = Gin; end

    Tcin = input('  Cell temperature Tc (degC) [Enter = 25]: ');
    if ~isempty(Tcin), cfg.Tc = Tcin; end

    Pin = input('  Desired PV active power output at this G/Tc (kW) [Enter = default array size]: ');
    if ~isempty(Pin), cfg.targetPpv_kW = Pin; end

    Sw = input('  Switch type - PWM / diode / thyristor [PWM]: ', 's');
    if ~isempty(Sw), cfg.switchType = Sw; end

    fprintf('\n');
    pv = pv_cell_model_core(t, cfg);
    fprintf('  Irradiance G = %.0f W/m^2, Temperature Tc = %.1f degC\n', pv.cfgUsed.G, pv.cfgUsed.Tc);
    fprintf('  Array sizing used: Nmod=%d (series), Np=%.3f (parallel, effective)\n', pv.cfgUsed.Nmod, pv.cfgUsed.Np);
    fprintf('  Pdc = %.1f W | Vdc = %.1f V | Idc = %.2f A\n', pv.Pdc, pv.Vdc, pv.Idc);
    fprintf('  Ppv = %.2f kW | Qpv = %.1f VAR | Injected-current THD = %.2f %%\n', pv.Ppv/1000, pv.Qpv, pv.THD_I);
    figure('Name','pv_cell_model demo: injected current waveform');
    plot(t, pv.Iabc, 'LineWidth', 1.1); grid on;
    xlabel('Time (s)'); ylabel('Current (A)'); legend('Ia','Ib','Ic');
    title(sprintf('PV Injected Current (G=%.0f W/m^2, Tc=%.1f degC, %s)', ...
        pv.cfgUsed.G, pv.cfgUsed.Tc, pv.cfgUsed.switchType));
    return
end
pv = pv_cell_model_core(t, cfg);
end

function pv = pv_cell_model_core(t, cfg)


d.G          = 1000;     % irradiance (W/m^2)
d.Tc         = 25;       % cell temperature (deg C)
d.Gref       = 1000;     d.Tref = 25;
d.Ns         = 60;       % cells in series per module
d.Np         = 1;        % parallel strings 
d.Nmod       = 20;       % modules in series per string
d.Voc_cell   = 0.721;    % V, per-cell open-circuit voltage
d.Isc_cell   = 8.21;     % A, per-cell short-circuit current
d.Ki         = 0.0032;   % A/degC, Isc temperature coefficient
d.Kv         = -0.0023;  % V/degC (per volt), Voc temperature coefficient
d.n          = 1.3;      % diode ideality factor
d.Rs         = 0.004;    % ohm, series resistance (per cell, lumped)
d.Rsh        = 20;       % ohm, shunt resistance (per cell, lumped)
d.f0         = 50;       % Hz, fundamental grid frequency
d.Vll_rated  = 12660;    % V, nominal line-to-line voltage at POI (bus 18)
d.PFcmd      = 0.98;     % inverter power-factor command
d.switchType = 'PWM';   
d.phaseA0    = 0;        % rad, phase of phase-A POI voltage (sin reference)

f = fieldnames(d);
if nargin < 2, cfg = struct(); end
for i = 1:numel(f)
    if ~isfield(cfg, f{i}), cfg.(f{i}) = d.(f{i}); end
end


% 1. DC SIDE -- single-diode model, solved for the MPP

q  = 1.602176634e-19;  k = 1.380649e-23;
Tc_K  = cfg.Tc + 273.15;
Vt    = cfg.n * cfg.Ns * k * Tc_K / q;

Iph = (cfg.Isc_cell + cfg.Ki*(cfg.Tc - cfg.Tref)) * (cfg.G/cfg.Gref);
Voc_string = cfg.Ns * (cfg.Voc_cell + cfg.Kv*(cfg.Tc - cfg.Tref));
Io  = Iph / (exp(Voc_string/Vt) - 1);

Vsweep = linspace(0.01, Voc_string, 400);
Isweep = zeros(size(Vsweep));
Iguess = Iph;
for j = 1:numel(Vsweep)
    V = Vsweep(j);
    I = Iguess;
    for it = 1:30
        expTerm = exp((V + I*cfg.Rs)/Vt);
        F  = Iph - Io*(expTerm - 1) - (V + I*cfg.Rs)/cfg.Rsh - I;
        dF = -Io*(cfg.Rs/Vt)*expTerm - cfg.Rs/cfg.Rsh - 1;
        Inew = I - F/dF;
        if abs(Inew - I) < 1e-8, I = Inew; break; end
        I = Inew;
    end
    Isweep(j) = max(I, 0);
    Iguess = I;
end

Psweep = Vsweep .* Isweep;
[Pdc_string, idxMPP] = max(Psweep); 
Vdc_string = Vsweep(idxMPP);
Idc_string = Isweep(idxMPP);

pv.Vdc = Vdc_string * cfg.Nmod;
eta_inv  = 0.97;

if isfield(cfg, 'targetPpv_kW') && ~isempty(cfg.targetPpv_kW)
    
    Pdc_target = (cfg.targetPpv_kW*1000) / eta_inv;
    Pdc_perNp  = Vdc_string * cfg.Nmod * Idc_string;
    cfg.Np     = Pdc_target / max(Pdc_perNp, eps);
end

pv.Idc = Idc_string * cfg.Np;
pv.Pdc = pv.Vdc * pv.Idc;


% 2. AC SIDE -- inverter power command

pv.Ppv   = eta_inv * pv.Pdc;                                   % W
pv.Qpv   = pv.Ppv * tan(acos(min(max(cfg.PFcmd,-1),1)));       % VAR 

Vph_rms  = cfg.Vll_rated / sqrt(3);
Sinv     = hypot(pv.Ppv, pv.Qpv);
I1_rms   = Sinv / (3 * Vph_rms);
I1_pk    = sqrt(2) * I1_rms;
phi1     = -atan2(pv.Qpv, pv.Ppv);            % current angle rel. to voltage

w0 = 2*pi*cfg.f0;
t  = t(:);
ang = w0*t + cfg.phaseA0;                     % phase-A voltage reference angle
pv.Vabc = Vph_rms*sqrt(2) * [ sin(ang), ...
                               sin(ang - 2*pi/3), ...
                               sin(ang + 2*pi/3) ];


% 3. HARMONICS

switch lower(cfg.switchType)
    case 'diode'
        harmOrders = [5  7  11 13];
        harmPct    = [23 11 7  5];
    case 'thyristor'
        harmOrders = [5  7  11 13 17 19];
        harmPct    = [28 18 9  7  4  3];
    otherwise            % 'PWM'
        harmOrders = [3   5   7   11  13];
        harmPct    = [1.0 1.5 1.0 0.8 0.6];
end
harmPhase = zeros(size(harmOrders));          
pv.harmonics = table(harmOrders(:), harmPct(:), harmPhase(:), ...
    'VariableNames', {'Order','PercentOfFundamental','Phase_rad'});

pv.Iabc = zeros(numel(t), 3);
for ph = 1:3
    phShift = -(ph-1)*2*pi/3;                          
    
    wave = I1_pk * sin(ang + phi1 + phShift);
    for hIdx = 1:numel(harmOrders)
        h     = harmOrders(hIdx);
        Ih_pk = (harmPct(hIdx)/100) * I1_pk;
        if mod(h,3) == 0
            hPhShift = 0;                              
        else
            hPhShift = phShift * h;                    
        end
        
        wave = wave + Ih_pk * sin(h*ang + h*phi1 + harmPhase(hIdx) + hPhShift);
    end
    pv.Iabc(:,ph) = wave;
end


pv.THD_I = 100 * sqrt(sum((harmPct/100).^2 * (I1_pk^2))) / I1_pk;

pv.cfgUsed = cfg;   

end