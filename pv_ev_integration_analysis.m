

%% 0. Choose the operation mode 
fprintf('\n=== Select operation mode ===\n');
fprintf('  1 : PV only integration\n');
fprintf('  2 : EV only integration\n');
fprintf('  3 : Both PV and EV integration (combined)\n');
opMode = round(ask_number('  Mode', 3, 1, 3));
pvOn   = (opMode == 1) || (opMode == 3);
evOn   = (opMode == 2) || (opMode == 3);
modeNames = {'PV only','EV only','PV + EV (combined)'};
fprintf('Running mode %d: %s\n', opMode, modeNames{opMode});

useSimulink  = false;    

% PV plants 
if pvOn
    
    pvBus = ask_bus_list('PV bus number(s), e.g. 18 or [18 33]', 18);
else
    pvBus = zeros(1,0);
end
nPV   = numel(pvBus);

G_in = zeros(1,nPV);  Tc_in = zeros(1,nPV);  PV_in = zeros(1,nPV);
for m = 1:nPV
    fprintf('\n--- PV plant %d of %d at Bus %d (Enter = default) ---\n', m, nPV, pvBus(m));
    G_in(m)  = ask_number('  1) Irradiance G (W/m^2)           ', 1000, 20, 1500);
    Tc_in(m) = ask_number('  2) Cell temperature Tc (degC)     ',   25, -20,  100);
    PV_in(m) = ask_number('  3) Desired PV power (kW)          ',  200, 1,  3000);
end

%  EV charging station 
filterOn = false;  nBr = 2;
if evOn
    fprintf('\n--- EV charging station (Enter = default) ---\n');
    evBus  = round(ask_number('  EV bus number                      ', 20,   2,   33));
    P_ev   = ask_number('  P_EV  active demand (kW)           ', 300,   1, 5000);
    Q_ev   = ask_number('  Q_EV  reactive demand (kVAr, +lag) ', 100, -5000, 5000);
    evSel  = round(ask_number('  Charger: 1=6-pulse diode, 2=12-pulse, 3=active front-end', 1, 1, 3));
    evTypes = {'diode6','pulse12','active'};
    evType  = evTypes{evSel};

    filterOn = round(ask_number('  Design & connect an R-L-C harmonic filter at the EV bus? (1=yes, 0=no)', 1, 0, 1)) == 1;
    if filterOn
        nBr = round(ask_number('  Number of filter branches (tuned to the dominant harmonics)', 2, 1, 3));
    end
end

%% 1. Baseline copy 
Vpu_pre = Vpu;  P_pu_pre = P_pu;  Q_pu_pre = Q_pu;
S_pu_pre = S_pu; PF_pre = PF;
V3_pre = V3;    I3_pre = I3;
kPV = zeros(1,nPV);
for m = 1:nPV
    kk = find(busTo == pvBus(m), 1);
    assert(~isempty(kk), 'Bus %d not found in busTo.', pvBus(m));
    kPV(m) = kk;                                    
end

%% 2. Build one PV plant per bus (each from pv_cell_model) 
pvList   = cell(1,nPV);
Ppv_kW   = zeros(1,nPV);   Qpv_kvar = zeros(1,nPV);
for m = 1:nPV
    k  = kPV(m);
    tt = t(idxSteady);
    c  = 2/numel(tt) * sum(V3_pre(idxSteady,1,k) .* exp(-1j*2*pi*f0*tt));
    phaseA0 = angle(c) + pi/2;                      % sin-reference phase of Va at THIS bus

    cfg = struct('G',G_in(m), 'Tc',Tc_in(m), 'targetPpv_kW',PV_in(m), ...
                 'f0',f0, 'Vll_rated',V_LL_base*Vpu_pre(k), 'phaseA0',phaseA0);
    pvm = pv_cell_model(t, cfg);                    % single-diode PV + inverter model
    pvList{m}   = pvm;
    Ppv_kW(m)   = pvm.Ppv/1000;
    Qpv_kvar(m) = pvm.Qpv/1000;

    fprintf('\n--- PV @ Bus %d ---\n', pvBus(m));
    fprintf('  G = %.0f W/m^2, Tc = %.1f degC, Nmod = %d, Np = %.2f, %s\n', ...
            pvm.cfgUsed.G, pvm.cfgUsed.Tc, pvm.cfgUsed.Nmod, pvm.cfgUsed.Np, pvm.cfgUsed.switchType);
    fprintf('  Pdc = %.1f kW | Vdc = %.0f V | Idc = %.1f A\n', pvm.Pdc/1e3, pvm.Vdc, pvm.Idc);
    fprintf('  Ppv = %.1f kW | Qpv = %.1f kVAr | Spv = %.1f kVA | I-THD = %.2f %%\n', ...
            Ppv_kW(m), Qpv_kvar(m), hypot(Ppv_kW(m),Qpv_kvar(m)), pvm.THD_I);

    P_check = sum(mean(V3_pre(idxSteady,:,k) .* pvm.Iabc(idxSteady,:), 1)) / 1e3;
    fprintf('  Sign/phase check: v.i at Bus %d = %.1f kW (expect ~ +%.0f kW)\n', ...
            pvBus(m), P_check, Ppv_kW(m));
    if P_check < 0.8*Ppv_kW(m)
        warning(['PV power check failed at Bus %d. Flip the sign of pv.Iabc and/or check ' ...
                 'that your abc phase sequence matches the model.'], pvBus(m));
    end
end
if nPV > 0
    fprintf('\nTotal PV: %.1f kW, %.1f kVAr  (%d plant(s))\n', sum(Ppv_kW), sum(Qpv_kvar), nPV);
end

%% 3. Network solve: pre vs post PV (all plants superimposed) 
Vbfs_pre  = ieee33_bfs([],    [],       [],         V_LL_base, S_base);   
Vbfs_post = ieee33_bfs(pvBus, Ppv_kW,   Qpv_kvar,   V_LL_base, S_base);

dV_bus   = Vbfs_post(busTo) - Vbfs_pre(busTo);      % change per bus (p.u.)
Vpu_post = Vpu_pre + dV_bus(:);                     

fprintf('\n--- Load-flow check (fundamental, constant-PQ loads) ---\n');
for m = 1:nPV
    fprintf('  Bus %2d : %.4f -> %.4f p.u.  (%+.4f)\n', pvBus(m), ...
            Vbfs_pre(pvBus(m)), Vbfs_post(pvBus(m)), Vbfs_post(pvBus(m))-Vbfs_pre(pvBus(m)));
end
[vmPre,  bPre ] = min(Vbfs_pre(2:end));   [vmPost, bPost] = min(Vbfs_post(2:end));
fprintf('  Weakest bus: %d (%.4f) -> %d (%.4f)\n', bPre+1, vmPre, bPost+1, vmPost);
if max(Vbfs_post) > 1.05
    warning('Overvoltage: max bus voltage %.4f p.u. (> 1.05). PV is too large for this feeder.', max(Vbfs_post));
end

%% 4. Post-PV branch currents
I3_post = I3_pre;
for m = 1:nPV
    pathm = branch_path_to_bus(busTo, pvBus(m));    
    for kk = pathm
        I3_post(:,:,kk) = I3_post(:,:,kk) - pvList{m}.Iabc;
    end
end

V3_post = V3_pre;
for k = 1:numBranches
    V3_post(:,:,k) = V3_pre(:,:,k) * (Vpu_post(k)/Vpu_pre(k));
end

%%  4b.full Simulink re-solve 

if useSimulink
    Tsim = 50e-6;  Tstop = 0.2;                    
    tpv  = (0:Tsim:Tstop)';
    ramp = min(tpv/0.04, 1);                       
    simIn2 = Simulink.SimulationInput(simName);
    simIn2 = simIn2.setVariable('siminCB', siminCB);
    for m = 1:nPV
        cfgm = pvList{m}.cfgUsed;                   
        pvS  = pv_cell_model(tpv, cfgm);
        simIn2 = simIn2.setVariable(sprintf('pvInjection_%d', pvBus(m)), ...
                                    timeseries(pvS.Iabc .* ramp, tpv));
    end
    simout2 = sim(simIn2);
    [Vpu_post, V3_post, I3_post, P_pu_post, Q_pu_post, S_pu_post, PF_post] = ...
        process_run(simout2.simOutputV_inst, simout2.simOutputI_inst, ...
                    numBranches, idxSteady, V_base, S_base);
end

%% 5. Post-PV P, Q, S, PF from waveforms
if ~useSimulink
    [~, ~, ~, P_pu_post, Q_pu_post, S_pu_post, PF_post] = ...
        process_run_from3(V3_post, I3_post, numBranches, idxSteady, V_base, S_base);
end
% keep the PV-only results
Vpu_pv = Vpu_post;  P_pu_pv = P_pu_post;  Q_pu_pv = Q_pu_post;  PF_pv = PF_post;
V3_pv  = V3_post;   I3_pv   = I3_post;

%%  5b. EV charging station (non-linear load) on top of PV 
Vpu_all = Vpu_pv;  V3_all = V3_pv;  I3_all = I3_pv;
if evOn
    kEV = find(busTo == evBus, 1);
    assert(~isempty(kEV), 'EV bus %d not found in busTo.', evBus);

    tt  = t(idxSteady);
    c   = 2/numel(tt) * sum(V3_pre(idxSteady,1,kEV) .* exp(-1j*2*pi*f0*tt));
    cfgEV = struct('P_kW',P_ev, 'Q_kvar',Q_ev, 'chargerType',evType, ...
                   'f0',f0, 'Vll_rated',V_LL_base*Vpu_pre(kEV), 'phaseA0',angle(c)+pi/2);
    ev = ev_charger_model(t, cfgEV);

    fprintf('\n--- EV charging station @ Bus %d (%s) ---\n', evBus, evType);
    fprintf('  P = %.1f kW | Q = %.1f kVAr | S = %.1f kVA | I1 = %.2f A rms\n', ...
            ev.Pev/1e3, ev.Qev/1e3, ev.Sev/1e3, ev.I1_rms);
    fprintf('  Current THD = %.2f %% | displacement PF = %.3f | true PF = %.3f\n', ...
            ev.THD_I, ev.PFdisp, ev.PFtrue);
    Pchk = sum(mean(V3_pre(idxSteady,:,kEV) .* ev.Iabc(idxSteady,:), 1)) / 1e3;
    fprintf('  Sign/phase check: v.i at Bus %d = %.1f kW (expect ~ +%.0f kW drawn)\n', evBus, Pchk, ev.Pev/1e3);
    if Pchk < 0.8*ev.Pev/1e3
        warning('EV power check failed: check the abc phase sequence.');
    end

    % Load flow + currents with the EV (grid-side net P/Q and current; no filter yet)
    [Vpu_all, V3_all, I3_all] = add_ev_to_network(ev, evBus, pvBus, Ppv_kW, Qpv_kvar, ...
        Vbfs_post, Vpu_pv, V3_pv, I3_pv, busTo, V_LL_base, S_base);
end
[~, ~, ~, P_pu_all, Q_pu_all, S_pu_all, PF_all] = ...
    process_run_from3(V3_all, I3_all, numBranches, idxSteady, V_base, S_base);

% comparative summary table 
if pvOn && evOn, allName = 'PV + EV'; elseif evOn, allName = 'EV only'; else, allName = 'PV only'; end
caseNames = {'Baseline'};  Vcase = Vpu_pre(:);  Pcase = P_pu_pre(1)*S_base/1e3;
if pvOn
    caseNames{end+1} = 'PV only';  Vcase(:,end+1) = Vpu_pv(:);  Pcase(end+1) = P_pu_pv(1)*S_base/1e3;
end
if evOn
    caseNames{end+1} = allName;    Vcase(:,end+1) = Vpu_all(:); Pcase(end+1) = P_pu_all(1)*S_base/1e3;
end
fprintf('\n=== Comparison summary (mode %d: %s) ===\n', opMode, modeNames{opMode});
fprintf('%-10s | %8s | %6s | %8s | %6s | %-14s | %s\n', 'Case','Vmin(pu)','@Bus','Vmax(pu)','@Bus','P from source','Buses < 0.95');
for cix = 1:numel(caseNames)
    [vmin, imin] = min(Vcase(:,cix));  [vmax, imax] = max(Vcase(:,cix));
    fprintf('%-10s | %8.4f | %6d | %8.4f | %6d | %9.1f kW   | %d\n', caseNames{cix}, ...
            vmin, busTo(imin), vmax, busTo(imax), Pcase(cix), sum(Vcase(:,cix) < 0.95));
end
if min(Vpu_all) < 0.95
    warning('Voltage below 0.95 p.u. (minimum %.4f p.u.) in case "%s".', min(Vpu_all), allName);
end

% PF < 0 means ACTIVE power flows toward the substation (reverse flow), not a bad PF.
rev = find(P_pu_all < 0);
if ~isempty(rev)
    fprintf('\nReverse active-power flow on %d branch(es) (PF sign flips negative):\n', numel(rev));
    for r = rev(:)'
        fprintf('  branch -> Bus %2d : P = %+7.1f kW, Q = %7.1f kVAr (magnitude), |PF| = %.3f\n', ...
                busTo(r), P_pu_all(r)*S_base/1e3, Q_pu_all(r)*S_base/1e3, abs(PF_all(r)));
    end
end

%%  6. Plots
lbl = strjoin(arrayfun(@(m) sprintf('Bus %d: %.0f kW', pvBus(m), Ppv_kW(m)), ...
              1:nPV, 'UniformOutput', false), ', ');

sets = struct('name',{'Baseline'}, 'V',{Vpu_pre}, 'P',{P_pu_pre}, 'PF',{PF_pre}, 'mk',{'-o'});
if pvOn
    sets(end+1) = struct('name',['PV only (' lbl ')'], 'V',Vpu_pv, 'P',P_pu_pv, 'PF',PF_pv, 'mk','-s');
end
if evOn
    if pvOn
        nm = sprintf('PV + EV (EV @ Bus %d: %.0f kW)', evBus, P_ev);
    else
        nm = sprintf('EV only (Bus %d: %.0f kW)', evBus, P_ev);
    end
    sets(end+1) = struct('name',nm, 'V',Vpu_all, 'P',P_pu_all, 'PF',PF_all, 'mk','-^');
end
nS = numel(sets);

figure('Name','Voltage Profile Comparison');
hold on;
for ii = 1:nS, plot(busTo, sets(ii).V, sets(ii).mk, 'LineWidth', 1.3); end
for m = 1:nPV, xline(pvBus(m),'--k',sprintf('PV %d',pvBus(m))); end
if evOn, xline(evBus,'--r',sprintf('EV %d',evBus)); end
yline(0.95,':r','0.95 p.u.'); grid on; hold off;
xlabel('Bus Number'); ylabel('Voltage (p.u.)');
legend({sets.name}, 'Location','best');
title(sprintf('Bus Voltage Profile Comparison - %s', modeNames{opMode}));

figure('Name','Voltage change vs baseline');
dVmat = cell2mat(arrayfun(@(ii) sets(ii).V - Vpu_pre, 2:nS, 'UniformOutput', false));
bar(busTo, dVmat); grid on;
xlabel('Bus Number'); ylabel('\DeltaV (p.u.)'); legend({sets(2:nS).name}, 'Location','best');
title('Voltage change relative to baseline');

figure('Name','Active Power and PF Comparison');
subplot(2,1,1); hold on;
for ii = 1:nS, plot(busTo, sets(ii).P, sets(ii).mk, 'LineWidth', 1.2); end
yline(0,':k'); grid on; hold off;
xlabel('Bus'); ylabel('P (p.u.)'); legend({sets.name}, 'Location','best');
title('Branch Active Power (negative = reverse flow)');
subplot(2,1,2); hold on;
for ii = 1:nS, plot(busTo, sets(ii).PF, sets(ii).mk, 'LineWidth', 1.2); end
yline(0,':k'); grid on; hold off;
xlabel('Bus'); ylabel('Power Factor'); legend({sets.name}, 'Location','best');
title('Branch Power Factor (PF < 0  <=>  P < 0, reverse flow)');

%% 7. THD / FFT for each PV bus (post PV + EV) 
for m = 1:nPV
    k = kPV(m);
    [thdF, fF, mF] = local_thd_fft(I3_all(idxSteady,1,k), Fs, f0);
    [thdP, fP, mP] = local_thd_fft(pvList{m}.Iabc(idxSteady,1), Fs, f0);
    fprintf('\nBus %d: feeder-current THD (phase A) = %.3f %% | PV inverter THD = %.2f %%\n', ...
            pvBus(m), thdF, thdP);
    figure('Name',sprintf('Current spectra, PV bus %d', pvBus(m)));
    subplot(2,1,1);
    stem_spec(fF, mF, 0.002, f0); grid on; xlim([0 25*f0]); xticks(0:5*f0:25*f0);
    xlabel('Frequency (Hz)'); ylabel('Magnitude (A)');
    title(sprintf('Feeder current into Bus %d (post), THD = %.2f%%', pvBus(m), thdF));
    subplot(2,1,2);
    stem_spec(fP, mP, 0.002, f0); grid on; xlim([0 25*f0]); xticks(0:5*f0:25*f0);
    xlabel('Frequency (Hz)'); ylabel('Magnitude (A)');
    title(sprintf('PV inverter current @ Bus %d, %s (THD = %.2f%%)', pvBus(m), pvList{m}.cfgUsed.switchType, thdP));
end

%%  8. EV station: harmonics, THD and IEEE 519 compliance 
if evOn
    H = 50;  hAll = 1:H;
    N1 = numel(idxSteady);
    % (a) current harmonics of the EV load itself 
    Aev   = harm_amps(ev.Iabc(idxSteady,1), Fs, f0, H);      % peak amplitudes
    IL    = Aev(1)/sqrt(2);                                   
    Ih_pc = 100 * Aev ./ max(Aev(1), eps);                    
    TDD   = 100 * sqrt(sum((Aev(2:end)/sqrt(2)).^2)) / IL;    

    % (b) short-circuit strength at the PCC from the feeder path impedance
    [~,~,Rb,Xb] = ieee33_data();
    pth  = branch_path_to_bus(busTo, evBus);
    Zpcc = sum(Rb(pth)) + 1j*sum(Xb(pth));                    % ohm (ideal source assumed)
    Isc  = (V_LL_base/sqrt(3)) / abs(Zpcc);               
    ratio = Isc/IL;
    [tddLim, limFun] = ieee519_current_limits(ratio);

    fprintf('\n=== IEEE 519 current check at Bus %d (PCC) ===\n', evBus);
    fprintf('  Isc = %.0f A, I_L = %.2f A  ->  Isc/I_L = %.0f\n', Isc, IL, ratio);
    fprintf('  TDD = %.2f %%  (limit %.1f %%)  ->  %s\n', TDD, tddLim, passfail(TDD <= tddLim));
    fprintf('  %-6s %-10s %-10s %s\n', 'Order','% of I_L','Limit %','');
    for h = 2:H
        if Ih_pc(h) > 0.05
            fprintf('  %-6d %-10.2f %-10.2f %s\n', h, Ih_pc(h), limFun(h), passfail(Ih_pc(h) <= limFun(h)));
        end
    end

    % (c) harmonic VOLTAGE distortion along the feeder (all harmonic sources)
    Ih_src = []; busSrc = []; Ih_pvOnly = []; busPV = [];
    for m = 1:nPV
        Am = harm_amps(pvList{m}.Iabc(idxSteady,1), Fs, f0, H) / sqrt(2);
        Ih_src = [Ih_src; Am]; busSrc = [busSrc pvBus(m)]; 
    end
    Ih_pvOnly = Ih_src;  busPV = busSrc;
    Ih_src = [Ih_src; Aev/sqrt(2)];  busSrc = [busSrc evBus];
    hList  = 2:H;
    V1_LN  = V_base * Vpu_all;
    if nPV > 0
        VTHD_pv = harmonic_voltage_profile(busPV, Ih_pvOnly(:,hList), hList, busTo, V_base*Vpu_pv);
    else
        VTHD_pv = zeros(numBranches,1);
    end
    VTHD_all = harmonic_voltage_profile(busSrc, Ih_src(:,hList),    hList, busTo, V1_LN);
    [vLimInd, vLimTHD] = ieee519_voltage_limits(V_LL_base);
    fprintf('\n=== IEEE 519 voltage check (%.2f kV class: THD <= %.1f %%, individual <= %.1f %%) ===\n', ...
            V_LL_base/1e3, vLimTHD, vLimInd);
    [vmx, ix] = max(VTHD_all);
    fprintf('  Max voltage THD = %.3f %% at Bus %d  ->  %s\n', vmx, busTo(ix), passfail(vmx <= vLimTHD));
    fprintf('  Voltage THD at EV bus %d = %.3f %%\n', evBus, VTHD_all(kEV));
    fprintf('  (worst-case: harmonic currents added in magnitude, ideal source, fundamental-frequency\n');
    fprintf('   network R + jhX, no capacitors/resonance, triplens treated like positive/negative sequence)\n');

    % plots 
    figure('Name','EV Charger Current Waveform');
    n3 = min(numel(t), round(3/(f0*Ts)));
    plot(t(1:n3), ev.Iabc(1:n3,:), 'LineWidth', 1.1); grid on;
    xlabel('Time (s)'); ylabel('Current (A)'); legend('Ia','Ib','Ic');
    title(sprintf('EV charger current (%s, THD = %.1f%%)', evType, ev.THD_I));

    figure('Name','EV Station Current Spectrum and IEEE 519 Limits');
    subplot(2,1,1);
    hp = 2:25;
    hS = stem_spec(hp, Ih_pc(hp), 0.002, []); hold on;
    hL = stairs(hp-0.5, arrayfun(limFun, hp), 'r', 'LineWidth', 1.5); hold off; grid on;
    xticks(hp); xlim([1 25.5]);
    xlabel('Harmonic order h'); ylabel('% of I_L');
    legend([hS hL], {'EV harmonic current','IEEE 519 individual limit'}, 'Location','northeast');
    title(sprintf('EV current spectrum @ Bus %d: TDD = %.2f%% (limit %.0f%%), Isc/I_L = %.0f', ...
                  evBus, TDD, tddLim, ratio));
    subplot(2,1,2);
    [thdFe, fFe, mFe] = local_thd_fft(I3_all(idxSteady,1,kEV), Fs, f0);
    stem_spec(fFe, mFe, 0.002, f0); grid on; xlim([0 25*f0]); xticks(0:5*f0:25*f0);
    xlabel('Frequency (Hz)'); ylabel('Magnitude (A)');
    title(sprintf('Feeder current into Bus %d (post PV + EV), THD = %.2f%%', evBus, thdFe));

    figure('Name','Harmonic Voltage Distortion Profile');
    hold on;  leg = {};
    if nPV > 0
        plot(busTo, VTHD_pv, '-s', 'LineWidth', 1.3);  leg{end+1} = 'PV harmonics only';
    end
    plot(busTo, VTHD_all, '-^', 'LineWidth', 1.3);
    if nPV > 0, leg{end+1} = 'PV + EV harmonics'; else, leg{end+1} = 'EV harmonics'; end
    yline(vLimTHD,'--r',sprintf('IEEE 519 THD limit %.0f%%', vLimTHD)); hold off; grid on;
    xlabel('Bus Number'); ylabel('Voltage THD (%)');
    legend(leg, 'Location','best');
    title('Estimated voltage THD along the feeder');
end

%% 9. R-L-C harmonic filter: optimal design, connection, verification 
if evOn && filterOn
    fprintf('\n=== Optimal R-L-C filter design at Bus %d ===\n', evBus);
    hDes    = 2:H;
    Vll_bus = V_LL_base * Vpu_pre(kEV);
    fin = struct('hOrders',hDes, 'Ih_rms',Aev(hDes)/sqrt(2), 'I1_rms',IL, ...
                 'P_kW',P_ev, 'Q_kvar',Q_ev, 'Vll',Vll_bus, 'f0',f0, ...
                 'Zsrc',[real(Zpcc) imag(Zpcc)], 'nBranch',nBr);
    flt = design_rlc_filter(fin);                   

    fprintf('  %-7s %-8s %-10s %-9s %-10s %-10s %-10s\n', 'Branch','Tuned h','Q (kVAr)','Q factor','R (ohm)','L (mH)','C (uF)');
    for b = 1:numel(flt.R)
        fprintf('  %-7d %-8.3f %-10.1f %-9.1f %-10.3f %-10.2f %-10.3f\n', b, flt.h_tuned(b), ...
                flt.Q_kvar(b), flt.Qfactor(b), flt.R_ohm(b), flt.L_mH(b), flt.C_uF(b));
    end
    fprintf('  (per-phase, wye-equivalent; targets: orders %s)\n', mat2str(flt.tuneOrders));
    fprintf('  Total filter %.1f kVAr (limit %.1f kVAr for PCC PF >= 0.95) | max capacitor voltage %.3f p.u. | constraints OK: %d\n', ...
            flt.Qtotal_kvar, flt.QtotMax_kvar, max(flt.Vc_pu), flt.constraintsOK);

    
    cfgF = cfgEV;  cfgF.filter = flt;  cfgF.Zsrc = fin.Zsrc;
    evF  = ev_charger_model(t, cfgF);
    [Vpu_f, V3_f, I3_f] = add_ev_to_network(evF, evBus, pvBus, Ppv_kW, Qpv_kvar, ...
        Vbfs_post, Vpu_pv, V3_pv, I3_pv, busTo, V_LL_base, S_base);
    [~, ~, ~, P_pu_f, ~, ~, PF_f] = ...
        process_run_from3(V3_f, I3_f, numBranches, idxSteady, V_base, S_base);

    % grid-side current harmonics 
    Ag      = harm_amps(evF.Igrid(idxSteady,1), Fs, f0, H);
    I1g     = Ag(1)/sqrt(2);
    Ih_pc_0 = Ih_pc;                                     % pre-filter
    Ih_pc_f = 100 * Ag ./ max(Aev(1), eps);              % post-filter
    TDD_0   = 100*sqrt(sum((Aev(2:end)/sqrt(2)).^2)) / IL;
    TDD_f   = 100*sqrt(sum((Ag(2:end)/sqrt(2)).^2)) / IL;
    ITHD_0  = 100*sqrt(sum((Aev(2:end)/sqrt(2)).^2)) / (Aev(1)/sqrt(2));
    ITHD_f  = 100*sqrt(sum((Ag(2:end)/sqrt(2)).^2)) / I1g;

    fprintf('\n=== IEEE 519 compliance at Bus %d: PRE-filter vs POST-filter (Isc/I_L = %.0f) ===\n', evBus, ratio);
    fprintf('  TDD  : %6.2f %% -> %6.2f %%   (limit %.1f %%)   %s -> %s\n', TDD_0, TDD_f, tddLim, ...
            passfail(TDD_0 <= tddLim), passfail(TDD_f <= tddLim));
    fprintf('  ITHD : %6.2f %% -> %6.2f %%   (relative to own fundamental)\n', ITHD_0, ITHD_f);
    fprintf('  %-6s %-10s %-10s %-9s %-6s %-6s\n', 'Order','Pre %I_L','Post %I_L','Limit %','Pre','Post');
    allOK_f = (TDD_f <= tddLim);
    for h = 2:H
        if max(Ih_pc_0(h), Ih_pc_f(h)) > 0.05
            fprintf('  %-6d %-10.2f %-10.2f %-9.2f %-6s %-6s\n', h, Ih_pc_0(h), Ih_pc_f(h), limFun(h), ...
                    passfail(Ih_pc_0(h) <= limFun(h)), passfail(Ih_pc_f(h) <= limFun(h)));
        end
        if Ih_pc_f(h) > limFun(h), allOK_f = false; end
    end
    if allOK_f
        fprintf('  ==> POST-filter: fully IEEE 519 compliant at the PCC.\n');
    else
        fprintf('  ==> POST-filter: improved, but some limits are still exceeded (see FAIL above).\n');
        fprintf('      Try one more branch (target the next dominant harmonic) or a larger reactive budget.\n');
    end

    % harmonic voltage distortion along the feeder with the filter
    Ih_src_f = [Ih_pvOnly; Ag/sqrt(2)];
    VTHD_f   = harmonic_voltage_profile(busSrc, Ih_src_f(:,hList), hList, busTo, V_base*Vpu_f);
    fprintf('\n  Voltage THD at EV bus %d : %.3f %% -> %.3f %%   (max on feeder %.3f %% -> %.3f %%, limit %.1f %%)\n', ...
            evBus, VTHD_all(kEV), VTHD_f(kEV), max(VTHD_all), max(VTHD_f), vLimTHD);
    if evF.Qgrid < 0, pfTxt = 'leading'; else, pfTxt = 'lagging'; end
    fprintf('  Grid-side PF at the PCC : %.3f -> %.3f (%s) | net Q %.1f -> %.1f kVAr | filter loss %.0f W\n', ...
            ev.PFgrid, evF.PFgrid, pfTxt, ev.Qgrid/1e3, evF.Qgrid/1e3, evF.filterLoss_W);
    fprintf('  Voltage at Bus %d : %.4f -> %.4f p.u. | minimum on feeder %.4f -> %.4f p.u.\n', ...
            evBus, Vpu_all(kEV), Vpu_f(kEV), min(Vpu_all), min(Vpu_f));

    % plots 
    figure('Name','Filter: Current Spectrum and THD, Pre vs Post');
    subplot(2,1,1);
    hp = 2:25;
    h1 = stem_spec(hp-0.12, Ih_pc_0(hp), 0.002, []); hold on;      % pre-filter
    h2 = stem_spec(hp+0.12, Ih_pc_f(hp), 0.002, []);               % post-filter
    h3 = stairs(hp-0.5, arrayfun(limFun, hp), 'r', 'LineWidth', 1.5); hold off; grid on;
    xticks(hp); xlim([1 25.5]);
    xlabel('Harmonic order h'); ylabel('% of I_L');
    legend([h1 h2 h3], {'Pre-filter','Post-filter','IEEE 519 limit'}, 'Location','northeast');
    title(sprintf('Grid-side current @ Bus %d: TDD %.1f%% -> %.1f%% (limit %.0f%%)', evBus, TDD_0, TDD_f, tddLim));
    subplot(2,1,2);
    bar([TDD_0 TDD_f; ITHD_0 ITHD_f]); grid on;
    set(gca,'XTickLabel',{'TDD','ITHD'}); ylabel('%');
    yline(tddLim,'--r',sprintf('TDD limit %.0f%%', tddLim));
    legend('Pre-filter','Post-filter','Location','northeast');
    title('Current distortion, pre vs post filter');

    figure('Name','Filter: Phase-A Current Waveform, Pre vs Post');
    n3 = min(numel(t), round(3/(f0*Ts)));
    plot(t(1:n3), evF.Iabc(1:n3,1), 'LineWidth', 1.2); hold on;
    plot(t(1:n3), evF.Igrid(1:n3,1), 'LineWidth', 1.2); hold off; grid on;
    xlabel('Time (s)'); ylabel('Current (A)');
    legend('Charger load current (pre-filter)','Grid current (post-filter)','Location','best');
    title(sprintf('EV station current, filter connected (ITHD %.1f%% -> %.1f%%)', ITHD_0, ITHD_f));

    figure('Name','Filter: Frequency Response');
    hh  = 0.5:0.05:30;
    Zfh = arrayfun(@(h) rlc_bank_Z(flt, h, f0), hh);
    Zsh = real(Zpcc) + 1j*hh*imag(Zpcc);
    Gh  = Zfh ./ (Zfh + Zsh);
    subplot(2,1,1);
    semilogy(hh, abs(Zfh), 'LineWidth', 1.3); hold on; semilogy(hh, abs(Zsh), 'LineWidth', 1.3); hold off; grid on;
    for tt_ = flt.tuneOrders, xline(tt_,':k'); end
    xlabel('Harmonic order h'); ylabel('|Z| (ohm)'); legend('Filter |Z_f|','Grid |Z_s|','Location','best');
    title('Filter impedance vs grid impedance');
    subplot(2,1,2);
    hG = plot(hh, abs(Gh), 'LineWidth', 1.3); hold on; yline(1,':r','no attenuation');
    ch = hAll(Aev > 0.001*Aev(1) & hAll > 1);
    hM = plot(ch, abs(interp1(hh, Gh, ch)), 'ko', 'MarkerFaceColor','k'); hold off; grid on;
    xlabel('Harmonic order h'); ylabel('|I_{grid}/I_{load}|'); ylim([0 1.5]);
    legend([hG hM], {'Grid current fraction','Charger harmonics'}, 'Location','best');
    title('Fraction of each harmonic that still reaches the grid (values > 1 = resonance)');

    figure('Name','Filter: Voltage THD Profile, Pre vs Post');
    plot(busTo, VTHD_all, '-^', 'LineWidth', 1.3); hold on;
    plot(busTo, VTHD_f, '-o', 'LineWidth', 1.3);
    yline(vLimTHD,'--r',sprintf('IEEE 519 THD limit %.0f%%', vLimTHD)); hold off; grid on;
    xlabel('Bus Number'); ylabel('Voltage THD (%)');
    legend('Without filter','With filter','Location','best');
    title('Estimated voltage THD along the feeder');

    figure('Name','Filter: Bus Voltage Profile');
    plot(busTo, Vpu_pre, '-o', 'LineWidth', 1.3); hold on;
    plot(busTo, Vpu_all, '-s', 'LineWidth', 1.3);
    plot(busTo, Vpu_f,   '-^', 'LineWidth', 1.3);
    xline(evBus,'--r',sprintf('EV %d',evBus)); yline(0.95,':r','0.95 p.u.'); hold off; grid on;
    xlabel('Bus Number'); ylabel('Voltage (p.u.)');
    legend('Baseline', [allName ' (no filter)'], [allName ' + filter'], 'Location','best');
    title('Bus voltage profile: effect of the filter''s reactive support');
end


%% LOCAL FUNCTIONS 
function [Vpu_x, V3_x, I3_x] = add_ev_to_network(evs, evBus, pvBus, Ppv_kW, Qpv_kvar, ...
        Vbfs_pv, Vpu_pv, V3_pv, I3_pv, busTo, V_LL_base, S_base)

    Vbfs_x = ieee33_bfs([pvBus evBus], [Ppv_kW, -evs.Pgrid/1e3], [Qpv_kvar, -evs.Qgrid/1e3], ...
                        V_LL_base, S_base);
    dV     = Vbfs_x(busTo) - Vbfs_pv(busTo);            
    Vpu_x  = Vpu_pv + dV(:);
    I3_x   = I3_pv;                                    
    for kk = branch_path_to_bus(busTo, evBus)
        I3_x(:,:,kk) = I3_x(:,:,kk) + evs.Igrid;
    end
    V3_x = V3_pv;
    for k = 1:size(V3_pv,3)
        V3_x(:,:,k) = V3_pv(:,:,k) * (Vpu_x(k)/Vpu_pv(k));
    end
end

function Z = rlc_bank_Z(flt, h, f0)
% Impedance of the parallel bank of series R-L-C branches at harmonic order h.
    w  = 2*pi*f0;
    Zb = flt.R(:).' + 1j*(h*w*flt.L(:).' - 1./(h*w*flt.C(:).'));
    Z  = 1 / sum(1 ./ Zb);
end

function hs = stem_spec(x, y, thr, f0)

    y = y(:).';  x = x(:).';
    y(y < thr*max(y)) = NaN;
    hs = stem(x, y, 'filled', 'LineWidth', 1.4, 'MarkerSize', 6);
    if ~isempty(f0)
        [~, ord] = sort(y, 'descend');            
        for i = ord(1:min(8, sum(~isnan(y))))
            text(x(i), y(i), sprintf(' h%d', round(x(i)/f0)), 'VerticalAlignment','bottom', 'FontSize', 8);
        end
    end
end

function A = harm_amps(sig, Fs, f0, H)
% Peak amplitude of harmonics 1..H from one fundamental cycle of data.
    L = numel(sig);  Y = fft(sig - mean(sig)) / L;  A = zeros(1,H);
    for h = 1:H
        idx = round(h*f0*L/Fs) + 1;
        if idx <= floor(L/2) + 1, A(h) = 2*abs(Y(idx)); end
    end
end

function [tddLim, limFun] = ieee519_current_limits(ratio)

    tbl = [ 4.0 2.0 1.5 0.6 0.3  5.0;      
            7.0 3.5 2.5 1.0 0.5  8.0;      
           10.0 4.5 4.0 1.5 0.7 12.0;      
           12.0 5.5 5.0 2.0 1.0 15.0;      
           15.0 7.0 6.0 2.5 1.4 20.0];     
    if     ratio < 20,   row = 1;
    elseif ratio < 50,   row = 2;
    elseif ratio < 100,  row = 3;
    elseif ratio < 1000, row = 4;
    else,                row = 5; end
    tddLim = tbl(row, 6);
    limFun = @(h) tbl(row, 1 + (h>=11) + (h>=17) + (h>=23) + (h>=35)) * (1 - 0.75*(mod(h,2)==0));
end

function [indiv, thd] = ieee519_voltage_limits(Vll)
% IEEE 519 Table 1: voltage distortion limits 
    if     Vll <= 1e3,   indiv = 5.0; thd = 8.0;
    elseif Vll <= 69e3,  indiv = 3.0; thd = 5.0;
    elseif Vll <= 161e3, indiv = 1.5; thd = 2.5;
    else,                indiv = 1.0; thd = 1.5; end
end

function [VTHD, Vh] = harmonic_voltage_profile(srcBus, Ih_rms, hList, busTo, V1_LN)

    [~,~,R,X] = ieee33_data();
    nB = numel(busTo);  Vh = zeros(nB, numel(hList));
    pathB = cell(nB,1);
    for b = 1:nB, pathB{b} = branch_path_to_bus(busTo, busTo(b)); end
    for s = 1:numel(srcBus)
        ps = branch_path_to_bus(busTo, srcBus(s));
        for b = 1:nB
            common = intersect(pathB{b}, ps);
            for ih = 1:numel(hList)
                Z = sum(R(common)) + 1j*hList(ih)*sum(X(common));
                Vh(b,ih) = Vh(b,ih) + Ih_rms(s,ih) * abs(Z);
            end
        end
    end
    VTHD = 100 * sqrt(sum(Vh.^2, 2)) ./ V1_LN(:);
end

function s = passfail(ok)
    if ok, s = 'PASS'; else, s = 'FAIL'; end
end

function b = ask_bus_list(promptText, defaultVal)

    while true
        r = input(sprintf('%s [%s]: ', promptText, num2str(defaultVal)));
        if isempty(r), b = defaultVal; return; end
        if isnumeric(r) && isvector(r) && all(r == round(r)) && all(r >= 2 & r <= 33)
            if numel(unique(r)) < numel(r)
                fprintf('     Repeated bus numbers removed.\n');
            end
            b = unique(r, 'stable'); return;
        end
        fprintf('     Enter integer bus numbers between 2 and 33, e.g. 18 or [18 33].\n');
    end
end

function x = ask_number(promptText, defaultVal, lo, hi)

    while true
        r = input(sprintf('%s [%g]: ', promptText, defaultVal));
        if isempty(r), x = defaultVal; return; end
        if isnumeric(r) && isscalar(r) && isfinite(r) && r >= lo && r <= hi
            x = r; return;
        end
        fprintf('     Please enter a number between %g and %g.\n', lo, hi);
    end
end

function [frm, too, R, X, PL, QL] = ieee33_data()
% Standard IEEE 33-bus (Baran & Wu) branch impedances
    frm = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17  2 19 20 21  3 23 24  6 26 27 28 29 30 31 32];
    too = [2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33];
    R = [0.0922 0.4930 0.3660 0.3811 0.8190 0.1872 0.7114 1.0300 1.0440 0.1966 0.3744 1.4680 ...
         0.5416 0.5910 0.7463 1.2890 0.7320 0.1640 1.5042 0.4095 0.7089 0.4512 0.8980 0.8960 ...
         0.2030 0.2842 1.0590 0.8042 0.5075 0.9744 0.3105 0.3410];
    X = [0.0470 0.2511 0.1864 0.1941 0.7070 0.6188 0.2351 0.7400 0.7400 0.0650 0.1238 1.1550 ...
         0.7129 0.5260 0.5450 1.7210 0.5740 0.1565 1.3554 0.4784 0.9373 0.3083 0.7091 0.7011 ...
         0.1034 0.1447 0.9337 0.7006 0.2585 0.9630 0.3619 0.5302];
    PL = [100 90 120 60 60 200 200 60 60 45 60 60 120 60 60 60 90 90 90 90 90 90 420 420 60 60 60 120 200 150 210 60];
    QL = [ 60 40  80 30 20 100 100 20 20 30 35 35  80 10 20 20 40 40 40 40 40 50 200 200 25 25 20  70 600  70 100 40];
    
end

function V = ieee33_bfs(pvBuses, Ppv_kW, Qpv_kvar, Vll, Sb)

    [frm, too, R, X, PL, QL] = ieee33_data();
    Zb = Vll^2 / Sb;
    Z  = (R + 1j*X) / Zb;
    S  = zeros(33,1);
    S(2:33) = (PL + 1j*QL).' * 1e3 / Sb;   
    for m = 1:numel(pvBuses)
        S(pvBuses(m)) = S(pvBuses(m)) - (Ppv_kW(m) + 1j*Qpv_kvar(m)) * 1e3 / Sb;
    end
    V = ones(33,1) * (1 + 0j);
    for it = 1:200
        I  = conj(S ./ V);
        Ib = zeros(33,1);
        for j = numel(too):-1:1                      
            Ib(too(j)) = Ib(too(j)) + I(too(j));
            Ib(frm(j)) = Ib(frm(j)) + Ib(too(j));
        end
        Vn = V;
        for j = 1:numel(too)                          
            Vn(too(j)) = Vn(frm(j)) - Z(j) * Ib(too(j));
        end
        if max(abs(Vn - V)) < 1e-10, V = Vn; break; end
        V = Vn;
    end
    V = abs(V);
end

function path = branch_path_to_bus(busTo, bus)

    parent = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 2 19 20 21 3 23 24 6 26 27 28 29 30 31 32];
    path = [];  b = bus;
    while b > 1
        path = [find(busTo == b, 1), path]; 
        b = parent(b-1);
    end
end

function [Vpu, V3, I3, P_pu, Q_pu, S_pu, PF] = process_run(Vts, Its, nB, idxSteady, V_base, S_base)

    n = numel(Vts.Time);
    V3r = reshape(Vts.Data, n, 3, nB);  I3 = reshape(Its.Data, n, 3, nB);
    V3 = zeros(size(V3r));
    for k = 1:nB
        V3(:,1,k) = (V3r(:,1,k) - V3r(:,3,k))/3;
        V3(:,2,k) = (V3r(:,2,k) - V3r(:,1,k))/3;
        V3(:,3,k) = (V3r(:,3,k) - V3r(:,2,k))/3;
    end
    [Vpu, ~, ~, P_pu, Q_pu, S_pu, PF] = process_run_from3(V3, I3, nB, idxSteady, V_base, S_base);
end

function [Vpu, V3, I3, P_pu, Q_pu, S_pu, PF] = process_run_from3(V3, I3, nB, idxSteady, V_base, S_base)
    P_pu = zeros(nB,1); Q_pu = P_pu; S_pu = P_pu; PF = P_pu; Vm = P_pu;
    for k = 1:nB
        Vk = V3(idxSteady,:,k);  Ik = I3(idxSteady,:,k);
        Vr = sqrt(mean(Vk.^2,'all'));  Ir = sqrt(mean(Ik.^2,'all'));
        Pt = sum(mean(Vk.*Ik,1));  St = 3*Vr*Ir;
        P_pu(k) = Pt/S_base;  S_pu(k) = St/S_base;
        PF(k)   = Pt/max(St,eps);
        Q_pu(k) = S_pu(k)*sin(acos(min(max(PF(k),-1),1)));
        Vm(k)   = Vr;
    end
    Vpu = Vm / V_base;
end

function [thdPercent, freqAxis, magSpectrum] = local_thd_fft(signal, Fs, f0)
    L = numel(signal);  win = signal - mean(signal);
    P2 = abs(fft(win)/L);  P1 = P2(1:floor(L/2)+1);  P1(2:end-1) = 2*P1(2:end-1);
    freqAxis = Fs*(0:floor(L/2))/L;  magSpectrum = P1;
    [~, fi] = min(abs(freqAxis - f0));  fundMag = P1(fi);
    hs = 0;
    for h = 2:25
        [~, hi] = min(abs(freqAxis - h*f0));
        if hi <= numel(P1), hs = hs + P1(hi)^2; end
    end
    thdPercent = 100*sqrt(hs)/fundMag;
end