clc; clear;


%% parameters 
S_base    = 100e6;                             % system base power (VA) - 100 MVA
V_LL_base = 12660;                             % nominal line-to-line base voltage (V)
V_base    = V_LL_base / sqrt(3);               % line-to-neutral base voltage (V)
I_base    = S_base / (sqrt(3) * V_LL_base);    % line current base (A)
f0        = 50;                                % fundamental source frequency (Hz)

simName   = 'IEEE33BusTestSystem_Instantaneous';   

%% IEEE 33-bus order
busTo = [2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 ...
         23 24 25 26 27 28 29 30 31 32 33];
numBranches = numel(busTo);

%% circuit breaker input vector
noBranches = 32;
noNoP = 5; % open circuit breaker
timeTS = (0:0.001:0.2)';
CBStats = [ones(length(timeTS),noBranches), zeros(length(timeTS),noNoP)];
siminCB = timeseries(CBStats,timeTS);

%% Run the simulation
load_system(simName);
simIn = Simulink.SimulationInput(simName);
simIn = simIn.setVariable('siminCB',siminCB);
simout = sim(simIn);

%% 1. Extract instantaneous voltage and current waveforms
Vts = simout.simOutputV_inst;
Its = simout.simOutputI_inst;

t    = Vts.Time;
Vraw = Vts.Data;                   % nSamples x 96 (32 branches x [A B C])
Iraw = Its.Data;

nSamples = numel(t);
Ts = mean(diff(t));
Fs = 1/Ts;
N_cycle = round((1/f0)/Ts);        % samples per fundamental cycle(50Hz=f0)(1/f0=20ms)

V3_raw = reshape(Vraw, nSamples, 3, numBranches);   % [time x (Vab,Vbc,Vca) x branch]
I3_raw = reshape(Iraw, nSamples, 3, numBranches);   % [time x (Ia,Ib,Ic) x branch]


V3 = zeros(size(V3_raw));
for k = 1:numBranches
    V_ab = V3_raw(:, 1, k);
    V_bc = V3_raw(:, 2, k);
    V_ca = V3_raw(:, 3, k);

    V3(:, 1, k) = (V_ab - V_ca) / 3;   % Van
    V3(:, 2, k) = (V_bc - V_ab) / 3;   % Vbn
    V3(:, 3, k) = (V_ca - V_bc) / 3;   % Vcn
end

I3 = I3_raw;   

idxSteady = (nSamples - N_cycle + 1):nSamples;  % final fundamental cycle only
idxPlot   = max(1, nSamples - 4*N_cycle + 1):nSamples;   % last 4 cycles: steady state only
% (the first 2-3 cycles contain the start-up transient)


%% 2. 3-phase steady-state P, Q, S, PF, V, I per bus

P_pu = zeros(numBranches,1);
Q_pu = zeros(numBranches,1);
S_pu = zeros(numBranches,1);
PF   = zeros(numBranches,1);
Vmag = zeros(numBranches,1);
Imag = zeros(numBranches,1);

for k = 1:numBranches
    Vk = V3(idxSteady, :, k);      % N_cycle x 3  [Va Vb Vc]
    Ik = I3(idxSteady, :, k);      % N_cycle x 3  [Ia Ib Ic]

    Vrms_LN     = sqrt(mean(Vk.^2, 'all'));    %  line-to-neutral RMS
    Irms_branch = sqrt(mean(Ik.^2, 'all'));    %  branch RMS current

    P_phase_avg = mean(Vk .* Ik, 1);
    P_total = sum(P_phase_avg);

    S_total = 3 * Vrms_LN * Irms_branch;

    P_pu(k) = P_total / S_base;
    S_pu(k) = S_total / S_base;

    PF(k)   = P_total / max(S_total, eps);                          
    Q_pu(k) = S_pu(k) * sin(acos(min(max(PF(k), -1), 1)));          

    Vmag(k) = Vrms_LN;
    Imag(k) = Irms_branch;
end

Vpu = Vmag / V_base;
Ipu = Imag / I_base;


if median(Vpu) < 0.5 || median(Vpu) > 1.5
    warning('Median bus voltage is %.3f p.u. - check the LN reconstruction / phase order.', median(Vpu));
end


%% 3. Excel File 1 - Bus Parameters (p.u.)

T1 = table(busTo', Vpu, P_pu, Q_pu, S_pu, PF, ...
    'VariableNames', {'Bus','Voltage_pu','Active_P_pu','Reactive_Q_pu','Apparent_S_pu','Power_Factor'});
writetable(T1, 'IEEE33_Bus_Parameters_PU.xlsx', 'Sheet', 'Bus_Parameters');

%% 4. Excel File 2 - Line/Branch Current (A and p.u.)
branchNames = arrayfun(@(k) sprintf('Branch_to_Bus%d', busTo(k)), 1:numBranches, 'UniformOutput', false)';
T2 = table(branchNames, Imag, Ipu, ...
    'VariableNames', {'Branch','Current_A','Current_pu'});
writetable(T2, 'IEEE33_Line_Currents_PU.xlsx', 'Sheet', 'Line_Currents');


%% 5. Weakest-bus analysis (lowest voltage magnitude, p.u.) 
[VpuSorted, sortIdx] = sort(Vpu, 'ascend');
weakest3Bus = busTo(sortIdx(1:3));
weakest3Vpu = VpuSorted(1:3);

fprintf('\n--- 3 Weakest Buses (lowest voltage, p.u.) ---\n');
for k = 1:3
    fprintf('  Bus %2d : %.4f p.u.\n', weakest3Bus(k), weakest3Vpu(k));
end

kMin = sortIdx(1);                 % branch index of the single weakest bus
weakestBusNum = busTo(kMin);

%% 6. Weakest-bus waveforms
figure('Name','Weakest Bus Voltage Waveforms');
plot(t(idxPlot), squeeze(V3(idxPlot,1,kMin)), 'LineWidth', 1.1); hold on;
plot(t(idxPlot), squeeze(V3(idxPlot,2,kMin)), 'LineWidth', 1.1);
plot(t(idxPlot), squeeze(V3(idxPlot,3,kMin)), 'LineWidth', 1.1); hold off; grid on;
xlabel('Time (s)'); ylabel('Voltage (V)');
legend('Va','Vb','Vc');
title(sprintf('3-Phase Voltage - Weakest Bus %d (last 4 cycles)', weakestBusNum));

figure('Name','Weakest Bus Current Waveforms');
plot(t(idxPlot), squeeze(I3(idxPlot,1,kMin)), 'LineWidth', 1.1); hold on;
plot(t(idxPlot), squeeze(I3(idxPlot,2,kMin)), 'LineWidth', 1.1);
plot(t(idxPlot), squeeze(I3(idxPlot,3,kMin)), 'LineWidth', 1.1); hold off; grid on;
xlabel('Time (s)'); ylabel('Current (A)');
legend('Ia','Ib','Ic');
title(sprintf('3-Phase Current - Weakest Bus %d Branch (last 4 cycles)', weakestBusNum));


%% 7. THD and FFT spectrum for Phase-A V and I at the weakest bus

Va_steady = V3(idxSteady, 1, kMin);
Ia_steady = I3(idxSteady, 1, kMin);

[thdV, freqV, magV] = local_thd_fft(Va_steady, Fs, f0);
[thdI, freqI, magI] = local_thd_fft(Ia_steady, Fs, f0);

fprintf('\n--- THD at Weakest Bus %d (Phase A) ---\n', weakestBusNum);
fprintf('  Voltage THD : %.3f %%\n', thdV);
fprintf('  Current THD : %.3f %%\n\n', thdI);

figure('Name','Weakest Bus Voltage FFT Spectrum');
stem_spec(freqV, magV, 0.002, f0); grid on; xlim([0, 25*f0]); xticks(0:5*f0:25*f0);   % 25th harmonic
xlabel('Frequency (Hz)'); ylabel('Magnitude (V)');
title(sprintf('Phase-A Voltage Spectrum - Weakest Bus %d', weakestBusNum));

figure('Name','Weakest Bus Current FFT Spectrum');
stem_spec(freqI, magI, 0.002, f0); grid on; xlim([0, 25*f0]); xticks(0:5*f0:25*f0);   % 25th harmonic
xlabel('Frequency (Hz)'); ylabel('Magnitude (A)');
title(sprintf('Phase-A Current Spectrum - Weakest Bus %d', weakestBusNum));

%%    KCL check of every branch against the expected load current

w0 = 2*pi*f0;
idxPrev = idxSteady - N_cycle;                 % the cycle just before the final one
assert(idxPrev(1) >= 1, 'Simulation too short for the settling check (need >= 2 cycles).');
phasorOf = @(x, idx) (2/numel(idx)) * sum(reshape(x(idx),[],1) .* exp(-1j*w0*reshape(t(idx),[],1)));    % peak phasor (cos ref.)

Vph = zeros(numBranches,3);   Iph = zeros(numBranches,3);                 % rms phasors
rmsLast = zeros(numBranches,1); rmsPrev = zeros(numBranches,1); dcPct = zeros(numBranches,1);
for k = 1:numBranches
    for ph = 1:3
        Vph(k,ph) = phasorOf(V3(:,ph,k), idxSteady) / sqrt(2);
        Iph(k,ph) = phasorOf(I3(:,ph,k), idxSteady) / sqrt(2);
    end
    rmsLast(k) = sqrt(mean(I3(idxSteady,:,k).^2, 'all'));
    rmsPrev(k) = sqrt(mean(I3(idxPrev,:,k).^2,  'all'));
    dcPct(k)   = 100 * max(abs(mean(I3(idxSteady,:,k), 1))) / max(rmsLast(k), eps);
end
settlePct = 100 * abs(rmsLast - rmsPrev) ./ max(rmsLast, eps);           

Imags    = abs(Iph);   Imean = mean(Imags, 2);
unbalPct = 100 * max(abs(Imags - Imean), [], 2) ./ max(Imean, eps);       % magnitude unbalance
wrap180  = @(a) mod(a + 180, 360) - 180;
angB     = rad2deg(angle(Iph(:,2) .* conj(Iph(:,1))));                    % expect -120 deg
angC     = rad2deg(angle(Iph(:,3) .* conj(Iph(:,1))));                    % expect +120 deg
angDev   = max(abs(wrap180(angB + 120)), abs(wrap180(angC - 120)));       % deviation from 120 deg


[frmT, tooT, PLt, QLt] = local_ieee33_topology();
kclErr = nan(numBranches,1);
if isequal(busTo(:), tooT(:))
    for kb = 1:numBranches
        b        = busTo(kb);
        kids     = find(frmT == b);                                       % branches leaving bus b
        Iload_m  = Iph(kb,1) - sum(Iph(kids,1));                          % measured load current (phase A)
        Sph      = (PLt(kb) + 1j*QLt(kb)) * 1e3 / 3;                      % per-phase load (VA)
        Iload_e  = conj(Sph / Vph(kb,1));                                 
        kclErr(kb) = 100 * abs(Iload_m - Iload_e) / abs(Iload_e);
    end
else
    warning('busTo differs from the standard IEEE 33 ordering: KCL check skipped.');
end

fprintf('\n=== Current diagnostics (final cycle) ===\n');
kFirst = 1;  k18 = find(busTo == 18, 1);  k33 = find(busTo == 33, 1);
fprintf('  Expected from a constant-PQ load flow : first branch ~210 A (%.3f p.u.), Bus-18 branch ~4.9 A, Bus-33 branch ~3.6 A\n', 210.4/I_base);
fprintf('  Simulated                             : first branch %.1f A (%.3f p.u.), Bus-18 branch %.1f A, Bus-33 branch %.1f A\n', ...
        Imag(kFirst), Ipu(kFirst), Imag(k18), Imag(k33));
fprintf('  (currents FALL along the feeder by design: each branch carries all downstream load)\n');

flag = (settlePct > 1) | (dcPct > 1) | (unbalPct > 2) | (angDev > 2) | (kclErr > 5);
if any(flag)
    fprintf('  Branches outside tolerance (settle>1%%, DC>1%%, unbal>2%%, angle>2deg, KCL>5%%):\n');
    fprintf('  %-16s %8s %7s %8s %9s %8s\n', 'Branch','Settle%','DC%','Unbal%','AngDev','KCL%');
    for k = find(flag).'
        fprintf('  %-16s %8.2f %7.2f %8.2f %9.2f %8.2f\n', sprintf('to Bus %d', busTo(k)), ...
                settlePct(k), dcPct(k), unbalPct(k), angDev(k), kclErr(k));
    end
    if any(settlePct > 1)
        fprintf('  -> "Settle" flags mean the start-up transient has not died out: increase the stop time\n');
        fprintf('     (or start from steady state in powergui) before trusting currents/THD.\n');
    end
    if any(kclErr > 5)
        fprintf('  -> KCL flags point to a mis-mapped or mis-oriented branch measurement at that bus.\n');
    end
else
    fprintf('  All branches settled, balanced (angles 120 deg) and KCL-consistent: no current anomalies.\n');
end

figure('Name','Branch Current Diagnostics');
subplot(3,1,1);
plot(busTo, Imag, '-o', 'LineWidth', 1.2); grid on;
xlabel('Bus (branch feeds this bus)'); ylabel('RMS current (A)');
title('Branch RMS current (falls along each feeder: each branch carries the downstream load)');
subplot(3,1,2);
stem(busTo, kclErr, 'filled'); hold on; yline(5, '--r', '5 %'); hold off; grid on;
xlabel('Bus'); ylabel('KCL residual (%)');
title('Measured bus load current vs expected (should be within a few %)');
subplot(3,1,3);
plot(busTo, unbalPct, '-o', busTo, angDev, '-s', 'LineWidth', 1.2); grid on;
xlabel('Bus'); ylabel('% / degrees'); legend('Magnitude unbalance (%)','Angle deviation from 120 (deg)','Location','best');
title('Three-phase balance of each branch current');


%% Local function: single-sided FFT magnitude spectrum + THD (up to 25th harmonic)
function [thdPercent, freqAxis, magSpectrum] = local_thd_fft(signal, Fs, f0)
    L = numel(signal);
    win = signal - mean(signal);           % remove DC offset
    Yf = fft(win);
    P2 = abs(Yf/L);
    P1 = P2(1:floor(L/2)+1);
    P1(2:end-1) = 2*P1(2:end-1);
    freqAxis = Fs*(0:floor(L/2))/L;
    magSpectrum = P1;

    [~, fundIdx] = min(abs(freqAxis - f0));
    fundMag = P1(fundIdx);

    maxHarm = 25;
    harmSumSq = 0;
    for h = 2:maxHarm
        [~, hIdx] = min(abs(freqAxis - h*f0));
        if hIdx <= numel(P1)
            harmSumSq = harmSumSq + P1(hIdx)^2;
        end
    end
    thdPercent = 100 * sqrt(harmSumSq) / fundMag;
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

function [frm, too, PL, QL] = local_ieee33_topology()
% Standard IEEE 33-bus (Baran & Wu) topology and loads (kW, kVAr) for buses 2..33.
    frm = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17  2 19 20 21  3 23 24  6 26 27 28 29 30 31 32];
    too = [2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33];
    PL  = [100 90 120 60 60 200 200 60 60 45 60 60 120 60 60 60 90 90 90 90 90 90 420 420 60 60 60 120 200 150 210 60];
    QL  = [ 60 40  80 30 20 100 100 20 20 30 35 35  80 10 20 20 40 40 40 40 40 50 200 200 25 25 20  70 600  70 100 40];
end