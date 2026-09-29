function flt = design_rlc_filter(in)

d.nBranch = 2;        d.tuneOrders = [];   d.pfMin = 0.95;   d.QtotMax_kvar = [];
d.QfacRange = [20 100]; d.detune = [0.90 0.995];
d.wTDD = 1;  d.wTHD = 0.5;  d.VcMax_pu = 1.10;  d.nStarts = 15;  d.seed = 1;
fn = fieldnames(d);
for i = 1:numel(fn)
    if ~isfield(in, fn{i}) || isempty(in.(fn{i})), in.(fn{i}) = d.(fn{i}); end
end
req = {'hOrders','Ih_rms','I1_rms','P_kW','Q_kvar','Vll','f0','Zsrc'};
for i = 1:numel(req)
    assert(isfield(in, req{i}), 'design_rlc_filter: missing input field "%s".', req{i});
end

hOrd = in.hOrders(:).';   Ih = in.Ih_rms(:).';
w    = 2*pi*in.f0;        Vph = in.Vll/sqrt(3);
Rs   = in.Zsrc(1);        Xs  = in.Zsrc(2);
Zs   = @(h) Rs + 1j*h*Xs;
I1   = in.I1_rms;
phi1 = -atan2(in.Q_kvar, in.P_kW);
Iev1 = I1*exp(1j*phi1);                       % EV fundamental current phasor (V phase = 0)
idxP = find(Ih > 1e-3*I1);                    % harmonics actually present in the EV current
hP   = hOrd(idxP);   IhP = Ih(idxP);

%  choose the tuned orders
if isempty(in.tuneOrders)
    [~, ord] = sort(Ih, 'descend');
    tuneOrders = sort(hOrd(ord(1:min(in.nBranch, numel(ord)))));
else
    tuneOrders = sort(in.tuneOrders(:).');
end
nB = numel(tuneOrders);

% total kVAr allowed: PCC power factor must stay >= pfMin (lagging or leading)
QtotMax = in.QtotMax_kvar;
if isempty(QtotMax)
    QtotMax = in.Q_kvar + tan(acos(in.pfMin))*in.P_kW;
end
QtotMax = max(QtotMax, 1);                    % kVAr

% optimization
rng(in.seed);
opt = optimset('MaxIter',3000,'MaxFunEvals',6000,'TolX',1e-8,'TolFun',1e-10,'Display','off');
bestJ = inf;  bestZ = zeros(1,3*nB);
for s = 1:in.nStarts
    if s == 1, z0 = zeros(1,3*nB); else, z0 = 1.5*randn(1,3*nB); end
    [z, Jz] = fminsearch(@cost, z0, opt);
    if Jz < bestJ, bestJ = Jz; bestZ = z; end
end

[Qk, ht, Qf] = decode(bestZ);
[R, L, C, Xc] = components(Qk, ht, Qf);
m = evaluate(R, L, C, Xc);

% output
flt.tuneOrders  = tuneOrders;
flt.h_tuned     = ht;
flt.Q_kvar      = Qk;
flt.Qtotal_kvar = sum(Qk);
flt.Qfactor     = Qf;
flt.R = R;  flt.L = L;  flt.C = C;
flt.R_ohm = R;  flt.L_mH = 1e3*L;  flt.C_uF = 1e6*C;
flt.TDD_pre   = m.TDDpre;
flt.TDD_post  = m.TDD;
flt.ITHD_post = m.THD;
flt.maxAmp    = m.maxAmp;
flt.Vc_pu     = m.Vc_pu;
flt.Ploss_W   = m.Ploss;
flt.cost      = bestJ;
flt.QtotMax_kvar = QtotMax;
flt.constraintsOK = (sum(Qk) <= QtotMax*1.01) && (m.maxAmp <= 0.99) && (max(m.Vc_pu) <= in.VcMax_pu*1.01);


    function [Qk, ht, Qf] = decode(z)
        sg = 1 ./ (1 + exp(-max(min(z(:), 30), -30)));
        sg = reshape(sg, 3, nB);
        Qk = 0.5 + (QtotMax - 0.5) * sg(1,:);                       % kVAr per branch
        ht = tuneOrders .* (in.detune(1) + diff(in.detune)*sg(2,:));% tuned order
        Qf = in.QfacRange(1) + diff(in.QfacRange)*sg(3,:);          % quality factor
    end

    function [R, L, C, Xc] = components(Qk, ht, Qf)
        % branch net fundamental var  Q = Vll^2 / (Xc - Xl),  Xl = Xc/ht^2
        Xc = in.Vll^2 ./ (Qk*1e3 .* (1 - 1./ht.^2));
        Xl = Xc ./ ht.^2;
        Xn = Xc ./ ht;                                              % sqrt(L/C) in ohm
        R  = Xn ./ Qf;
        C  = 1 ./ (w*Xc);
        L  = Xl / w;
    end

    function m = evaluate(R, L, C, Xc)
        Zbr = @(h) R + 1j*(h*w*L - 1./(h*w*C));                     % branch impedances
        nH  = numel(hP);
        Zfh = zeros(1,nH);  G = zeros(1,nH);
        for kk = 1:nH
            zb = Zbr(hP(kk));
            Zfh(kk) = 1 / sum(1 ./ zb);                             % parallel bank impedance
            G(kk)   = Zfh(kk) / (Zfh(kk) + Zs(hP(kk)));             % grid-current fraction
        end
        Ig   = IhP .* abs(G);                                       % grid-side harmonic current
        z1   = Zbr(1);
        If1  = Vph / (1 / sum(1 ./ z1));                            % fundamental filter current
        Ig1  = Iev1 + If1;                                          % grid fundamental current
        m.TDDpre = 100*sqrt(sum(Ih.^2))/I1;
        m.TDD    = 100*sqrt(sum(Ig.^2))/I1;                         % relative to I_L
        m.THD    = 100*sqrt(sum(Ig.^2))/abs(Ig1);                   % relative to grid fundamental
        m.maxAmp = max([abs(G), 0]);                                % worst harmonic amplification
        Vc = zeros(1, nB);                                          % capacitor rms voltage per branch
        for b = 1:nB
            acc = abs(Vph / z1(b) * Xc(b))^2;                       % fundamental part
            for kk = 1:nH
                zb  = Zbr(hP(kk));
                Ibh = IhP(kk) * (1 - G(kk)) * (Zfh(kk) / zb(b));    % this branch's harmonic current
                acc = acc + abs(Ibh * Xc(b) / hP(kk))^2;
            end
            Vc(b) = sqrt(acc) / Vph;
        end
        m.Vc_pu = Vc;
        m.Ploss = 3 * real(Vph * conj(If1));                        % fundamental filter loss (W)
    end

    function J = cost(z)
        [Qk_, ht_, Qf_] = decode(z);
        [R_, L_, C_, Xc_] = components(Qk_, ht_, Qf_);
        m_ = evaluate(R_, L_, C_, Xc_);
        pQ = max(0, sum(Qk_) - QtotMax)^2 / QtotMax^2;              % PF / kVAr limit
        pA = max(0, m_.maxAmp - 0.98)^2;                            % no harmonic amplification
        pV = sum(max(0, m_.Vc_pu - in.VcMax_pu).^2);                % capacitor stress
        J  = in.wTDD*m_.TDD + in.wTHD*m_.THD + 1e4*pQ + 1e3*pA + 1e3*pV;
        if ~isfinite(J), J = 1e9; end
    end
end
