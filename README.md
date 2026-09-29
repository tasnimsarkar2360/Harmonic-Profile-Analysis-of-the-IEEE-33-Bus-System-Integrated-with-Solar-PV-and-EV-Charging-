# Harmonic & Voltage Profile Analysis of the IEEE 33-Bus System with Solar PV and EV Charging Stations

**EEE 306 Final Project, Group 4**
Department of EEE, Bangladesh University of Engineering and Technology (BUET)

## Overview
Distribution feeders were built for one-way power flow from the substation to passive loads.
Rooftop solar PV and DC fast-charging EV stations change that: they reverse power flow, shift
bus voltages and inject harmonic currents. This project models the standard IEEE 33-bus radial
feeder (12.66 kV, 100 MVA base, 50 Hz) in MATLAB/Simulink. It measures how PV and EV integration
affect the voltage profile and harmonic distortion (THD/TDD). It then designs an optimal passive
R-L-C shunt filter that brings the EV bus back within IEEE 519-2014 limits.

## Key Results
| Metric | Before | After |
|---|---|---|
| EV-only current TDD at the PCC (6-pulse, 300 kW) | 27.14 % ❌ (limit 15 %) | **11.74 % ✅** with the R-L-C filter |
| Weakest bus (Bus 18) voltage | 0.9119 p.u. | **0.9313 p.u.** with 200 kW PV |
| Maximum voltage THD on the feeder | — | ≤ 0.59 % (limit 5 %) |

## Processing Pipeline
1. **Baseline**: a Simulink time-domain run, followed by a steady-state P/Q/S/PF, voltage and current extraction, weakest-bus identification and FFT/THD analysis.
2. **PV model**: a single-diode PV array solved at MPP, feeding a grid-tie inverter (PWM, diode or thyristor harmonic profile).
3. **EV model**: an EV charger as a non-linear load (6-pulse diode, 12-pulse or active front-end).
4. **Filter design**: a multi-branch tuned R-L-C shunt filter, optimized for minimum grid-side TDD/THD.
5. **Integration & compliance**: PV-only, EV-only or combined PV+EV cases, checked against the IEEE 519 current and voltage limits before and after filtering.

## Repository Structure
| File | Purpose |
|---|---|
| `BaseCode.m` | Runs the Simulink model and converts line-to-line voltages to line-to-neutral. Computes per-branch P, Q, S, PF, voltage and current in p.u., finds the 3 weakest buses, and plots waveforms and FFT spectra. Runs diagnostics: settling, DC offset, phase balance and a KCL check against the Baran & Wu loads. Exports the two Excel files. |
| `pv_cell_model.m` | Modular PV array and inverter model. The single-diode DC model is solved by Newton-Raphson and swept to MPP. It auto-sizes parallel strings to a target kW, generates a 3-phase injected current with characteristic harmonics, and can also run in interactive mode. |
| `ev_charger_model.m` | EV charging station as a harmonic-producing load, with 6-pulse / 12-pulse / active-front-end spectra. It can optionally connect an R-L-C filter in parallel and return grid-side current, THD, PF and filter loss. Includes a standalone demo mode. |
| `design_rlc_filter.m` | Optimal single- or multi-branch tuned R-L-C filter. For each branch it optimizes the kVAr, detuned order and quality factor, using base-MATLAB `fminsearch` with multi-start (no toolbox needed). Constraints: PCC power factor ≥ 0.95, no harmonic amplification (resonance), and a capacitor voltage limit. |
| `pv_ev_integration_analysis.m` | The main study script, run interactively (mode 1 = PV only, 2 = EV only, 3 = PV + EV). It runs a backward/forward-sweep load flow anchored to the Simulink baseline, handles multiple PV buses, and checks IEEE 519 TDD and individual-harmonic limits using the Isc/I_L ratio. It also estimates the harmonic voltage THD profile along the feeder, reports reverse power flow, and compares pre- and post-filter results with plots. |
| `IEEE33_Bus_Parameters_PU.xlsx` | Baseline output: bus voltage, P, Q, S (p.u.) and power factor. |
| `IEEE33_Line_Currents_PU.xlsx` | Baseline output: branch currents in A and p.u. |
| `EEE306_Group4.pptx` | Final project presentation. |

## Requirements
- MATLAB (R2020b or later recommended) with Simulink and Simscape Electrical (Specialized Power Systems)
- The Simulink model `IEEE33BusTestSystem_Instantaneous.slx`. It needs "To Workspace" blocks that log `simOutputV_inst` and `simOutputI_inst` before the RMS blocks.
- No Optimization Toolbox is needed; the filter optimizer uses only base MATLAB.

## How to Run
```matlab
% 1. Baseline simulation, exports the Excel files and plots
BaseCode

% 2. PV / EV integration study (uses the workspace from step 1)
pv_ev_integration_analysis
%    -> choose mode 1/2/3, PV bus(es), irradiance, temperature, kW,
%       EV bus, P/Q, charger type, and whether to design a filter

% Standalone demos
pv_cell_model          % interactive PV model
ev_charger_model       % default 300 kW 6-pulse charger
```

## Default Study Case
- PV: Bus 18 (the weakest bus), 200 kW, G = 1000 W/m², Tc = 25 °C, PWM inverter
- EV: Bus 20, 300 kW + 100 kVAr, 6-pulse diode charger
- Filter: 2 branches tuned just below the 5th and 7th harmonics

## Assumptions & Limitations
- Balanced three-phase operation; fundamental-frequency constant-PQ loads
- Harmonic spectra are typical textbook values for each converter type (IEEE 519 defines limits, not spectra); replace them with datasheet or measured spectra if available
- Harmonic voltage estimates use a worst-case magnitude sum over network R + jhX, with an ideal source and no shunt capacitors
- Transient switching, protection coordination and unbalanced single-phase loading are out of scope

## Team (Group 4)
Safwan Zaman Chowdhury (2206070) · Md. Montashim Billaha Chisty (2206081) · Tasnim Sarkar (2206084) ·
Naimul Hasan Ayon (2206085) · Abhishek Pandit (2206090) · Md. Shadman Sifat (2206092) ·
Md. Redwoan Islam (2206097)

## References
- IEEE Std 519-2014, *Recommended Practice and Requirements for Harmonic Control in Electric Power Systems*
- M. E. Baran and F. F. Wu, "Network reconfiguration in distribution systems for loss reduction and load balancing," *IEEE Trans. Power Delivery*, 1989
