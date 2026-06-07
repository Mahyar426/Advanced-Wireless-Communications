<div align="center">

![header](https://readme-typing-svg.demolab.com?font=Fira+Code&size=26&pause=1000&color=00D9FF&center=true&vCenter=true&width=800&lines=📡+Advanced+Wireless+Communications;DS-CDMA+·+MIMO+Capacity+·+SINR+Beamforming.;From+equations+to+interactive+GUIs+in+MATLAB.)

```
██╗    ██╗██╗██████╗ ███████╗██╗     ███████╗███████╗███████╗
██║    ██║██║██╔══██╗██╔════╝██║     ██╔════╝██╔════╝██╔════╝
██║ █╗ ██║██║██████╔╝█████╗  ██║     █████╗  ███████╗███████╗
██║███╗██║██║██╔══██╗██╔══╝  ██║     ██╔══╝  ╚════██║╚════██║
╚███╔███╔╝██║██║  ██║███████╗███████╗███████╗███████║███████║
 ╚══╝╚══╝ ╚═╝╚═╝  ╚═╝╚══════╝╚══════╝╚══════╝╚══════╝╚══════╝
   multi-user · mimo · beamforming · from scratch
```

[![MATLAB](https://img.shields.io/badge/MATLAB-R2023b%2B-orange?style=flat-square&logo=mathworks)](https://www.mathworks.com/)
[![App Designer](https://img.shields.io/badge/UI-App_Designer-00D9FF?style=flat-square)](https://www.mathworks.com/products/matlab/app-designer.html)
[![Domain](https://img.shields.io/badge/Domain-Physical_Layer_%7C_Wireless-blueviolet?style=flat-square)](#)
[![Status](https://img.shields.io/badge/Status-Complete-brightgreen?style=flat-square)](#)

</div>

---

Three simulation suites built from scratch in MATLAB, covering core physical-layer algorithms in modern wireless systems — complete with interactive App Designer GUIs and full BER/capacity/radiation-pattern visualization.

| Module | Topic | Key Algorithms |
|--------|-------|----------------|
| [BER Analysis](#1-multi-user-ber-analysis--ds-cdma) | Multi-user DS-CDMA detection | SUMF · Decorrelator · MMSE · ML-Optimal |
| [MIMO Capacity](#2-mimo-channel-capacity--water-filling) | MIMO channel capacity | Water-filling · Ergodic/Outage capacity · Monte Carlo |
| [Precoding](#3-sinr-constrained-beamforming--precoding) | SINR-constrained beamforming | Fixed-point iteration · Newton's method · Radiation patterns |

---

## 1. Multi-User BER Analysis — DS-CDMA

**The core challenge:** in a spread-spectrum system with multiple simultaneous users, how well can different receivers separate their signals? This simulation benchmarks four receiver architectures head-to-head across SNR.

### What it does
- Generates random spreading sequences and simulates an N-user DS-CDMA channel with correlated noise
- Correlated noise is generated via spectral decomposition of **R** — not just AWGN
- Computes the **BER vs. SNR curve** for all four receiver types simultaneously:
  - **SUMF** — Single-User Matched Filter (treats interference as noise): `sign(y)`
  - **Decorrelator** — inverts the cross-correlation matrix: `sign(R⁻¹ · y)`
  - **MMSE** — balances noise vs. interference: `F = A⁻¹ · (R + σ²A⁻²)⁻¹`
  - **ML-Optimal** — exhaustive Euclidean distance search over all 2ᴷ combinations using `R^(1/2)` whitening
- Derives **Near-Far resistance** from `diag(R⁻¹)` — quantifies robustness to amplitude disparities
- Full simulation exposed via interactive **MATLAB App Designer GUI** — tune users, SNR range, sequence length, and seeds live

---

## 2. MIMO Channel Capacity — Water-Filling

**The core challenge:** how much information can a multi-antenna channel carry, and how much does knowing the channel at the transmitter actually help?

### What it does
- Simulates Rayleigh-fading MIMO (Nₜ × Nᵣ) over thousands of Monte Carlo realizations: `H ~ CN(0, 1/2)` per element
- Computes two capacity benchmarks per realization:
  - **Without CSIT** — uniform power: `C = log₂ det(I + SNR/Nₜ · HH†)`
  - **With CSIT** — optimal water-filling over SVD eigenvalues
- Custom **bisection water-filling solver**: converges to `|Σγᵢ − Pₜ| < 10⁻⁶`
- Extracts **ergodic capacity** (mean) and **outage capacity** (configurable percentile via empirical CDF interpolation)
- Interactive GUI to sweep antenna counts, SNR, and percentile thresholds

---

## 3. SINR-Constrained Beamforming — Precoding

**The core challenge:** design transmit beamforming vectors for multiple users simultaneously, guaranteeing each user meets a minimum SINR target while minimizing total transmit power.

### What it does
- Models a ULA transmitter with configurable inter-element spacing and per-user AoA
- Channel model: `H(i,k) = (λ/4πd) · exp(−j·2π·spacing·k·cos(αᵢ))` — near-field ULA steering vectors
- Solves SINR-constrained precoding via two-stage iterative optimization:
  1. **Fixed-point iteration** — Lagrange multipliers via self-consistent equations (relative norm convergence)
  2. **Newton's method** — refinement using analytically derived Jacobian `∂f/∂λ` (faster than numerical differentiation)
- Computes optimal **beamforming matrix W** and per-user power allocation via linear system solve
- Sweeps 0° → 180° over 2000 angular points to render **polar radiation patterns** per user beam
- Reports total transmit power in dBW and per-stage iteration counts

---

## Repository Structure

```
├── A1/
│   ├── RX_BER_Calculator.m       # Core BER engine (SUMF/DECO/MMSE/OPT)
│   ├── plot_BER_AppDesigner.m    # BER plotting helper
│   └── assignment1.mlapp         # Interactive GUI
│
├── A2/
│   ├── mainFunctions.m           # Capacity + water-filling engine
│   └── assignment2.mlapp         # Interactive GUI
│
└── A3/
    ├── calculate_precoding.m     # Full precoding solver
    └── app1.mlapp                # Interactive GUI
```

---

## Skills Demonstrated

- Physical-layer algorithm implementation from mathematical formulations — no toolbox shortcuts
- Numerical optimization: bisection solvers, fixed-point iteration, Newton's method with analytic Jacobians
- Linear algebra: SVD, matrix factorizations, pseudoinverses
- Monte Carlo simulation design and statistical analysis (ergodic & outage metrics)
- MATLAB App Designer for interactive, parameter-driven simulation GUIs

---

*M.Sc. Communications Engineering — Politecnico di Torino*
