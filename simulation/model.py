"""
Mechanistic phage-host population dynamics model for SusI simulation studies.

Purpose
-------
Generate synthetic turbidity (OD) curves whose underlying biology is known
exactly, so that the behaviour of SusI, ti/tc, VI and SupI can be traced back
to defined physiological parameters (adsorption rate constant, burst size,
latent period, phage loss rate, frequency of insensitive cells).

State variables
---------------
    A         primary (fast) growth resource
    B         secondary (slow) growth resource
    S         susceptible cells
    T         reversibly phage-tolerant cells (non-heritable)
    I_1..I_n  infected cells, n sequential stages (Erlang latent period)
    R         heritably phage-insensitive cells
    V         free phage
    D         non-viable light-scattering debris

Equations
---------
    qA      = A / (KmA + A)              primary-resource activity, 0..1
    qB      = B / (KmB + B)              secondary-resource activity, 0..1
    mu      = muA * qA + muB * qB        specific growth rate
    k_eff   = (n / tau) * max(qA, q_floor)       lysis rate
    beta_eff= beta * qA                          burst size

    dA/dt   = -(muA * qA * Nv) / YA
    dB/dt   = -(muB * qB * Nv) / YB
    dS/dt   = mu * S - delta * S * V - k_tol * S + k_rev * T
    dT/dt   = mu_T_rel * mu * T - eta * delta * T * V + k_tol * S - k_rev * T
    dI_1/dt = delta * (S + eta * T) * V - k_eff * I_1
    dI_j/dt = k_eff * (I_{j-1} - I_j)            j = 2..n
    dR/dt   = mu_R_rel * mu * R
    dV/dt   = beta_eff * k_eff * I_n - delta * (S + eta*T + I_tot) * V - omega * V
    dD/dt   = kappa_release * k_eff * I_n - omega_D * D

    Nv = S + T + R                       resource-consuming cells
    N  = S + T + R + I_tot               viable biomass
    OD = (N + kappa_debris * D) / cells_per_OD + OD_blank

Why two resources
-----------------
A single resource cannot reproduce two features of the experimental control
curves at once. A late tc (8-12 h, as reported) requires a slow approach to
plateau, but a single Monod resource achieves that only by lowering the
exponential growth rate, which then makes regrowth of an insensitive
subpopulation unrealistically slow: at mu = 0.45 /h, 160 founder cells per well
need more than 24 h to become detectable, so ti saturates at the algorithm cap
for every resistance frequency below about 1e-4 and SusI appears insensitive to
resistance frequency. That would be a modelling artefact, and it would give a
wrong answer to the reviewer question about resistance frequency.

Splitting growth into a fast primary phase and a slow saturating secondary
phase gives a realistic exponential rate (muA ~ 0.9 /h) together with the long
post-exponential creep that places tc near 10 h.

Why phage replication is tied to host physiological state
---------------------------------------------------------
In a closed batch culture with a growth-independent burst size, a replicating
lytic phage always wins completely: the free phage titre overshoots the level
at which killing merely balances division by orders of magnitude, so the
culture is always driven to the detection floor and the suppression depth
saturates. Under such a model the only thing that could differ between phages
is how long suppression lasts, which would make SusI redundant with ti/tc by
construction. That conclusion would be an artefact of the model.

Phage replication requires an actively growing host, so burst size and lysis
rate fall as the primary resource is exhausted. This creates a race between
phage amplification and culture saturation:

* high adsorption rate and large burst clear the culture while the host is
  still growing, giving deep (sharp) lysis;
* low adsorption rate, small burst or high phage loss means the phage is
  overtaken by resource exhaustion, a substantial fraction of cells is never
  lysed, and turbidity stays high (blunt lysis).

The controlling quantity is the phage's effective reproductive number in the
culture,

    R_phage = beta * (delta * N0) / (delta * N0 + omega)

that is, burst size times the probability that a free phage adsorbs to a host
before being lost. Suppression depth saturates above R_phage of roughly 5.

Why the tolerant state T and the debris term D are present
----------------------------------------------------------
* T is a reversible, non-heritable tolerant state (entered at k_tol, left at
  k_rev, adsorbing phage at relative rate eta). This is the adaptive,
  non-mutational form of phage insensitivity raised in review, and is distinct
  from R, which is heritable. Disabled by default (k_tol = 0).
* D is debris released on lysis, contributing light scattering at relative
  weight kappa_debris. It makes the objection that OD cannot distinguish
  viable cells from debris quantifiable rather than only arguable. Disabled by
  default (kappa_debris = 0).

Cryptic growth
--------------
`recycle` returns a fraction of lysed biomass to the primary resource pool.
Lysed cells release their contents, and surviving or resistant cells can grow
on them; this is standard cryptic growth. It matters here because when lysis
occurs at high culture density most of the resource has already been consumed,
so without recycling an insensitive subpopulation has almost nothing left to
grow on and regrowth after ti is barely visible. With recycling, regrowth after
a high-density lysis event is pronounced, as it is in the experimental curves.

Discrete founders
-----------------
The insensitive subpopulation is seeded as a continuous concentration f * N0.
In a 200 uL well holding 1.6e7 cells, f = 1e-7 is about 1.6 founders per well
and f = 1e-9 is 0.016, so most wells would contain none and regrowth would
require de novo mutation. `Params.founders_per_well()` reports this number;
below roughly f = 1e-6 the discrete stochastic regime is the realistic one and
well-to-well variability is expected. `mutation_rate` enables de novo S -> R.

All rates are per hour; concentrations are per mL.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
import numpy as np
from scipy.integrate import solve_ivp


# --------------------------------------------------------------------------
# Parameters
# --------------------------------------------------------------------------

@dataclass
class Params:
    # --- host growth (calibrated: tc ~ 9.5 h, plateau ~ 1.14 OD) ---
    muA: float = 0.90          # 1/h   primary-phase maximum specific growth rate
    KmA: float = 0.02          # -     primary-resource half-saturation
    YA: float = 6.0e8          # CFU/mL produced from the primary resource
    muB: float = 0.07          # 1/h   secondary-phase maximum specific growth rate
    KmB: float = 0.35          # -     secondary-resource half-saturation
    YB: float = 3.0e8          # CFU/mL produced from the secondary resource
    A0: float = 1.0            # -     initial primary resource
    B0: float = 1.0            # -     initial secondary resource
    q_floor: float = 0.05      # -     minimum relative lysis rate at exhaustion
    recycle: float = 0.0       # -     fraction of lysed biomass returned to the
                               #       primary resource pool (cryptic growth)

    mu_R_rel: float = 1.0      # -     growth of R relative to S (< 1 = fitness cost)
    mu_T_rel: float = 1.0      # -     growth of T relative to S
    N0: float = 8.0e7          # CFU/mL inoculum (OD ~ 0.1 at cells_per_OD = 8e8)
    f: float = 1.0e-6          # -     frequency of pre-existing insensitive cells
    mutation_rate: float = 0.0  # per replication, de novo S -> R

    # --- reversible, non-heritable tolerance ---
    k_tol: float = 0.0         # 1/h   S -> T entry rate (0 disables tolerance)
    k_rev: float = 0.5         # 1/h   T -> S return rate
    eta: float = 0.02          # -     adsorption to T relative to S

    # --- phage ---
    delta: float = 1.0e-7      # mL/(PFU*h) adsorption rate constant
    beta: float = 100.0        # -     burst size at full host activity
    tau: float = 0.5           # h     mean latent period
    n_stages: int = 10         # -     Erlang shape of the latent period
    omega: float = 0.05        # 1/h   lumped free-phage loss rate
    moi: float = 1.0           # -     multiplicity of infection (V0 = moi * N0)

    # --- optics / debris ---
    cells_per_OD: float = 8.0e8
    OD_blank: float = 0.02
    kappa_release: float = 1.0
    kappa_debris: float = 0.0
    omega_D: float = 0.0

    # --- measurement ---
    t_end: float = 48.0        # h
    dt_obs: float = 0.25       # h  (15 min, as in the study)
    noise_sd: float = 0.004    # OD
    seed: int | None = 0
    well_volume: float = 0.2   # mL

    def label(self) -> str:
        return (f"delta={self.delta:.1e} beta={self.beta:g} tau={self.tau:g} "
                f"omega={self.omega:g} f={self.f:.1e} moi={self.moi:g}")

    def r_phage(self) -> float:
        """Effective reproductive number of the phage at inoculation."""
        aN = self.delta * self.N0
        return self.beta * aN / (aN + self.omega)

    def founders_per_well(self) -> float:
        return self.f * self.N0 * self.well_volume

    def tolerant_fraction(self) -> float:
        tot = self.k_tol + self.k_rev
        return self.k_tol / tot if tot > 0 else 0.0


# --------------------------------------------------------------------------
# Right-hand side
# --------------------------------------------------------------------------
# state layout: [A, B, S, T, I_1..I_n, R, V, D]

def _rhs(t, y, p: Params):
    n = p.n_stages
    A = max(y[0], 0.0)
    B = max(y[1], 0.0)
    S = max(y[2], 0.0)
    T = max(y[3], 0.0)
    I = np.clip(y[4:4 + n], 0.0, None)
    R = max(y[4 + n], 0.0)
    V = max(y[5 + n], 0.0)
    D = max(y[6 + n], 0.0)

    I_tot = I.sum()
    Nv = S + T + R

    qA = A / (p.KmA + A) if A > 0 else 0.0
    qB = B / (p.KmB + B) if B > 0 else 0.0
    mu = p.muA * qA + p.muB * qB

    k = (n / p.tau) * max(qA, p.q_floor)
    beta_eff = p.beta * qA

    adsorb_S = p.delta * S * V
    adsorb_T = p.eta * p.delta * T * V
    growth_S = mu * S
    mutate = p.mutation_rate * growth_S

    lysis_flux = k * I[-1]
    dA = -(p.muA * qA * Nv) / p.YA + p.recycle * lysis_flux / p.YA
    dB = -(p.muB * qB * Nv) / p.YB
    dS = growth_S - adsorb_S - mutate - p.k_tol * S + p.k_rev * T
    dT = p.mu_T_rel * mu * T - adsorb_T + p.k_tol * S - p.k_rev * T

    dI = np.empty(n)
    dI[0] = adsorb_S + adsorb_T - k * I[0]
    if n > 1:
        dI[1:] = k * (I[:-1] - I[1:])

    dR = p.mu_R_rel * mu * R + mutate
    dV = beta_eff * lysis_flux - p.delta * (S + p.eta * T + I_tot) * V - p.omega * V
    dD = p.kappa_release * lysis_flux - p.omega_D * D

    return np.concatenate(([dA], [dB], [dS], [dT], dI, [dR], [dV], [dD]))


def simulate(p: Params, with_phage: bool = True):
    """Integrate the model. Returns (t_obs, od_obs, od_true, states)."""
    n = p.n_stages
    R0 = p.f * p.N0
    phi = p.tolerant_fraction()
    T0 = (p.N0 - R0) * phi
    S0 = p.N0 - R0 - T0
    V0 = p.moi * p.N0 if with_phage else 0.0
    y0 = np.concatenate(([p.A0], [p.B0], [S0], [T0], np.zeros(n),
                         [R0], [V0], [0.0]))

    t_obs = np.arange(0.0, p.t_end + 1e-9, p.dt_obs)
    sol = solve_ivp(_rhs, (0.0, p.t_end), y0, args=(p,), t_eval=t_obs,
                    method="LSODA", rtol=1e-6, atol=1e-2, max_step=0.5)
    if not sol.success:
        raise RuntimeError(f"integration failed: {sol.message}")

    Y = sol.y
    S = np.clip(Y[2], 0, None)
    T = np.clip(Y[3], 0, None)
    I_tot = np.clip(Y[4:4 + n], 0, None).sum(axis=0)
    R = np.clip(Y[4 + n], 0, None)
    D = np.clip(Y[6 + n], 0, None)

    biomass = S + T + R + I_tot + p.kappa_debris * D
    od_true = biomass / p.cells_per_OD + p.OD_blank

    rng = np.random.default_rng(p.seed)
    od_obs = od_true + rng.normal(0.0, p.noise_sd, size=od_true.shape)

    states = dict(t=t_obs, A=np.clip(Y[0], 0, None), B=np.clip(Y[1], 0, None),
                  S=S, T=T, I=I_tot, R=R, V=np.clip(Y[5 + n], 0, None), D=D)
    return t_obs, od_obs, od_true, states


def simulate_pair(p: Params):
    """Return (t, od_phage, od_control), both with measurement noise."""
    t, od_p, _, _ = simulate(p, with_phage=True)
    pc = Params(**asdict(p))
    pc.seed = None if p.seed is None else p.seed + 10_000
    _, od_c, _, _ = simulate(pc, with_phage=False)
    return t, od_p, od_c
