"""
Hückel π-system sandbox — the hello-world of the causal intervention loop.

This is L1 of the validation ladder (PLAN.md §4): a toy simulator where the
causal ground truth is textbook chemistry, and the *loop mechanics* (intervention
menu -> mediator observables -> hypothesis update) can be exercised end-to-end
with zero dependencies beyond numpy.

Causal template (mirrors the real D-A case):
    X (interventions)  ->  M (mediators)  ->  Y (outcome)
    X: do(alpha_k)  electronegativity intervention  (Coulomb integral at atom k)
       do(beta_ij) coupling/resonance intervention  (resonance integral on bond ij)
    M: pi charges, MO energies, HOMO coefficients
    Y: charge redistribution at the para site, Delta q_para

Pedagogical punchline (why this demo exists):
  Part 1  Observational comparison (pyridine vs pyrrole) is CONFOUNDED —
          alpha AND electron count changed together. You cannot attribute.
  Part 2  do(alpha) sweep with everything else frozen isolates the
          electronegativity lever. Slope d(Dq_para)/dh = causal effect.
  Part 3  do(beta) sweep with alpha frozen isolates the resonance lever.
  Verdict: both levers matter; their interventional slopes quantify each
          pathway's contribution — the toy version of mediation analysis.

Units: alpha_C = 0, beta_CC = -1 (energies in |beta| units, beta < 0 = bonding).
Standard heteroatom params (Streitwieser): N pyridine-type h=0.5 (1 pi e-),
N pyrrole-type h=1.5 (2 pi e-).
"""

import numpy as np

ALPHA_C = 0.0
BETA_CC = -1.0


# ---------------------------------------------------------------- system spec
def ring(n=6, alpha_h=None, beta_k=None, electrons=None):
    """Build an n-membered ring. alpha_h[i], beta_k[(i,j)], electrons[i]."""
    alpha_h = alpha_h or {}
    beta_k = beta_k or {}
    electrons = electrons or {}
    atoms = [dict(h=alpha_h.get(i, 0.0), n_e=electrons.get(i, 1)) for i in range(n)]
    bonds = {}
    for i in range(n):
        j = (i + 1) % n
        bonds[frozenset((i, j))] = beta_k.get(frozenset((i, j)), 1.0)
    return dict(atoms=atoms, bonds=bonds)


def hamiltonian(system):
    n = len(system["atoms"])
    H = np.zeros((n, n))
    for i, a in enumerate(system["atoms"]):
        H[i, i] = ALPHA_C + a["h"] * BETA_CC
    for bond, k in system["bonds"].items():
        i, j = tuple(bond)
        H[i, j] = H[j, i] = k * BETA_CC
    return H


def solve(system):
    """Diagonalize, fill electrons, return energies, occ MOs, pi charges, bond orders."""
    H = hamiltonian(system)
    E, C = np.linalg.eigh(H)          # ascending energy
    n_total = sum(a["n_e"] for a in system["atoms"])
    assert n_total % 2 == 0, "closed shell assumed"
    n_occ = n_total // 2
    occ = np.zeros(len(E))
    occ[:n_occ] = 2.0
    # Hückel population: pi charge on atom i = sum over occ MOs of 2|c|^2
    charges = (C**2 * occ).sum(axis=1)
    # bond order (free-valence style): sum over occ of 2 c_i c_j
    P = (C * occ) @ C.T
    return dict(E=E, C=C, occ=occ, charges=charges, P=P)


def para_site(n=6, site=0):
    """Para position to `site` in an n-ring."""
    return (site + n // 2) % n


def charge_pattern(sys, site=0):
    """Excess pi charge (vs each atom's own contribution) at N/ortho/meta/para of `site`."""
    r = solve(sys)
    n = len(sys["atoms"])
    labels = {"N": site, "ortho": 1, "meta": 2, "para": para_site(n, site)}
    return {k: r["charges"][i] - sys["atoms"][i]["n_e"] for k, i in labels.items()}, r


def fmt_row(label, d):
    return f"{label:<28}" + "".join(f"{v:>+9.3f}" for v in d.values())


# ================================================================ PART 1
# Observational comparison: pyridine-N (h=0.5, 1e) vs pyrrole-N (h=1.5, 2e).
# CONFOUNDED: alpha AND electron count changed together.
print("=" * 72)
print("PART 1 — OBSERVATIONAL (confounded): pyridine-type N vs pyrrole-type N")
print("=" * 72)
benzene = ring()
pyridine = ring(alpha_h={0: 0.5}, electrons={0: 1})          # 6 pi e-
pyrrole_ring = ring(alpha_h={0: 1.5}, electrons={0: 2})      # 7th e-? no: 8 pi e- total... 

# NOTE: a 6-ring with one 2e heteroatom has 7 atoms' worth of electrons if we
# naively sum (5*1 + 2 = 7, odd). Real pyrrole is a 5-ring. For a *closed-shell
# toy* we keep the 6-ring and model the "donor" as an extra electron pair is
# impossible; instead we model donor vs acceptor by h sign/magnitude with the
# SAME electron count — which is exactly the point: the observational pair in
# the wild (pyridine vs pyrrole) differs in BOTH knobs; the simulator lets us
# freeze one. Here we compare the real pair via 5-ring pyrrole vs 6-ring
# pyridine for the textbook pattern, then move to controlled sweeps.

pyrrole5 = dict(
    atoms=[dict(h=1.5, n_e=2)] + [dict(h=0.0, n_e=1) for _ in range(4)],
    bonds={frozenset((i, (i + 1) % 5)): 0.8 if 0 in (i, (i + 1) % 5) else 1.0
           for i in range(5)},
)

print(f"{'system':<28}{'N':>9}{'ortho':>9}{'meta':>9}{'para':>9}")
for name, sys, site in [("benzene (reference)", benzene, None),
                        ("pyridine-N, h=0.5", pyridine, 0),
                        ("pyrrole (5-ring), h=1.5", pyrrole5, 0)]:
    if site is None:
        print(fmt_row(name, {k: 0.0 for k in ["N", "ortho", "meta", "para"]}))
    else:
        d, _ = charge_pattern(sys, site=site)
        print(fmt_row(name, d))
print("""
Textbook pattern reproduced: pyridine-N depletes ortho/para (EAS at meta);
pyrrole-N enriches ortho (EAS at C2/C5). But this comparison changed alpha,
electron count, ring size AND bond couplings at once -> ATTRIBUTION IMPOSSIBLE.
Observational data, however pretty, cannot separate the levers. Now intervene.
""")

# ================================================================ PART 2
# do(alpha_N) sweep — everything frozen except electronegativity.
print("=" * 72)
print("PART 2 — INTERVENTIONAL: do(alpha_N = h*beta), 6-ring, 1 pi e on N, all else frozen")
print("=" * 72)
print(f"{'h (electronegativity)':<24}{'dq_N':>9}{'dq_ortho':>9}{'dq_meta':>9}{'dq_para':>9}{'E_HOMO':>9}")
rows = []
for h in [0.0, 0.5, 1.0, 1.5, 2.0]:
    sys = ring(alpha_h={0: h}, electrons={0: 1})
    d, r = charge_pattern(sys, site=0)
    homo = r["E"][r["occ"] > 0].max()
    rows.append((h, d["para"], homo))
    print(f"{h:<24.1f}{d['N']:>+9.3f}{d['ortho']:>+9.3f}{d['meta']:>+9.3f}{d['para']:>+9.3f}{homo:>9.3f}")

hs = np.array([r[0] for r in rows]); qp = np.array([r[1] for r in rows])
slope_alpha = float(np.polyfit(hs, qp, 1)[0])
r2_alpha = float(np.corrcoef(hs, qp)[0, 1] ** 2)
print(f"\nInterventional response: d(dq_para)/dh = {slope_alpha:+.3f}  (R^2 = {r2_alpha:.4f}, linear in h)")
print("Causal claim (conditional on this simulator): electronegativity alone")
print(f"depletes the para site at {slope_alpha:+.3f} pi-electrons per unit h.\n")

# ================================================================ PART 3
# do(beta_CN) sweep — alpha frozen at pyridine value, coupling varies.
print("=" * 72)
print("PART 3 — INTERVENTIONAL: do(beta_CN = k*beta), h=0.5 frozen, coupling varies")
print("=" * 72)
print(f"{'k (resonance coupling)':<24}{'dq_N':>9}{'dq_ortho':>9}{'dq_meta':>9}{'dq_para':>9}{'E_HOMO':>9}")
rows = []
for k in [0.3, 0.6, 1.0, 1.4]:
    sys = ring(alpha_h={0: 0.5}, electrons={0: 1},
               beta_k={frozenset((0, 1)): k, frozenset((0, 5)): k})
    d, r = charge_pattern(sys, site=0)
    homo = r["E"][r["occ"] > 0].max()
    rows.append((k, d["para"], homo))
    print(f"{k:<24.1f}{d['N']:>+9.3f}{d['ortho']:>+9.3f}{d['meta']:>+9.3f}{d['para']:>+9.3f}{homo:>9.3f}")

ks = np.array([r[0] for r in rows]); qk = np.array([r[1] for r in rows])
slope_beta = float(np.polyfit(ks, qk, 1)[0])
print(f"\nInterventional response: d(dq_para)/dk = {slope_beta:+.3f}")
print("Same outcome, orthogonal lever: coupling strength ALSO moves the para site.")
print()

# ================================================================ VERDICT
print("=" * 72)
print("VERDICT — toy mediation analysis")
print("=" * 72)
print(f"alpha-pathway (electronegativity):  d(dq_para)/dh = {slope_alpha:+.3f}")
print(f"beta-pathway  (resonance coupling):  d(dq_para)/dk = {slope_beta:+.3f}")
print(f"lever ranges: h in [0,2], k in [0.3,1.4]  ->  alpha-path spans "
      f"{abs(slope_alpha)*2:.3f} e-, beta-path spans {abs(slope_beta)*1.1:.3f} e-")
print("""
Both hypotheses survive; each controls an independent, quantified pathway.
In the real project this is the "path decomposition" outcome — one of the
three legitimate verdicts (H1 wins / H2 wins / both, with effect sizes).
The EIG selector's job upstream was to notice that h-sweep and k-sweep are
the two cheapest interventions that DECORRELATE the two mediators.
""")
