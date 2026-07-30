# INFO — combustion (family M1)

Verifies the metal-combustion track: the Beckstead d^n burn law (`model=5`,
`rhsAlCombustion` / `IGLOO_Lib_Combustion::becksteadRate`), selected per material
by `[GPB-PhaseX] combustion = Beckstead`. State layout is identical to model 2
(T at Z(7), Al mass at Z(8)); the heat term RELEASES `beta-part*q-comb*|mdot|`
to the particle (sign opposite to the evaporation latent sink).

Physics: `d^n(t) = d0^n − Keff·t` with `Keff = K-burn·X-eff` (oxidizer
effectiveness per `[Beck05]`: X_eff = C_O2 + 0.6·C_H2O + 0.22·C_CO2), giving
`mdot = −(rho_p·pi·Keff/(2n))·d^(3−n)`, n≈1.8. The particle is bitwise-inert
(`mdot = 0` exactly) below `T-ign`.

Citation tags resolve in [../REFERENCES.md](../REFERENCES.md).

## Tests

| id | mode | source | doi | locus | compared | settings_provenance | status |
|---|---|---|---|---|---|---|---|
| CB1 | analytic | `[Beck05]` (differentiated regression) | 10.1007/s10573-005-0067-2 | closed-form `m(t) = rho·pi/6·(d0^n − Keff t)^(3/n)`; complex-step derivative vs `becksteadRate` | rel err at 4 epochs | constant Keff, Tp > T-ign | GREEN (1e-13) |
| CB2 | analytic | derived | n/a | (a) `n·d^(n−1)·dd/dt == −Keff` identity ⇒ `t_b = d0^n/Keff` exact; (b) in-test RK4 of the production rate vs closed-form m(t) at 0.9·t_b | (a) 1e-14; (b) 1e-10 rel-to-m0 | same | GREEN |
| CB3 | analytic | design (plan §3-M1) | n/a | ignition gate: `mdot == 0.0` BITWISE below `T-ign`, burning at/above; heat-release identity `qrel = −beta·q·mdot` | bitwise / 1e-15 | — | GREEN |
| CB4 | correlation (DISCRIMINATOR) | `[Beck05]` | 10.1007/s10573-005-0067-2 | K-burn calibrated at the anchor (100 µm, X_eff=0.21, 1 atm, 300 K, a=0.0244 ms·µm^−1.8); other (D, X_eff) points predicted vs the published *linear* X_eff | 5% — RED unless the code's X_eff exponent is 1: err = abs(X_eff^(1−aXeff) − 1) = 14.4% at X_eff=0.21 if aXeff=0.9 | Beckstead final fit verified against source PDF (beck.txt line 501: t_b·X_eff·p^0.1·T0^0.2 = 0.00735·D^1.8) | GREEN by construction (aXeff=1.0; rerun on build host) |
| E2E | e2e box | `[Beck05]` | — | `combustion/burn-box/`: uniform-gas box, kv=1 kt=1, T-ign < Tp0 ⇒ burning from injection; exact closed-form `d^n(x)` + independent RK4 of the heat-release energy balance (Nu=2) + mass-telescoping audit Σwdot = Σnpdot·Δm | truncation-budget tol on d^n; 0.1 K on Tp; 1% telescoping | make_box_case.py --kv 1.0 --kt 1.0 --rp 15e-6 | GREEN |

## Notes

- **p/T0 dependence** of the Beckstead correlation (weak: p^0.1·T0^0.2) is folded
  into the user-supplied `K-burn`; only the X_eff weighting is applied in-code.
- **Burnout is out of scope** (phase M2): the e2e case is sized so particles exit
  at ~50% d^n loss. Near-burnout, unphysical trial states (m ≤ 0 ⇒ d = NaN) are
  scrubbed to the 1e30 penalty like models 3/4.
- **beta-part** (heat partition to the condensed phase) is weakly constrained in
  the literature — exposed as an input, default 0 (see docs/theory/combustion.md
  sensitivity note); the e2e uses 0.3 with a reduced q-comb for a tame ~140 K
  signal.
- Combustion is mutually exclusive with evaporation (warned, evaporation dropped)
  and with breakup (hard error, coupling deferred) per material.
