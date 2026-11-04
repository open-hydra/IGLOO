# INFO — combustion (metal combustion)

Verifies the metal-combustion track: the Beckstead d^n burn law (`model=5`,
`rhsAlCombustion` / `IGLOO_Lib_Combustion::becksteadRate`), selected per material
by the `combustion=Beckstead` token on the material line of `INPUT/phase.txt`. State layout is identical to model 2
(T at Z(7), Al mass at Z(8)); the heat term RELEASES `beta-part*q-comb*|mdot|`
to the particle (sign opposite to the evaporation latent sink).

Physics: `d^n(t) = d0^n − Keff·t` with `Keff = K-burn·X-eff` (oxidizer
effectiveness per `[Beck05]`: X_eff = C_O2 + 0.6·C_H2O + 0.22·C_CO2), giving
`mdot = −(rho_p·pi·Keff/(2n))·d^(3−n)`, n≈1.8. The particle is bitwise-inert
(`mdot = 0` exactly) below `T-ign`.

The oxide cap (`model=8`, `rhsAlCap`): with `xi-cap > 0` and `product=<name>` (a 0-group line of the phase file
whose table zone gives the cap its density, heat capacity and absolute enthalpy datum) the parcel carries the cap
Z(9) beside the metal Z(8). The cap retains `ξν` of the metal burned (`ν = M(Al₂O₃)/(2 M(Al))`), the burn law on
the metal's diameter is scaled by the exposed-surface factor `F(f_ox)` of an oxide lobe `[Gal21]` and by the
stretch `c = ∫₀¹ ds/F` that keeps Beckstead's burn time, and drag, heat and radiation act on the equal-volume
sphere of metal and cap. With `xi-ox0 = x0 > 0` the parcel is injected as metal and oxide at the equal volume,
its cap starting from `x0` of its mass, and the stretch is re-evaluated from that cap,
`c = (d_Al0/d_p0)^n ∫₀¹ ds/F(f_ox(s; x0))`, so the particle burns out in Beckstead's time for its injected diameter.

Citation tags resolve in [../REFERENCES.md](../REFERENCES.md).

## Tests

| id | mode | source | doi | locus | compared | settings_provenance | status |
|---|---|---|---|---|---|---|---|
| CB1 | analytic | `[Beck05]` (differentiated regression) | 10.1007/s10573-005-0067-2 | closed-form `m(t) = rho·pi/6·(d0^n − Keff t)^(3/n)`; complex-step derivative vs `becksteadRate` | rel err at 4 epochs | constant Keff, Tp > T-ign | GREEN (1e-13) |
| CB2 | analytic | derived | n/a | (a) `n·d^(n−1)·dd/dt == −Keff` identity ⇒ `t_b = d0^n/Keff` exact; (b) in-test RK4 of the production rate vs closed-form m(t) at 0.9·t_b | (a) 1e-14; (b) 1e-10 rel-to-m0 | same | GREEN |
| CB3 | analytic | design | n/a | ignition gate: `mdot == 0.0` BITWISE below `T-ign`, burning at/above; heat-release identity `qrel = −beta·q·mdot` | bitwise / 1e-15 | — | GREEN |
| CB4 | correlation (DISCRIMINATOR) | `[Beck05]` | 10.1007/s10573-005-0067-2 | K-burn calibrated at the anchor (100 µm, X_eff=0.21, 1 atm, 300 K) with a = 0.0244, Eq. 7's coefficient for n = 1.5 (0.00735 for n = 1.8), paired here with n = 1.8: a cancels in the verdict, K-burn being calibrated with it; other (D, X_eff) points predicted vs the published *linear* X_eff | 5% — RED unless the code's X_eff exponent is 1: err = abs(X_eff^(1−aXeff) − 1) = 14.4% at X_eff=0.21 if aXeff=0.9 | Beckstead final fit verified against source PDF (beck.txt line 501: t_b·X_eff·p^0.1·T0^0.2 = 0.00735·D^1.8) | GREEN by construction (aXeff=1.0) |
| E2E | e2e box | `[Beck05]` | — | `combustion/burn-box/`: uniform-gas box, kv=1 kt=1, T-ign < Tp0 ⇒ burning from injection; exact closed-form `d^n(x)` + independent RK4 of the heat-release energy balance (Nu=2) + mass-telescoping audit Σwdot = Σnpdot·Δm | truncation-budget tol on d^n; 0.1 K on Tp; 1% telescoping | make_box_case.py --kv 1.0 --kt 1.0 --rp 15e-6 | GREEN |
| E-AIR | e2e box | `[Beck05]` Eq. 7 + Fig. 1 (digitized, drawn, not gated) | 10.1007/s10573-005-0067-2 | `combustion/beckstead-air/`: 25 particles of 10–1000 µm burning out in air at 1 atm and 300 K (kv=1 kt=1, T-ign below Tp0, no heat release); Eq. 7 coded from the paper (K_eff = 1.416955485e-6 m^1.8/s): the input, the burnout witness (bc 0, d_p ≤ d_floor), each exit time against (D0^n − d_p^n)/K_eff and the floor's time, d^n = D0^n − K_eff t on every row, T_p and u_p unchanged | 1e-6; 1e-5; the row's print budget; half an F12.6 unit | make_box_case.py --kv 1.0 --kt 1.0 --tg 300 --ug 0.04 --rp-sweep (25 radii) | GREEN |
| CB5 | analytic | design | n/a | the oxide cap (ODE model 8): one `rhsAlCap` call, `F(9) = −ξν F(8)` | 1e-15 rel.; RED `F(9) = 0`: 1.0 | ξ 0.2, ν typed from the molar masses | GREEN (1.7e-16) |
| CB6 | analytic | `[Beck05]`, `[Gal21]` | 10.1016/j.combustflame.2021.01.039 | `rhsAlCap` at f_ox = 0.1, 0.3, 0.5, 0.9: `F(8)` = `becksteadRate`'s ṁ on the metal's diameter × F(f_ox)·c = 1.294539, 1.185062, 1.084161, 0.785019 | 3e-6 rel.; RED F·c skipped 0.274 | c = `capStretch(0.2, 1.8)` packed | GREEN (4.1e-7) |
| CB7 | analytic | — | n/a | `nuOx` = (2·26.9815385 + 3·15.999)/(2·26.9815385) = 1.889441 | 1e-6; RED M(O) = 16: 5.6e-5 | — | GREEN (0) |
| CB8 | analytic | design | n/a | `computeSource` of a capped parcel at 900 K: E = ṅ(m_Al(c_Al T + hOff_Al) + m_ox(c_res T + hOff_res) + m v²/2) | 1e-14 rel.; RED the metal's c_p and datum for the whole parcel: 1.12 | absolute data of `burn-box-cap`'s table | GREEN (0) |
| CB9 | analytic | `[Gal21]` | 10.1016/j.combustflame.2021.01.039 | `lobeF` at ten f_ox against Gallier's exposed-surface factor at 2570 K (typed to 10 digits); decreasing over 1001 points | 3e-6 abs.; RED weight linear in f 1.4e-5, Gallier's densities dropped 5.3e-2 | 513 nodes uniform in the lobe angle | GREEN (9.5e-7) |
| CB10 | analytic | `[Gal21]` | — | `capStretch(0.2, 1.8)` = ∫₀¹ ds/F = 1.397686 | 1e-5 rel.; RED c = ∫F ds: 0.457 | midpoint, 4096 cells | GREEN (9.2e-7) |
| CB11 | analytic | `[Gal21]` | — | `capStretchX0(0.2, 1.8, 0.1, 2950, 3950)` = (d_Al0/d_p0)^n ∫₀¹ ds/F(f_ox(s; x0)) = 1.388383 (the initial cap r0 = x0/(1 − x0)) | 1e-5 rel.; RED the integrand without the initial cap 4.0e-2, r0 = x0 3.6e-3, the factor dropped 4.9e-2 | midpoint, 4096 cells | GREEN (7.6e-7) |
| E-CAP | e2e box | `[Beck05]`, `[Gal21]` | — | `combustion/burn-box-cap/`: burn-box with `xi-cap=0.2 product=OX` on an absolute two-zone table; every trajectory row against an oracle on Gallier's exact lobe (RK4 of the burn law and the mixture heat balance), the eulerian ρ_p/n_p identity with the cap's mass, the mass telescoping (1 − ξν) of the metal burned | m, d 4e-6 rel.; T 0.1 K; E3 1e-4; telescoping 1e-5 | burn-box's `bc.txt`/`solfile.tec`; `out-file = ALL`, `mollify = off` | GREEN |
| E-RES | e2e box | — (definitions) | — | `combustion/burn-box-residue/`: residue=on, product OX freezing; one residue per parcel burned out, born at its parent's exit record with the cap's mass; the mass and energy books with the birth flux withdrawn; nucleation and relaxation; siblings with d-residue and two groups (IDs 1–50) and with `[IGLOO-ODE] fail-after` (the residue after a solver failure) | m 3e-5, d 5e-6; books 1e-5; T the F12.6 budget | burn-box's `bc.txt`/`solfile.tec`, burn-box-cap's table; `out-file = S`, `out-time = on`, `mollify = off` | GREEN |
| E-OX0 | e2e box | `[Beck05]`, `[Gal21]` | — | `combustion/burn-box-ox0/`: burn-box-residue's capped line with `xi-ox0=0.1` and an inert product; the whole particle's mass and number rate at injection, every trajectory row against an oracle on Gallier's exact lobe with the cap from its initial oxide and the re-evaluated stretch, the burnout inside Beckstead's time for the injected diameter, the residue's mass and flow, the mass and energy books, the scatter quanta | m, d 4e-6 rel.; T 0.1 K; injection m 2e-6; residue m 3e-5; books 1e-5; scatter 1e-3 | burn-box's `bc.txt`/`solfile.tec`, burn-box-cap's table; `out-file = S`, `out-time = on`, `mollify = off` | GREEN |

## Notes

- **p/T0 dependence** of the Beckstead correlation (weak: p^0.1·T0^0.2) is folded
  into the user-supplied `K-burn`; only the X_eff weighting is applied in-code.
- **Burnout** fires at `m ≤ 1e-15 kg` (the metal's mass under model 8): `burn-box` and `burn-box-cap` are sized
  so that the parcels exit at ~50 % d^n loss; `burn-box-residue` burns them out inside the box, as `beckstead-air`
  does at 10–1000 µm. Near burnout,
  unphysical trial states (m ≤ 0 ⇒ d = NaN) are scrubbed to the 1e30 penalty like models 3/4.
- **The burnout residue** (`residue=on`): the cap becomes a parcel of the product at burnout, at a solver failure
  (the parcel reverts to its last accepted step) or at a non-finite state; with `residue=off` the cap goes to the
  gas where the parcel is consumed; a parcel that leaves the domain carries its cap out. `d-residue` sizes the
  residues, the mass flow kept.
- **beta-part** (heat partition to the condensed phase) is weakly constrained in
  the literature — exposed as an input, default 0 (see docs/theory/combustion.md
  sensitivity note); the e2e uses 0.3 with a reduced q-comb for a tame ~140 K
  signal.
- Combustion is mutually exclusive with evaporation (warned, evaporation disabled)
  and with breakup (hard error: `combustion=Beckstead with breakup is not supported`) per material.
