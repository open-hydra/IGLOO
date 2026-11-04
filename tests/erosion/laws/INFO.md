# INFO — erosion/laws (the wall-erosion laws of `[IGLOO-Erosion]`)

Production `okaVolume`, `ngRatio` and `erosionRatio` (`IGLOO_Lib_Erosion`, `src/lib/Lib_Erosion.f90`), called as
`depositErosion` calls them, against the closed forms of their papers and the anchors the papers state.

## References

- `[Oka05a]` Y.I. Oka, K. Okamura, T. Yoshida, Practical estimation of erosion damage caused by solid particle
  impact, Part 1: effects of impact parameters on a predictive equation, *Wear* 259 (2005) 95–101: eqs. 1–2,
  Tables 2–3, Fig. 2.
- `[Oka05b]` Y.I. Oka, T. Yoshida, Part 2: mechanical properties of materials directly associated with erosion
  damage, *Wear* 259 (2005) 102–109: eqs. 1–4, Table 1.
- `[NG68]` J.H. Neilson, A. Gilchrist, Erosion by a stream of solid particles, *Wear* 11 (1968) 111–122: eqs. 3–7,
  Table I, Figs. 5 and 11 (constants in ft·lbf/g and ft/s; 1 ft·lbf/g = 1355.818 J/kg, 1 ft = 0.3048 m).

## Tests (`test_erosion_laws.f90`, ctest `test_erosion_laws`, labels `unit;erosion`)

| id | compared | tolerance | measured |
|---|---|---|---|
| O1 | Oka's g(α) for silica (Table 3: n1 = 0.71 Hv^0.14, n2 = 2.4 Hv^−0.94) on aluminium ALR (Hv 0.40, Table 2): the angle of its maximum, which the paper puts at 20° for silica on aluminium; g(90°) = 1, g(0) = 0 | [19°, 21°]; 4 ulp; exact | 20.29°; 0; 0 |
| O2 | the same on the stainless steel 304R (Hv 2.90): the maximum of g and its angle, the curve of Fig. 2b | [1.25, 1.40]; [35°, 45°] | 1.301 at 40.51° |
| O3 | E90 = K (a Hv)^(k1 b) (v/v′)^k2 (D/D′)^k3 (Part 2 eq. 4) against its closed form at four (v, D, a, b), and with a = b = 1 Part 1's K Hv^k1 (v/v′)^k2 (D/D′)^k3; silica's K = 65 mm³/kg at v′ = 104 m/s, D′ = 326 µm, Hv = 1 (Table 1) | 1e-14 relative; 1e-12 | within |
| O4 | `erosionRatio` of an Oka material = 1e-9 m³/mm³ × wall density × `okaVolume`, nine angles | 1e-15 relative | within |
| N1 | Neilson-Gilchrist, aluminium and 210 µm alumina at 354, 424, 494 ft/s (Table I): E(90°) = V²/(2ε); the branches of eqs. 3 and 4 meet at α0 = π/(2n); the maximum of E (0.001° grid) lies below α0 and satisfies eq. 7, φ/ε = sin(nα) − n cos(nα)/(2 tan α); at 424 ft/s the maximum and its angle (Fig. 5) | 1e-14; 1e-9; 2e-3 on φ/ε; [1.5e-3, 1.7e-3], [16.5°, 18°] | eq. 7: 0.4520 for 0.4523, 0.4268 for 0.4264, 0.3926 for 0.3927; 1.634e-3 at 17.19° |
| N2 | glass spheres on glass (Fig. 11: V = 260 ft/s, K = 149.3 ft/s, ε = 46 ft·lbf/g, no cutting): no erosion below asin(K/V) = 35.05°, "no erosion occurs at angles of attack less than 35°"; erosion at 35.1°; E(90°) = (V − K)²/(2ε) | below 1e-20 of E(90°); 1e-14 | within |
| N3 | cutting alone (ε → ∞, n = 4): V² cos²α sin(nα)/(2φ) below α0, V² cos²α/(2φ) above, every 5° | 1e-14 relative | within |
| N4 | `erosionRatio` of a Neilson-Gilchrist material = `ngRatio` | 1e-15 relative | within |

The figure readings (O2, the Fig. 5 band of N1) are loose by design; the closed forms and the paper's own relations
(eq. 7, the 35° threshold, the 20° maximum) carry the gate.

## RED (measured: each mutation in `Lib_Erosion.f90`, the tree rebuilt)

| mutation | what fails |
|---|---|
| Oka's g without Hv: (1 + (1 − sin α))^n2 | O1 (maximum at 11.4°), O2 (75.0°, 1.001) |
| E90 without b: (a Hv)^k1 | O3 (up to 0.31 relative) |
| E90 without (D/D′)^k3 | O3 (up to 0.25) |
| 1e-6 for 1e-9 m³/mm³ | O4 (999 relative) |
| Neilson-Gilchrist: sin(nα) above α0 too | N3 (up to 1.98) |
| no threshold: (V sin α − K)² below K | N2 (erosion at 0.1°) |
| α0 = π/n | N3 |
| φ and ε swapped in the deformation term | N1, N2, N3 |
