# REFERENCES — IGLOO verification suite master bibliography

Deduplicated bibliography keyed by short citation tags. Every `[tag]` cited in
any per-family `INFO.md` must resolve here. The T9 CI gate enforces this.

> **Status note (T0 seed):** The entries below were reproduced from the
> verification-phase1-plan §2.5 worked-examples table. DOIs are best-guess
> from the published metadata commonly associated with each reference and
> **must be re-verified against the original sources** before any external
> citation (paper, presentation, archival report). Mark as `verified` in the
> right column once an entry has been cross-checked against the source PDF.

## Bibliography

| Tag | Full citation | DOI / ISBN / id | Verified? |
|---|---|---|---|
| `[CGW78]` | Clift, R.; Grace, J. R.; Weber, M. E. *Bubbles, Drops, and Particles.* Academic Press, 1978. (Dover reprint 2005.) | ISBN 0-12-176950-X (Academic); ISBN 0-486-44580-1 (Dover) | no |
| `[SN33]` | Schiller, L.; Naumann, A. "Über die grundlegenden Berechnungen bei der Schwerkraftaufbereitung." *Z. Ver. Dtsch. Ing.*, 77, 1933, pp. 318–320. | (historical; no DOI) | no |
| `[RM52]` | Ranz, W. E.; Marshall, W. R. "Evaporation from drops, Parts I & II." *Chem. Eng. Prog.*, 48 (1952), Part I pp. 141–146, Part II pp. 173–180. | (no DOI assigned) | no |
| `[Lef89]` | Lefebvre, A. H. *Atomization and Sprays.* Hemisphere/Taylor & Francis, 1989. | ISBN 0-89116-697-5 | no |
| `[ORA87]` | O'Rourke, P. J.; Amsden, A. A. "The TAB Method for Numerical Calculation of Spray Droplet Breakup." *SAE Technical Paper* 872089, 1987. | DOI:10.4271/872089 | yes (2026-07-15, source PDF; Ck/Cd/Cf/Cb/K + eq. 22 child balance) |
| `[RD86]` | Reitz, R. D.; Diwakar, R. "Effect of Drop Breakup on Fuel Sprays." *SAE Technical Paper* 860469, 1986. | DOI:10.4271/860469 | yes (2026-07-15, source PDF; forms/criteria/stable sizes; defaults follow [RD87] lineage) |
| `[PE87]` | Pilch, M.; Erdman, C. A. "Use of breakup time data and velocity history data to predict the maximum size of stable fragments for acceleration-induced breakup of a liquid drop." *Int. J. Multiphase Flow*, 13(6), 1987, pp. 741–757. | DOI:10.1016/0301-9322(87)90063-2 | yes (2026-07-15, source PDF; all constants incl. 0.75/3 vel-poly, 4.5/1.2/1.64) |
| `[Reitz87]` | Reitz, R. D. "Modeling atomization processes in high-pressure vaporizing sprays." *Atomisation and Spray Technology*, 3, 1987, pp. 309–337. | (no DOI) | yes (2026-07-15, source PDF; KH Λ/Ω fits, Z/T defs — A13) |
| `[BR99]` | Beale, J. C.; Reitz, R. D. "Modeling spray atomization with the Kelvin-Helmholtz/Rayleigh-Taylor hybrid model." *Atomization and Sprays*, 9, 1999, pp. 623–650. | DOI:10.1615/AtomizSpr.v9.i6.40 | yes (2026-07-15, source PDF `KHRT-hybrid.pdf`; RT eqs 8–11 — A14. `Reitz-KHRT1999.pdf` is a 0-char scan) |
| `[Tan97]` | Tanner, F. X. "Liquid Jet Atomization and Droplet Breakup Modeling of Non-Evaporating Diesel Fuel Sprays." *SAE Technical Paper* 970050, 1997. | DOI:10.4271/970050 | yes (2026-07-15, source PDF; ETAB eqs 1–11 — A15) |
| `[Tan98]` | Tanner, F. X.; Weisser, G. "Simulation of Liquid Jet Atomization for Fuel Sprays by Means of a Cascade Drop Breakup Model." *SAE Technical Paper* 980808, 1998. | DOI:10.4271/980808 | yes (2026-07-15, source PDF; k1=k2=0.2222, We_t=80 table) |
| `[RD87]` | Reitz, R. D.; Diwakar, R. "Structure of High-Pressure Fuel Sprays." *SAE Technical Paper* 870598, 1987. | DOI:10.4271/870598 | yes (2026-07-23, source PDF `Reitz-StructureHighPressureFuel-1987.pdf`; bag D=π Eq.7, stripping **C=20 curve-fit** Eqs 6/8, criteria We_r>6 / We_r/√Re>0.5 — B-VAL-4) |
| `[MHB98]` | Miller, R. S.; Harstad, K.; Bellan, J. "Evaluation of equilibrium and non-equilibrium evaporation models for many-droplet gas-liquid flow simulations." *Int. J. Multiphase Flow*, 24(6), 1998, pp. 1025–1055. | DOI:10.1016/S0301-9322(98)00028-7 | yes (2026-07-13, source PDF; **2026-07-29 Appendix A read**: publishes all property correlations + sources — air & water Harpole (1981), decane Abramzon & Sirignano (1989), heptane Park & Aggarwal (1995), benzene/hexane Reid et al. (1987); species with no `Γ`/`Sc` (incl. **water**) run at `Le=1`; properties frozen at `T_R=T_WB` eq. (27); `p_sat` = single-anchor C-C eq. (10); M7 energy eq. carries `f₂=β/(e^β−1)` eq. (19)) |
| `[WL92]` | Wong, S. C.; Lin, A. C. "Internal temperature distributions of droplets vaporizing in high-temperature convective flows." *J. Fluid Mech.*, 237, 1992, pp. 671–687. | DOI:10.1017/S0022112092003574 | not needed as a source: the decane size+temperature data enter only through `[MHB98]` Fig. 4, which is **pixel-measured** from the paper's own plot (E-VAL-2b `test_mhb98_decane`); the exp circles are not gated |
| `[Harp81]` | Harpole, G. M. "Droplet evaporation in high temperature environments." *J. Heat Transfer*, 103, 1981. | (via `[MHB98]` App. A) | not read directly — cited as the source of MHB98's **air and water** property correlations, which are reproduced verbatim in their Appendix A and used by `mhb98-water` |
| `[AS89]` | Abramzon, B.; Sirignano, W. A. "Droplet vaporization model for spray combustion calculations." *Int. J. Heat Mass Transfer*, 32(9), 1989, pp. 1605–1618. | DOI:10.1016/0017-9310(89)90043-4 | not read directly — cited as the source of MHB98's **decane** property set (App. A), used by `test_mhb98_decane`. NOTE their Γ_V gives Le=4.40 while MHB98's own Fig. 4 requires Le≈1 (see FINDINGS 2026-07-29) |
| `[Beck05]` | Beckstead, M. W. "Correlating Aluminum Burning Times." *Combust. Explos. Shock Waves*, 41(5), 2005, pp. 533–546. | DOI:10.1007/s10573-005-0067-2 | yes (2026-07-13, source PDF; X_eff exponent = 1.0, final fit) |
| `[Shim06]` | Shimada, T. et al. "Computational Fluid Dynamics of Multiphase Flows in Solid Rocket Motors." *JAXA Special Publication* JAXA-SP-05-035E, 2006. **Proximate provenance of the IGLOO drag (eqs. 12–24) and Nu (eqs. 45–50) catalogs.** | ISSN JAXA-SP-05-035E | yes (2026-07-15, source PDF in papers/; A11/A12 found against it) |
| `[SP8039]` | "Solid Rocket Motor Performance Analysis and Prediction." *NASA Space Vehicle Design Criteria*, NASA-SP-8039, 1971. (Primary for the "JAXA1–4" Nu correlations and the Crowe drag entry, per [Shim06] ref [12].) | NASA-SP-8039 | yes (2026-07-16, source PDF p. 26 Tables I/II; confirms A7/A11/A12; its own C-H entry has a sign typo) |
| `[Put61]` | Putnam, A. "Integrable form of droplet drag coefficient." *ARS Journal*, 31, 1961, p. 1467. | (no DOI) | closed by decision 2026-07-16 (unavailable; [Shim06] eq. 17 = 0.4392 adopted as authoritative) |
| `[Hen76]` | Henderson, C. B. "Drag Coefficient of Spheres in Continuum and Rarefied Flows." *AIAA Journal*, 14(6), 1976. | DOI:10.2514/3.61409 | closed by decision 2026-07-16 (unavailable; blend 4/3 kept — [Shim06] eq. 22's 3/4 is internally inconsistent) |
| `[WY66]` | Wen, C. Y.; Yu, Y. H. *Chem. Engr. Prog. Symp. Series*, 62, 1966, p. 100. | (no DOI) | closed by decision 2026-07-16 (unavailable; [Shim06] eq. 16 = 0.43 adopted as authoritative) |

## Adding a new reference

1. Append a row above with a new short tag (e.g. `[Sh06]` for Shashank et al. 2006).
2. Include full citation + DOI (or ISBN for books, or `(no DOI)` for sources without one).
3. Set `Verified?` to `no`; flip to `yes` once cross-checked against the source PDF.
4. Reference the tag from any `INFO.md` that uses the source.
| `[TC2012]` | Tonini, S.; Cossali, G. E. "An analytical model of liquid drop evaporation in gaseous environment." *Int. J. Thermal Sciences*, 57, 2012, pp. 45–53. | DOI:10.1016/j.ijthermalsci.2012.01.017 | yes (2026-07-13, source PDF; m̂=m_ev/(4πRd·Dv·ρ∞), Stefan-Fuchs) |
| `[ATC24]` | Antonov, D. V.; Tonini, S.; Cossali, G. E.; Al Qubeissi, M.; Sazhin, S. S. "Three approaches to modelling the heating and evaporation of drops." *Int. J. Multiphase Flow*, 179, 2024, 104922 (open access). | DOI:10.1016/j.ijmultiphaseflow.2024.104922 | yes (2026-07-13, OA PDF; "Eq 9" = TC Stefan-Fuchs, OCR caveats) |
| `[Vie15]` | Vié, A.; Doisneau, F.; Massot, M. "On the Anisotropic Gaussian velocity closure for inertial-particle laden flows." *Commun. Comput. Phys.*, 17(1), 2015, pp. 1–46. | DOI:10.4208/cicp.021213.140514a | yes (2026-07-20, source PDF papers/Vie.pdf; §5.1 Eqs 5.1–5.4 compressive-field Lagrangian solution; Eq. 5.4 oscillatory branch has a sin sign typo — fails V_p0=0, ODE authoritative, `+sin` confirmed by scipy to 1e-13) |
| `[God53]` | Godsave, G. A. E. "Studies of the combustion of drops in a fuel spray — the burning of single drops of fuel." *4th Symp. (Int.) on Combustion*, 4(1), 1953, pp. 818–830. | DOI:10.1016/S0082-0784(53)80107-4 | no (best-guess DOI; classical d²-law origin, cited by d2law) |
| `[Spa53]` | Spalding, D. B. "The combustion of liquid fuels." *4th Symp. (Int.) on Combustion*, 4(1), 1953, pp. 847–864. | DOI:10.1016/S0082-0784(53)80110-4 | no (best-guess DOI; classical d²-law origin, cited by d2law) |
