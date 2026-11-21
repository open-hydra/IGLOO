# sym-plane-skate — parcels launched along a symmetry plane

**Purpose.** A census control of the integration loop at a symmetry plane (code 300): parcels that graze the
plane, bounce off it or approach it end where they should, at the outlet, and none of the loop's give-up guards
(`stuck in cell`, `no net progress`, the inner-loop and outer caps) or a solver failure ends one of them. The
same streams run at gas order 2 in [`sym-plane-skate-ord2`](../sym-plane-skate-ord2/INFO.md).

**Fixture.** `drag-stokes`'s box (`INPUT/solfile.tec` symlinked: 60×5×5 cells, uniform gas along +x at
U = 10 m/s, w = 0, T = 600 K) with face 5 (z = 0) a symmetry plane and face 1 a wall:
`tools/make_box_case.py --sym-face 5 --no-inlet` writes `INPUT/bc.txt`; `phase.txt` and `properties.dat` from
`common/`. Twelve assigned streams (`[IGLOO-BC]` x, y, z) start at x = 0.004 m, z = 1e-4 m, three per
(speed, incidence) on y = 0.015, 0.025, 0.035 m, heading toward the plane: 0.1 m/s at 1° and 10° (IDs 1–6),
50 m/s at 1° and 10° (IDs 7–12); d = 11.89 µm, ρ_p = 2950 kg/m³, Schiller-Naumann drag, SDIRK4 at its default
tolerances. The gas carries no normal velocity and there is no body force, so the plane alone can raise a parcel
above its start height.

What the streams do: the 50 m/s, 1° streams meet the plane at |v_n|/|v| = 0.016, below the grazing threshold
(0.02): `bcDef` projects them onto it at x = 0.0099 m and they slide to the outlet at the grazing standoff
(z = 1e-6 m). The 50 m/s, 10° streams reflect at x = 0.0046 m and leave at z = 5.9e-3 m. Drag stops the 0.1 m/s
streams 2e-6 m (1°) and 1.7e-5 m (10°) below their start height, short of the plane.

**Gates** (`check.py`).
- C1 no give-up guard and no solver failure: no `stuck in cell`, `no net progress`, `Run_ODESolver err=`,
  `EXIT OF SDIRK4`, `Inner loop >`, `hit outer maxIter`, `non-finite` or `initialized out of domain` line; the run
  ends normally (`Stop condition`, empty `run_err.txt`).
- C2 twelve exit records, IDs 1–12 once each, every one at the outlet x = 0.15 m.
- C3 the plane was reached: IDs 7–9 leave sliding on it (z ≤ 1e-5 m), IDs 10–12 leave reflected (z > 1e-4 m),
  IDs 1–6 leave between the plane and their start height.

**RED.** The plane made an exit (`bcDef`'s 300 branch skipped): IDs 7–12 leave at z = 0, x = 0.0099 and
0.0046 m (C2, C3), at both gas orders. The residency guard tripped by a single segment in one cell: 12
`stuck in cell` lines at gas order 1 (C1, C2); gas order 2 resets the count at every crossing and stays green.

**Not gated here.** A parcel pressed onto the plane (a body force toward it), whose reflections repeat inside
one cell; the drift of the gas toward the plane, gated by `test_ghost_bc` and `wall-approach`.
