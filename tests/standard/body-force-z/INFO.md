# body-force-z — body-force drift across k-faces

**Purpose.** Every other box case moves its parcels in x and y only; here the drift is along z, through the
k-faces of the mesh, and the gas temperature varies along z so that the cell holding a parcel is observable: at
gas-order 1 a parcel relaxes toward the temperature of the k-layer that holds it.

**Fixture** (`python3 ../../tools/make_box_case.py INPUT --kv 1.0 --zgradient 2000`; `solfile.tec` committed,
`bc.txt` is `body-force`'s, symlinked, `phase.txt` and `properties.dat` from `common/`): `body-force`'s box and
parcels (25 at the inlet cell centres, kV = 1, kT = 1, d = 11.89 µm, Stokes, Ranz-Marshall) with the gas
temperature T_c(k) = 600 + 2000 z_c(k) K at the cell centres (610 … 690 K by layer). `input.ini` is
`body-force`'s with `body-accel = 0.0 0.0 -600.0` and `out-time = on`: τ = 1.287e-3 s, v_∞ = g_z τ = −0.772 m/s,
a drift of 10.6 mm over the 15 ms crossing, so every parcel crosses one k-face (cells 10 mm thick) and the bottom
row (z_c = 5 mm) leaves through z = 0.

**Oracle** (`check.py`). With Stokes drag the motion decouples from the temperature: u ≡ u_g, v ≡ 0,
v_z(t) = v_∞(1 − e^{−t/τ}), z(t) = z₀ + v_∞ t − v_∞ τ(1 − e^{−t/τ}). The temperature is a classical RK4 of
dT/dt = K(t)(T_c(k(t)) − T), K = 6 Nu k_g/(ρ_p c_p d²), Nu = 2 + 0.6 Re^½ Pr^⅓, Re = ρ_g |v_z| d/μ, with k(t)
the layer of z(t), switched at the k-face planes (bisection on the closed form to 1e-14 s); checked at half the
step (6.9e-12 K).

**Gates.**
- B0 each first row is the injection state (x = 0, the inlet cell centre, (u_g, 0, 0), T = T_c(k₀), t = 0).
- B1 v_z on the closed form (`body-force`'s tolerance), u = u_g, v = 0, d = d_p on every row; B2 z on the closed
  form within 2e-6 m.
- K1 every parcel that reaches x = 0.15 drifts more than a cell height (it crosses a k-face).
- K2 the five bottom-row parcels leave at z = 0 at the closed-form x where z(x) = 0, within the print (5e-7 m
  + 2e-8 m); measured 4.1e-7 m. The other twenty leave at x = 0.15.
- K3 T on every row within tol_T = 5e-7 + |dT/dt|·5e-9 t + 1e-6 K (the T print, the time print, the integrator)
  plus, for each k-face crossed before the row, K·|T_c(new) − T|·ε/|v|, ε = 1e-8 m: a segment that crosses a face
  is accepted once a step moves less than ε (`Lib_Integration`), so it ends up to ε/|v| ≈ 1e-9 s past the face,
  still in the old cell, where the gas temperature differs by a layer step (≈ 1.5e-5 K per crossing). Measured:
  every crossing row 4.39e-10 s past the plane; moving the oracle's switches to those times leaves 0.33 of the
  print tolerance, so the term is the whole of the excess. K3 worst 0.385 of its tolerance.

**RED.** A parcel that keeps sampling its previous cell after a crossing (gas-order 1: the gas properties not
refreshed for the new cell): K3 fails on all 20 crossing parcels, by up to 19.9 K (1.2e6 × the tolerance);
B1–K2 hold. The k index of the crossing itself is not observable: exchanging the k − 1 / k + 1 steps of the
interior-face step in `checkBoundary`, or of the first guess of the cell search (`searchInBlock`, tier 0), leaves
every number unchanged — the search's ±1 box around the old cell finds the cell that holds the parcel either
way. The gate pins the transport across k-faces and the cell of the gas sample, not that step.
