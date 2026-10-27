#!/usr/bin/env python3
"""The scatter cloud's weight quantum as the solver sizes it, for oracles.

    dNscat = Ndot tau_ref / (fsample-traj nStreams),  tau_ref = L_ref / V_ref

Ndot is the streams' number flow at injection (each stream's mass flow over its injection mass), L_ref the diagonal
of the mesh's bounding box, V_ref the mean injection speed, fsample-traj 100 by default. logged() reads the value the
run printed (`>> scatter cloud [<material>]: dNscat= ...`, ES10.3); half_ulp() is half a unit in its last place.
"""
import math
import re


def logged(text, mat="A"):
    m = re.search(r"scatter cloud \[%s\]: dNscat=\s*([0-9.Ee+-]+)" % re.escape(mat), text)
    return float(m.group(1)) if m else None


def expected(ndot, l_ref, v_ref, n_streams, traj_sample=100):
    return ndot * (l_ref / v_ref) / (traj_sample * n_streams)


def half_ulp(v):
    """Half a unit in the last place of an ES10.3 print (d.ddd E+xx)."""
    return 0.0 if v == 0.0 else 0.5e-3 * 10.0 ** math.floor(math.log10(abs(v)))
