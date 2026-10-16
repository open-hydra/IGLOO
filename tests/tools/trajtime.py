#!/usr/bin/env python3
"""The parcel time along a trajectory written without the out-time column: dt = dx / u by the trapezoid rule."""


def cumulative_time(xs, us):
    """[0, t1, ...] with t[k] = t[k-1] + (x[k] - x[k-1]) / (0.5 * (u[k-1] + u[k]))."""
    t = [0.0]
    for k in range(1, len(xs)):
        t.append(t[-1] + (xs[k] - xs[k - 1]) / (0.5 * (us[k - 1] + us[k])))
    return t
