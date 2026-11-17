# INFO — combustion/burnout-cold (a parcel at its removal size that does not burn, e2e)

## What it verifies
A burning model removes a parcel at the removal mass only while it burns: at or above `T-ign` under `Beckstead`,
`T-burn` under `direct`, the temperature the burning branches test (`solout`, `src/lib/Lib_Integration.f90`). A
parcel at the removal mass that does not burn stays an inert parcel. The evaporating models 2 and 4 remove at the
mass alone and are not concerned.

## Case
`burn-box`'s stream (its `bc.txt`, `phase.txt` and `solfile.tec` symlinked, `common/properties.dat`, its settings):
30 µm aluminium parcels, `combustion=Beckstead`, in the 600 K uniform gas of the shared box. Its own `models.txt`:

    A combustion Beckstead K-burn=4.5e-07 n-burn=1.8 X-eff=0.5 T-ign=2000 beta-part=0.3 q-comb=1000000 d-burnout=6.0e-05

`d-burnout` 60 µm is twice the injected diameter, so every parcel is at its removal mass from injection; `T-ign`
2000 K is far above any temperature a parcel reaches in the 600 K gas, so none ignites.

## Pass criterion
`check.py` reads `OUTPUT/source.tec`: the mass source `wdot(A)` must be exactly 0 in every cell, each segment's
mass in equal to its mass out. Measured: 0 in all 1500 cells.

## RED
A removal at the mass alone (the burnout test without the burning temperature) hands every parcel to the gas at its
first step: `wdot(A)` sums to 1.545455e-02 kg/s over the 1500 cells, largest 2.293091e-04 kg/s; `check.py` prints `[FAIL] a cold parcel is removed: the stream loses mass`.

Full test matrix in the family overview [../INFO.md](../INFO.md).
