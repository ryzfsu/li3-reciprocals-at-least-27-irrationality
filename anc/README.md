# Ancillary replay bundle

The scripts in this directory read the TeX source one level above the
directory and verify the exact endpoint data printed in the paper. The JSON
files are the complete seed arrays and the rational certificate input.

Run from this directory. The checks use Python 3. `check_seed.py` additionally
requires `python-flint`; `check_energy.py`, both replay scripts, and
`check_contour.py` require `mpmath`.

```sh
python3 check_energy.py
python3 check_seed.py
python3 check_profile.py
python3 check_tau.py
python3 check_contour.py
```

`check_energy.py` parses the two tables in Appendices A and B and recomputes the positive
pair and penalty bounds, the negative two-ray pair kernel, every interval box,
and both infinite-tail inequalities. It uses the paper convention
`L=-integral integral log|z-w|`: each directed lower bound is asserted to be
at most the reference upper value, above its printed target, and the pair-energy
sign is checked. It writes the replay receipts to `energy_receipts.json` and
leaves the input `certificate_data.json` unchanged. `check_seed.py` rebuilds both seed
determinants from the displayed definitions, verifies the degree-bounded
identities at 307 and 553 integer points, and checks the factorizations and
residue tables. `check_profile.py` recomputes the conservative arithmetic profile
and its ten low-layer increments. `check_tau.py` uses rational logarithm
enclosures for the four endpoint margins and the two ray slopes.
`check_contour.py` is a finite diagnostic of the residue mechanism for an
undifferentiated `F`, together with a separate triple-zero, sine-bound, and
contour-versus-series consistency check. It is not an `F''` certificate and
does not replace the uniform analytic proof.

The JSON files contain the exact coefficient arrays, density counts, and
interval boxes used by these replays. `results.txt` records source hashes and
the observed outputs.
