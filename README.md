# Trilogarithm at reciprocal integers of absolute value at least 27

This directory contains an unconditional proof that
\(\operatorname{Li}_3(1/a)\) is irrational for every integer \(a\ge 27\)
and every integer \(a\le -27\).  The source is split into `main.tex`, `body.tex`,
`rank-lemma.tex`, `tower-proof.tex` and `appendices.tex`; the compiled manuscript is produced from
this directory.

The `anc/` directory is a compact replay bundle. Run its checks
from that directory with Python 3:

```sh
python3 check_energy.py
python3 check_seed.py
python3 check_profile.py
python3 check_tau.py
python3 check_contour.py
```

The energy replay recomputes both 256-cell endpoint certificates from the
tables in Appendices A and B and the interval data, including exact masses,
caps, pair kernels, and tails. It writes the resulting receipts to
`anc/energy_receipts.json` and leaves `certificate_data.json` unchanged. It uses the paper convention
`L=-integral integral log|z-w|` and asserts lower-versus-upper consistency,
the printed targets, and the pair-energy sign. The seed replay rebuilds the two coefficient
towers and then checks their modular residues. `check_seed.py` requires
`python-flint`; the energy and contour replays require `mpmath`.
The tau replay checks the four endpoint inequalities and ray-slope
thresholds. The contour replay is a finite consistency diagnostic; the
uniform analytic estimate is proved in the manuscript. `results.txt` records
the source hashes and the observed outputs.
