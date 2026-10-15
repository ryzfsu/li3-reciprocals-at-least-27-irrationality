#!/usr/bin/env python3
"""Recompute both directed endpoint certificates from the public tables."""
from fractions import Fraction
import json
import re
from pathlib import Path

from replay_positive import replay as replay_positive
from replay_negative import replay as replay_negative

HERE = Path(__file__).resolve().parent


def table(tex, label, zeros=()):
    """Rows of a printed density table; indices in `zeros` are the rows the
    caption states to be zero and omits."""
    marker = r"\label{" + label + r"}"
    at = tex.index(marker)
    start = tex.rfind(r"\begin{longtable}", 0, at)
    end = tex.index(r"\end{longtable}", at)
    out = {}
    for line in tex[start:end].splitlines():
        line = line.split("%", 1)[0].strip()
        if not line.endswith(r"\\"):
            continue
        vals = []
        for field in line[:-2].split("&"):
            field = field.replace(r"\,", "").replace(" ", "")
            if re.fullmatch(r"-?\d+", field):
                vals.append(int(field))
            else:
                vals.append(None)
        for pos in range(0, len(vals), 3):
            row = vals[pos:pos + 3]
            if len(row) == 3 and all(x is not None for x in row):
                i, left, right = row
                if 0 <= i < 256:
                    assert i not in out
                    out[i] = (left, right)
    for i in zeros:
        assert i not in out
        out[i] = (0, 0)
    assert set(out) == set(range(256))
    return [out[i] for i in range(256)]


def main():
    tex = (HERE.parent / "appendices.tex").read_text()
    data = json.loads((HERE / "certificate_data.json").read_text())
    positive = table(tex, "tab:positive-density", range(162, 253))
    negative = table(tex, "tab:negative-density", range(172, 256))
    den = int(data["positive"]["denominator"])
    positive_receipts = {}
    levels = {
        27: Fraction(1364590351, 268435456),
        28: Fraction(687124203, 134217728),
    }
    for q in (27, 28):
        rows = [pair[0 if q == 27 else 1] for pair in positive]
        assert rows == data["positive"]["counts"][str(q)]
        assert sum(rows) == 128 * den
        assert all(0 <= x <= 3 * den for x in rows)
        receipt = replay_positive({
            "q": Fraction(q),
            "dens": [Fraction(x, den) for x in rows],
            "ell": levels[q],
        })
        positive_receipts[str(q)] = receipt
        print(f"positive q={q}: table, mass, cap, pair, penalty, tail=PASS")

    negative_receipts = {}
    for q in (27, 28):
        rows = [pair[0 if q == 27 else 1] for pair in negative]
        record = data["negative"][str(q)]
        assert rows == record["counts"]
        assert sum(rows) == int(record["denominator"])
        assert all(x >= 0 for x in rows)
        negative_receipts[str(q)] = replay_negative(q, HERE)
        print(f"negative q={q}: table, mass, boxes, pair, tail=PASS")

    receipts = {
        "positive": {"directed_bounds": positive_receipts},
        "negative": {"replay_receipt": negative_receipts},
    }
    (HERE / "energy_receipts.json").write_text(
        json.dumps(receipts, indent=2) + "\n"
    )
    print("energy VERIFY_OK")


if __name__ == "__main__":
    main()
