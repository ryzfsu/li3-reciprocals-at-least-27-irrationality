#!/usr/bin/env python3
"""Directed rational checks for the four endpoint arithmetic margins."""

from fractions import Fraction as F


S_PLUS = F(106301329053850833773025543473,
           23199292058634066375000000000)
S_MINUS = F(8108653419015955202861020329389,
            1747888832138152054871520000000)


def log_interval(x, terms=64):
    """Rational lower and upper bounds for log(x), x in [1,2]."""
    x = F(x)
    assert F(1) <= x <= F(2)
    z = (x - 1) / (x + 1)
    lo = 2 * sum(z ** (2 * k + 1) / (2 * k + 1) for k in range(terms))
    rem = 2 * z ** (2 * terms + 1) / ((2 * terms + 1) * (1 - z * z))
    return lo, lo + rem


def log_integer(n):
    """Reduce an integer to [1,2] and enclose its logarithm."""
    n = int(n)
    assert n >= 1
    powers = 0
    x = F(n)
    while x > 2:
        x /= 2
        powers += 1
    l2 = log_interval(F(2))
    lx = log_interval(x)
    return powers * l2[0] + lx[0], powers * l2[1] + lx[1]


def main():
    energies = {
        ("+", 27): F("4.9972574688"),
        ("+", 28): F("5.0224702721"),
        ("-", 27): F("4.9720000000"),
        ("-", 28): F("4.9900000000"),
    }
    margins = {}
    for sign, q in (("+", 27), ("+", 28), ("-", 27), ("-", 28)):
        lq = log_integer(q)
        other = log_integer(q - 1 if sign == "+" else q + 1)
        if sign == "+":
            tau_hi = S_PLUS + F(1, 4) * other[1] - F(1, 8) * lq[0]
        else:
            tau_hi = S_MINUS + F(225, 1024) * other[1] - F(1, 8) * lq[0]
        margin = energies[(sign, q)] - tau_hi
        assert margin > 0
        margins[f"{sign}{q}"] = margin
        print(f"tau({sign}{q}) upper and energy margin: PASS")
    assert margins["+27"] > F("0.0126")
    assert margins["-27"] > F("0.0126")
    assert F(1, 2) - F(27, 200) == F(73, 200)
    assert F(5, 16) - F(97, 1024) == F(223, 1024)
    print("slope_plus_net=73/200 slope_minus_net=223/1024 PASS")


if __name__ == "__main__":
    main()
