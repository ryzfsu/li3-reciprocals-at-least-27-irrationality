# Timestamps

[OpenTimestamps](https://opentimestamps.org) proofs for the assets of release `v1`.
Each proof commits to the SHA-256 of the file below.

| file | SHA-256 |
|---|---|
| `li3-reciprocals-at-least-27-irrationality-v1.pdf` | `c2e397777c32b215719a6c9cc9c3e08afbc27a7084179423865550879a1100d4` |
| `li3-reciprocals-at-least-27-irrationality-v1.tar.gz` | `9a2739ee2dcbf845ce2513cfb09d412abd9b2026ace22f5a1bd5d45ae0134d19` |

To verify, download the release asset into this directory and run

    pip install opentimestamps-client
    ots verify timestamps/<file>.ots
