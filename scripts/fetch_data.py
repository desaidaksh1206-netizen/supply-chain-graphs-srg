"""Step 1: download the BEA tables (via EPA's useeior) into data/raw and verify checksums."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from supplygraphs import data  # noqa: E402

if __name__ == "__main__":
    digests = data.fetch(force="--force" in sys.argv)
    for name, d in digests.items():
        print(f"{d}  {name}")
    print(f"{len(digests)} files verified against data/SHA256SUMS.json")
