"""Step 2: build the direct-requirements matrices and select the 400 industries."""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from supplygraphs import data  # noqa: E402

if __name__ == "__main__":
    summary = data.build()
    print(json.dumps({k: v for k, v in summary.items() if not k.startswith("summary_")}, indent=2))
