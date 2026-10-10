"""Write feature analysis results."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import List

from .models import FeatureResult


def write_results(results: List[FeatureResult], out_path: str) -> None:
    fieldnames = list(FeatureResult.__dataclass_fields__.keys())
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for result in results:
            row = result.__dict__.copy()
            row["breadth"] = f"{row['breadth']:.4f}"
            row["mean_depth"] = f"{row['mean_depth']:.4f}"
            row["mean_identity"] = f"{row['mean_identity']:.4f}"
            writer.writerow(row)
