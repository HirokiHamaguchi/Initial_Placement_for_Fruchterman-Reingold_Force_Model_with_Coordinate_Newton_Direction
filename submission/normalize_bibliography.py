"""Apply the agreed final bibliography omissions to the Zotero export."""

from collections import Counter
from pathlib import Path
import re


BIB = Path(__file__).resolve().parents[1] / "doc/main/FruchtermanReingoldByRandomSubspace.bib"


def main() -> None:
    original = BIB.read_bytes()
    text = original.decode("utf-8")
    counts: Counter[str] = Counter()
    records = re.split(r"(?m)(?=^@)", text)
    for index, record in enumerate(records):
        arxiv = re.search(r"(?im)^\s*archiveprefix\s*=\s*\{arxiv\}", record) is not None
        lines = []
        for line in record.splitlines(keepends=True):
            field = re.match(r"^\s*(\w+)\s*=", line)
            name = field.group(1).lower() if field else ""
            if name in {"month", "day", "address", "location"} or (arxiv and name == "eprint"):
                # This export stores fields on one line. Refuse to truncate multiline values.
                if not re.fullmatch(r"\s*\w+\s*=\s*(?:\{[^\r\n]*\}|\w+)\s*,?\s*", line):
                    raise ValueError(f"Unsupported multiline field: {line!r}")
                counts[name] += 1
            else:
                lines.append(line)
        records[index] = "".join(lines)
    result = "".join(records).encode("utf-8")
    # Preserve all other bytes, including DOI, year, publisher, and access dates.
    if result != original:
        BIB.write_bytes(result)
    print(f"Removed fields: {dict(counts)}")


if __name__ == "__main__":
    main()
