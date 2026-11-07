#!/usr/bin/env python3
"""Parse a concatenated JEFF-4.0 ENDF-6 fission-yield file into CONFLUX XML.

The JEFF-4.0 release stores all materials in a single ENDF-6 text file. This
script separates MF=8/MT=454 (independent) and MF=8/MT=459 (cumulative)
records into the XML layout read by :class:`conflux.FPYEngine.FissionIstp`.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path
import xml.etree.ElementTree as ET


def endf_float(field: str) -> float:
    """Read an ENDF real field, with or without an explicit ``E`` exponent."""
    value = field.strip()
    if not value:
        return 0.0
    if "e" not in value.lower():
        for index in range(len(value) - 1, 0, -1):
            if value[index] in "+-" and value[index - 1].isdigit():
                value = f"{value[:index]}e{value[index:]}"
                break
    return float(value)


def fields(line: str) -> list[float]:
    """Return the six 11-column ENDF data fields from one record."""
    return [endf_float(line[index:index + 11]) for index in range(0, 66, 11)]


def element_symbols() -> dict[int, str]:
    """Load atomic-number to element-symbol mappings distributed with CONFLUX."""
    symbols: dict[int, str] = {}
    path = Path(__file__).parents[1] / "betaDB" / "Z_to_element.csv"
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            symbols[int(row["Z"])] = row["Element"]
    return symbols


def parse_sections(source: Path) -> dict[tuple[int, int], dict]:
    """Read every MF=8 fission-yield section and its ENDF LIST records."""
    sections: dict[tuple[int, int], dict] = {}
    active: dict[tuple[int, int], dict] = {}

    with source.open() as handle:
        for line in handle:
            if len(line) < 80:
                continue
            try:
                mat = int(line[66:70])
                mf = int(line[70:72])
                mt = int(line[72:75])
                ns = int(line[75:80])
            except ValueError:
                continue
            if mf != 8 or mt not in (454, 459):
                continue

            key = (mat, mt)
            data = fields(line)
            if ns == 1:
                section = {
                    "fission_za": int(data[0]), "awr": data[1], "le": int(data[2]),
                    "lists": [], "pending": None,
                }
                sections[key] = section
                active[key] = section
                continue

            section = active.get(key)
            if section is None:
                continue
            pending = section["pending"]
            if pending is None:
                nni, nfp = int(data[4]), int(data[5])
                if nni and nfp:
                    section["pending"] = {
                        "energy": data[0], "interpolation": int(data[2]), "nni": nni,
                        "nfp": nfp, "values": [],
                    }
                continue

            pending["values"].extend(data)
            if len(pending["values"]) >= pending["nni"]:
                if pending["nni"] != 4 * pending["nfp"]:
                    raise ValueError(
                        f"MAT {mat}, MT {mt}: NNI={pending['nni']} is not "
                        f"4*NFP={4 * pending['nfp']}"
                    )
                pending["values"] = pending["values"][:pending["nni"]]
                section["lists"].append(pending)
                section["pending"] = None

    incomplete = [key for key, section in sections.items() if section["pending"] is not None]
    if incomplete:
        raise ValueError(f"Incomplete ENDF LIST record(s): {incomplete}")
    return sections


def write_material_xml(mat: int, sections: dict[int, dict], output_dir: Path,
                       symbols: dict[int, str]) -> Path:
    """Write independent and cumulative yield sections for one material."""
    representative = next(iter(sections.values()))
    fission_za = representative["fission_za"]
    z, a = divmod(fission_za, 1000)
    symbol = symbols.get(z, f"Z{z}")
    filename = output_dir / f"nfpy_{mat:04d}_{z}-{symbol}-{a}.xml"
    root = ET.Element(f"nfpy_{mat:04d}_{z}-{symbol}-{a}")

    for mt in (454, 459):
        section = sections.get(mt)
        if section is None:
            continue
        head = ET.SubElement(
            root, "HEAD", FissionZA=str(section["fission_za"]),
            AWR=f"{section['awr']:.12g}", LE=str(section["le"]),
            MT="IFP" if mt == 454 else "CFP",
        )
        for record in section["lists"]:
            node = ET.SubElement(
                head, "LIST", Ei=f"{record['energy']:.12g}",
                Ii=str(record["interpolation"]), NFPi=str(record["nfp"]),
            )
            values = record["values"]
            for offset in range(0, len(values), 4):
                za, fps, yield_value, yield_error = values[offset:offset + 4]
                ET.SubElement(
                    node, "CONT", ZA=str(int(za)), FPS=f"{fps:.12g}",
                    Y=f"{yield_value:.12g}", DY=f"{yield_error:.12g}",
                )

    ET.indent(root, space="\t")
    ET.ElementTree(root).write(filename, encoding="utf-8", xml_declaration=True)
    return filename


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, help="JEFF-4.0 concatenated ENDF-6 yield file")
    parser.add_argument(
        "--output-dir", type=Path,
        default=Path(__file__).parent / "JEFF40",
        help="directory for CONFLUX XML files (default: %(default)s)",
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    sections = parse_sections(args.source)
    materials: dict[int, dict[int, dict]] = defaultdict(dict)
    for (mat, mt), section in sections.items():
        materials[mat][mt] = section

    symbols = element_symbols()
    output = [write_material_xml(mat, material, args.output_dir, symbols)
              for mat, material in sorted(materials.items())]
    print(f"Parsed {len(sections)} MF=8 sections for {len(output)} materials into {args.output_dir}")
    for filename in output:
        print(filename.name)


if __name__ == "__main__":
    main()
