#!/usr/bin/env python3
"""Configure genuine Arial faces and validate editable figure exports.

No font files are distributed or copied. ``check`` creates a per-run Fontconfig
configuration and shell environment for figure rendering. ``verify`` checks actual
PDF font names/embedding and SVG text families, rather than trusting labels in
the plotting code. Both commands retain validation results when a check fails.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET


FACES = {
    "regular": (False, False), "bold": (True, False),
    "italic": (False, True), "bolditalic": (True, True),
}
FC_FIELDS = ("family", "style", "weight", "slant", "file", "index", "postscriptname")
FC_FORMAT = "\x1f".join("%{" + name + "}" for name in FC_FIELDS) + "\n"
PDF_ARIAL_NAMES = {
    "Arial", "ArialMT", "Arial-Regular", "Arial-RegularMT", "Arial-Bold",
    "Arial-BoldMT", "Arial-Italic", "Arial-ItalicMT", "Arial-BoldItalic",
    "Arial-BoldItalicMT",
}


class FontError(RuntimeError):
    pass


def run(args, env=None):
    result = subprocess.run(args, env=env, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, check=False)
    if result.returncode:
        raise FontError(f"Command failed ({result.returncode}): {shlex.join(map(str, args))}\n"
                        + result.stderr.strip())
    return result.stdout


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def save_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".tmp")
    temp.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n")
    temp.replace(path)


def fc_records(output):
    records = []
    for line in output.splitlines():
        values = line.split("\x1f")
        if len(values) != len(FC_FIELDS):
            raise FontError("Fontconfig returned an unrecognised metadata record")
        records.append(dict(zip(FC_FIELDS, values)))
    return records


def actual_arial(record):
    # Family aliases used by fc-match are not proof of the font's own identity.
    return any(family.strip().casefold() == "arial"
               for family in record["family"].split(","))


def face_of(record):
    try:
        # Arial's four static faces have numeric Fontconfig weight/slant fields.
        weight = float(record["weight"])
        slant = float(record["slant"])
    except ValueError:
        return None
    bold = weight >= 180
    italic = slant >= 90
    if not bold and not 40 <= weight <= 100:
        return None
    if italic and slant not in (100, 110):
        return None
    if not italic and slant != 0:
        return None
    return next(name for name, pair in FACES.items() if pair == (bold, italic))


def check_fonts(font_dir, out):
    out = Path(out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    report = {"state": "FAIL", "family": "Arial", "faces": {}, "errors": [],
              "fonts_copied": False, "font_directory": None,
              "scope": "Font identity, four actual faces and Fontconfig resolution"}
    for program in ("fc-scan", "fc-match", "fc-cache"):
        if shutil.which(program) is None:
            report["errors"].append("Required Fontconfig command is unavailable: " + program)
    config = out / "fonts.conf"
    env_file = out / "FONT_ENV.sh"
    report_file = out / "FONT_PREFLIGHT.json"
    # A failed rerun must not leave a usable environment from an earlier success.
    if env_file.exists():
        env_file.unlink()
    env = os.environ.copy()
    try:
        if report["errors"]:
            raise FontError("Install Fontconfig before running the Arial preflight")
        base_config = Path("/etc/fonts/fonts.conf")
        if not base_config.is_file():
            raise FontError("System Fontconfig configuration is missing: /etc/fonts/fonts.conf")
        root = ET.Element("fontconfig")
        ET.SubElement(root, "include", {"ignore_missing": "no"}).text = str(base_config)
        if font_dir:
            directory = Path(font_dir).expanduser().resolve(strict=True)
            if not directory.is_dir():
                raise FontError("--font-dir must be a directory containing genuine Arial files")
            report["font_directory"] = str(directory)
            candidates = sorted(p for p in directory.rglob("*")
                                if p.is_file() and p.suffix.lower() in (".ttf", ".otf", ".ttc"))
            if not candidates:
                raise FontError("No TTF, OTF or TTC font files found in " + str(directory))
            inventory = []
            for candidate in candidates:
                inventory.extend(fc_records(run(["fc-scan", "--format", FC_FORMAT, str(candidate)])))
            report["supplied_font_metadata"] = inventory
            present = {face_of(record) for record in inventory if actual_arial(record)}
            missing = sorted(set(FACES) - present)
            if missing:
                observed = sorted({record["family"] for record in inventory})
                raise FontError("The supplied directory does not contain all four genuine Arial faces. "
                                f"Missing: {', '.join(missing)}. Actual families: {', '.join(observed)}. "
                                "Use licensed Arial regular/bold/italic/bolditalic files; aliases and "
                                "renamed Nimbus/Liberation fonts are not accepted.")
            ET.SubElement(root, "dir").text = str(directory)
        ET.SubElement(root, "cachedir").text = str(out / "fontconfig-cache")
        config.write_text('<?xml version="1.0"?>\n<!DOCTYPE fontconfig SYSTEM "urn:fontconfig:fonts.dtd">\n'
                          + ET.tostring(root, encoding="unicode") + "\n")
        env["FONTCONFIG_FILE"] = str(config)
        # Do not inherit a possibly unrelated user FONTCONFIG_PATH override.
        env.pop("FONTCONFIG_PATH", None)
        run(["fc-cache", "-f"], env=env)
        for face, (bold, italic) in FACES.items():
            pattern = f"Arial:weight={'bold' if bold else 'regular'}:slant={'italic' if italic else 'roman'}"
            records = fc_records(run(["fc-match", "--format", FC_FORMAT, pattern], env=env))
            if len(records) != 1:
                raise FontError(f"Expected one Fontconfig match for Arial {face}")
            record = records[0]
            report["faces"][face] = record
            if not actual_arial(record) or face_of(record) != face:
                report["errors"].append(
                    f"Arial {face} resolves to {record['family']} / {record['style']} "
                    f"({record['file']}); a genuine Arial {face} face is required.")
                continue
            path = Path(record["file"])
            if not path.is_file():
                report["errors"].append(f"Resolved font file is missing: {path}")
                continue
            # Independently inspect the file instead of accepting match aliases.
            scans = fc_records(run(["fc-scan", "--format", FC_FORMAT, str(path)], env=env))
            selected = [item for item in scans if item["index"] == record["index"]]
            if len(selected) != 1 or not actual_arial(selected[0]) or face_of(selected[0]) != face:
                report["errors"].append(f"Actual file metadata does not confirm Arial {face}: {path}")
                continue
            record["sha256"] = sha256(path)
        locations = {(item["file"], item["index"]) for item in report["faces"].values()}
        if len(locations) != 4:
            report["errors"].append("Four distinct actual font faces were not resolved; synthetic faces are rejected")
        if not report["errors"]:
            report["state"] = "PASS"
            report["environment"] = {"PBP_FIGURE_FONT": "Arial", "FONTCONFIG_FILE": str(config)}
            env_file.write_text("# Generated only after genuine Arial preflight passes.\n"
                                "unset FONTCONFIG_PATH\nexport PBP_FIGURE_FONT=Arial\n"
                                + "export FONTCONFIG_FILE=" + shlex.quote(str(config)) + "\n")
            report["environment_file"] = str(env_file)
    except (FontError, OSError) as exc:
        report["errors"].append(str(exc))
    report["fontconfig_file"] = str(config) if config.exists() else None
    save_json(report_file, report)
    return report


def pdf_fonts(path):
    if shutil.which("pdffonts"):
        output = run(["pdffonts", str(path)])
        lines = output.splitlines()
        header = next((line for line in lines if line.startswith("name ")), None)
        if header is None:
            raise FontError("pdffonts returned no recognisable header")
        starts = [header.index(name) for name in ("name", "type", "encoding", "emb", "sub", "uni", "object ID")]
        records = []
        for line in lines[lines.index(header) + 2:]:
            if not line.strip():
                continue
            values = [line[start:end].strip() for start, end in zip(starts, starts[1:] + [len(line)])]
            records.append({"name": values[0], "type": values[1], "embedded": values[3] == "yes"})
        return "pdffonts", records
    try:
        import fitz
    except ImportError as exc:
        raise FontError("PDF audit requires Poppler pdffonts or PyMuPDF; neither is available") from exc
    records = {}
    with fitz.open(path) as document:
        for page in document:
            for item in page.get_fonts(full=True):
                xref, _extension, kind, name = item[:4]
                embedded = bool(document.extract_font(xref)[3]) if xref else False
                records[(xref, name)] = {"name": name, "type": kind, "embedded": embedded}
    return "PyMuPDF", list(records.values())


def pdf_text_count(path):
    if shutil.which("pdftotext"):
        return len(run(["pdftotext", str(path), "-"]).strip())
    try:
        import fitz
    except ImportError as exc:
        raise FontError("PDF editable-text audit requires pdftotext or PyMuPDF") from exc
    with fitz.open(path) as document:
        return sum(len(page.get_text().strip()) for page in document)


def audit_pdf(path):
    engine, fonts = pdf_fonts(path)
    errors = []
    if not fonts:
        errors.append("PDF contains no inspectable fonts")
    for font in fonts:
        base = re.sub(r"^[A-Z]{6}\+", "", font["name"])
        font["base_name"] = base
        if base not in PDF_ARIAL_NAMES:
            errors.append("Non-Arial PDF font: " + font["name"])
        if re.sub(r"\s+", "", font["type"]).casefold() == "type3":
            errors.append("Type 3 fonts are not accepted: " + font["name"])
        if not font["embedded"]:
            errors.append("PDF font is not embedded: " + font["name"])
    text_count = pdf_text_count(path)
    if not text_count:
        errors.append("PDF contains no extractable editable text")
    return {"file": path.name, "format": "PDF", "engine": engine, "fonts": fonts,
            "text_characters": text_count, "errors": errors, "state": "FAIL" if errors else "PASS"}


def declarations(style):
    return {key.strip().lower(): value.strip().replace("!important", "").strip()
            for key, value in re.findall(r"([\w-]+)\s*:\s*([^;}]+)", style)}


def arial_family(value):
    families = [part.strip().strip("\"'").casefold() for part in value.split(",")]
    return families == ["arial"]


def local_name(element):
    return element.tag.rsplit("}", 1)[-1]


def selector_matches(selector, element, parents):
    # Current svglite exports use inline styles. Support common simple class/id/
    # descendant selectors too; fail unresolved text rather than guessing a font.
    if re.search(r"[:>+~\[]", selector):
        return False
    parts = selector.split()
    if not parts:
        return False
    def matches(token, node):
        match = re.fullmatch(r"([\w*-]+)?((?:[.#][\w-]+)*)", token)
        if not match:
            return False
        tag, constraints = match.groups()
        if tag and tag != "*" and tag != local_name(node):
            return False
        for kind, value in re.findall(r"([.#])([\w-]+)", constraints):
            if kind == "#" and node.get("id") != value:
                return False
            if kind == "." and value not in node.get("class", "").split():
                return False
        return True
    if not matches(parts[-1], element):
        return False
    current = parents.get(element)
    for part in reversed(parts[:-1]):
        while current is not None and not matches(part, current):
            current = parents.get(current)
        if current is None:
            return False
        current = parents.get(current)
    return True


def audit_svg(path):
    root = ET.parse(path).getroot()
    parents = {child: parent for parent in root.iter() for child in parent}
    errors, families, rules = [], set(), []
    for node in root.iter():
        styles = [node.get("style", "")]
        if local_name(node) == "style":
            css = "".join(node.itertext())
            if re.search(r"@font-face\b", css, re.I):
                errors.append("Embedded @font-face rules require an independent font-identity audit")
            for selector, body in re.findall(r"([^{}]+)\{([^{}]*)\}", css):
                props = declarations(body)
                for item in selector.split(","):
                    rules.append((item.strip(), props))
                styles.append(body)
        for style in styles:
            props = declarations(style)
            if "font" in props:
                errors.append("SVG font shorthand is not accepted; declare font-family explicitly")
            if "font-family" in props:
                families.add(props["font-family"])
                if not arial_family(props["font-family"]):
                    errors.append("Non-Arial SVG font-family: " + props["font-family"])
        if node.get("font-family"):
            families.add(node.get("font-family"))
            if not arial_family(node.get("font-family")):
                errors.append("Non-Arial SVG font-family: " + node.get("font-family"))
    def effective_family(node):
        value = node.get("font-family")
        matched = [(selector, props) for selector, props in rules
                   if "font-family" in props and selector_matches(selector, node, parents)]
        matched.sort(key=lambda item: (item[0].count("#"), item[0].count("."),
                                      len(re.findall(r"(?:^|\s)[A-Za-z]", item[0]))))
        for _selector, props in matched:
            value = props["font-family"]
        value = declarations(node.get("style", "")).get("font-family", value)
        if value is None and node in parents:
            return effective_family(parents[node])
        return value
    text_nodes = [node for node in root.iter() if local_name(node) == "text"
                  and "".join(node.itertext()).strip()]
    if not text_nodes:
        errors.append("SVG contains no editable text elements")
    checked = 0
    for text in text_nodes:
        for node in text.iter():
            if local_name(node) not in ("text", "tspan", "textPath"):
                continue
            family = effective_family(node)
            if family is None:
                errors.append("SVG text has no explicitly resolved Arial family: "
                              + "".join(node.itertext()).strip()[:60])
            elif not arial_family(family):
                errors.append("SVG text resolves to non-Arial family: " + family)
            checked += 1
    return {"file": path.name, "format": "SVG", "font_families": sorted(families),
            "text_elements": len(text_nodes), "checked_text_spans": checked,
            "errors": sorted(set(errors)), "state": "FAIL" if errors else "PASS"}


def verify_figures(figures, out, expected_count=9):
    figures = Path(figures).resolve()
    report = {"state": "FAIL", "figures_directory": str(figures), "expected_figures": expected_count,
              "files": [], "errors": [],
              "scope": "Actual embedded PDF fonts and editable SVG text; raster pixels are not font-auditable"}
    if not figures.is_dir():
        report["errors"].append("Figures directory is missing: " + str(figures))
    else:
        pdfs, svgs = sorted(figures.glob("*.pdf")), sorted(figures.glob("*.svg"))
        for kind, files in (("PDF", pdfs), ("SVG", svgs)):
            if len(files) != expected_count:
                report["errors"].append(f"Expected {expected_count} {kind} figures, found {len(files)}")
        if {path.stem for path in pdfs} != {path.stem for path in svgs}:
            report["errors"].append("PDF and SVG figure names do not match")
        for path, function in [(p, audit_pdf) for p in pdfs] + [(p, audit_svg) for p in svgs]:
            try:
                item = function(path)
            except (FontError, OSError, ET.ParseError, ValueError) as exc:
                item = {"file": path.name, "state": "FAIL", "errors": [str(exc)]}
            item["sha256"] = sha256(path)
            report["files"].append(item)
            report["errors"].extend(path.name + ": " + error for error in item["errors"])
        if not report["errors"]:
            report["state"] = "PASS"
    save_json(out, report)
    return report


# API used by the figure export workflow.
def check(font_dir, out):
    return check_fonts(font_dir, out)


def verify(figures, out, expected_count=9):
    return verify_figures(figures, out, expected_count)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    check = commands.add_parser("check", help="Require genuine Arial regular/bold/italic/bolditalic")
    check.add_argument("--font-dir", type=Path, help="Directory with user-licensed Arial font files")
    check.add_argument("--out", type=Path, required=True, help="Font validation/configuration directory")
    verify = commands.add_parser("verify", help="Audit PDF embedding and editable SVG Arial text")
    verify.add_argument("--figures", type=Path, required=True)
    verify.add_argument("--out", type=Path, required=True, help="Output FONT_AUDIT.json path")
    verify.add_argument("--expected-count", type=int, default=9)
    args = parser.parse_args()
    if args.command == "check":
        report = check_fonts(args.font_dir, args.out)
    else:
        if args.expected_count < 1:
            parser.error("--expected-count must be positive")
        report = verify_figures(args.figures, args.out, args.expected_count)
    print(json.dumps({"state": report["state"], "errors": report["errors"],
                      **({"environment": report["environment"]} if "environment" in report else {})}, indent=2))
    return 0 if report["state"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
