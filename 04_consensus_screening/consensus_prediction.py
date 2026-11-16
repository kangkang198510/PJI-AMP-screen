"""Count one vote per tool family; consume predictions, do not run AMP models."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import platform


def read_table(path):
    with open(path, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError("Missing or duplicate CSV headers")
        rows = list(reader)
        if any(None in r or any(v is None for v in r.values()) for r in rows):
            raise ValueError("CSV row width differs from header")
    if not rows:
        raise ValueError("Empty input")
    return rows


def load_mapping(path, key, value):
    rows = read_table(path)
    result = {}
    for row in rows:
        if key not in row or value not in row or not row[key].strip() or not row[value].strip():
            raise ValueError(f"Mapping needs nonempty {key},{value}")
        if row[key] in result:
            raise ValueError(f"Duplicate mapping key: {row[key]}")
        result[row[key]] = row[value]
    return result


def consensus(rows, models, min_votes, family_map, family_rules=None):
    if not models or len(models) != len(set(models)):
        raise ValueError("Models must be an explicit, unique list")
    if any(m not in family_map or not family_map[m].strip() for m in models):
        raise ValueError("Every selected model needs an explicit tool-family mapping")
    families = list(dict.fromkeys(family_map[m] for m in models))
    if not 1 <= min_votes <= len(families):
        raise ValueError("min-votes must be between 1 and the selected FAMILY count")
    family_rules = family_rules or {}
    members = {f: [m for m in models if family_map[m] == f] for f in families}
    for f, group in members.items():
        if len(group) > 1 and family_rules.get(f) not in {"all", "any", "majority"}:
            raise ValueError(f"Multiple outputs in {f}: explicit all/any/majority rule required; no default")
    required = {"peptide_id", "sequence", "model", "prediction"}
    if not all(required <= r.keys() for r in rows):
        raise ValueError("Required columns: peptide_id,sequence,model,prediction")
    records = {}
    seen = set()
    for r in rows:
        if r["model"] not in models:
            continue
        pid, seq, model = r["peptide_id"].strip(), r["sequence"].strip().upper(), r["model"]
        if not pid or not seq or any(a not in "ACDEFGHIKLMNPQRSTVWY" for a in seq):
            raise ValueError("Empty ID or noncanonical amino-acid sequence; review modifications explicitly")
        if r["prediction"] not in {"0", "1"}:
            raise ValueError("prediction must be explicitly supplied as 0 or 1, never missing")
        if (pid, model) in seen:
            raise ValueError(f"Duplicate peptide/model pair: {pid}, {model}")
        seen.add((pid, model))
        item = records.setdefault(pid, {"sequence": seq, "votes": {}})
        if item["sequence"] != seq:
            raise ValueError(f"Inconsistent sequence for {pid}")
        item["votes"][model] = int(r["prediction"])
    if not records:
        raise ValueError("No predictions match the selected models")
    result = []
    for pid, item in sorted(records.items()):
        if set(item["votes"]) != set(models):
            raise ValueError(f"Incomplete model coverage for {pid}; missing is not negative")
        positives = [m for m in models if item["votes"][m]]
        positive_families = []
        for f, group in members.items():
            k = sum(item["votes"][m] for m in group)
            rule = family_rules.get(f)
            passed = (k == 1) if len(group) == 1 else (
                k == len(group) if rule == "all" else k > 0 if rule == "any" else k > len(group) / 2)
            if passed:
                positive_families.append(f)
        result.append(dict(peptide_id=pid, sequence=item["sequence"], n_models=len(models),
                           n_families=len(families), votes=len(positive_families),
                           positive_models=";".join(positives), positive_families=";".join(positive_families),
                           consensus_pass=int(len(positive_families) >= min_votes)))
    return sorted(result, key=lambda r: (-r["votes"], r["peptide_id"]))


def write_table(path, rows):
    path = Path(path)
    if path.exists() or Path(str(path) + ".audit.json").exists():
        raise ValueError("Output or audit already exists; use a new output path")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def audit(output, inputs, settings):
    info = {"python": platform.python_version(), "seed": None,
            "randomness": "none", "settings": settings,
            "inputs": [{"name": Path(p).name, "sha256": hashlib.sha256(Path(p).read_bytes()).hexdigest()}
                       for p in inputs]}
    with Path(str(output) + ".audit.json").open("x", encoding="utf-8") as handle:
        json.dump(info, handle, indent=2)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--models", nargs="+", required=True)
    parser.add_argument("--family-map", required=True, help="CSV model,family; one vote per family")
    parser.add_argument("--family-rules", help="Optional CSV family,rule (all/any/majority); required when using multiple variants in a family")
    parser.add_argument("--min-votes", type=int, required=True,
                        help="Minimum positive TOOL FAMILIES; explicit reconstruction setting")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    try:
        mapping = load_mapping(args.family_map, "model", "family")
        rules = load_mapping(args.family_rules, "family", "rule") if args.family_rules else {}
        inputs = [args.input, args.family_map] + ([args.family_rules] if args.family_rules else [])
        outputs = [Path(args.output).resolve(), Path(str(args.output) + ".audit.json").resolve()]
        if any(Path(p).resolve() in outputs for p in inputs):
            raise ValueError("Outputs must not overwrite an input")
        result = consensus(read_table(args.input), args.models, args.min_votes, mapping, rules)
        write_table(args.output, result)
        audit(args.output, inputs, {"models": args.models, "min_family_votes": args.min_votes,
                                   "family_map": mapping, "family_rules": rules, "voting_unit": "tool_family"})
    except (ValueError, OSError) as exc:
        parser.exit(2, f"Error: {exc}\n")


if __name__ == "__main__":
    main()
