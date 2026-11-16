"""Join explicit author review to consensus; no inferred synthesis or toxicity decision."""
import argparse
from pathlib import Path
from consensus_prediction import read_table, write_table, audit


def select(consensus_rows, reviews):
    required = {"peptide_id", "decision", "reason"}
    if not all(required <= r.keys() for r in reviews):
        raise ValueError("Review requires peptide_id,decision,reason")
    ids = [r["peptide_id"] for r in consensus_rows]
    if len(set(ids)) != len(ids):
        raise ValueError("Duplicate consensus peptide IDs")
    review = {}
    for r in reviews:
        pid = r["peptide_id"]
        if pid in review or pid not in ids:
            raise ValueError("Duplicate or unknown review ID")
        if r["decision"] not in {"include", "exclude", "pending"} or not r["reason"].strip():
            raise ValueError("Each review needs include/exclude/pending and a reason")
        review[pid] = r
    result = []
    for r in consensus_rows:
        if r.get("consensus_pass") not in {"0", "1"}:
            raise ValueError("Invalid consensus_pass")
        decision = review.get(r["peptide_id"], {"decision": "pending", "reason": "Author review required"})
        if decision["decision"] == "include" and r["consensus_pass"] != "1":
            raise ValueError("Cannot include a peptide failing the chosen consensus rule; change rule explicitly")
        result.append({**r, "review_decision": decision["decision"], "review_reason": decision["reason"],
                       "selected": int(r["consensus_pass"] == "1" and decision["decision"] == "include")})
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--review", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    try:
        outputs = [Path(args.output).resolve(), Path(str(args.output) + ".audit.json").resolve()]
        if any(Path(p).resolve() in outputs for p in [args.input, args.review]):
            raise ValueError("Outputs must not overwrite an input")
        result = select(read_table(args.input), read_table(args.review))
        write_table(args.output, result)
        audit(args.output, [args.input, args.review], {"selection": "consensus plus explicit author review"})
    except (ValueError, OSError, KeyError) as exc:
        parser.exit(2, f"Error: {exc}\n")


if __name__ == "__main__":
    main()
