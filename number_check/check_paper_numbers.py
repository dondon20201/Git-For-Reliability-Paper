"""Check that every number stated in the paper exists in a committed results file.

Usage:
    python check_paper_numbers.py PAPER.tex RESULTS_DIR [RESULTS_DIR ...]

PAPER.tex is the LaTeX source of the paper (main.tex) or the supplement
(supplement.tex). Each RESULTS_DIR is searched recursively for CSV and JSON
files; every numeric value in them goes into a pool. Typical use, from the
repository root:

    python number_check/check_paper_numbers.py main.tex expected data

A number in the paper passes when some value in the pool rounds to it at the
precision the paper prints (".54" matches 0.5446; "27" matches 27.31). Years,
citation years and the numbers listed in ALLOWED below are skipped. The script
prints every number that has no match and exits with status 1 if there is any.

What this does and does not show: a pass means each printed number could have
come from a committed results file. It does not show that the number came from
the RIGHT file or the right row. That link is given for the main results by
the claim-status table and by analysis_manifest.csv.
"""
import csv, json, os, re, sys

# Numbers in the text that are not results: interval levels, the reliability
# benchmarks, counts that are part of a name, and thousands separators.
ALLOWED = {"95", "90", ".7", ".70", ".8", ".80", ".5", ".50", "0.7", "0.8", "18", "13"}


def pool_from(dirs):
    vals = []
    def add(x):
        try:
            v = float(x)
            if v == v and abs(v) != float("inf"):
                vals.append(abs(v))
        except (TypeError, ValueError):
            pass
    def walk_json(o):
        if isinstance(o, dict):
            for v in o.values(): walk_json(v)
        elif isinstance(o, list):
            for v in o: walk_json(v)
        else:
            add(o)
    for d in dirs:
        for root, _, files in os.walk(d):
            for f in files:
                p = os.path.join(root, f)
                if f.endswith(".csv"):
                    with open(p, newline="", encoding="utf-8") as fh:
                        for row in csv.reader(fh):
                            for cell in row:
                                add(cell)
                                for m in re.findall(r"-?\d*\.?\d+", cell):
                                    add(m)
                elif f.endswith(".json"):
                    with open(p, encoding="utf-8") as fh:
                        walk_json(json.load(fh))
    return sorted(set(vals))


def paper_numbers(path):
    tex = open(path, encoding="utf-8").read()
    if "\\begin{document}" in tex:
        tex = tex.split("\\begin{document}", 1)[1]                 # skip the preamble
    if "\\begin{abstract}" in tex:
        tex = tex.split("\\begin{abstract}", 1)[1]
    tex = re.split(r"\\section\*?\{References\}", tex)[0]
    tex = re.sub(r"(?<!\\)%.*", "", tex)                             # comments
    tex = re.sub(r"\\(renewcommand|setlength|addtolength|setstretch)\b.*", "", tex)  # layout settings
    tex = re.sub(r"\\begin\{(longtable|tabular\*?|tabularx)\}(\{[^{}]*\})?\{([^{}]|\{[^{}]*\})*\}", " ", tex)  # column specs
    tex = re.sub(r"\d*\.?\d+\s*(in|pt|cm|mm|em|ex)\b", " ", tex)      # lengths
    tex = re.sub(r"\\(ref|label|cite\w*|includegraphics|url|href)(\[[^\]]*\])?\{[^}]*\}", " ", tex)
    tex = re.sub(r"\\[a-zA-Z]+\*?(\[[^\]]*\])?", " ", tex)          # commands
    tex = re.sub(r"(\d),(\d{3})", r"\1\2", tex)                     # 510,750 -> 510750
    out = []
    for i, line in enumerate(tex.split("\n"), 1):
        for m in re.finditer(r"(?<![\w.])-?\d*\.?\d+", line):
            out.append((m.group(0).lstrip("-"), line.strip()))
    return out


def matches(n, pool):
    x = float(n)
    d = len(n.split(".")[1]) if "." in n else 0
    tol = 0.5 * 10 ** -d + 1e-12
    import bisect
    i = bisect.bisect_left(pool, x - tol)
    return i < len(pool) and pool[i] <= x + tol


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    pool = pool_from(sys.argv[2:])
    nums = paper_numbers(sys.argv[1])
    missing = []
    for n, ctx in nums:
        if n in ALLOWED or re.fullmatch(r"(19|20)\d\d", n):
            continue
        if not matches(n, pool):
            missing.append((n, ctx))
    checked = len(nums)
    print(f"{checked} numbers read from {sys.argv[1]}; {len(pool)} values in the results pool")
    if missing:
        print(f"{len(missing)} number(s) with no match in any results file:")
        for n, ctx in missing:
            print(f"  {n:>10}  ...{ctx[:110]}")
        sys.exit(1)
    print("Every number in the paper matches a value in a committed results file.")


if __name__ == "__main__":
    main()
