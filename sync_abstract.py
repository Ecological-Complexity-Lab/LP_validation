#!/usr/bin/env python3
"""Copy the paper's abstract from main.tex into README.md and docs/index.html.

The abstract has been rewritten several times, and each copy drifted
separately. This makes the sync a command rather than a manual pass.

    python3 sync_abstract.py            # rewrite both copies
    python3 sync_abstract.py --check    # report drift, change nothing, exit 1 if stale

The landing page drops the sentence pointing readers to the online guide,
because it is that guide. The README keeps it, with the URL as a link.
"""
import argparse, re, sys, textwrap
from pathlib import Path

TEX = Path("/Users/shai/Library/CloudStorage/Dropbox/Apps/Overleaf/"
           "LP evaluation and validation/main.tex")
README, INDEX = Path("README.md"), Path("docs/index.html")
GUIDE_SENTENCE = r"\s*We provide a worked empirical analysis and an online guide at [^.]*\."
GUIDE_URL = "http://lpguide.ecomplab.com"


def balanced(s, start):
    """Return the contents of the {...} group whose opening brace is at `start`."""
    depth = 0
    for i in range(start, len(s)):
        if s[i] == "{":
            depth += 1
        elif s[i] == "}":
            depth -= 1
            if depth == 0:
                return s[start + 1:i]
    raise ValueError("unbalanced braces in the abstract")


def read_abstract(tex_path):
    tex = tex_path.read_text(encoding="utf-8")
    m = re.search(r"\\section\*\{Abstract\}\s*\n(.*?)\n\s*\n", tex, re.S)
    if not m:
        sys.exit("could not find \\section*{Abstract} in " + str(tex_path))
    a = m.group(1).strip()
    if a.startswith(r"\textbf{"):
        a = balanced(a, len(r"\textbf"))
    a = re.sub(r"\\(?:emph|textit|textbf)\{([^{}]*)\}", r"\1", a)
    a = re.sub(r"\\url\{([^}]*)\}", r"\1", a)
    a = re.sub(r"\\cite[tp]?\{[^}]*\}", "", a)
    a = a.replace(r"\%", "%").replace("~", " ")
    return re.sub(r"\s+", " ", a).strip()


def replace_span(text, pattern, group, new, label):
    m = re.search(pattern, text, re.S)
    if not m:
        sys.exit(f"could not locate the abstract in {label}")
    return text[:m.start(group)] + new + text[m.end(group):], m.group(group)


def norm(t):
    t = re.sub(r"<[^>]+>", " ", t)
    t = re.sub(r"\[([^\]]*)\]\([^)]*\)", r"\1", t)      # unwrap markdown links
    return re.sub(r"\s+", " ", t).strip()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tex", type=Path, default=TEX)
    ap.add_argument("--check", action="store_true",
                    help="report drift and exit 1 if stale; write nothing")
    args = ap.parse_args()

    if not args.tex.is_file():
        sys.exit(f"cannot read {args.tex}\n"
                 "If it lives in Dropbox/Overleaf, make sure the folder is available "
                 "locally, or pass --tex with another path.")

    abstract = read_abstract(args.tex)
    for_site = re.sub(GUIDE_SENTENCE, "", abstract)
    for_readme = abstract.replace(GUIDE_URL, f"[{GUIDE_URL}]({GUIDE_URL})")

    targets = [
        (README, r"(# Abstract:\s*\n)(.*?)(\n\s*\n)", 2, for_readme, abstract),
        (INDEX, r"(<div class=\"fold-body\">\s*<p>)(.*?)(</p>)", 2,
         textwrap.fill(for_site, 110, initial_indent="      ",
                       subsequent_indent="      ").strip(), for_site),
    ]

    stale = []
    for path, pat, grp, new, expect in targets:
        if not path.is_file():
            sys.exit(f"missing {path} (run from the repository root)")
        src = path.read_text(encoding="utf-8")
        updated, current = replace_span(src, pat, grp, new, str(path))
        if norm(current) == norm(expect):
            print(f"  {path}: up to date")
            continue
        stale.append(path)
        if args.check:
            print(f"  {path}: STALE")
        else:
            path.write_text(updated, encoding="utf-8")
            print(f"  {path}: updated")

    if args.check and stale:
        print(f"\n{len(stale)} file(s) stale. Run without --check to update.")
        return 1
    if not stale:
        print("\nBoth copies match the paper.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
