#!/usr/bin/env python3
"""Tangle: write every Quarto chapter's Rust code cell into its target .rs
file, at the location marked by `// @tangle:start <id>` / `// @tangle:end
<id>`. This is the ongoing sync direction (qmd chapters are authoritative);
run automatically as a Quarto `pre-render` hook, and can be run standalone:

    python3 book/tools/tangle.py

It never touches anything outside the marked regions, so hand-written prose
and file scaffolding (imports, struct/impl declarations, attributes) are
left exactly as they are in the .rs files.
"""
import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common


def extract_tangle_blocks(qmd_text):
    """Return a list of (tangle_id, filename, code) for every tangled fenced
    code block in a chapter file."""
    blocks = []
    lines = qmd_text.splitlines()
    i = 0
    while i < len(lines):
        m = common.FENCE_OPEN_RE.match(lines[i])
        if m:
            attrs = dict(common.ATTR_RE.findall(m.group("attrs")))
            tid = attrs.get("tangle")
            filename = attrs.get("filename")
            i += 1
            body = []
            while i < len(lines) and lines[i] != "```":
                body.append(lines[i])
                i += 1
            if tid and filename:
                blocks.append((tid, filename, "\n".join(body)))
        i += 1
    return blocks


def collect_all_blocks():
    by_file = {}
    for qmd_path in sorted(glob.glob(os.path.join(common.CHAPTERS_DIR, "**", "*.qmd"), recursive=True)):
        text = open(qmd_path, encoding="utf-8").read()
        for tid, filename, code in extract_tangle_blocks(text):
            by_file.setdefault(filename, {})[tid] = (code, qmd_path)
    return by_file


def apply_to_file(relpath, blocks):
    path = os.path.join(common.REPO_DIR, relpath)
    text = open(path, encoding="utf-8").read()
    lines = text.split("\n")
    out = []
    i = 0
    seen = set()
    while i < len(lines):
        m = common.TANGLE_START_RE.match(lines[i])
        if m:
            indent, tid = m.group(1), m.group(2)
            out.append(lines[i])
            i += 1
            # skip existing content up to (not including) the end marker
            while i < len(lines) and not common.TANGLE_END_RE.match(lines[i]):
                i += 1
            if tid in blocks:
                code, qmd_path = blocks[tid]
                out.extend(code.split("\n"))
                seen.add(tid)
            else:
                sys.stderr.write(f"warning: no chapter provides tangle id {tid!r} (in {relpath}); leaving unchanged\n")
            continue
        out.append(lines[i])
        i += 1
    new_text = "\n".join(out)
    changed = new_text != text
    if changed:
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(new_text)
    missing = set(blocks) - seen
    for tid in missing:
        _, qmd_path = blocks[tid]
        sys.stderr.write(
            f"warning: chapter {qmd_path} targets tangle id {tid!r} but no such marker exists in {relpath}\n"
        )
    return changed


def main():
    by_file = collect_all_blocks()
    any_changed = False
    for relpath, blocks in sorted(by_file.items()):
        changed = apply_to_file(relpath, blocks)
        any_changed = any_changed or changed
        print(f"tangle: {relpath}: {len(blocks)} functions{' (updated)' if changed else ' (unchanged)'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
