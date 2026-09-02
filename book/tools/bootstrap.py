#!/usr/bin/env python3
"""One-time (re-runnable) bootstrap: scans the current, verified .rs sources,
wraps every `fn` item in `// @tangle:start/end <id>` markers, and generates
one Quarto chapter per function GROUP (see common.CHAPTER_GROUPS) under
book/chapters/, seeded with the exact current source text. Related/small
functions (e.g. a struct's methods, a family of sort helpers) share one
chapter with one fenced code block per function.

After this runs, book/tools/tangle.py becomes the ongoing sync direction
(qmd -> rs). Safe to re-run: if a source file already has markers, it is
left untouched; existing chapter files have only their fenced code blocks
refreshed, leaving hand-edited prose untouched.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common
import rustscan


def insert_markers(path, modname):
    text = open(path, encoding="utf-8").read()
    if "@tangle:start" in text:
        print(f"  {path}: markers already present, skipping insertion")
        return text
    fns = rustscan.find_functions(text)
    # Insert from the end of the file backwards so earlier offsets stay valid.
    for f in reversed(fns):
        indent = " " * f["indent"]
        start_marker = f"{indent}// @tangle:start {common.tangle_id(modname, f['qualifier'], f['name'])}\n"
        end_marker = f"{indent}// @tangle:end {common.tangle_id(modname, f['qualifier'], f['name'])}\n"
        text = text[: f["end"] + 1] + "\n" + end_marker.rstrip("\n") + text[f["end"] + 1 :]
        text = text[: f["start"]] + start_marker + text[f["start"] :]
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)
    print(f"  {path}: inserted {len(fns)} marker pairs")
    return text


def render_chapter(title, filename, modname, entries, multi):
    """entries: list of (label, tid, code). `multi` controls whether each
    function gets its own '### label' subheading (only useful when the
    chapter bundles more than one function)."""
    parts = [f'---\ntitle: "{title}"\n---\n\n']
    parts.append(f"Part of `{filename}` (`{modname}` binary/module).\n\n")
    parts.append(
        "<!-- TODO: describe what "
        + ("these functions do" if multi else "this function does")
        + ", and, where relevant, which Fortran subroutine/statement "
        + ("they port" if multi else "it ports")
        + ". -->\n\n"
    )
    for label, tid, code in entries:
        if multi:
            parts.append(f"### `{label}`\n\n")
        parts.append(f'```{{.rust filename="{filename}" tangle="{tid}"}}\n{code}\n```\n\n')
    return "".join(parts).rstrip("\n") + "\n"


def replace_code_cell(existing_text, tid, filename, code):
    """Re-run support: replace just the fenced code block for this tangle id,
    keeping any hand-edited prose (including subheadings) around it. If no
    fence with this tid exists yet (e.g. a function newly added to a group),
    append a fresh one at the end."""
    lines = existing_text.splitlines(keepends=True)
    out = []
    i = 0
    replaced = False
    while i < len(lines):
        m = common.FENCE_OPEN_RE.match(lines[i].rstrip("\n"))
        if m:
            attrs = dict(common.ATTR_RE.findall(m.group("attrs")))
            if attrs.get("tangle") == tid:
                out.append(lines[i])
                out.append(code + "\n")
                i += 1
                while i < len(lines) and lines[i].rstrip("\n") != "```":
                    i += 1
                if i < len(lines):
                    out.append(lines[i])
                    i += 1
                replaced = True
                continue
        out.append(lines[i])
        i += 1
    if not replaced:
        out.append(
            f'\n```{{.rust filename="{filename}" tangle="{tid}"}}\n{code}\n```\n'
        )
    return "".join(out)


def write_group_chapter(modname, slug, title, funcs, filename, text_before, fn_index):
    path = common.chapter_path(modname, slug)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    multi = len(funcs) > 1
    entries = []
    for qualifier, name in funcs:
        f = fn_index[(qualifier, name)]
        code = text_before[f["start"] : f["end"] + 1]
        tid = common.tangle_id(modname, qualifier, name)
        label = f"{qualifier}::{name}" if qualifier else name
        entries.append((label, tid, code))

    if os.path.exists(path):
        print(f"    {path}: exists, leaving prose untouched, refreshing {len(entries)} code cell(s)")
        text = open(path, encoding="utf-8").read()
        for _label, tid, code in entries:
            text = replace_code_cell(text, tid, filename, code)
        new_text = text
    else:
        new_text = render_chapter(title, filename, modname, entries, multi)

    with open(path, "w", encoding="utf-8") as fh:
        fh.write(new_text)
    return path


def write_quarto_yml(modules):
    """Regenerate _quarto.yml's book/chapters list from the current set of
    generated chapters, in SOURCE_FILES order. Everything outside the
    `chapters:` list (title, format options, pre-render hook, ...) is
    hand-maintained and left untouched."""
    path = os.path.join(common.BOOK_DIR, "_quarto.yml")
    lines = ["  chapters:", "    - index.qmd"]
    for modname, relpath, part_title in common.SOURCE_FILES:
        lines.append(f'    - part: "{part_title}"')
        lines.append("      chapters:")
        for chapter_path in modules[modname]:
            rel = os.path.relpath(chapter_path, common.BOOK_DIR)
            lines.append(f"        - {rel}")
    new_block = "\n".join(lines) + "\n"

    if not os.path.exists(path):
        print(f"  {path}: does not exist yet, skipping (create it first)")
        return
    text = open(path, encoding="utf-8").read()
    start_marker = "  chapters:"
    idx = text.find(start_marker)
    if idx == -1:
        print(f"  {path}: no 'chapters:' key found, leaving untouched")
        return
    # The chapters list runs until the next top-level (non-indented, non-blank)
    # key, i.e. a line that doesn't start with whitespace.
    after = text[idx:]
    rest_lines = after.split("\n")
    end = 1
    while end < len(rest_lines):
        line = rest_lines[end]
        if line and not line[0].isspace():
            break
        end += 1
    tail = "\n".join(rest_lines[end:])
    new_text = text[:idx] + new_block + tail
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(new_text)
    print(f"  {path}: regenerated chapter list ({sum(len(v) for v in modules.values())} chapters)")


def main():
    modules = {}
    for modname, relpath, _part_title in common.SOURCE_FILES:
        path = os.path.join(common.REPO_DIR, relpath)
        print(f"Processing {relpath} ...")
        text_before = open(path, encoding="utf-8").read()
        fns_before = rustscan.find_functions(text_before)
        insert_markers(path, modname)

        fn_index = {(f["qualifier"], f["name"]): f for f in fns_before}
        groups = common.CHAPTER_GROUPS[modname]
        grouped_keys = {key for _slug, _title, funcs in groups for key in funcs}
        actual_keys = set(fn_index)
        if grouped_keys != actual_keys:
            missing = actual_keys - grouped_keys
            extra = grouped_keys - actual_keys
            raise SystemExit(
                f"common.CHAPTER_GROUPS['{modname}'] is out of sync with {relpath}: "
                f"missing from groups={sorted(missing)}, groups reference nonexistent fns={sorted(extra)}"
            )

        chapters = []
        for slug, title, funcs in groups:
            p = write_group_chapter(modname, slug, title, funcs, relpath, text_before, fn_index)
            chapters.append(p)
        modules[modname] = chapters

    write_quarto_yml(modules)

    # Round-trip self-check: markers + chapters should reconstruct byte-identical
    # source when tangled back (checked by the caller running tangle.py + diff).
    print("\nDone. Run tangle.py and diff against a backup to confirm a lossless round trip.")
    return modules


if __name__ == "__main__":
    main()
