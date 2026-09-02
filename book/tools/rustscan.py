"""Minimal Rust source scanner used to locate `fn` items for tangling.

Handles exactly the syntax present in this codebase: line comments (//),
string literals ("..."), simple/escaped char literals ('x', '\\n', b'x'),
and loop labels/lifetimes ('main:, 'a) which look like char literals but
are not. No block comments or raw strings are used in this project.
"""
import re

FN_LINE_RE = re.compile(
    r'^(?P<indent> *)(?P<prefix>(?:pub(?:\([^)]*\))?\s+)?(?:async\s+)?fn\s+)(?P<name>[A-Za-z_][A-Za-z0-9_]*)'
)
IMPL_LINE_RE = re.compile(
    r'^impl(?:<[^>]*>)?\s+(?:[A-Za-z_][A-Za-z0-9_:]*\s+for\s+)?(?P<type>[A-Za-z_][A-Za-z0-9_]*)'
)


def _try_skip_token(text, i):
    """If text[i] starts a comment/string/char-literal, return the index
    just past it; otherwise return None."""
    n = len(text)
    c = text[i]
    if c == '/' and i + 1 < n and text[i + 1] == '/':
        j = text.find('\n', i)
        return n if j == -1 else j
    if c == '"':
        j = i + 1
        while j < n:
            if text[j] == '\\':
                j += 2
                continue
            if text[j] == '"':
                return j + 1
            j += 1
        return n
    if c == "'":
        if i + 1 < n and text[i + 1] == '\\':
            j = i + 2
            if j < n and text[j] == 'u' and j + 1 < n and text[j + 1] == '{':
                k = text.find('}', j)
                j = (k + 1) if k != -1 else (j + 1)
            else:
                j += 1
            if j < n and text[j] == "'":
                return j + 1
            return i + 1  # not a char literal after all; skip just the quote
        if i + 2 < n and text[i + 2] == "'":
            return i + 3
        return i + 1  # lifetime/label apostrophe: skip just the quote
    return None


def _find_body_span(text, from_idx):
    """Given an index at or before a fn's `{`, return (body_start, body_end)
    indices (body_end = index of the matching closing brace, inclusive)."""
    n = len(text)
    i = from_idx
    while i < n:
        skip_to = _try_skip_token(text, i)
        if skip_to is not None:
            i = skip_to
            continue
        if text[i] == '{':
            break
        i += 1
    else:
        raise RuntimeError('no function body found')
    body_start = i
    depth = 1
    i += 1
    while i < n and depth > 0:
        skip_to = _try_skip_token(text, i)
        if skip_to is not None:
            i = skip_to
            continue
        if text[i] == '{':
            depth += 1
        elif text[i] == '}':
            depth -= 1
        i += 1
    body_end = i - 1
    return body_start, body_end


def find_impl_blocks(text):
    """Return a list of (start, end, type_name) char spans for `impl` blocks
    at column 0 (end is inclusive index of the closing brace)."""
    blocks = []
    for m in re.finditer(r'^impl.*\{', text, re.MULTILINE):
        im = IMPL_LINE_RE.match(m.group(0))
        type_name = im.group('type') if im else '?'
        brace_idx = m.end() - 1
        _, body_end = _find_body_span(text, brace_idx)
        blocks.append((m.start(), body_end, type_name))
    return blocks


def find_functions(text):
    """Return a list of dicts describing every `fn` item in the file, in
    source order: name, qualified (impl-qualified) name, indent, line_start
    (char offset of the start of the fn's line), end (char offset of the
    closing brace, inclusive)."""
    impls = find_impl_blocks(text)
    results = []
    line_start = 0
    for line in text.splitlines(keepends=True):
        m = FN_LINE_RE.match(line)
        if m:
            fn_col_in_line = m.start('name')
            name_idx = line_start + fn_col_in_line
            indent = len(m.group('indent'))
            fn_kw_idx = line_start + len(m.group('indent'))
            body_start, body_end = _find_body_span(text, fn_kw_idx)
            qualifier = None
            for (s, e, type_name) in impls:
                if s < line_start < e:
                    qualifier = type_name
                    break
            results.append({
                'name': m.group('name'),
                'qualifier': qualifier,
                'indent': indent,
                'start': line_start,
                'end': body_end,
            })
        line_start += len(line)
    return results
