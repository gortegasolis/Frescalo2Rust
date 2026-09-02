"""Shared constants/helpers for the bootstrap and tangle scripts."""
import os
import re

BOOK_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO_DIR = os.path.dirname(BOOK_DIR)
CHAPTERS_DIR = os.path.join(BOOK_DIR, "chapters")

# (module name, path relative to repo root, book part title)
# Order here is both processing order and book part order: shared support
# code first, then the three binaries in pipeline order (sampdist ->
# neighsim -> frescalo), matching the README.
SOURCE_FILES = [
    ("lib", "src/lib.rs", "src/lib.rs — shared support code"),
    ("sampdist", "src/bin/sampdist.rs", "src/bin/sampdist.rs — SAMPDIST"),
    ("neighsim", "src/bin/neighsim.rs", "src/bin/neighsim.rs — NEIGHSIM"),
    ("frescalo", "src/bin/frescalo.rs", "src/bin/frescalo.rs — FRESCALO"),
]

# Chapter groups: functions that are small, closely related, or naturally
# read together (e.g. a struct's methods, a family of sort/format helpers)
# share one chapter instead of getting one chapter each. Each entry is
# (slug, title, [(qualifier_or_None, name), ...]) in the order they appear
# in the chapter. The tangle markers/ids in the .rs files are unaffected --
# grouping is purely a book-organization concern.
CHAPTER_GROUPS = {
    "lib": [
        ("name_helpers", "Name helpers: `blank_name` and `name_to_string`",
         [(None, "blank_name"), (None, "name_to_string")]),
        ("rec30", "`Rec30`: fixed-width name/site record",
         [(None, "make_rec30"), (None, "rec30_field")]),
        ("data_reader", "`DataReader`: token-stream file reader",
         [("DataReader", "new"), ("DataReader", "rewind"), ("DataReader", "getd")]),
        ("getnum", "`getnum`: numeric token parsing",
         [(None, "getnum")]),
        ("binfnd", "`binfnd`: binary search",
         [(None, "binfnd")]),
        ("addwrd", "`addwrd` / `addwrd_guarded`: word-list insertion",
         [(None, "addwrd"), (None, "addwrd_guarded")]),
        ("sort_wrappers", "`sort30` / `sort_real`: name and real-array sort wrappers",
         [(None, "sort30"), (None, "sort_real")]),
        ("sort_helpers", "`isort` / `sort2`: index and parallel-array sort helpers",
         [(None, "isort"), (None, "sort2")]),
        ("fmt", "`ifmt` / `ffmt`: Fortran `Iw` / `Fw.d` formatting",
         [(None, "ifmt"), (None, "ffmt")]),
        ("rec_ctor", "`Rec`: constructors",
         [("Rec", "new"), ("Rec", "raw")]),
        ("rec_accessors", "`Rec`: field accessors and line writer",
         [("Rec", "s"), ("Rec", "name"), ("Rec", "x"), ("Rec", "i"), ("Rec", "f"), ("Rec", "writeln")]),
        ("cout", "`cout`: console output helper",
         [(None, "cout")]),
        ("list_directed_input", "`ld_i` / `ld_f` / `ld_line`: list-directed input helpers",
         [(None, "ld_i"), (None, "ld_f"), (None, "ld_line")]),
        ("stdin_fields", "`read_stdin_line` / `read_a20`: stdin line and name-field readers",
         [(None, "read_stdin_line"), (None, "read_a20")]),
        ("numeric_fields", "`read_f8_4` / `read_int_listdirected`: numeric field readers",
         [(None, "read_f8_4"), (None, "read_int_listdirected")]),
        ("file_io", "`filin` / `filout` / `hold`: file-open and modal-hold helpers",
         [(None, "filin"), (None, "filout"), (None, "hold")]),
        ("arr2_ctor", "`Arr2`: constructor and getter",
         [("Arr2", "new"), ("Arr2", "at")]),
        ("arr2_mutators", "`Arr2`: mutators",
         [("Arr2", "set"), ("Arr2", "add")]),
    ],
    "sampdist": [
        ("main", "`main`", [(None, "main")]),
    ],
    "neighsim": [
        ("main", "`main`", [(None, "main")]),
    ],
    "frescalo": [
        ("main", "`main`: CLI entry point and orchestration", [(None, "main")]),
        ("writeln", "`writeln`: output line formatting helper", [(None, "writeln")]),
        ("trend_writers", "`write_trend` / `write_zero_trend`: trend-file writers",
         [(None, "write_trend"), (None, "write_zero_trend")]),
        ("core", "`fresca` / `tfcalc`: rescaling and time-factor computation",
         [(None, "fresca"), (None, "tfcalc")]),
    ],
}

TANGLE_START_RE = re.compile(r'^([ \t]*)// @tangle:start (\S+)[ \t]*$')
TANGLE_END_RE = re.compile(r'^[ \t]*// @tangle:end (\S+)[ \t]*$')

# Pandoc/Quarto fenced code attribute block, e.g.:
#   ```{.rust filename="src/lib.rs" tangle="lib__blank_name"}
FENCE_OPEN_RE = re.compile(
    r'^```\{\.rust(?P<attrs>[^}]*)\}\s*$'
)
ATTR_RE = re.compile(r'(\w+)="([^"]*)"')


def tangle_id(modname, qualifier, name):
    return f"{modname}__{qualifier}__{name}" if qualifier else f"{modname}__{name}"


def chapter_path(modname, slug):
    return os.path.join(CHAPTERS_DIR, modname, f"{slug}.qmd")
