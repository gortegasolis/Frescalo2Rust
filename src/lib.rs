//! Shared support code for the Rust port of Mark Hill's FRESCALO suite
//! (Frescalo_1.f, neighsim_1.f, sampdist_1.f, January-June 2011).
//!
//! Fidelity conventions:
//! * Fortran `real` is single precision, so all floating-point work uses `f32`.
//! * Fortran `character*10` names are represented as blank-padded `[u8; 10]`
//!   byte arrays; comparisons are byte-wise, exactly like Fortran character
//!   comparison of equal-length strings.
//! * Fortran `character*30` packed records are `[u8; 30]`.
//! * All arrays keep the Fortran 1-based indexing: element 0 is unused.
//! * Formatted output replicates Fortran I/F edit descriptors, including
//!   blank padding, the trailing "." of F w.0, dropping of a leading zero
//!   when the field would otherwise overflow, and asterisks on overflow.

use std::fs::{File, OpenOptions};
use std::io::{self, BufRead, BufReader, BufWriter, Write};

pub type Name = [u8; 10];
pub type Rec30 = [u8; 30];

pub const BLANK: u8 = b' ';

// @tangle:start lib__blank_name
pub fn blank_name() -> Name {
    [BLANK; 10]
}
// @tangle:end lib__blank_name

// @tangle:start lib__name_to_string
pub fn name_to_string(n: &Name) -> String {
    String::from_utf8_lossy(n).trim_end().to_string()
}
// @tangle:end lib__name_to_string

// @tangle:start lib__make_rec30
pub fn make_rec30(a: &Name, b: &Name, c: &Name) -> Rec30 {
    let mut r = [BLANK; 30];
    r[0..10].copy_from_slice(a);
    r[10..20].copy_from_slice(b);
    r[20..30].copy_from_slice(c);
    r
}
// @tangle:end lib__make_rec30

/// field index: 0, 1 or 2 (Fortran d(1), d(2), d(3))
// @tangle:start lib__rec30_field
pub fn rec30_field(r: &Rec30, i: usize) -> Name {
    let mut n = [BLANK; 10];
    n.copy_from_slice(&r[i * 10..i * 10 + 10]);
    n
}
// @tangle:end lib__rec30_field

// ---------------------------------------------------------------------------
// getd: read one record and split into the first three blank-delimited words.
// Faithful port of the Fortran subroutine getd, including these quirks:
//   * words 1 and 2 are truncated to 10 characters, word 3 to 9 characters;
//   * the scan for the next word starts two positions after the end of the
//     previous word (k2+2), matching the original exactly;
//   * a blank line, or a line with fewer words than expected, returns with
//     the output fields UNCHANGED (the caller then re-processes stale data);
//   * end of file returns `false` (Fortran iend = 1).
// ---------------------------------------------------------------------------
pub struct DataReader {
    reader: BufReader<File>,
}

impl DataReader {
    // @tangle:start lib__DataReader__new
    pub fn new(f: File) -> Self {
        DataReader {
            reader: BufReader::new(f),
        }
    }
    // @tangle:end lib__DataReader__new

    /// Fortran rewind: seeking a BufReader discards its buffer.
    // @tangle:start lib__DataReader__rewind
    pub fn rewind(&mut self) {
        use std::io::Seek;
        let _ = self.reader.seek(std::io::SeekFrom::Start(0));
    }
    // @tangle:end lib__DataReader__rewind

    /// Returns false at end of file (iend = 1).
    // @tangle:start lib__DataReader__getd
    pub fn getd(&mut self, w1: &mut Name, w2: &mut Name, w3: &mut Name) -> bool {
        let mut line = String::new();
        match self.reader.read_line(&mut line) {
            Ok(0) => return false, // EOF -> iend = 1
            Ok(_) => {}
            Err(_) => return false,
        }
        if line.ends_with('\n') {
            line.pop();
        }
        // Fortran read with format a80: take at most 80 bytes, blank-padded.
        // Windows CRLF files: strip the trailing CR, as the records seen by
        // the original Windows executables never contained it.
        let bytes = line.as_bytes();
        let mut len = bytes.len();
        if len > 0 && bytes[len - 1] == b'\r' {
            len -= 1;
        }
        let mut b = [BLANK; 80];
        let n = len.min(80);
        b[..n].copy_from_slice(&bytes[..n]);

        // k1 = first non-blank
        let mut k = 0usize;
        while k < 80 && b[k] == BLANK {
            k += 1;
        }
        if k >= 80 {
            return true; // blank record: return with fields unchanged
        }
        let k1 = k;
        // k2 = last character of word 1
        let mut k = k1;
        while k < 80 && b[k] != BLANK {
            k += 1;
        }
        if k >= 80 {
            return true; // no terminating blank: return unchanged
        }
        let k2 = k - 1;
        // k3 = first character of word 2 (scan starts at k2+2)
        let mut k = k2 + 2;
        while k < 80 && b[k] == BLANK {
            k += 1;
        }
        if k >= 80 {
            return true;
        }
        let k3 = k;
        // k4 = last character of word 2
        let mut k = k3;
        while k < 80 && b[k] != BLANK {
            k += 1;
        }
        if k >= 80 {
            return true;
        }
        let k4 = k - 1;
        // k5 = first character of word 3 (scan starts at k4+2)
        let mut k = k4 + 2;
        while k < 80 && b[k] == BLANK {
            k += 1;
        }
        let k5 = if k >= 80 { 79 } else { k }; // Fortran: if(k.gt.79) k=80
        // k6 = last character of word 3
        let mut k = k5;
        while k < 80 && b[k] != BLANK {
            k += 1;
        }
        let mut k6 = if k >= 80 { 79 } else { k - 1 };
        if k5 == 79 {
            k6 = 79;
        }

        let mut nw1 = [BLANK; 10];
        for k in k1..=k2 {
            if k < k1 + 10 {
                nw1[k - k1] = b[k];
            }
        }
        let mut nw2 = [BLANK; 10];
        for k in k3..=k4 {
            if k < k3 + 10 {
                nw2[k - k3] = b[k];
            }
        }
        let mut nw3 = [BLANK; 10];
        for k in k5..=k6 {
            if k < k5 + 9 {
                nw3[k - k5] = b[k];
            }
        }
        *w1 = nw1;
        *w2 = nw2;
        *w3 = nw3;
        true
    }
    // @tangle:end lib__DataReader__getd
}

// ---------------------------------------------------------------------------
// getnum: parse a real number from a character*10 word (Fortran getnum).
// If the word has no decimal point before its first blank, a point is appended
// at the first blank position (this defeats the implied-decimals rule of the
// F10.4 descriptor).  Blanks are then ignored (BN mode) and the rest parsed;
// any parse error yields 0 (Fortran err=200 branch).
// ---------------------------------------------------------------------------
// @tangle:start lib__getnum
pub fn getnum(weight: &Name) -> f32 {
    let mut w = *weight;
    let mut idot = false;
    let mut kblank: Option<usize> = None;
    for k in 0..10 {
        if w[k] == b'.' {
            idot = true;
        }
        if w[k] == BLANK {
            kblank = Some(k);
            break;
        }
    }
    if !idot {
        if let Some(k) = kblank {
            w[k] = b'.';
        }
        // (If all 10 characters are non-blank with no dot, the original
        //  writes past w(10); we simply leave the word unchanged.)
    }
    let s: Vec<u8> = w.iter().filter(|&&c| c != BLANK).copied().collect();
    if s.is_empty() {
        return 0.0;
    }
    let s = String::from_utf8_lossy(&s);
    s.parse::<f32>().unwrap_or(0.0)
}
// @tangle:end lib__getnum

// ---------------------------------------------------------------------------
// binfnd: binary search in a sorted 1-based list of names.
// Returns the 1-based index, or 0 if not found.
// ---------------------------------------------------------------------------
// @tangle:start lib__binfnd
pub fn binfnd(ma: &[Name], n: usize, na: &Name) -> usize {
    let mut imin = 1usize;
    let mut iamin = ma[imin];
    let mut imax = n;
    let mut iamax = ma[imax];
    loop {
        if imax - imin <= 1 {
            if iamin == *na {
                return imin;
            }
            if iamax == *na {
                return imax;
            }
            return 0;
        }
        let imid = (imax + imin) / 2;
        let iamid = ma[imid];
        if *na <= iamid {
            imax = imid;
            iamax = iamid;
        } else {
            imin = imid;
            iamin = iamid;
        }
    }
}
// @tangle:end lib__binfnd

// ---------------------------------------------------------------------------
// addwrd: append samp to the sorted list sa (1-based, length *m) if new.
// ---------------------------------------------------------------------------
// @tangle:start lib__addwrd
pub fn addwrd(sa: &mut [Name], m: &mut usize, samp: &Name) {
    let i = if *m == 0 { 0 } else { binfnd(sa, *m, samp) };
    if i != 0 {
        return;
    }
    *m += 1;
    sa[*m] = *samp;
    sa[1..=*m].sort_unstable();
}
// @tangle:end lib__addwrd

/// neighsim variant of addwrd with the extra overflow guard.
// @tangle:start lib__addwrd_guarded
pub fn addwrd_guarded(sa: &mut [Name], mm: usize, m: &mut usize, samp: &Name) {
    let i = if *m == 0 { 0 } else { binfnd(sa, *m, samp) };
    if i != 0 {
        return;
    }
    *m += 1;
    if *m > mm {
        return;
    }
    sa[*m] = *samp;
    sa[1..=*m].sort_unstable();
}
// @tangle:end lib__addwrd_guarded

// ---------------------------------------------------------------------------
// Sorting.  The Fortran routines are deterministic heapsorts:
//   * sort10/sort30/isort/sort order plain keys ascending; equal keys are
//     interchangeable, so Rust's sort_unstable gives identical results.
//   * sort2 orders by (dict, type) lexicographically ascending - including
//     its tie-breaking rules - so an ascending sort on the pair is exact.
// Fortran real comparison treats -0.0 == 0.0; partial_cmp does the same.
// ---------------------------------------------------------------------------
// @tangle:start lib__sort30
pub fn sort30(dict: &mut [Rec30], n: usize) {
    dict[1..=n].sort_unstable();
}
// @tangle:end lib__sort30

// @tangle:start lib__sort_real
pub fn sort_real(dict: &mut [f32], n: usize) {
    dict[1..=n].sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
}
// @tangle:end lib__sort_real

// @tangle:start lib__isort
pub fn isort(dict: &mut [i32], n: usize) {
    dict[1..=n].sort_unstable();
}
// @tangle:end lib__isort

// @tangle:start lib__sort2
pub fn sort2(dict: &mut [f32], typ: &mut [i32], n: usize) {
    let mut pairs: Vec<(f32, i32)> = (1..=n).map(|i| (dict[i], typ[i])).collect();
    pairs.sort_unstable_by(|a, b| {
        a.0.partial_cmp(&b.0)
            .unwrap()
            .then_with(|| a.1.cmp(&b.1))
    });
    for (k, (d, t)) in pairs.into_iter().enumerate() {
        dict[1 + k] = d;
        typ[1 + k] = t;
    }
}
// @tangle:end lib__sort2

// ---------------------------------------------------------------------------
// Fortran-style formatted output.
// ---------------------------------------------------------------------------

/// Integer edit descriptor Iw: right-justified, asterisks on overflow.
// @tangle:start lib__ifmt
pub fn ifmt(v: i64, w: usize) -> String {
    let s = v.to_string();
    if s.len() > w {
        "*".repeat(w)
    } else {
        format!("{:>width$}", s, width = w)
    }
}
// @tangle:end lib__ifmt

/// Real edit descriptor Fw.d.
/// * rounds to d decimals (correct rounding of the exact binary value);
/// * Fw.0 keeps the trailing decimal point, as Fortran does ("  123.");
/// * if the field overflows, a leading zero is dropped ("0.74" -> ".74"),
///   and asterisks are printed if it still does not fit;
/// * Infinity/NaN follow gfortran conventions.
// @tangle:start lib__ffmt
pub fn ffmt(v: f32, w: usize, d: usize) -> String {
    if v.is_nan() {
        return if w >= 3 {
            format!("{:>width$}", "NaN", width = w)
        } else {
            "*".repeat(w)
        };
    }
    if v.is_infinite() {
        let (full, short) = if v > 0.0 {
            ("Infinity", "Inf")
        } else {
            ("-Infinity", "-Inf")
        };
        if w >= full.len() {
            return format!("{:>width$}", full, width = w);
        }
        if w >= short.len() {
            return format!("{:>width$}", short, width = w);
        }
        return "*".repeat(w);
    }
    let mut s = format!("{:.*}", d, v);
    if d == 0 {
        s.push('.');
    }
    if s.len() > w {
        if let Some(rest) = s.strip_prefix("0.") {
            s = format!(".{}", rest);
        } else if let Some(rest) = s.strip_prefix("-0.") {
            s = format!("-.{}", rest);
        }
    }
    if s.len() > w {
        return "*".repeat(w);
    }
    format!("{:>width$}", s, width = w)
}
// @tangle:end lib__ffmt

/// Record builder for formatted output lines (byte-oriented so that a10
/// name fields are reproduced exactly).
pub struct Rec {
    pub buf: Vec<u8>,
}

impl Rec {
    // @tangle:start lib__Rec__new
    pub fn new() -> Self {
        Rec { buf: Vec::new() }
    }
    // @tangle:end lib__Rec__new
    // @tangle:start lib__Rec__raw
    pub fn raw(&mut self, s: &[u8]) -> &mut Self {
        self.buf.extend_from_slice(s);
        self
    }
    // @tangle:end lib__Rec__raw
    // @tangle:start lib__Rec__s
    pub fn s(&mut self, s: &str) -> &mut Self {
        self.buf.extend_from_slice(s.as_bytes());
        self
    }
    // @tangle:end lib__Rec__s
    // @tangle:start lib__Rec__name
    pub fn name(&mut self, n: &Name) -> &mut Self {
        self.buf.extend_from_slice(n);
        self
    }
    // @tangle:end lib__Rec__name
    // @tangle:start lib__Rec__x
    pub fn x(&mut self, n: usize) -> &mut Self {
        self.buf.extend(std::iter::repeat(BLANK).take(n));
        self
    }
    // @tangle:end lib__Rec__x
    // @tangle:start lib__Rec__i
    pub fn i(&mut self, v: i64, w: usize) -> &mut Self {
        self.s(&ifmt(v, w));
        self
    }
    // @tangle:end lib__Rec__i
    // @tangle:start lib__Rec__f
    pub fn f(&mut self, v: f32, w: usize, d: usize) -> &mut Self {
        self.s(&ffmt(v, w, d));
        self
    }
    // @tangle:end lib__Rec__f
    // @tangle:start lib__Rec__writeln
    pub fn writeln<W: Write>(&self, w: &mut W) {
        w.write_all(&self.buf).unwrap();
        w.write_all(b"\n").unwrap();
    }
    // @tangle:end lib__Rec__writeln
}

// ---------------------------------------------------------------------------
// Console (stdout) output.  Prompts are flushed so that they appear before
// input is read.  List-directed output (write(*,*)) is approximated: leading
// blank, character items concatenated, integers in a field of width 12.
// ---------------------------------------------------------------------------
// @tangle:start lib__cout
pub fn cout(s: &str) {
    let mut o = io::stdout();
    let _ = o.write_all(s.as_bytes());
    let _ = o.write_all(b"\n");
    let _ = o.flush();
}
// @tangle:end lib__cout

// @tangle:start lib__ld_i
pub fn ld_i(v: i64) -> String {
    format!("{:12}", v)
}
// @tangle:end lib__ld_i

// @tangle:start lib__ld_f
pub fn ld_f(v: f32) -> String {
    format!(" {:>12}", format!("{}", v))
}
// @tangle:end lib__ld_f

/// Approximate a Fortran list-directed write: a leading blank followed by the
/// (already concatenated) items.
// @tangle:start lib__ld_line
pub fn ld_line(body: &str) {
    cout(&format!(" {}", body));
}
// @tangle:end lib__ld_line

// ---------------------------------------------------------------------------
// Console input.  All console reads in the originals are formatted or
// list-directed reads from unit *; hitting EOF without an end= branch is a
// runtime error in Fortran, so we exit with an error message.
// ---------------------------------------------------------------------------
// @tangle:start lib__read_stdin_line
fn read_stdin_line(stdin: &mut impl BufRead) -> String {
    let mut s = String::new();
    match stdin.read_line(&mut s) {
        Ok(0) => {
            eprintln!("\nFortran runtime error: End of file");
            std::process::exit(2);
        }
        Ok(_) => {
            if s.ends_with('\n') {
                s.pop();
                if s.ends_with('\r') {
                    s.pop();
                }
            }
            s
        }
        Err(e) => {
            eprintln!("\nFortran runtime error: {}", e);
            std::process::exit(2);
        }
    }
}
// @tangle:end lib__read_stdin_line

/// read(*,1000) filein  with format a20: first 20 characters of the record,
/// trailing blanks trimmed for OPEN.
// @tangle:start lib__read_a20
pub fn read_a20(stdin: &mut impl BufRead) -> String {
    let line = read_stdin_line(stdin);
    let bytes = line.as_bytes();
    let n = bytes.len().min(20);
    String::from_utf8_lossy(&bytes[..n]).trim_end().to_string()
}
// @tangle:end lib__read_a20

/// read(*,2010,err=...) x  with format f8.4 (BN blank mode):
/// blanks are ignored; an all-blank field reads as 0; without a decimal
/// point the field has 4 implied decimals.  Returns None on the err= branch.
// @tangle:start lib__read_f8_4
pub fn read_f8_4(stdin: &mut impl BufRead) -> Option<f32> {
    let line = read_stdin_line(stdin);
    let bytes = line.as_bytes();
    let n = bytes.len().min(8);
    let field: Vec<u8> = bytes[..n].iter().filter(|&&c| c != BLANK).copied().collect();
    if field.is_empty() {
        return Some(0.0);
    }
    let s = String::from_utf8_lossy(&field).into_owned();
    if s.contains('.') {
        s.parse::<f32>().ok()
    } else {
        match s.parse::<i64>() {
            Ok(v) => Some(v as f32 / 10000.0),
            Err(_) => None,
        }
    }
}
// @tangle:end lib__read_f8_4

/// read(*,*) neigh : list-directed integer read.
// @tangle:start lib__read_int_listdirected
pub fn read_int_listdirected(stdin: &mut impl BufRead) -> i32 {
    let line = read_stdin_line(stdin);
    match line.split_whitespace().next() {
        Some(tok) => match tok.parse::<i32>() {
            Ok(v) => v,
            Err(_) => {
                eprintln!("\nFortran runtime error: Bad integer for item 1 in list input");
                std::process::exit(2);
            }
        },
        None => {
            eprintln!("\nFortran runtime error: End of file");
            std::process::exit(2);
        }
    }
}
// @tangle:end lib__read_int_listdirected

// ---------------------------------------------------------------------------
// filin / filout / hold
// ---------------------------------------------------------------------------
// @tangle:start lib__filin
pub fn filin(stdin: &mut impl BufRead) -> (String, File) {
    loop {
        let name = read_a20(stdin);
        match File::open(&name) {
            Ok(f) => return (name, f),
            Err(_) => {
                cout("");
                cout("  *** ERROR *** File does not exist");
                cout(" Type another name");
            }
        }
    }
}
// @tangle:end lib__filin

// @tangle:start lib__filout
pub fn filout(stdin: &mut impl BufRead) -> (String, File) {
    loop {
        let name = read_a20(stdin);
        match OpenOptions::new().write(true).create_new(true).open(&name) {
            Ok(f) => return (name, f),
            Err(_) => {
                cout("");
                cout("  *** ERROR *** File already exists");
                cout(" Type another name");
            }
        }
    }
}
// @tangle:end lib__filout

/// Fortran hold: prompt and wait for <RETURN>, then stop.
// @tangle:start lib__hold
pub fn hold(stdin: &mut impl BufRead) -> ! {
    let mut o = io::stdout();
    let _ = o.write_all(b"\n\nPress <RETURN> to exit\n\n\n----------------------\n");
    let _ = o.flush();
    let mut s = String::new();
    let _ = stdin.read_line(&mut s);
    std::process::exit(0);
}
// @tangle:end lib__hold

// ---------------------------------------------------------------------------
// 2-D array helper with Fortran 1-based indexing (row 0 / column 0 unused).
// ---------------------------------------------------------------------------
pub struct Arr2<T: Copy> {
    pub data: Vec<T>,
    pub stride: usize,
}

impl<T: Copy> Arr2<T> {
    /// Allocate rows 0..=nrows and columns 0..=ncols.
    // @tangle:start lib__Arr2__new
    pub fn new(nrows: usize, ncols: usize, init: T) -> Self {
        Arr2 {
            data: vec![init; (nrows + 1) * (ncols + 1)],
            stride: ncols + 1,
        }
    }
    // @tangle:end lib__Arr2__new
    #[inline(always)]
    // @tangle:start lib__Arr2__at
    pub fn at(&self, i: usize, j: usize) -> T {
        self.data[i * self.stride + j]
    }
    // @tangle:end lib__Arr2__at
    #[inline(always)]
    // @tangle:start lib__Arr2__set
    pub fn set(&mut self, i: usize, j: usize, v: T) {
        self.data[i * self.stride + j] = v;
    }
    // @tangle:end lib__Arr2__set
    #[inline(always)]
    // @tangle:start lib__Arr2__add
    pub fn add(&mut self, i: usize, j: usize, v: T)
    where
        T: std::ops::AddAssign,
    {
        self.data[i * self.stride + j] += v;
    }
    // @tangle:end lib__Arr2__add
}

pub type LogWriter = BufWriter<File>;
