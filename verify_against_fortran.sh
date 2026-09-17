#!/usr/bin/env bash
# Verify the Rust port against a native gfortran build of the original
# Fortran sources, using the example datasets in Original_frescalo/.
#
# Usage: ./verify_against_fortran.sh
#
# Builds both implementations, runs the full pipeline
# (sampdist -> neighsim -> frescalo) with each, and diffs every output file.
# All output files must be byte-identical (modulo the line-ending of the
# Windows reference files, which is a platform artifact).
#
# Note: the original Fortran sources rely on argument aliasing that gfortran
# miscompiles at -O2 (coordinates read as 0), so the reference build uses -O0.
# Input files are CR-stripped for the Fortran runs: the original Windows
# executables never saw the CR (text-mode I/O), and the Rust port strips CR
# itself to emulate that behaviour.

set -euo pipefail
cd "$(dirname "$0")"

WORK=verify_work
rm -rf "$WORK"
mkdir -p "$WORK/f" "$WORK/r"

echo "== Building Rust port =="
cargo build --release --quiet

echo "== Building original Fortran (gfortran -O0) =="
gfortran -std=legacy -O0 -o "$WORK/f/sampdist_f" Original_frescalo/sampdist_1.f
gfortran -std=legacy -O0 -o "$WORK/f/neighsim_f" Original_frescalo/neighsim_1.f
gfortran -std=legacy -O0 -o "$WORK/f/frescalo_f" Original_frescalo/Frescalo_1.f

for d in f r; do
    for f in Samp_locations.txt Training_vasc.txt Test.txt NotBench.txt; do
        tr -d '\r' < "Original_frescalo/$f" > "$WORK/$d/$f"
    done
done

echo "== Running sampdist =="
( cd "$WORK/f" && printf 'Samp_locations.txt\ndist.txt\n200\n\n' | ./sampdist_f >/dev/null )
( cd "$WORK/r" && printf 'Samp_locations.txt\ndist.txt\n200\n\n' | ../../target/release/sampdist >/dev/null )

echo "== Running neighsim =="
( cd "$WORK/f" && printf 'Training_vasc.txt\ndist.txt\nsim.txt\nweights.txt\n100\n\n' | ./neighsim_f >/dev/null )
( cd "$WORK/r" && printf 'Training_vasc.txt\ndist.txt\nsim.txt\nweights.txt\n100\n\n' | ../../target/release/neighsim >/dev/null )

echo "== Running frescalo (defaults, no exclusions) =="
( cd "$WORK/f" && printf 'log.txt\nTest.txt\nweights.txt\n\nsamples.txt\nfrequencies.txt\ntrends.txt\n\n\n\n' | ./frescalo_f >/dev/null )
( cd "$WORK/r" && printf 'log.txt\nTest.txt\nweights.txt\n\nsamples.txt\nfrequencies.txt\ntrends.txt\n\n\n\n' | ../../target/release/frescalo >/dev/null )

echo "== Running frescalo (NotBench.txt exclusions, phi=0.80, blim=0.15) =="
mkdir -p "$WORK/f2" "$WORK/r2"
for d in f2 r2; do
    for f in Test.txt NotBench.txt; do tr -d '\r' < "Original_frescalo/$f" > "$WORK/$d/$f"; done
    cp "$WORK/f/weights.txt" "$WORK/$d/weights.txt"
done
( cd "$WORK/f2" && printf 'log.txt\nTest.txt\nweights.txt\nNotBench.txt\nsamples.txt\nfrequencies.txt\ntrends.txt\n0.80\n0.15\n\n' | ../f/frescalo_f >/dev/null )
( cd "$WORK/r2" && printf 'log.txt\nTest.txt\nweights.txt\nNotBench.txt\nsamples.txt\nfrequencies.txt\ntrends.txt\n0.80\n0.15\n\n' | ../../target/release/frescalo >/dev/null )

status=0
for pair in "f r dist.txt" "f r sim.txt" "f r weights.txt" \
            "f r log.txt" "f r samples.txt" "f r frequencies.txt" "f r trends.txt" \
            "f2 r2 log.txt" "f2 r2 samples.txt" "f2 r2 frequencies.txt" "f2 r2 trends.txt"; do
    set -- $pair
    if diff -q "$WORK/$1/$3" "$WORK/$2/$3" >/dev/null; then
        echo "IDENTICAL: $1/$3"
    else
        echo "DIFFERS:   $1/$3"
        status=1
    fi
done

if [ "$status" -eq 0 ]; then
    echo "ALL OUTPUTS BYTE-IDENTICAL"
else
    echo "SOME OUTPUTS DIFFER"
fi
exit $status
