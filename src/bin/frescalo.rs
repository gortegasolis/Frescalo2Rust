//! Rust port of Frescalo_1.f:
//! FRESCALO - Trend_analysis using local frequencies
//! written by Mark Hill, January-June 2011
//!
//! Calculates sampling effort multipliers and rescaled species probabilities
//! for each location (subroutine fresca), and time factors for each species
//! (subroutine tfcalc), following Hill (2011) Methods in Ecology and
//! Evolution 2: 502-512.

use frescalo::*;
use std::io::{self, BufRead, BufWriter, Write};

// Input limits (Fortran parameter values)
const MM: usize = 4000; // max number of samples
const NN: usize = 2000; // max number of species
const NNT: usize = 100; // max number of time periods
const NNDAT: usize = 2_000_000; // size of main data matrix
const NNWGT: usize = 500_000; // size of neighbourhood weight data

const FMAX: f32 = 0.99999; // < 1 so that log(1-f) can be taken
const FMIN: f32 = 1.0E-10; // ensures frequencies do not sum to zero
const TOL: f32 = 0.0003; // convergence limit for rescaling
const IREPMX: i32 = 100; // maximum number of iterations
const PHIDEF: f32 = 0.74; // default phi, the local frequency
const BLMDEF: f32 = 0.2703; // default limit for benchmark species

// unused_assignments: the jx/iitx marker updates inside the zero-padding
// blocks are redundant but kept to mirror the Fortran control flow exactly.
#[allow(unused_assignments)]
// @tangle:start frescalo__main
fn main() {
    let mut stdin = io::stdin().lock();

    // Allocations (1-based indexing; element 0 unused), mirroring the
    // original's statically dimensioned arrays.
    let mut sa = vec![blank_name(); MM + 2];
    let mut sp = vec![blank_name(); NN + 2];
    let mut tim = vec![blank_name(); NNT + 2];
    let mut bnchx = vec![blank_name(); NN + 2];
    let mut numwgt = vec![0i32; MM + 2];
    let mut idat = vec![0i32; MM + 2];
    let mut iitot = vec![0i32; MM + 2];
    let mut jrank = vec![0i32; NN + 2];
    let mut iocc = vec![0i32; MM + 2];
    let mut jocc = vec![0i32; NN + 2];
    let mut lendat: Arr2<i32> = Arr2::new(NN + 1, NNT + 1, 0);
    let mut idata: Arr2<i32> = Arr2::new(MM + 1, NN + 1, 0);
    let mut ibench: Arr2<i32> = Arr2::new(MM + 1, NN + 1, 0);
    let mut sampef: Arr2<f32> = Arr2::new(MM + 1, NNT + 1, 1.0E-7);
    let mut ffij: Arr2<f32> = Arr2::new(MM + 1, NN + 1, 0.0);
    let mut smpint = vec![0.0f32; MM + 2];
    let mut fff = vec![0.0f32; MM + 2];
    let mut wgttot = vec![0.0f32; MM + 2];
    let mut wgtt2 = vec![0.0f32; MM + 2];
    let mut f = vec![0.0f32; NN + 2];
    let mut ff = vec![0.0f32; NN + 2];
    let mut ffff = vec![0.0f32; MM + 2];
    let mut abtot = vec![0.0f32; MM + 2];
    let mut bwght = vec![1.0f32; NN + 2]; // 0.001 for unsuitable bench species
    let mut ddijt = vec![[b' '; 30]; NNDAT + 2];
    let mut wneigh = vec![[b' '; 30]; NNWGT + 2];

    // Banner (format 1998; the leading // produces two blank records)
    cout("");
    cout("");
    cout(" FRESCALO - Trend_analysis using local frequencies");
    cout(" written by Mark Hill, January-June 2011");
    cout("");
    cout(" Input limits:");
    cout(&format!("    Number of samples = {}", ifmt(MM as i64, 7)));
    cout(&format!("    Number of species = {}", ifmt(NN as i64, 7)));
    cout(&format!("    Number of time periods = {}", ifmt(NNT as i64, 7)));
    cout(&format!("    Number of observations = {}", ifmt(NNDAT as i64, 7)));
    cout(&format!(
        "    Number of neighbourhood weights = {}",
        ifmt(NNWGT as i64, 7)
    ));
    cout("");
    cout(" Type name of log file ...");

    let (fileou, flog) = filout(&mut stdin);
    let mut unit10 = BufWriter::new(flog);

    // Log banner (format 1999)
    writeln(&mut unit10, " Log file for FRESCALO");
    writeln(&mut unit10, "");
    writeln(&mut unit10, " Input limits:");
    writeln(&mut unit10, &format!("    Number of samples = {}", ifmt(MM as i64, 7)));
    writeln(&mut unit10, &format!("    Number of species = {}", ifmt(NN as i64, 7)));
    writeln(&mut unit10, &format!("    Number of time periods = {}", ifmt(NNT as i64, 7)));
    writeln(&mut unit10, &format!("    Number of observations = {}", ifmt(NNDAT as i64, 7)));
    writeln(
        &mut unit10,
        &format!("    Number of neighbourhood weights = {}", ifmt(NNWGT as i64, 7)),
    );
    writeln(&mut unit10, "");
    writeln(&mut unit10, " Type name of log file ...");
    writeln(&mut unit10, &format!("{:<20}", fileou)); // 1000 format(a20)

    cout(" Type occurrence input file [sample species time] ...");
    writeln(&mut unit10, " Type occurrence input file [sample species time] ...");
    let (filein, focc) = filin(&mut stdin);
    writeln(&mut unit10, &format!("{:<20}", filein));
    let mut occ_reader = DataReader::new(focc);

    cout(" Type neighbourhood weight input file ...");
    writeln(&mut unit10, " Type neighbourhood weight input file ...");
    let (filein, fwgt) = filin(&mut stdin);
    writeln(&mut unit10, &format!("{:<20}", filein));
    let mut wgt_reader = DataReader::new(fwgt);

    cout(" Type file with species to exclude from benchmarks ");
    cout("     or press <RETURN> if no exclusions...");
    writeln(&mut unit10, " Type file with species to exclude from benchmarks ");
    writeln(&mut unit10, "     or press <RETURN> if no exclusions...");
    // label 20: on open error the original silently re-reads, without
    // reprinting the prompt or any error message.
    let mut nbnchx: usize = 0;
    loop {
        let filein = read_a20(&mut stdin);
        if filein.is_empty() {
            writeln(&mut unit10, "<No exclusions>");
            break;
        }
        match std::fs::File::open(&filein) {
            Ok(fx) => {
                writeln(&mut unit10, &format!("{:<20}", filein));
                // read(2,1001,end=25) with format a10: first 10 bytes of each
                // line, blank-padded.
                let br = io::BufReader::new(fx);
                let mut lines = br.lines();
                for jj in 1..=NN {
                    match lines.next() {
                        Some(Ok(line)) => {
                            let bytes = line.as_bytes();
                            let mut len = bytes.len();
                            if len > 0 && bytes[len - 1] == b'\r' {
                                len -= 1;
                            }
                            let mut nm = blank_name();
                            let c = len.min(10);
                            nm[..c].copy_from_slice(&bytes[..c]);
                            nbnchx = jj;
                            bnchx[nbnchx] = nm;
                        }
                        _ => break,
                    }
                }
                break;
            }
            Err(_) => continue,
        }
    }

    cout(" Type name of sample stats output file...");
    writeln(&mut unit10, " Type name of sample stats output file...");
    let (fileou, fsam) = filout(&mut stdin);
    writeln(&mut unit10, &format!("{:<20}", fileou));
    let mut unit7 = BufWriter::new(fsam);

    cout(" Type name of rescaled frequency file...");
    writeln(&mut unit10, " Type name of rescaled frequency file...");
    let (fileou, ffq) = filout(&mut stdin);
    writeln(&mut unit10, &format!("{:<20}", fileou));
    let mut unit8 = BufWriter::new(ffq);

    cout(" Type name of trend output file...");
    writeln(&mut unit10, " Type name of trend output file...");
    let (fileou, ftr) = filout(&mut stdin);
    writeln(&mut unit10, &format!("{:<20}", fileou));
    let mut unit9 = BufWriter::new(ftr);

    // Target value of phi (label 40)
    let mut phibig = PHIDEF;
    loop {
        let p = format!(
            " Type target value of local frequency phi (default={})...",
            ffmt(PHIDEF, 4, 2)
        );
        cout(&p);
        writeln(&mut unit10, &p);
        if let Some(v) = read_f8_4(&mut stdin) {
            phibig = v;
        }
        if phibig == 0.0 {
            phibig = PHIDEF;
        }
        if phibig > 0.95 || phibig < 0.50 {
            cout(" ***ERROR*** Outside range 0.50 to 0.95");
            continue;
        }
        break;
    }
    // 2011 format has a trailing slash: blank line follows
    let p = format!(" Target value is {}", ffmt(phibig, 4, 2));
    cout(&p);
    cout("");
    writeln(&mut unit10, &p);
    writeln(&mut unit10, "");

    // Benchmark limit (label 46)
    let mut blim = BLMDEF;
    loop {
        let p = format!(" Type value of Benchmark Limit (default={})...", ffmt(BLMDEF, 4, 2));
        cout(&p);
        writeln(&mut unit10, &p);
        if let Some(v) = read_f8_4(&mut stdin) {
            blim = v;
        }
        if blim == 0.0 {
            blim = BLMDEF;
        }
        if blim > 0.5 || blim < 0.08 {
            cout(" ***ERROR*** Outside range 0.08 to 0.5");
            continue;
        }
        break;
    }
    // 2014 format has a trailing slash: blank line follows
    let p = format!(" Benchmark limit is {}", ffmt(blim, 4, 2));
    cout(&p);
    cout("");
    writeln(&mut unit10, &p);
    writeln(&mut unit10, "");

    // ------------------------------------------------------------------
    // Read in data
    // ------------------------------------------------------------------
    let mut m: usize = 0;
    let mut n: usize = 0;
    let mut nt: usize = 0;
    let mut ndtji: usize = 0;
    let mut nwgt: usize = 0;

    let mut samp = blank_name();
    let mut samp1 = blank_name();
    let mut spec = blank_name();
    let mut time = blank_name();
    let mut weight = blank_name();

    // Read hectad weights (samp, samp1, weight); stored as (samp1,samp,weight)
    // so that all targets of which samp1 is a neighbour are together.
    loop {
        if nwgt == 0 {
            ld_line("Reading in smoothing weights from samples");
        }
        if !wgt_reader.getd(&mut samp, &mut samp1, &mut weight) {
            break;
        }
        nwgt += 1;
        if nwgt % 20000 == 0 {
            ld_line(&format!(
                "Weights {}{}{}",
                name_to_string(&samp),
                name_to_string(&samp1),
                ld_i(nwgt as i64)
            ));
        }
        if nwgt > NNWGT {
            let p = format!(
                " No of neighbourhood weights > maximum which is{}",
                ifmt(NNWGT as i64, 5)
            );
            cout(&p);
            writeln(&mut unit10, &p);
            hold(&mut stdin);
        }
        addwrd(&mut sa, &mut m, &samp);
        addwrd(&mut sa, &mut m, &samp1);
        wneigh[nwgt] = make_rec30(&samp1, &samp, &weight);
    }

    cout(" Sorting local frequency weights ...");
    sort30(&mut wneigh, nwgt);

    // Read the occurrence data [sample species time]
    loop {
        if !occ_reader.getd(&mut samp, &mut spec, &mut time) {
            break;
        }
        let i = binfnd(&sa, m, &samp);
        if i == 0 {
            if samp != samp1 {
                // 2996: location in species data but not in neighbourhood weights
                let mut r = Rec::new();
                r.s("*** ").name(&samp).s(
                    " location ignored - in species data but not listed in neighbourhood weights",
                );
                r.writeln(&mut unit10);
            }
            samp1 = samp;
            continue;
        }
        addwrd(&mut sp, &mut n, &spec);
        addwrd(&mut tim, &mut nt, &time);
        if nt > NNT {
            let p = format!(
                " Number of time periods exceeds maximum which is{}",
                ifmt(NNT as i64, 5)
            );
            cout(&p);
            writeln(&mut unit10, &p);
            hold(&mut stdin);
        }
        ndtji += 1;
        if ndtji % 20000 == 0 {
            ld_line(&format!(
                "{}{}{}{}",
                name_to_string(&samp),
                name_to_string(&spec),
                name_to_string(&time),
                ld_i(ndtji as i64)
            ));
        }
        ddijt[ndtji] = make_rec30(&samp, &spec, &time);
    }

    // 2505 / 2506: report actual numbers
    let rep = |unit10: &mut LogWriter| {
        writeln(unit10, "");
        writeln(unit10, " Actual numbers in data");
        writeln(unit10, &format!("    Number of samples      {}", ifmt(m as i64, 8)));
        writeln(unit10, &format!("    Number of species      {}", ifmt(n as i64, 8)));
        writeln(unit10, &format!("    Number of time periods {}", ifmt(nt as i64, 8)));
        writeln(unit10, &format!("    Number of observations {}", ifmt(ndtji as i64, 8)));
        writeln(unit10, &format!("    Neighbourhood weights  {}", ifmt(nwgt as i64, 8)));
        writeln(unit10, &format!("    Benchmark exclusions   {}", ifmt(nbnchx as i64, 8)));
        writeln(unit10, "");
    };
    cout("");
    cout(" Actual numbers in data");
    cout(&format!("    Number of samples      {}", ifmt(m as i64, 8)));
    cout(&format!("    Number of species      {}", ifmt(n as i64, 8)));
    cout(&format!("    Number of time periods {}", ifmt(nt as i64, 8)));
    cout(&format!("    Number of observations {}", ifmt(ndtji as i64, 8)));
    cout(&format!("    Neighbourhood weights  {}", ifmt(nwgt as i64, 8)));
    cout(&format!("    Benchmark exclusions   {}", ifmt(nbnchx as i64, 8)));
    cout("");
    rep(&mut unit10);
    if nbnchx > 0 {
        writeln(&mut unit10, " Benchmark exclusions");
        for ib in 1..=nbnchx {
            let mut r = Rec::new();
            r.x(4).name(&bnchx[ib]);
            r.writeln(&mut unit10);
        }
    }

    cout(" Sorting main data ...");
    sort30(&mut ddijt, ndtji);

    // First calculate sum of weights for each sample
    for iwgt in 1..=nwgt {
        let dtji = wneigh[iwgt];
        let f_neigh = rec30_field(&dtji, 0); // samp1: the neighbour
        let f_targ = rec30_field(&dtji, 1); // samp: the target
        let f_wgt = rec30_field(&dtji, 2);
        let i = binfnd(&sa, m, &f_targ);
        let ii = binfnd(&sa, m, &f_neigh);
        let wgt = getnum(&f_wgt);
        wgttot[i] += wgt;
        wgtt2[i] += wgt * wgt;
        numwgt[ii] += 1;
    }

    // Now calculate the number of data items for each sample in main data
    let mut spec1 = blank_name();
    for idtji in 1..=ndtji {
        let dtji = ddijt[idtji];
        let f_samp = rec30_field(&dtji, 0);
        let f_spec = rec30_field(&dtji, 1);
        let i = binfnd(&sa, m, &f_samp);
        idat[i] += 1;
        if f_spec != spec1 {
            iitot[i] += 1;
            spec1 = f_spec;
        }
    }

    // Now calculate frequencies
    let mut idtji = 0usize;
    let mut iwgt = 0usize;
    let mut ii_leftover = 0usize; // tracks Fortran loop variable ii afterwards
    for ii in 1..=m {
        ii_leftover = ii + 1;
        if ii % 100 == 0 {
            ld_line(&format!("frequencies {}{}", name_to_string(&sa[ii]), ld_i(ii as i64)));
        }
        for j in 1..=n {
            jocc[j] = 0;
            idata.set(ii, j, 0);
        }
        for _ in 1..=idat[ii] {
            idtji += 1;
            let dtji = ddijt[idtji];
            let f_spec = rec30_field(&dtji, 1);
            let j = binfnd(&sp, n, &f_spec);
            jocc[j] = 1;
            idata.set(ii, j, 1);
        }
        for _ in 1..=numwgt[ii] {
            iwgt += 1;
            let dtji = wneigh[iwgt];
            let f_targ = rec30_field(&dtji, 1);
            let f_wgt = rec30_field(&dtji, 2);
            let i = binfnd(&sa, m, &f_targ);
            let wgt = getnum(&f_wgt);
            let denom = wgttot[i] + 1.0E-10;
            for j in 1..=n {
                let contrib = (jocc[j] as f32) * wgt / denom;
                ffij.add(i, j, contrib);
            }
        }
    }
    // Local frequencies are now in ffij(i,j).

    // Downweight species that are unsuitable as benchmarks
    for ib in 1..=nbnchx {
        let j = binfnd(&sp, n, &bnchx[ib]);
        if j != 0 {
            bwght[j] = 0.001;
        }
    }

    for i in 1..=m {
        // NB: the original mistakenly uses ii (leftover loop variable) here;
        // we reproduce the condition but avoid the out-of-bounds name lookup.
        if ii_leftover % 100 == 0 {
            let name = sa.get(ii_leftover).map(name_to_string).unwrap_or_default();
            ld_line(&format!("rescaling {}{}", name, ld_i(ii_leftover as i64)));
        }
        for j in 1..=n {
            f[j] = ffij.at(i, j);
            jocc[j] = idata.at(i, j);
        }
        let samp_i = sa[i];
        let itot = iitot[i];
        // wn2 is the Effective Number of weights in the neighbourhood
        let wn2 = wgttot[i] * wgttot[i] / (wgtt2[i] + 1.0E-12);
        let (phi1, spnum) = fresca(
            i,
            n,
            itot,
            &jocc,
            &mut f,
            &mut ff,
            &mut jrank,
            &samp_i,
            &sp,
            phibig,
            FMAX,
            FMIN,
            wn2,
            TOL,
            IREPMX,
            &mut unit7,
            &mut unit8,
        );
        ffff[i] = phi1;
        abtot[i] = 1.0E-7;
        for j in 1..=n {
            if ffij.at(i, j) != 0.0 {
                ffij.set(i, j, ff[j]);
            }
            let jj = jrank[j] as usize;
            let rank1 = (j as f32) / spnum;
            // The case j=1 is included because with small samples rank1 may
            // be greater than the limit for j=1
            if rank1 < blim || j == 1 {
                ibench.set(i, jj, 1);
                abtot[i] += bwght[jj];
            } else {
                ibench.set(i, jj, 0);
            }
        }
    }

    // The next task is to reorder the main data matrix
    for idtji in 1..=ndtji {
        let dtji = ddijt[idtji];
        let f_samp = rec30_field(&dtji, 0);
        let f_spec = rec30_field(&dtji, 1);
        let f_time = rec30_field(&dtji, 2);
        if idtji % 20000 == 0 {
            ld_line(&format!(
                "Re-ordering {}{}{}{}",
                name_to_string(&f_samp),
                name_to_string(&f_spec),
                name_to_string(&f_time),
                ld_i(idtji as i64)
            ));
        }
        ddijt[idtji] = make_rec30(&f_spec, &f_time, &f_samp);
    }
    cout(" Now doing second sort of reordered main data ...");
    sort30(&mut ddijt, ndtji);
    // From now onwards the sorted form is called ddjti.

    for idtji in 1..=ndtji {
        if idtji % 20000 == 0 {
            ld_line(&format!(
                "Main data to calc sampling effort{}{}",
                String::from_utf8_lossy(&ddijt[idtji]),
                ld_i(idtji as i64)
            ));
        }
        let dtji = ddijt[idtji];
        let f_spec = rec30_field(&dtji, 0);
        let f_time = rec30_field(&dtji, 1);
        let f_samp = rec30_field(&dtji, 2);
        let i = binfnd(&sa, m, &f_samp);
        let j = binfnd(&sp, n, &f_spec);
        let iit = binfnd(&tim, nt, &f_time);
        lendat.add(j, iit, 1);
        let contrib = (ibench.at(i, j) as f32) * bwght[j] / abtot[i];
        sampef.add(i, iit, contrib);
    }

    let mut idtji = 0usize;
    let mut jx = 1usize;
    let mut iitx = 1usize;
    // These are markers to ensure that we print out time periods at which
    // species were not recorded.

    // 2060 header
    writeln(
        &mut unit9,
        "Species__  Time______ TFactor St_Dev _Count ___spt ___est N>0.00 N>0.98",
    );

    // label 130 loop
    'main: loop {
        for i in 1..=m {
            iocc[i] = 0;
        }
        idtji += 1;
        if idtji > ndtji {
            break 'main;
        }
        let dtji = ddijt[idtji];
        let f_spec = rec30_field(&dtji, 0);
        let f_time = rec30_field(&dtji, 1);
        let mut f_samp = rec30_field(&dtji, 2);
        let j = binfnd(&sp, n, &f_spec);
        let iit = binfnd(&tim, nt, &f_time);
        for i in 1..=m {
            smpint[i] = sampef.at(i, iit);
            fff[i] = ffij.at(i, j);
        }
        let i = binfnd(&sa, m, &f_samp);
        iocc[i] = 1;
        for _ in 1..lendat.at(j, iit) {
            idtji += 1;
            if idtji > ndtji {
                break 'main;
            }
            let dtji = ddijt[idtji];
            f_samp = rec30_field(&dtji, 2);
            let i = binfnd(&sa, m, &f_samp);
            iocc[i] = 1;
        }
        let (tf, sd, spt, jtot, est, ic1, ic2) = tfcalc(&iocc, &smpint, &fff, m);
        // first pad out the output with times when the species was not recorded
        if sp[jx] < sp[j] {
            for iiit in iitx..=nt {
                write_zero_trend(&mut unit9, &sp[jx], &tim[iiit]);
            }
            for jj in (jx + 1)..j {
                for iiit in 1..=nt {
                    write_zero_trend(&mut unit9, &sp[jj], &tim[iiit]);
                }
            }
            jx = j;
            iitx = 1;
        }
        if tim[iitx] < tim[iit] {
            for iiit in iitx..iit {
                write_zero_trend(&mut unit9, &sp[j], &tim[iiit]);
            }
            iitx = iit;
        }
        // 2050 format(a10,1x,a10,f8.3,f7.3,i7,2f7.1,2i7)
        write_trend(&mut unit9, &sp[j], &tim[iit], tf, sd, jtot, spt, est, ic1, ic2);
        if j % 10 == 0 {
            let mut r = Rec::new();
            r.name(&sp[j])
                .x(1)
                .name(&tim[iit])
                .f(tf, 8, 3)
                .f(sd, 7, 3)
                .i(j as i64, 7);
            cout(&String::from_utf8_lossy(&r.buf));
        }
        jx = j;
        iitx = iit + 1;
    }

    // label 140
    for iiit in iitx..=nt {
        write_zero_trend(&mut unit9, &sp[jx], &tim[iiit]);
    }

    unit7.flush().unwrap();
    unit8.flush().unwrap();
    unit9.flush().unwrap();

    // Finally test whether given value of phi appears to be unrealistically low
    // (2513: leading // produces two blank records)
    sort_real(&mut ffff, m);
    let i985 = (0.985f32 * m as f32) as usize;
    let phi985 = ffff[i985];
    for _ in 0..2 {
        writeln(&mut unit10, "");
        cout("");
    }
    writeln(&mut unit10, &format!(" 98.5 percentile of input phi {}", ffmt(phi985, 5, 2)));
    writeln(&mut unit10, &format!(" Target value of phi          {}", ffmt(phibig, 5, 2)));
    cout(&format!(" 98.5 percentile of input phi {}", ffmt(phi985, 5, 2)));
    cout(&format!(" Target value of phi          {}", ffmt(phibig, 5, 2)));
    if phibig < phi985 {
        // 2514: /// then text, /, text, /
        for _ in 0..3 {
            cout("");
            writeln(&mut unit10, "");
        }
        cout(" *** BEWARE *** ");
        writeln(&mut unit10, " *** BEWARE *** ");
        cout("");
        writeln(&mut unit10, "");
        cout(" Target value of phi may be too small");
        writeln(&mut unit10, " Target value of phi may be too small");
        cout("");
        writeln(&mut unit10, "");
    }
    // 2503: leading and trailing slashes
    writeln(&mut unit10, "");
    writeln(&mut unit10, " Calculation reached completion");
    writeln(&mut unit10, "");
    cout("");
    cout(" Calculation reached completion");
    cout("");
    unit10.flush().unwrap();

    hold(&mut stdin);
}
// @tangle:end frescalo__main

// @tangle:start frescalo__writeln
fn writeln<W: Write>(w: &mut W, s: &str) {
    w.write_all(s.as_bytes()).unwrap();
    w.write_all(b"\n").unwrap();
}
// @tangle:end frescalo__writeln

/// 2050 format with a real data row.
#[allow(clippy::too_many_arguments)]
// @tangle:start frescalo__write_trend
fn write_trend<W: Write>(
    w: &mut W,
    sp: &Name,
    tim: &Name,
    tf: f32,
    sd: f32,
    jtot: i32,
    spt: f32,
    est: f32,
    ic1: i32,
    ic2: i32,
) {
    let mut r = Rec::new();
    r.name(sp)
        .x(1)
        .name(tim)
        .f(tf, 8, 3)
        .f(sd, 7, 3)
        .i(jtot as i64, 7)
        .f(spt, 7, 1)
        .f(est, 7, 1)
        .i(ic1 as i64, 7)
        .i(ic2 as i64, 7);
    r.writeln(w);
}
// @tangle:end frescalo__write_trend

/// 2050 format with a zero row (species not recorded at that time).
// @tangle:start frescalo__write_zero_trend
fn write_zero_trend<W: Write>(w: &mut W, sp: &Name, tim: &Name) {
    write_trend(w, sp, tim, 0.0, 0.0, 0, 0.0, 0.0, 0, 0);
}
// @tangle:end frescalo__write_zero_trend

// ---------------------------------------------------------------------------
// fresca: calculate sampling effort and probabilities of species occurrence
// based on frequencies in the neighbourhood.
//
// The algorithm seeks a value of alpha, the Sampling Effort Multiplier, such
// that when frequencies are rescaled as ff(j) = 1-exp(-f(j)*alpha), the
// frequency-weighted mean frequency phi achieves its target value phibig.
//   m     - serial number of location
//   n     - total number of species
//   itot  - number of species at location m
//   jocc  - recorded/not recorded (1/0)
//   f     - input species frequencies (clamped in place)
//   ff    - workspace, returned holding rescaled frequencies
//   jrank - workspace for sorting to calculate rank order
// ---------------------------------------------------------------------------
#[allow(clippy::too_many_arguments)]
// @tangle:start frescalo__fresca
fn fresca(
    m: usize,
    n: usize,
    itot: i32,
    jocc: &[i32],
    f: &mut [f32],
    ff: &mut [f32],
    jrank: &mut [i32],
    samp1: &Name,
    splist: &[Name],
    phibig: f32,
    fmax: f32,
    fmin: f32,
    wn2: f32,
    tol: f32,
    irepmx: i32,
    unit7: &mut LogWriter,
    unit8: &mut LogWriter,
) -> (f32, f32) {
    let mut alpha: f32 = 1.0;
    let mut phi: f32 = 0.0;
    let mut phi1: f32 = 0.0;
    let mut spnum: f32 = 0.0;
    let mut spnum1: f32 = 0.0;
    let mut ir = irepmx;
    let mut converged = false;
    for iter in 1..=irepmx {
        ir = iter;
        for j in 1..=n {
            if f[j] > fmax {
                f[j] = fmax;
            }
            if f[j] < fmin {
                f[j] = fmin;
            }
            ff[j] = -(1.0 - f[j]).ln();
        }
        let mut tot: f32 = 0.0;
        let mut tot2: f32 = 0.0;
        for j in 1..=n {
            // new frequency after recording intensity multiplied by alpha
            let ffij = 1.0 - (-ff[j] * alpha).exp();
            tot += ffij;
            tot2 += ffij * ffij;
        }
        phi = tot2 / tot;
        spnum = tot;
        if iter < 20 {
            // successive approximation based on linear relation
            alpha = alpha * (1.86f32 * ((1.0 - phi).ln() - (1.0 - phibig).ln())).exp();
        } else {
            // crude successive approximation - slower with big datasets
            alpha = alpha * phibig / phi;
        }
        if iter == 1 {
            phi1 = phi;
            spnum1 = tot;
        }
        if (phi - phibig).abs() < tol {
            converged = true;
            break;
        }
    }
    // Fortran `DO ir=1,irepmx ... enddo` leaves ir = irepmx+1 after normal
    // (non-EXIT) completion: the loop variable is incremented and tested
    // *before* the loop is abandoned. Confirmed against gfortran with an
    // isolated repro of this exact do/goto pattern. This only affects the
    // reported Iter count on non-convergence, a documented failure signal.
    if !converged {
        ir = irepmx + 1;
    }

    if m == 1 {
        // 2001 header
        writeln(
            unit7,
            "Location  Loc_no  No_spp Phi_in  Alpha  Wgt_n2 Phi_out  Spnum_in Spnum_out Iter",
        );
    }
    let alph = if alpha > 999.99 { 999.99 } else { alpha };
    // 2002 format(a10,2i7,f7.3,1x,f6.2,1x,f7.2,1x,f7.3,2f10.1,i5)
    let mut r = Rec::new();
    r.name(samp1)
        .i(m as i64, 7)
        .i(itot as i64, 7)
        .f(phi1, 7, 3)
        .x(1)
        .f(alph, 6, 2)
        .x(1)
        .f(wn2, 7, 2)
        .x(1)
        .f(phi, 7, 3)
        .f(spnum1, 10, 1)
        .f(spnum, 10, 1)
        .i(ir as i64, 5);
    r.writeln(unit7);

    for j in 1..=n {
        ff[j] = -f[j] + (j as f32) * 1.0E-12;
        jrank[j] = j as i32;
    }
    sort2(ff, jrank, n);
    for j in 1..=n {
        let jj = jrank[j] as usize;
        let mut fij = f[jj];
        if fij > fmax {
            fij = fmax;
        }
        let ffij = 1.0 - (alpha * (1.0 - fij).ln()).exp();
        let sdfij = (fij * (1.0 - fij) / wn2).sqrt();
        let mut ffff = fij + sdfij;
        if 1.0 - ffff < 1.0E-12 {
            ffff = 1.0 - 1.0E-12;
        }
        let fffff = fij - sdfij;
        let ffsd = 1.0 - (alpha * (1.0 - ffff).ln()).exp();
        let fffsd = 1.0 - (alpha * (1.0 - fffff).ln()).exp();
        // sdij is an estimate of the standard error
        let sdij = 0.5 * (ffsd - fffsd);
        if j == 1 && m == 1 {
            // 2003 header
            writeln(
                unit8,
                "Location   Species    Pres  Freq__  Freq_1 SD_Frq1  Rank  Rank_1",
            );
        }
        ff[jj] = ffij;
        if f[jj] > 0.00005 {
            // 2004 format(a10,1x,a10,1x,i4,3f8.4,1x,i5,f8.3)
            let mut r = Rec::new();
            r.name(samp1)
                .x(1)
                .name(&splist[jj])
                .x(1)
                .i(jocc[jj] as i64, 4)
                .f(f[jj], 8, 4)
                .f(ffij, 8, 4)
                .f(sdij, 8, 4)
                .x(1)
                .i(j as i64, 5)
                .f((j as f32) / spnum, 8, 3);
            r.writeln(unit8);
        }
    }
    (phi1, spnum)
}
// @tangle:end frescalo__fresca

// ---------------------------------------------------------------------------
// tfcalc: calculate a time factor tf for a species at a time, given the
// observed species total for that time; sd is its standard error.
//   iocc(i)   - 1 if the species is found at location i at the time, else 0
//   smpint(i) - sampling intensity at location i and time t
//   fff(i)    - smoothed time-independent frequency of the species at i
//   ic1       - number of samples with nonzero probability of occurrence
//   ic2       - number of cases where smpint*fff > 0.98
// Weights downweight cases where smpint < 0.1 (no systematic sampling).
// Returns (tf, sd, sptot, jtot, esttot, ic1, ic2).
// ---------------------------------------------------------------------------
// @tangle:start frescalo__tfcalc
fn tfcalc(iocc: &[i32], smpint: &[f32], fff: &[f32], m: usize) -> (f32, f32, f32, i32, f32, i32, i32) {
    const KMAX: i32 = 100;
    let mut tf: f32 = 1.0;
    let mut esttot: f32 = 0.0;
    let mut estvar: f32 = 0.0;
    let mut sptot: f32 = 0.0;
    let mut jtot: i32 = 0;
    let mut ic1: i32 = 0;
    let mut ic2: i32 = 0;
    for _ in 1..=KMAX {
        esttot = 0.0;
        estvar = 0.0;
        sptot = 0.0;
        jtot = 0;
        ic1 = 0;
        ic2 = 0;
        for i in 1..=m {
            let mut wgt: f32 = 1.0;
            if smpint[i] < 0.0995 {
                wgt = 10.0 * smpint[i] + 0.005;
            }
            // probability of finding the species, times the sampling intensity
            let mut pfac = smpint[i] * fff[i];
            if pfac > 0.0 {
                ic1 += 1;
            }
            if pfac > 0.98 {
                pfac = 0.98;
                ic2 += 1;
            }
            // fudge to allow taking logs; otherwise there is a danger pfac=1
            let plog = -(1.0 - pfac).ln();
            let estval = 1.0 - (-plog * tf).exp();
            esttot += wgt * estval;
            estvar += wgt * wgt * estval * (1.0 - estval);
            sptot += wgt * (iocc[i] as f32);
            jtot += iocc[i];
        }
        if (sptot - esttot).abs() < 0.0005 {
            break;
        }
        tf = tf * sptot / (esttot + 0.0000001);
    }

    // sptot1 is precisely 1 standard deviation bigger than sptot;
    // recalculate with this target value to obtain the standard error of tf.
    let sptot1 = sptot + estvar.sqrt();
    let mut tf1 = tf;
    let mut esttt1: f32;
    for _ in 1..=KMAX {
        esttt1 = 0.0;
        for i in 1..=m {
            let mut wgt: f32 = 1.0;
            if smpint[i] < 0.0995 {
                wgt = 10.0 * smpint[i] + 0.005;
            }
            let mut pfac = smpint[i] * fff[i];
            if pfac > 0.98 {
                pfac = 0.98;
            }
            let plog = -(1.0 - pfac).ln();
            let estval = 1.0 - (-plog * tf1).exp();
            esttt1 += wgt * estval;
        }
        if (sptot1 - esttt1).abs() < 0.0005 {
            break;
        }
        tf1 = tf1 * sptot1 / (esttt1 + 0.0000001);
    }

    let sd = tf1 - tf;
    (tf, sd, sptot, jtot, esttot, ic1, ic2)
}
// @tangle:end frescalo__tfcalc
