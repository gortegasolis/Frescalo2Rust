//! Rust port of neighsim_1.f:
//! NEIGHSIM - Neighbourhood similarity based on training-set species and
//! physical proximity
//! written by Mark Hill, January-June 2011
//!
//! Calculates floristic similarity between samples from a training set,
//! combined with physical proximity ranks, and writes neighbourhood weights
//! for use in FRESCALO.

use frescalo::*;
use std::io::{self, BufWriter, Write};

const MM: usize = 4000;
const NN: usize = 10000;
const NNDAT: usize = 5_000_000;
const BIG: f32 = 1_000_000.0;
const SMALL: f32 = 0.00005;

// @tangle:start neighsim__main
fn main() {
    let mut stdin = io::stdin().lock();

    let mut sa = vec![blank_name(); MM + 2];
    let mut sp = vec![blank_name(); NN + 2];
    let mut itot = vec![0i32; MM + 2];
    let mut iocc = vec![0i32; MM + 2];
    let mut index = vec![0i32; MM + 2];
    let mut jtot = vec![0i32; NN + 2];
    let mut sim = vec![0.0f32; MM + 2];
    let mut simil: Arr2<f32> = Arr2::new(MM + 1, MM + 1, 0.0);
    let mut iseqq: Arr2<i32> = Arr2::new(MM + 1, MM + 1, 0);
    let mut ddji = vec![[b' '; 30]; NNDAT + 2];
    let mut distii = vec![[b' '; 30]; NNDAT + 2];

    let mut m: usize = 0;
    let mut n: usize = 0;
    let mut ndist: usize = 0;
    let mut ndata: usize = 0;

    // Set up files for reading and writing
    cout("");
    cout(" NEIGHSIM - Neighbourhood similarity based on training-set species and physical proximity");
    cout(" written by Mark Hill, January-June 2011");
    cout("");
    // 1x,'NOTE PARAMETER LIMITS: Sites',i5,' Species',i6,1x,'Training-set number of records',i9
    cout(&format!(
        " NOTE PARAMETER LIMITS: Sites{} Species{} Training-set number of records{}",
        ifmt(MM as i64, 5),
        ifmt(NN as i64, 6),
        ifmt(NNDAT as i64, 9)
    ));
    cout("");
    cout(" Type name of input file with Training-set species data [sample species] ....");
    let (_f1, ftr) = filin(&mut stdin);
    let mut train_reader = DataReader::new(ftr);
    cout(" Type name of input file with physical distances ....");
    let (_f2, fdi) = filin(&mut stdin);
    let mut dist_reader = DataReader::new(fdi);
    cout(" Type name of Training-set similarity output file ...");
    let (_f3, fsim) = filout(&mut stdin);
    let mut unit9 = BufWriter::new(fsim);
    cout(" Type name of weights output file for use in Frescalo ...");
    let (_f4, fwgt) = filout(&mut stdin);
    let mut unit8 = BufWriter::new(fwgt);
    cout(" Type number of neighbours to include ...");
    let neigh = read_int_listdirected(&mut stdin);

    let mut samp = blank_name();
    let mut samp1 = blank_name();
    let mut spec = blank_name();
    let mut any = blank_name();

    // Read in distance data
    loop {
        if !dist_reader.getd(&mut samp, &mut samp1, &mut any) {
            break;
        }
        ndist += 1;
        if ndist > NNDAT {
            ld_line(&format!(
                " Too many data items in physical distance file - limit is{}",
                ld_i(NNDAT as i64)
            ));
            hold(&mut stdin);
        }
        if ndist % 20000 == 0 {
            ld_line(&format!(
                "{} Dist {}{}{}",
                ld_i(ndist as i64),
                name_to_string(&samp),
                name_to_string(&samp1),
                name_to_string(&any)
            ));
        }
        distii[ndist] = make_rec30(&samp, &samp1, &any);
    }

    // Read in training-set species data
    loop {
        if !train_reader.getd(&mut samp, &mut spec, &mut any) {
            break;
        }
        addwrd_guarded(&mut sa, MM, &mut m, &samp);
        if m > MM {
            ld_line(&format!(" Too many samples - limit is{}", ld_i(MM as i64)));
            hold(&mut stdin);
        }
        addwrd_guarded(&mut sp, NN, &mut n, &spec);
        if n > NN {
            ld_line(&format!(" Too many species - limit is{}", ld_i(NN as i64)));
            hold(&mut stdin);
        }
        ndata += 1;
        if ndata > NNDAT {
            ld_line(&format!(" Too many data items - limit is{}", ld_i(NNDAT as i64)));
            hold(&mut stdin);
        }
        if ndata % 20000 == 0 {
            ld_line(&format!(
                "{} Spdata {}{}{}",
                ld_i(ndata as i64),
                name_to_string(&samp),
                name_to_string(&spec),
                name_to_string(&any)
            ));
        }
        ddji[ndata] = make_rec30(&spec, &samp, &any);
    }

    cout(" Sorting main data ...");
    sort30(&mut ddji, ndata);
    cout(" Sort completed");

    for idata in 1..=ndata {
        if idata % 20000 == 0 {
            ld_line(&format!(" Calculating totals{}", ld_i(idata as i64)));
        }
        let dji = ddji[idata];
        let f_spec = rec30_field(&dji, 0);
        let f_samp = rec30_field(&dji, 1);
        let i = binfnd(&sa, m, &f_samp);
        let j = binfnd(&sp, n, &f_spec);
        itot[i] += 1;
        jtot[j] += 1;
    }

    // Now start calculating similarity
    let mut idata = 0usize;
    for j in 1..=n {
        for i in 1..=m {
            iocc[i] = 0;
        }
        for _ in 1..=jtot[j] {
            idata += 1;
            if idata % 20000 == 0 {
                ld_line(&format!(" Similarities{}", ld_i(idata as i64)));
            }
            let dji = ddji[idata];
            let f_spec = rec30_field(&dji, 0);
            let f_samp = rec30_field(&dji, 1);
            if f_spec != sp[j] {
                ld_line(&format!(
                    "Unequal species{}{}",
                    name_to_string(&f_spec),
                    name_to_string(&sp[j])
                ));
                hold(&mut stdin);
            }
            let i = binfnd(&sa, m, &f_samp);
            iocc[i] = -(i as i32);
        }
        isort(&mut iocc, m);
        // mcc is the length of nonzero items in iocc
        let mut mcc = m;
        for ic in 1..=m {
            if iocc[ic] == 0 {
                mcc = ic - 1;
                break;
            }
        }
        for icc1 in 1..=mcc {
            let i1 = (-iocc[icc1]) as usize;
            for icc2 in 1..=mcc {
                let i2 = (-iocc[icc2]) as usize;
                // number of species in common; later divided to get similarity
                simil.add(i1, i2, 1.0);
            }
        }
    }

    // Multiply those within preferred region by big
    let mut neigh1: i32 = 0;
    for idist in 1..=ndist {
        if idist % 20000 == 0 {
            // NB: the original prints samp/samp1 from the *previous* record here.
            ld_line(&format!(
                "{} Dist {}{}",
                ld_i(idist as i64),
                name_to_string(&samp),
                name_to_string(&samp1)
            ));
        }
        let dji = distii[idist];
        samp = rec30_field(&dji, 0);
        samp1 = rec30_field(&dji, 1);
        let i1 = binfnd(&sa, m, &samp);
        let i2 = binfnd(&sa, m, &samp1);
        if i1 != 0 && i2 != 0 {
            simil.set(i1, i2, simil.at(i1, i2) * BIG);
            let anyf = rec30_field(&dji, 2);
            let seq = getnum(&anyf);
            let iseq = seq as i32; // ifix: truncation towards zero
            iseqq.set(i1, i2, iseq);
            if neigh1 < iseq {
                neigh1 = iseq;
            }
        }
    }
    if neigh1 == 0 {
        // list-directed write to the weights file: leading blank
        unit8.write_all(b" Unrecognized sample names in distance data\n").unwrap();
    }

    for i1 in 1..=m {
        if i1 % 100 == 0 {
            ld_line(&format!("Writing output  {}{}", name_to_string(&sa[i1]), ld_i(i1 as i64)));
        }
        for i2 in 1..=m {
            sim[i2] = simil.at(i1, i2) * 2.0 / (itot[i1] + itot[i2]) as f32;
            index[i2] = i2 as i32;
        }
        sort2(&mut sim, &mut index, m);
        'neigh: for is2 in 1..=(neigh.max(0) as usize) {
            let i2 = m as isize - is2 as isize + 1;
            if i2 < 1 {
                break 'neigh;
            }
            let i2 = i2 as usize;
            let iis2 = index[i2] as usize;
            // 2030 format(2a10,i5,f6.3)
            let mut r = Rec::new();
            r.name(&sa[i1])
                .name(&sa[iis2])
                .i(is2 as i64, 5)
                .f(sim[i2] / BIG, 6, 3);
            r.writeln(&mut unit9);
            if neigh1 == 0 {
                continue 'neigh;
            }
            let t = (is2 as f32 - 1.0) / neigh as f32;
            let amult1 = (1.0 - t * t).powi(4);
            let t2 = (iseqq.at(i1, iis2) - 1) as f32 / neigh1 as f32;
            let amult2 = (1.0 - t2 * t2).powi(4);
            let amult = amult1 * amult2;
            if amult > SMALL {
                // 2031 format(2a10,3f7.4,2i6)
                let mut r = Rec::new();
                r.name(&sa[i1])
                    .name(&sa[iis2])
                    .f(amult, 7, 4)
                    .f(amult1, 7, 4)
                    .f(amult2, 7, 4)
                    .i(neigh as i64, 6)
                    .i(neigh1 as i64, 6);
                r.writeln(&mut unit8);
            }
        }
    }
    ld_line(&format!("neigh,neigh1{}{}", ld_i(neigh as i64), ld_i(neigh1 as i64)));

    unit8.flush().unwrap();
    unit9.flush().unwrap();
    hold(&mut stdin);
}
// @tangle:end neighsim__main
