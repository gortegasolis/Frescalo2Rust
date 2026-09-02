//! Rust port of sampdist_1.f:
//! SAMPDIST - Sample distances in neighbourhoods
//! written by Mark Hill, January-June 2011
//!
//! Calculates the Euclidean distances between sample locations (treated as
//! easting/northing) and writes, for each location, the nearest `neigh`
//! locations (including itself) in increasing distance order.

use frescalo::*;
use std::io::{self, BufWriter, Write};

const MM: usize = 400000;

// @tangle:start sampdist__main
fn main() {
    let mut stdin = io::stdin().lock();

    let mut sa = vec![blank_name(); MM + 2];
    let mut aeast = vec![0.0f32; MM + 2];
    let mut anorth = vec![0.0f32; MM + 2];
    let mut dist = vec![0.0f32; MM + 2];
    let mut index = vec![0i32; MM + 2];

    let mut m: usize = 0;

    // Set up files for reading and writing
    cout("");
    cout(" SAMPDIST - Sample distances in neighbourhoods");
    cout(" written by Mark Hill, January-June 2011");
    cout("");
    cout(" Type name of file with locations ....");
    let (_filein, fin) = filin(&mut stdin);
    let mut reader = DataReader::new(fin);
    cout(" Type name of output file with neighbourhood distances ...");
    let (_fileou, fout) = filout(&mut stdin);
    let mut unit9 = BufWriter::new(fout);
    cout(" Type number of neighbours to include ...");
    let neigh = read_int_listdirected(&mut stdin);

    // Set up index of samples
    let mut samp = blank_name();
    let mut east = blank_name();
    let mut north = blank_name();
    loop {
        if !reader.getd(&mut samp, &mut east, &mut north) {
            break;
        }
        addwrd(&mut sa, &mut m, &samp);
        if m % 100 == 0 {
            ld_line(&format!("{}  Sample  {}", ld_i(m as i64), name_to_string(&samp)));
        }
    }

    reader.rewind();

    loop {
        if !reader.getd(&mut samp, &mut east, &mut north) {
            break;
        }
        let i = binfnd(&sa, m, &samp);
        aeast[i] = getnum(&east);
        anorth[i] = getnum(&north);
        if i % 100 == 0 {
            ld_line(&format!(
                "{}   {}{}{}",
                ld_i(i as i64),
                name_to_string(&sa[i]),
                ld_f(aeast[i]),
                ld_f(anorth[i])
            ));
        }
    }

    // Now start calculating distances
    for i1 in 1..=m {
        if i1 % 100 == 0 {
            ld_line(&format!("Calculating distances   {}", ld_i(m as i64)));
        }
        for i2 in 1..=m {
            let de = aeast[i1] - aeast[i2];
            let dn = anorth[i1] - anorth[i2];
            dist[i2] = (de * de + dn * dn).sqrt();
            index[i2] = i2 as i32;
        }
        sort2(&mut dist, &mut index, m);
        // (the original would read past the sorted arrays if neigh > m)
        let nout = (neigh.max(0) as usize).min(m);
        for is2 in 1..=nout {
            let iis2 = index[is2] as usize;
            // 2030 format(2a10,i5,1x,f6.0)
            let mut r = Rec::new();
            r.name(&sa[i1])
                .name(&sa[iis2])
                .i(is2 as i64, 5)
                .x(1)
                .f(dist[is2], 6, 0);
            r.writeln(&mut unit9);
        }
    }

    unit9.flush().unwrap();
    hold(&mut stdin);
}
// @tangle:end sampdist__main
