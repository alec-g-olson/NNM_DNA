//! Regression guard: for perfect Watson-Crick duplexes (no bulges, no
//! mismatches, no loops), `aligned_duplex_thermo` must agree exactly with
//! `duplex_thermo`. The two code paths share no logic — `duplex_thermo`
//! iterates `pair_at` positions, `aligned_duplex_thermo` walks a Feature
//! list — so this is a strong end-to-end check that the new pathway didn't
//! introduce drift in the WC fast path.

use n3dna::{aligned_duplex_thermo, duplex_thermo, AlignedDuplex, Duplex, Sequence};

fn assert_close(a: f64, b: f64, label: &str) {
    let diff = (a - b).abs();
    assert!(
        diff < 1e-9,
        "{label}: aligned = {a}, duplex = {b}, diff = {diff}"
    );
}

fn check(seq_str: &str) {
    let top = Sequence::parse(seq_str).unwrap();
    let bottom = top.reverse_complement();

    let from_duplex = duplex_thermo(&Duplex::perfect(top.clone())).unwrap();
    let aligned = AlignedDuplex::align(top, bottom).unwrap();
    let from_aligned = aligned_duplex_thermo(&aligned).unwrap();

    assert_close(from_aligned.dh, from_duplex.dh, &format!("ΔH for {seq_str}"));
    assert_close(from_aligned.ds, from_duplex.ds, &format!("ΔS for {seq_str}"));
    assert_eq!(
        from_aligned.self_comp, from_duplex.self_comp,
        "self_comp mismatch for {seq_str}"
    );
}

#[test]
fn cgttga_matches() {
    // Non-self-complementary 6-mer with one of each base
    check("CGTTGA");
}

#[test]
fn cgcgcg_matches() {
    // Self-complementary palindrome — symmetry correction must agree
    check("CGCGCG");
}

#[test]
fn aaaa_matches() {
    // Edge case: all AT initiation, no GC anywhere
    check("AAAA");
}

#[test]
fn long_random_like_matches() {
    // Longer, no special structure
    check("ACGTACGTACGTACGT");
}

#[test]
fn dimer_matches() {
    // Smallest possible duplex
    check("AT");
}

#[test]
fn all_g_c_matches() {
    // All GC, both initiations are G·C
    check("GCGCGC");
}
