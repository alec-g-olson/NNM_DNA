//! Reaction-conditions and equilibrium tests for `ReactionConditions`,
//! `Thermo::equilibrium`, and `Thermo::tm_at_conditions`.
//!
//! Reaction model: A + B ⇌ AB for non-self-complementary, 2A ⇌ AA for
//! self-complementary. Tm defined as the temperature where 50% of the
//! limiting strand is in the duplex form.

use n3dna::{
    duplex_thermo, Duplex, ReactionConditions, SaltCorrection, Sequence, Thermo,
};

fn close(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() < tol
}

fn thermo_of(seq: &str) -> Thermo {
    duplex_thermo(&Duplex::perfect(Sequence::parse(seq).unwrap())).unwrap()
}

#[test]
fn default_conditions_match_primer3_defaults() {
    let c = ReactionConditions::default();
    assert!(close(c.strand_a_conc_mol, 50e-9, 1e-15));
    assert!(close(c.strand_b_conc_mol, 50e-9, 1e-15));
    assert!(close(c.mv_conc_mol, 0.050, 1e-9));
    assert!(close(c.dv_conc_mol, 0.0015, 1e-9));
    assert!(close(c.dntp_conc_mol, 0.0006, 1e-9));
    assert_eq!(c.salt_correction, SaltCorrection::SantaLucia);
}

#[test]
fn free_mg_subtracts_dntp_and_clamps_to_zero() {
    let c = ReactionConditions {
        dv_conc_mol: 0.0015,
        dntp_conc_mol: 0.0006,
        ..Default::default()
    };
    assert!(close(c.free_mg_conc(), 0.0009, 1e-9));

    let chelated = ReactionConditions {
        dv_conc_mol: 0.0005,
        dntp_conc_mol: 0.002, // dNTPs exceed Mg → all chelated
        ..Default::default()
    };
    assert_eq!(chelated.free_mg_conc(), 0.0);
    assert!(!chelated.has_divalent());

    let no_mg = ReactionConditions {
        dv_conc_mol: 0.0,
        dntp_conc_mol: 0.0,
        ..Default::default()
    };
    assert!(!no_mg.has_divalent());
}

#[test]
fn equilibrium_at_low_t_mostly_bound() {
    // Cold → nearly complete duplex formation. Note: even at 250 K, finite
    // micromolar concentrations cap saturation around ~99.6% for a 6-mer
    // (this is real biophysics, not a bug — K is large but ⟨x⟩ ≈ C·(1 −
    // 1/√(4KC)), so saturation approaches 1 only asymptotically).
    let t = thermo_of("CGTTGA");
    let c = ReactionConditions {
        mv_conc_mol: 1.0,
        dv_conc_mol: 0.0,
        dntp_conc_mol: 0.0,
        strand_a_conc_mol: 1e-6,
        strand_b_conc_mol: 1e-6,
        ..Default::default()
    };
    let eq = t.equilibrium(250.0, &c);
    assert!(eq.frac_a_bound > 0.99, "frac_a_bound = {}", eq.frac_a_bound);
    assert!(eq.ab_bound_mol > 0.99e-6);
}

#[test]
fn equilibrium_at_high_t_all_free() {
    // Very hot → no duplex.
    let t = thermo_of("CGTTGA");
    let c = ReactionConditions {
        mv_conc_mol: 1.0,
        strand_a_conc_mol: 1e-6,
        strand_b_conc_mol: 1e-6,
        ..Default::default()
    };
    let eq = t.equilibrium(380.0, &c);
    assert!(eq.frac_a_bound < 0.001, "frac_a_bound = {}", eq.frac_a_bound);
    assert!(eq.ab_bound_mol < 1e-9);
}

#[test]
fn equilibrium_tm_equal_strands_matches_closed_form() {
    // At 1 M Na⁺ with equal strand concentrations and no Mg, the
    // equilibrium-derived Tm should match the closed-form `tm(C_T)` from
    // the existing API. C_T = strand_a + strand_b.
    let t = thermo_of("CGTTGA");
    let c = ReactionConditions {
        mv_conc_mol: 1.0,
        dv_conc_mol: 0.0,
        dntp_conc_mol: 0.0,
        strand_a_conc_mol: 50e-9,
        strand_b_conc_mol: 50e-9,
        salt_correction: SaltCorrection::SantaLucia,
        ..Default::default()
    };
    let tm_equilibrium = t.tm_at_conditions(&c);
    let tm_closed = t.tm(100e-9); // C_T = 100 nM
    assert!(
        close(tm_equilibrium, tm_closed, 1e-3),
        "equilibrium = {tm_equilibrium}, closed = {tm_closed}"
    );
}

#[test]
fn excess_strand_a_raises_b_binding() {
    // Same Tm experiment: at temperature near the limiting Tm, the
    // less-abundant strand should be MORE bound than the more-abundant one
    // (because the excess provides extra binding partners).
    //
    // Note: must use a NON-palindromic sequence — for self-complementary
    // duplexes both "strands" are the same species, and frac_a == frac_b
    // by definition. AAAGCGTACG ≠ rev_comp(AAAGCGTACG) = CGTACGCTTT.
    let t = thermo_of("AAAGCGTACG");
    assert!(!t.self_comp, "test sequence must NOT be self-complementary");
    let c_equal = ReactionConditions {
        mv_conc_mol: 1.0,
        dv_conc_mol: 0.0,
        dntp_conc_mol: 0.0,
        strand_a_conc_mol: 1e-7,
        strand_b_conc_mol: 1e-7,
        ..Default::default()
    };
    let c_excess_a = ReactionConditions {
        strand_a_conc_mol: 1e-5, // 100x excess
        strand_b_conc_mol: 1e-7,
        ..c_equal
    };

    // At intermediate temperature, with 100x excess A, frac_b_bound should
    // be much higher than frac_a_bound.
    let t_test = 310.0;
    let eq_excess = t.equilibrium(t_test, &c_excess_a);
    assert!(
        eq_excess.frac_b_bound > eq_excess.frac_a_bound,
        "frac_b = {}, frac_a = {}",
        eq_excess.frac_b_bound,
        eq_excess.frac_a_bound,
    );

    // Tm with excess strand A should be HIGHER than with equal strands
    // (more excess → easier for B to find a partner → duplex stable to higher T).
    let tm_equal = t.tm_at_conditions(&c_equal);
    let tm_excess = t.tm_at_conditions(&c_excess_a);
    assert!(
        tm_excess > tm_equal + 1.0,
        "tm_excess = {tm_excess}, tm_equal = {tm_equal}"
    );
}

#[test]
fn mg_raises_tm_relative_to_no_mg() {
    // Adding Mg²⁺ stabilizes the duplex → higher Tm.
    let t = thermo_of("CGTACGTACG");
    let no_mg = ReactionConditions {
        mv_conc_mol: 0.050,
        dv_conc_mol: 0.0,
        dntp_conc_mol: 0.0,
        ..Default::default()
    };
    let with_mg = ReactionConditions {
        dv_conc_mol: 0.005, // 5 mM Mg²⁺ — well above dNTP
        dntp_conc_mol: 0.0,
        ..no_mg
    };
    let tm_no_mg = t.tm_at_conditions(&no_mg);
    let tm_mg = t.tm_at_conditions(&with_mg);
    assert!(
        tm_mg > tm_no_mg,
        "tm_mg = {tm_mg}, tm_no_mg = {tm_no_mg}"
    );
}

#[test]
fn dntp_chelation_reduces_mg_effect() {
    // Same Mg²⁺ concentration; more dNTPs → less effective Mg → lower Tm.
    let t = thermo_of("CGTACGTACG");
    let base = ReactionConditions {
        mv_conc_mol: 0.050,
        dv_conc_mol: 0.003, // 3 mM Mg
        ..Default::default()
    };
    let low_dntp = ReactionConditions {
        dntp_conc_mol: 0.0005,
        ..base
    };
    let high_dntp = ReactionConditions {
        dntp_conc_mol: 0.0025, // closer to Mg → most chelated
        ..base
    };
    let tm_low = t.tm_at_conditions(&low_dntp);
    let tm_high = t.tm_at_conditions(&high_dntp);
    assert!(
        tm_low > tm_high,
        "low dNTP Tm = {tm_low}, high dNTP Tm = {tm_high}"
    );
}

#[test]
fn self_complementary_equilibrium_uses_2a_ab() {
    // CGCGCG is a palindrome. For self-complementary, A_tot = a + b (same species).
    // At Tm, [AA] = A_tot/4 and fraction-bound = 0.5.
    let t = thermo_of("CGCGCG");
    assert!(t.self_comp);

    let c = ReactionConditions {
        mv_conc_mol: 1.0,
        dv_conc_mol: 0.0,
        dntp_conc_mol: 0.0,
        strand_a_conc_mol: 50e-9,
        strand_b_conc_mol: 50e-9,
        ..Default::default()
    };
    let tm_c = t.tm_at_conditions(&c);
    let tm_k = tm_c + 273.15;
    let eq = t.equilibrium(tm_k, &c);
    // At Tm, half the strands are bound. A_tot = 100 nM, so [AA] = 25 nM.
    let a_tot = 100e-9;
    assert!(close(eq.ab_bound_mol, a_tot / 4.0, 1e-9), "ab = {}", eq.ab_bound_mol);
    assert!(close(eq.frac_a_bound, 0.5, 1e-3));

    // And compare to the closed-form Tm for self-comp (1 M Na⁺ baseline).
    let tm_closed = t.tm(100e-9);
    assert!(
        close(tm_c, tm_closed, 1e-3),
        "equilibrium = {tm_c}, closed = {tm_closed}"
    );
}

#[test]
fn tm_at_50c_default_conditions_runs_cleanly() {
    // Smoke test: default conditions should produce a sane Tm for a
    // typical PCR-primer-like sequence.
    let t = thermo_of("AGTCGTACGTAGCATGC");
    let c = ReactionConditions::default();
    let tm = t.tm_at_conditions(&c);
    // Loose bounds — just confirming we get a physical value
    assert!(tm > 30.0 && tm < 90.0, "Tm = {tm}");
}
