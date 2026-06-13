//! Thermodynamic summation over a duplex and downstream quantities (ΔG, Tm).

use crate::align::{AlignError, AlignedDuplex};
use crate::base::Base;
use crate::duplex::Duplex;
use crate::feature::{extract_features, Feature};
use crate::params::{nn_params, INIT_AT, INIT_GC, SYMMETRY_CORRECTION};
use crate::params_loop::{bulge_length, internal_loop_length};
use crate::params_mismatch::mismatch_params;
use std::fmt;

/// Ideal-gas constant in cal/(mol·K).
pub const R: f64 = 1.987;

/// Aggregated thermodynamic parameters for a duplex.
///
/// `dh` is in kcal/mol; `ds` is in cal/(mol·K). `self_comp` records whether
/// the duplex is self-complementary, which affects the Tm formula's
/// concentration term.
///
/// `n_bp` and `gc_fraction` are populated to support the salt-correction
/// formulas in [`Thermo::tm_with_salt`]:
/// - `n_bp` is the number of base pairs in the duplex (count of paired
///   columns, ignoring bulges and loops). The SantaLucia salt correction uses
///   `n_bp - 1` as the per-strand phosphate count.
/// - `gc_fraction` is the fraction of paired columns that are G·C or C·G.
///   The Owczarzy 2004 correction is GC-dependent.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Thermo {
    pub dh: f64,
    pub ds: f64,
    pub self_comp: bool,
    pub n_bp: usize,
    pub gc_fraction: f64,
}

/// Salt-correction method for Tm prediction. Names match primer3-py's
/// `salt_corrections_method` parameter; `SantaLucia` is the default.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SaltCorrection {
    /// SantaLucia (1998), *PNAS* 95:1460. Adjusts ΔS° before recomputing Tm:
    ///   ΔS°(salt) = ΔS°(1 M) + 0.368 · (n_bp − 1) · ln([Na⁺]).
    /// Same author as the unified NN parameters this crate is built on.
    #[default]
    SantaLucia,
    /// Schildkraut & Lifson (1965), *Biopolymers* 3:195. Applies an additive
    /// offset to Tm:  Tm(salt) = Tm(1 M) + 16.6 · log₁₀([Na⁺]).
    Schildkraut,
    /// Owczarzy et al. (2004), *Biochemistry* 43:3537. Applies a sequence-
    /// GC-dependent correction on 1/Tm (Kelvin):
    ///   1/Tm(salt) = 1/Tm(1 M) + (4.29·fGC − 3.95)·10⁻⁵·ln[Na⁺]
    ///                          + 9.40·10⁻⁶·(ln[Na⁺])².
    Owczarzy,
}

impl Thermo {
    /// ΔG° at temperature `t_k` (Kelvin), in kcal/mol.
    pub fn delta_g(self, t_k: f64) -> f64 {
        self.dh - t_k * self.ds / 1000.0
    }

    /// Convenience: ΔG° at 37 °C, in kcal/mol.
    pub fn delta_g_37(self) -> f64 {
        self.delta_g(310.15)
    }

    /// Melting temperature in °C for total strand concentration `ct` (mol/L)
    /// at the reference 1 M [Na⁺] of the unified NN parameter set. For Tm at
    /// other salt concentrations use [`Thermo::tm_with_salt`].
    pub fn tm(self, ct: f64) -> f64 {
        let x = if self.self_comp { 1.0 } else { 4.0 };
        let dh_cal = self.dh * 1000.0;
        dh_cal / (self.ds + R * (ct / x).ln()) - 273.15
    }

    /// Melting temperature in °C with monovalent-cation salt correction.
    ///
    /// `ct` is total strand concentration (mol/L); `mv_conc_mol` is the
    /// monovalent cation concentration (mol/L; e.g. `0.050` for 50 mM Na⁺);
    /// `method` selects which correction formula to apply. At `mv_conc_mol`
    /// = 1.0 (1 M) all three methods reduce to [`Thermo::tm`].
    pub fn tm_with_salt(self, ct: f64, mv_conc_mol: f64, method: SaltCorrection) -> f64 {
        match method {
            SaltCorrection::SantaLucia => {
                let n = (self.n_bp as f64 - 1.0).max(0.0);
                let ds_corrected = self.ds + 0.368 * n * mv_conc_mol.ln();
                let x = if self.self_comp { 1.0 } else { 4.0 };
                let dh_cal = self.dh * 1000.0;
                dh_cal / (ds_corrected + R * (ct / x).ln()) - 273.15
            }
            SaltCorrection::Schildkraut => self.tm(ct) + 16.6 * mv_conc_mol.log10(),
            SaltCorrection::Owczarzy => {
                let tm_1m_k = self.tm(ct) + 273.15;
                let ln_na = mv_conc_mol.ln();
                let inv_tm = 1.0 / tm_1m_k
                    + (4.29 * self.gc_fraction - 3.95) * 1e-5 * ln_na
                    + 9.40e-6 * ln_na * ln_na;
                1.0 / inv_tm - 273.15
            }
        }
    }
}

/// Full reaction conditions for a hybridization experiment.
///
/// Defaults are primer3-py's PCR-like defaults: each strand 50 nM, 50 mM Na⁺,
/// 1.5 mM Mg²⁺, 0.6 mM dNTP, SantaLucia salt correction. Strand concentrations
/// are independent — the two oligos can be at different concentrations
/// (e.g. probe-target hybridization or PCR with primer-template stoichiometry).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ReactionConditions {
    /// Total concentration of strand A (the "top" strand), in mol/L.
    pub strand_a_conc_mol: f64,
    /// Total concentration of strand B (the "bottom" strand), in mol/L.
    pub strand_b_conc_mol: f64,
    /// Monovalent cation concentration (Na⁺/K⁺), in mol/L.
    pub mv_conc_mol: f64,
    /// Divalent cation concentration (Mg²⁺), in mol/L.
    pub dv_conc_mol: f64,
    /// dNTP concentration, in mol/L. Each dNTP chelates one Mg²⁺, reducing
    /// the effective divalent concentration available for duplex stabilization.
    pub dntp_conc_mol: f64,
    /// Which salt correction formula to apply when `dv_conc_mol == 0`. When
    /// Mg²⁺ is present, the Owczarzy 2008 ratio-regime algorithm is used
    /// regardless of this setting.
    pub salt_correction: SaltCorrection,
}

impl Default for ReactionConditions {
    fn default() -> Self {
        Self {
            strand_a_conc_mol: 50e-9,
            strand_b_conc_mol: 50e-9,
            mv_conc_mol: 0.050,
            dv_conc_mol: 0.0015,
            dntp_conc_mol: 0.0006,
            salt_correction: SaltCorrection::SantaLucia,
        }
    }
}

impl ReactionConditions {
    /// Free Mg²⁺ concentration after dNTP chelation: max(0, [Mg²⁺] − [dNTP]).
    pub fn free_mg_conc(&self) -> f64 {
        (self.dv_conc_mol - self.dntp_conc_mol).max(0.0)
    }

    /// True if Mg²⁺ is significant enough to trigger the Owczarzy 2008
    /// divalent-aware salt correction (vs the monovalent menu).
    pub fn has_divalent(&self) -> bool {
        self.free_mg_conc() > 1e-12
    }
}

/// Result of solving the A + B ⇌ AB equilibrium at a given temperature.
///
/// All concentrations are in mol/L. Bound fractions are dimensionless;
/// `frac_a_bound = ab_bound_mol / A_total` and similarly for B.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Equilibrium {
    pub a_free_mol: f64,
    pub b_free_mol: f64,
    pub ab_bound_mol: f64,
    pub frac_a_bound: f64,
    pub frac_b_bound: f64,
}

impl Thermo {
    /// Solve the strand-pairing equilibrium at temperature `t_k` (Kelvin)
    /// with the given reaction conditions.
    ///
    /// For **non-self-complementary** duplexes this is the standard
    /// `A + B ⇌ AB` quadratic with conservation `[A] + [AB] = A_tot` and
    /// `[B] + [AB] = B_tot`.
    ///
    /// For **self-complementary** duplexes the two strands are the same
    /// chemical species, so the equilibrium is `2A ⇌ AA` with conservation
    /// `[A] + 2[AA] = A_tot`. `A_tot` is taken as `strand_a + strand_b`
    /// (both inputs treated as concentrations of the single species).
    ///
    /// Uses the salt-corrected ΔS° (SantaLucia 1998) for `K(T)`. Other salt
    /// correction methods adjust Tm post-hoc and aren't reflected in this
    /// `K(T)` — use [`Thermo::tm_at_conditions`] for those.
    pub fn equilibrium(self, t_k: f64, conditions: &ReactionConditions) -> Equilibrium {
        let mv_eff = effective_mv_conc(conditions);
        let n = (self.n_bp as f64 - 1.0).max(0.0);
        let ds_corr = self.ds + 0.368 * n * mv_eff.ln();
        // K = exp(-ΔG/(R·T)) = exp(ΔS/R - ΔH·1000/(R·T))
        let k = (ds_corr / R - self.dh * 1000.0 / (R * t_k)).exp();

        if self.self_comp {
            // 2A ⇌ AA  (single species; strand_a + strand_b = total of that species)
            //   [A] + 2[AA] = A_tot
            //   K = [AA] / [A]²
            //   4K·x² − (4K·A + 1)·x + K·A² = 0,   x = [AA]
            //   x = ((4KA + 1) − √(8KA + 1)) / (8K)
            let a_tot = conditions.strand_a_conc_mol + conditions.strand_b_conc_mol;
            if k <= 0.0 || a_tot <= 0.0 {
                return Equilibrium {
                    a_free_mol: a_tot,
                    b_free_mol: a_tot,
                    ab_bound_mol: 0.0,
                    frac_a_bound: 0.0,
                    frac_b_bound: 0.0,
                };
            }
            let b = 4.0 * k * a_tot + 1.0;
            let disc = (8.0 * k * a_tot + 1.0).max(0.0);
            let x = ((b - disc.sqrt()) / (8.0 * k))
                .max(0.0)
                .min(a_tot / 2.0);
            let a_free = a_tot - 2.0 * x;
            let frac_bound = if a_tot > 0.0 { (2.0 * x) / a_tot } else { 0.0 };
            Equilibrium {
                a_free_mol: a_free,
                b_free_mol: a_free,
                ab_bound_mol: x,
                frac_a_bound: frac_bound,
                frac_b_bound: frac_bound,
            }
        } else {
            // A + B ⇌ AB
            //   [A] + [AB] = A_tot,   [B] + [AB] = B_tot
            //   K = [AB] / ([A][B])
            //   K·x² − (1 + K(A+B))·x + K·A·B = 0
            //   x = ((1 + K(A+B)) − √((1+K(A+B))² − 4K²·A·B)) / (2K)
            let a_tot = conditions.strand_a_conc_mol;
            let b_tot = conditions.strand_b_conc_mol;
            if k <= 0.0 || a_tot <= 0.0 || b_tot <= 0.0 {
                return Equilibrium {
                    a_free_mol: a_tot,
                    b_free_mol: b_tot,
                    ab_bound_mol: 0.0,
                    frac_a_bound: 0.0,
                    frac_b_bound: 0.0,
                };
            }
            let s = a_tot + b_tot;
            let p = a_tot * b_tot;
            let b = 1.0 + k * s;
            let disc = (b * b - 4.0 * k * k * p).max(0.0);
            let x = ((b - disc.sqrt()) / (2.0 * k))
                .max(0.0)
                .min(a_tot.min(b_tot));
            Equilibrium {
                a_free_mol: a_tot - x,
                b_free_mol: b_tot - x,
                ab_bound_mol: x,
                frac_a_bound: if a_tot > 0.0 { x / a_tot } else { 0.0 },
                frac_b_bound: if b_tot > 0.0 { x / b_tot } else { 0.0 },
            }
        }
    }

    /// Melting temperature in °C with full reaction-conditions handling.
    ///
    /// Tm is defined as the temperature where 50% of the *limiting* strand is
    /// in the AB duplex form. Found by binary search over the equilibrium
    /// state (which respects unequal strand concentrations, salt-corrected
    /// `K(T)`, and self-complementarity).
    ///
    /// Then if the user chose Schildkraut or Owczarzy 2004 salt correction
    /// (and no Mg²⁺), the corresponding Tm offset is applied post-search.
    /// If Mg²⁺ is present, the Owczarzy 2008 ratio-regime divalent correction
    /// is applied.
    pub fn tm_at_conditions(self, conditions: &ReactionConditions) -> f64 {
        // 1. Find Tm via equilibrium binary search (with SantaLucia ΔS).
        let limiting_tot = conditions
            .strand_a_conc_mol
            .min(conditions.strand_b_conc_mol);
        let target_ab = if self.self_comp {
            // 2A ⇌ AA. A_tot = strand_a + strand_b (same species).
            // Half-of-strands-bound corresponds to [AA] = A_tot / 4.
            (conditions.strand_a_conc_mol + conditions.strand_b_conc_mol) / 4.0
        } else {
            limiting_tot / 2.0
        };

        let tm_k = binary_search_tm(self, conditions, target_ab);
        let tm_c = tm_k - 273.15;

        // 2. Apply Tm-offset salt corrections.
        if conditions.has_divalent() {
            apply_owczarzy_2008(self, conditions, tm_c)
        } else {
            apply_monovalent_offset(self, conditions, tm_c)
        }
    }
}

fn binary_search_tm(thermo: Thermo, conditions: &ReactionConditions, target_ab: f64) -> f64 {
    // [AB] is monotonically decreasing in T. Search T in (200 K, 400 K).
    let mut lo = 200.0;
    let mut hi = 400.0;
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        let eq = thermo.equilibrium(mid, conditions);
        if eq.ab_bound_mol > target_ab {
            lo = mid;
        } else {
            hi = mid;
        }
        if hi - lo < 1e-6 {
            break;
        }
    }
    0.5 * (lo + hi)
}

/// Effective [Na⁺] used by the Tm-offset methods (Schildkraut, Owczarzy 2004).
/// Uses the equivalent-monovalent approximation when Mg²⁺ < dNTP (i.e. all Mg
/// is chelated). The actual divalent regime is handled by Owczarzy 2008
/// elsewhere.
fn effective_mv_conc(conditions: &ReactionConditions) -> f64 {
    let mg_free = conditions.free_mg_conc();
    conditions.mv_conc_mol + 120.0 * mg_free.sqrt()
}

fn apply_monovalent_offset(thermo: Thermo, conditions: &ReactionConditions, tm_c: f64) -> f64 {
    match conditions.salt_correction {
        // SantaLucia: already incorporated in K(T) inside equilibrium().
        SaltCorrection::SantaLucia => tm_c,
        SaltCorrection::Schildkraut => {
            tm_c + 16.6 * conditions.mv_conc_mol.log10()
                 - 16.6 * (1.0_f64).log10()  // 0; included for clarity
        }
        SaltCorrection::Owczarzy => {
            // Owczarzy 2004 1/Tm formula. Subtract the SantaLucia ΔS
            // contribution we baked into K(T) so we don't double-correct.
            //
            // The cleanest fix is to re-derive Tm at 1 M (= no salt correction
            // in K(T)) and then apply the Owczarzy 1/Tm formula. We do that
            // by running a second binary search with mv_conc temporarily 1 M.
            let conds_1m = ReactionConditions {
                mv_conc_mol: 1.0,
                dv_conc_mol: 0.0,
                dntp_conc_mol: 0.0,
                salt_correction: SaltCorrection::SantaLucia,
                ..*conditions
            };
            let limiting_tot = conds_1m
                .strand_a_conc_mol
                .min(conds_1m.strand_b_conc_mol);
            let target_ab = if thermo.self_comp {
                (conds_1m.strand_a_conc_mol + conds_1m.strand_b_conc_mol) / 4.0
            } else {
                limiting_tot / 2.0
            };
            let tm_1m_k = binary_search_tm(thermo, &conds_1m, target_ab);
            let ln_na = conditions.mv_conc_mol.ln();
            let inv_tm = 1.0 / tm_1m_k
                + (4.29 * thermo.gc_fraction - 3.95) * 1e-5 * ln_na
                + 9.40e-6 * ln_na * ln_na;
            1.0 / inv_tm - 273.15
        }
    }
}

/// Owczarzy 2008 ratio-regime correction for duplexes in Mg²⁺-containing
/// buffer. Three regimes based on R = √[Mg²⁺_free] / [Na⁺]:
///   - R < 0.22: monovalent-only (uses Owczarzy 2004 formula).
///   - R > 6.0:  divalent-only formula.
///   - 0.22 ≤ R ≤ 6.0: divalent formula with [Na⁺]-modified coefficients.
fn apply_owczarzy_2008(thermo: Thermo, conditions: &ReactionConditions, tm_c: f64) -> f64 {
    // Compute Tm at 1 M Na⁺ via the equilibrium framework (no salt correction
    // inside K(T)). This is the reference Tm₁ₘ that the Owczarzy 2008 formulas
    // adjust from.
    let conds_1m = ReactionConditions {
        mv_conc_mol: 1.0,
        dv_conc_mol: 0.0,
        dntp_conc_mol: 0.0,
        salt_correction: SaltCorrection::SantaLucia,
        ..*conditions
    };
    let limiting_tot = conds_1m
        .strand_a_conc_mol
        .min(conds_1m.strand_b_conc_mol);
    let target_ab = if thermo.self_comp {
        (conds_1m.strand_a_conc_mol + conds_1m.strand_b_conc_mol) / 4.0
    } else {
        limiting_tot / 2.0
    };
    let tm_1m_k = binary_search_tm(thermo, &conds_1m, target_ab);

    let mg = conditions.free_mg_conc();
    let mv = conditions.mv_conc_mol;
    let ratio = mg.sqrt() / mv;

    let ln_mg = mg.ln();
    let ln_mv = mv.ln();
    let fgc = thermo.gc_fraction;
    let n = (thermo.n_bp as f64 - 1.0).max(1.0);

    let inv_tm = if ratio < 0.22 {
        // Monovalent-only regime (Owczarzy 2004 formula).
        1.0 / tm_1m_k
            + (4.29 * fgc - 3.95) * 1e-5 * ln_mv
            + 9.40e-6 * ln_mv * ln_mv
    } else {
        // Divalent regime. Coefficients (Owczarzy 2008):
        let mut a = 3.92e-5;
        let b = -9.11e-6;
        let c = 6.26e-5;
        let mut d = 1.42e-5;
        let e = -4.82e-4;
        let f = 5.25e-4;
        let g = 8.31e-5;

        if ratio <= 6.0 {
            // Mixed regime: modify a and d as functions of ln[Na⁺].
            a *= 0.843 - 0.352 * mv.sqrt() * ln_mv;
            d *= 1.279 - 4.03e-3 * ln_mv - 8.03e-3 * ln_mv * ln_mv;
        }

        1.0 / tm_1m_k
            + a
            + b * ln_mg
            + fgc * (c + d * ln_mg)
            + (1.0 / (2.0 * n)) * (e + f * ln_mg + g * ln_mg * ln_mg)
    };

    let _ = tm_c; // tm_c was the SantaLucia equilibrium answer; not used here.
    1.0 / inv_tm - 273.15
}

#[derive(Debug, Clone, PartialEq)]
pub enum ThermoError {
    /// A nearest-neighbor stack contains two non-Watson-Crick pairs (tandem mismatch).
    TandemMismatch {
        position: usize,
        top: (Base, Base),
        bottom: (Base, Base),
    },
    /// A nearest-neighbor stack has a single mismatch of a type that is not
    /// yet implemented (currently anything other than G·T).
    UnsupportedMismatch {
        position: usize,
        top: (Base, Base),
        bottom: (Base, Base),
    },
    /// A bulge or internal loop has a length for which no parameter is defined
    /// (currently only `n == 0`; all positive lengths resolve via lookup or
    /// Jacobson-Stockmayer extrapolation).
    UnknownLoopLength(usize),
    /// The alignment structure was rejected during feature extraction. In
    /// practice this is reached only when an `AlignedDuplex` was constructed
    /// outside the standard validation path.
    Align(AlignError),
    /// A feature that requires flanking Watson-Crick context (Mismatch,
    /// length-1 Bulge) was found at a duplex terminus, where there is no
    /// flanking stem on one side.
    UnflankedFeature,
}

impl From<AlignError> for ThermoError {
    fn from(e: AlignError) -> Self {
        ThermoError::Align(e)
    }
}

impl fmt::Display for ThermoError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ThermoError::TandemMismatch { position, top, bottom } => write!(
                f,
                "tandem mismatch at NN stack position {position}: 5'-{}{}-3' / 3'-{}{}-5' (not supported)",
                top.0, top.1, bottom.0, bottom.1
            ),
            ThermoError::UnsupportedMismatch { position, top, bottom } => write!(
                f,
                "unsupported mismatch at NN stack position {position}: 5'-{}{}-3' / 3'-{}{}-5' (only G·T mismatches implemented)",
                top.0, top.1, bottom.0, bottom.1
            ),
            ThermoError::UnknownLoopLength(n) => write!(
                f,
                "no loop parameter available for length {n}"
            ),
            ThermoError::Align(e) => write!(f, "alignment error: {e}"),
            ThermoError::UnflankedFeature => f.write_str(
                "feature requires flanking WC context but is adjacent to a terminus",
            ),
        }
    }
}

impl std::error::Error for ThermoError {}

/// Compute nearest-neighbor thermodynamics for a duplex.
///
/// Iterates the NN stacks left-to-right. For each stack:
/// - both pairs Watson-Crick → look up in the 1998 unified table,
/// - one pair WC + one G·T mismatch → look up in the Allawi 1997 table,
/// - any other configuration → `ThermoError`.
///
/// Terminal pairs are guaranteed Watson-Crick by `Duplex::new`, so initiation
/// is applied unconditionally per end.
pub fn duplex_thermo(dx: &Duplex) -> Result<Thermo, ThermoError> {
    let mut dh = 0.0;
    let mut ds = 0.0;

    for i in 0..(dx.len() - 1) {
        let (t1, b1) = dx.pair_at(i);
        let (t2, b2) = dx.pair_at(i + 1);

        let p1_wc = b1 == t1.complement();
        let p2_wc = b2 == t2.complement();

        let (h, s) = if p1_wc && p2_wc {
            nn_params(t1, t2)
        } else if !p1_wc && !p2_wc {
            return Err(ThermoError::TandemMismatch {
                position: i,
                top: (t1, t2),
                bottom: (b1, b2),
            });
        } else {
            mismatch_params(t1, t2, b1, b2).ok_or(ThermoError::UnsupportedMismatch {
                position: i,
                top: (t1, t2),
                bottom: (b1, b2),
            })?
        };
        dh += h;
        ds += s;
    }

    for i in [0, dx.len() - 1] {
        let (t, _) = dx.pair_at(i);
        let (h, s) = if t.is_at() { INIT_AT } else { INIT_GC };
        dh += h;
        ds += s;
    }

    let self_comp = dx.is_self_complementary();
    if self_comp {
        dh += SYMMETRY_CORRECTION.0;
        ds += SYMMETRY_CORRECTION.1;
    }

    // Salt-correction metadata: every column is a paired column in this branch
    // (Duplex requires equal-length strands and WC termini).
    let n_bp = dx.len();
    let gc_count = (0..dx.len())
        .filter(|i| {
            let (t, _) = dx.pair_at(*i);
            matches!(t, Base::G | Base::C)
        })
        .count();
    let gc_fraction = gc_count as f64 / n_bp as f64;

    Ok(Thermo {
        dh,
        ds,
        self_comp,
        n_bp,
        gc_fraction,
    })
}

const T_REF: f64 = 310.15;

/// Compute nearest-neighbor thermodynamics for an aligned duplex that may
/// contain bulges, internal loops, and/or mismatches.
///
/// Walks the structural features emitted by
/// [`crate::feature::extract_features`]:
///
/// - **WcStem**: sum NN stacks over consecutive pairs *inside* the stem only.
///   Stack boundaries between adjacent `WcStem`s are not crossed — by
///   construction something (mismatch/bulge/loop) sits between them.
/// - **Mismatch**: two mismatch NN stacks, flanking the mismatched column.
/// - **Bulge**: length penalty from [`bulge_length`]. For length 1 only,
///   also add the "intervening NN" stack (paper's flipped-out approximation)
///   and the +0.5 kcal/mol AT flanking bonus per AT-pair flank.
/// - **InternalLoop**: length penalty only (no intervening NN, no AT bonus).
///
/// Initiation is applied to the first and last WcStem's terminal pair.
/// Self-complementarity is permissive: applied whenever `top == bottom`.
pub fn aligned_duplex_thermo(aligned: &AlignedDuplex) -> Result<Thermo, ThermoError> {
    let features = extract_features(aligned)?;
    if features.is_empty() {
        return Err(ThermoError::Align(AlignError::EmptyStrand));
    }

    let mut dh = 0.0;
    let mut ds = 0.0;

    // Initiation
    let first_top = first_stem_first_top(&features)?;
    let last_top = last_stem_last_top(&features)?;
    for t in [first_top, last_top] {
        let (h, s) = if t.is_at() { INIT_AT } else { INIT_GC };
        dh += h;
        ds += s;
    }

    for k in 0..features.len() {
        match &features[k] {
            Feature::WcStem { tops, .. } => {
                for w in tops.windows(2) {
                    let (h, s) = nn_params(w[0], w[1]);
                    dh += h;
                    ds += s;
                }
            }
            Feature::Mismatch {
                top: mm_top,
                bottom: mm_bot,
                ..
            } => {
                let (prev_top, prev_bot) = prev_stem_last(&features, k)?;
                let (next_top, next_bot) = next_stem_first(&features, k)?;
                let (h1, s1) = mismatch_params(prev_top, *mm_top, prev_bot, *mm_bot)
                    .ok_or(ThermoError::UnsupportedMismatch {
                        position: 0,
                        top: (prev_top, *mm_top),
                        bottom: (prev_bot, *mm_bot),
                    })?;
                let (h2, s2) = mismatch_params(*mm_top, next_top, *mm_bot, next_bot)
                    .ok_or(ThermoError::UnsupportedMismatch {
                        position: 0,
                        top: (*mm_top, next_top),
                        bottom: (*mm_bot, next_bot),
                    })?;
                dh += h1 + h2;
                ds += s1 + s2;
            }
            Feature::Bulge { bases, .. } => {
                let l = bases.len();
                let (h_loop, s_loop) = bulge_length(l)?;
                dh += h_loop;
                ds += s_loop;
                if l == 1 {
                    let (prev_top, _) = prev_stem_last(&features, k)?;
                    let (next_top, _) = next_stem_first(&features, k)?;
                    // Eq. 13 intervening NN: as if the bulged base weren't there.
                    let (h_nn, s_nn) = nn_params(prev_top, next_top);
                    dh += h_nn;
                    ds += s_nn;
                    // +0.5 kcal/mol ΔG°37 per AT-pair flanking the bulge.
                    let at_count = (prev_top.is_at() as u32) + (next_top.is_at() as u32);
                    if at_count > 0 {
                        ds += -0.5 * at_count as f64 * 1000.0 / T_REF;
                    }
                }
            }
            Feature::InternalLoop {
                top_bases,
                bottom_bases,
                ..
            } => {
                let l = top_bases.len() + bottom_bases.len();
                let (h_loop, s_loop) = internal_loop_length(l)?;
                dh += h_loop;
                ds += s_loop;
            }
        }
    }

    let self_comp = aligned.is_self_complementary();
    if self_comp {
        dh += SYMMETRY_CORRECTION.0;
        ds += SYMMETRY_CORRECTION.1;
    }

    // Salt-correction metadata: count paired columns (WC + mismatch). Bulges
    // and internal-loop columns are excluded — they're not phosphate-paired
    // and shouldn't contribute to the SantaLucia electrostatic correction.
    let (n_bp, gc_count) = features.iter().fold((0usize, 0usize), |(n, g), f| match f {
        Feature::WcStem { tops, .. } => {
            let gc = tops.iter().filter(|b| matches!(b, Base::G | Base::C)).count();
            (n + tops.len(), g + gc)
        }
        Feature::Mismatch { .. } => (n + 1, g),
        Feature::Bulge { .. } | Feature::InternalLoop { .. } => (n, g),
    });
    let gc_fraction = if n_bp == 0 { 0.0 } else { gc_count as f64 / n_bp as f64 };

    Ok(Thermo {
        dh,
        ds,
        self_comp,
        n_bp,
        gc_fraction,
    })
}

fn first_stem_first_top(features: &[Feature]) -> Result<Base, ThermoError> {
    match features.first() {
        Some(Feature::WcStem { tops, .. }) => Ok(*tops.first().unwrap()),
        _ => Err(ThermoError::Align(AlignError::TerminalBulge)),
    }
}

fn last_stem_last_top(features: &[Feature]) -> Result<Base, ThermoError> {
    match features.last() {
        Some(Feature::WcStem { tops, .. }) => Ok(*tops.last().unwrap()),
        _ => Err(ThermoError::Align(AlignError::TerminalBulge)),
    }
}

fn prev_stem_last(features: &[Feature], k: usize) -> Result<(Base, Base), ThermoError> {
    if k == 0 {
        return Err(ThermoError::UnflankedFeature);
    }
    match &features[k - 1] {
        Feature::WcStem { tops, bottoms, .. } => {
            Ok((*tops.last().unwrap(), *bottoms.last().unwrap()))
        }
        _ => Err(ThermoError::UnflankedFeature),
    }
}

fn next_stem_first(features: &[Feature], k: usize) -> Result<(Base, Base), ThermoError> {
    if k + 1 >= features.len() {
        return Err(ThermoError::UnflankedFeature);
    }
    match &features[k + 1] {
        Feature::WcStem { tops, bottoms, .. } => {
            Ok((*tops.first().unwrap(), *bottoms.first().unwrap()))
        }
        _ => Err(ThermoError::UnflankedFeature),
    }
}
