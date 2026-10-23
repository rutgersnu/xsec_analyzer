# ANNIE νNCQE single-bin cross section in xsec_analyzer

This directory adds the ANNIE neutral-current quasielastic (νNCQE) analysis of the PRD draft (Oct 2026) to the ANNIE port of xsec_analyzer.

Input: ToolAnalysis MC files that carry the reweighting branches. They go through the νNCQE selection, and the framework computes:
- the flux-averaged cross section;
- its interaction-model (GENIE multisim), flux and MC-statistics uncertainties;
- the unisims;
- the data-statistics and data-driven background terms.

The selection and the extraction are line-by-line ports of the Python analysis that produced the paper numbers (`NCQE_XS_analysis_v3`, branch `adopted-extraction`). `validate/reference_python.py` checks the two against each other.

## Quick start

```bash
make ncqe                                # bin/annie_ncqe_prep, bin/univmake, bin/ncqe_xsec (needs ROOT; tested with 6.40)

# frozen CV sample of the paper (3922 files, no unisim branches)
ncqe/run_ncqe.sh cv /path/to/MC_GENIE-O16 4.41733e20 ncqe/systcalc_cv.conf

# 44-knob reweighting production (two of its 3884 files are empty and are skipped)
ncqe/run_ncqe.sh knob44 /path/to/empirical_all_weights_no_downscaling 4.37208e20 ncqe/systcalc_44knob.conf --unisims
```

- **Outputs** go to `$NCQE_OUT/<tag>` (default `~/builds/ncqe_xsa/<tag>`).
- **Run time** is about 3 min for the 3922-file CV sample with 18 parallel prep jobs.
- **Subsets:** `--nfiles N` processes only the first N files (and `--list FILE` takes the files from a list). The POT argument must then be the simulated POT of those files.

## Pipeline

| Step | Program | What it does |
|---|---|---|
| 1 | `bin/annie_ncqe_prep MC_<n>.root out.root [--unisims]` | νNCQE selection per cluster (`src/NCQESelection.hh`), true category per trigger, weight vectors; writes a slim `phaseIITriggerTree` |
| 2 | `hadd` | One merged MC file. The framework scales every MC file to the data POT and sums them, so the files must be merged first |
| 3 | `bin/univmake files.txt ncqe/bin_config_ncqe.txt univ.root file_properties.txt` | Universe histograms (unchanged framework code) |
| 4 | `bin/ncqe_xsec univ.root systcalc.conf ncqe/ncqe_inputs.conf file_properties.txt out.json` | Extraction and uncertainties through `NCQESystematicsCalculator` |

`run_ncqe.sh` does all four steps. It also writes the file list, the file-properties table and the two one-entry placeholder files.

The placeholder beam-on file carries the data POT and trigger count. The framework requires a beam-on and a beam-off entry per run, but the data enter only through the constants in `ncqe_inputs.conf`.

It also deletes the univmake output before each run: the framework reuses a cached `total_*` directory without checking it.

## Selection (`src/NCQESelection.hh`, `src/annie_ncqe_prep.cpp`)

Per tank cluster, in this order:
1. **Prompt cluster:** `clusterNumber == 0`, time before the end of the prompt buffer (1835.35 ns), at least 10 hits.
2. **No veto activity:**
   - `TankMRDCoinc != 1` and `NoVeto != 0`;
   - no MRD or FMV hit in [0, 300] ns.
3. **In the spill:** 26.795–1548.775 ns.
4. **Single prompt cluster:** no other cluster in the prompt window.
5. **Charge balance and energy:** q_b < 0.6; vertex found; reco FV (r < 0.7 m, |y| < 1.0 m); 5 ≤ E_reco ≤ 12 MeV with E_reco = (PE + 2.90)/5.90.
6. **C_V direction cut:** PE-weighted pair directions of the hits within 13 ns of the first hit, accepted inside a rotated parabola.

Stored flags:

| Flag | Content |
|---|---|
| `sel_ncqe_nobunch` | all of the above |
| `sel_ncqe_inwin` | also within ±3.5 ns of a predicted νNCQE bunch centre (CCinc centre − 0.785 ns) |
| `sel_ncqe` | also within the first 26 bunches (t ≤ 507.19 ns); the paper's selection |

A trigger is kept if it has a selected cluster (any bunch) or is a true FV νNCQE interaction (the efficiency denominator). No trigger had more than one selected cluster in either sample.

**Floating point.** The port reproduces NumPy bit for bit:
- `std::pow` where NumPy squares a scalar;
- the BLAS fused-multiply-add order in `linalg.norm`;
- `-ffp-contract=off`, so the compiler fuses nothing else.

### True categories (`EventCategory.hh`, values 12–18)

The true vertex is converted to metres and centred on the tank. External means outside the tank (r > 1.524 m or |y| > 1.98 m); FV means r < 0.7 m and |y| < 1.0 m.

| Value | Category | Definition |
|---|---|---|
| 12 | νNCQE (FV) | signal: NC, QE, ¹⁶O, final-state ν_μ or ν_e, true vertex in the FV |
| 13 | νNCQE (OOFV) | the same with the vertex in the tank but outside the FV; signal-like |
| 14 | NC other | other NC in the FV |
| 15 | CC | CC in the FV |
| 16 | ν̄NCQE | antineutrino NCQE on ¹⁶O in the FV |
| 17 | OOFV other | everything else in the tank outside the FV |
| 18 | External | vertex outside the tank |

### Weights written by the prep

| Branch | Content |
|---|---|
| `weight_TunedCentralValue_UBGenie` | tuned CV (missing → 0, as the Python CV yields drop such events) |
| `weight_All_UBGenie` | `weight_All0_UBGenie`, the 44- (or 42-) knob GENIE multisim, 100 universes, used as stored |
| `weight_flux_all` | universe-by-universe product of the 13 flux knobs; the framework multiplies it by the tuned CV |
| `weight_<slot>_UBGenie` (8, with `--unisims`) | unisim slots, 2 universes each |
| `DirtMu`, `MRDEff` | `DirtMu` = α_dirt = 0.27 for External, 1 otherwise; `MRDEff` = 1. UniverseMaker multiplies every fill by `DirtMu × MRDEff`, so the full dirt sample is used with weight α. No downsampling. |

Missing universe entries are padded with NaN. The framework's `safe_weight()` maps NaN, inf, negative values and values above 30 to 1, the same rule as the Python `safe()`.

The unisim slots reuse PhaseIITreeMaker's fixed branch names (`LoadReweightGenieEventConfig_NCQE`):

| Slot | Dial | Universes 0 / 1 |
|---|---|---|
| `XSecShape_CCMEC` | NC 2p2h normalisation | +100% / −100% |
| `DecayAngMEC` | M_A^NCEL | ±1σ |
| `Theta_Delta2Npi` | η^NCEL | ±1σ |
| `NormCCCOH` | M_A^CCQE | ±1σ |
| `ThetaDelta2NRad` | Coulomb CCQE | ±1σ |
| `RPA_CCQE` | RPA strength | ±1σ |
| `AxFFCCQEshape`, `VecFFCCQEshape` | form-factor shapes | 1 / 0 |

## Extraction and uncertainties (`src/NCQESystematicsCalculator.hh`)

The observable evaluated in every universe is the cross section itself (paper Eqs. 25–28):

```
f       = (N_obs − D_CIT − D_sky − M_bkg) / (M_FV + M_OOFVν)
M_bkg   = NC other + CC + ν̄NCQE + OOFV other + External
⟨σ⟩     = f · T / (POT_data · Φ · N_T)      T = true FV νNCQE (POT-scaled)
```

Conventions, the same as the paper:
- **External** is held at its CV yield in every universe, because it is constrained by data.
- **GENIE universes and unisims** recompute T, so the signal model enters only through the efficiency.
- **Flux universes** hold T, and therefore ⟨σ_GENIE⟩ and the integrated flux, at the CV.
- **Each systematic** is the mean of (σ_u − σ_CV)² over the universes, i.e. the RMS about the CV (systcalc average flag 1).
- **MC statistics** (`MCstat`) use the delta method, with independent selected yields per category plus the true FV signal events that were not selected.
- **Data statistics** (`BNBstat`) are √N_obs propagated through Eq. 26 (T2K form). The CIT statistics are configured as an `EXTstat` term.
- **The skyshine term** comes from `d_sky_err`.

`ncqe_inputs.conf` holds N_obs, D_CIT ± error, D_sky ± error, Φ and N_T. With `scale_oofv 0`, OOFV νNCQE is subtracted as background instead, which is the configuration used before 2026-10-09; use it for cross-checks only.

`ncqe_xsec` prints and writes to JSON:
- the CV yields per category and the extraction;
- every covariance in the systcalc file as a percentage of σ;
- per-universe distributions (RMS, mean, standard deviation, 16th/84th percentiles), or the two shifts of each unisim;
- the MC-statistics breakdown.

## Validation (2026-10-16)

`validate/reference_python.py` runs the v3 Python selection and conventions on the same files with the same POT:
- `select` compares the selected events one by one;
- `reference` computes the numbers;
- `compare` checks them against `ncqe_xsec`'s JSON.

### Frozen CV sample (`MC_GENIE-O16`, 3922 files, 4.41733e20 POT)

Same events as the paper CV; its own multisim weights are the earlier 42-knob set.

| Quantity | xsec_analyzer | Reference |
|---|---|---|
| Selected clusters (any bunch / in window / CV) | 9059 / 4211 / 1234 | v3 selection: identical event by event |
| CV yields (FV, OOFVν, NCother, CC, ν̄, OOFV other, Ext) | 152.608 / 140.376 / 47.763 / 6.995 / 2.912 / 113.007 / 78.231 | paper Table V |
| f, ⟨σ⟩ (10⁻³⁸ cm²) | 1.4603, 2.9955 | paper Eqs. 29–30 |
| Interaction (42-knob weights) / flux / MC stat. | 6.99 / 9.48 / 5.18 % | systematics ladder step L7 |
| Data stat. / CIT stat. / skyshine | 6.49 / 2.32 / 3.01 % | paper Table IV and Sec. VII |
| Earlier configuration (`scale_oofv 0`, off-beam CIT 127.44, skyshine 14.1): f, interaction, flux, MC stat. | 1.5804, 15.25 / 15.45 / 8.58 % | ladder step L3 |

### 44-knob production (3882 readable files, 4.37208e20 POT), self-contained

| Quantity | xsec_analyzer |
|---|---|
| Selected clusters (any bunch / in window / CV) | 8971 / 4104 / 1293, identical event by event to the v3 selection |
| f, ⟨σ⟩ | 1.3034, 2.675 |
| Interaction (44 knobs, MicroBooNE 1σ) | 10.31 % |
| Flux | 10.61 % |
| MC statistics | 5.58 % |
| NC 2p2h unisim | ∓8.57 % |
| Other unisims | ≤ 0.9 % |

The production's yields and f match `syst_reevaluation/final_44knob/final_syst_log.txt`.

For both samples, every number agrees with the Python reference on the same files to ≤ 4×10⁻¹⁰ (relative), the precision of the JSON output. That covers the CV yields, f, σ, each of the 100 GENIE and 100 flux universe shifts, the 16 unisim shifts, MC statistics, and the data-driven terms.

The outputs and comparison logs are in `validate/results_20261016/`.

**Why these are not the paper's 8.6 / 10.0 %.** The production was reprocessed with re-randomised spill timing, so its first-26-bunch subsample is a different third of the events from the frozen CV. Its CV has more background (M_bkg 278.6 vs 248.9) and a lower f (1.30 vs 1.46), which raises every background-driven fractional uncertainty.

The paper's interaction and flux values instead transfer per-category universe/CV ratios measured on the production onto the frozen CV yields. That "ratio transfer" between two samples is not a single-sample framework operation.

When an MC sample is reprocessed with consistent timing and all 44 knobs, this chain gives the corresponding numbers directly.

## Framework changes outside `ncqe/`

- **`src/FilePropertiesManager.hh`:** optional simulated POT as the fourth column of MC lines, read through `mc_pot()`. `SystematicsCalculator.hh` uses it in place of the hardcoded 3.5e19, which remains the fallback, so existing CC tables are unchanged.
- **`src/EventCategory.hh`:** νNCQE categories 12–18 appended, with labels and colours. CC values are unchanged.
- **`src/UniverseMaker.hh`, `src/EventCategory.hh`, `src/FilePropertiesManager.hh`:** missing standard headers added.
- **`Makefile`:**
  - the `ncqe` target is new;
  - the `bin/univmake` rule now puts the ROOT libraries after the source, because linkers using `--as-needed` (the Ubuntu default) drop libraries listed first;
  - `bin/univmake` gains header dependencies.
