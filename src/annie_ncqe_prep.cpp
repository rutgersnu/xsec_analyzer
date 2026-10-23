/**
 * annie_ncqe_prep
 *
 * Prepare ANNIE ToolAnalysis MC files (phaseIITriggerTree + phaseIITankClusterTree) for the
 * nuNCQE single-bin cross-section analysis in xsec_analyzer:
 *
 * * Apply the nuNCQE cluster selection (src/NCQESelection.hh, a port of NCQE_XS_analysis_v3)
 *   and store per-trigger selection flags
 * * Assign the true event category (EventCategory.hh, kNCQE_*)
 * * Write the weight vectors used by UniverseMaker: TunedCV, the GENIE multisim (All0), the
 *   product of the 13 flux knobs, and optionally the 8 unisim slots
 * * Set DirtMu = alpha_dirt for External (dirt) events, so that the full dirt sample is used with
 *   weight alpha instead of a downsampled subset
 * * Keep only triggers that pass the selection (any bunch) or are true FV nuNCQE (the efficiency
 *   denominator); all other triggers contribute nothing to the analysis
 *
 * Usage:
 *   annie_ncqe_prep INPUT.root OUTPUT.root [--unisims] [--keep-all] [--file-number N]
 *   annie_ncqe_prep --placeholder OUTPUT.root
 *
 * The placeholder mode writes the one-entry beam-on / beam-off stand-in files that the framework's
 * file-properties bookkeeping requires (no MC truth, no weights, nothing selected).
 *
 * S. Doran 2026-10, based on A. Mastbaum's annie_stv_prep.
 */

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <regex>
#include <string>
#include <vector>

#include "TBranch.h"
#include "TFile.h"
#include "TTree.h"

#include "EventCategory.hh"
#include "NCQESelection.hh"

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

// Flux knobs in the order of NCQE_XS_analysis_v3/lib/reweight.py (FLUX_KNOBS); the product is
// formed in this order so that it matches the Python value bit for bit
const std::vector<std::string> kFluxKnobs = {
  "horncurrent_FluxUnisim", "expskin_FluxUnisim",
  "piplus_PrimaryHadronSWCentralSplineVariation", "piminus_PrimaryHadronSWCentralSplineVariation",
  "kplus_PrimaryHadronFeynmanScaling", "kzero_PrimaryHadronSanfordWang",
  "kminus_PrimaryHadronNormalization", "pioninexsec_FluxUnisim", "pionqexsec_FluxUnisim",
  "piontotxsec_FluxUnisim", "nucleoninexsec_FluxUnisim", "nucleonqexsec_FluxUnisim",
  "nucleontotxsec_FluxUnisim"
};

// Unisim slots of the 44-knob reprocessing. PhaseIITreeMaker writes fixed branch names, so the
// reweighting config (LoadReweightGenieEventConfig_NCQE) re-used them for these dials:
//   XSecShape_CCMEC = NormNCMEC (+100% / -100%), DecayAngMEC = MaNCEL, Theta_Delta2Npi = EtaNCEL,
//   NormCCCOH = MaCCQE, ThetaDelta2NRad = CoulombCCQE, RPA_CCQE, AxFFCCQEshape, VecFFCCQEshape
const std::vector<std::string> kUnisimSlots = {
  "XSecShape_CCMEC", "DecayAngMEC", "Theta_Delta2Npi", "NormCCCOH", "ThetaDelta2NRad",
  "RPA_CCQE", "AxFFCCQEshape", "VecFFCCQEshape"
};

constexpr size_t kNumGenieUniv = 100;
constexpr size_t kNumFluxUniv = 100;
constexpr size_t kNumUnisimUniv = 2;

// Python category (lib/selection._categories) + OOFV split -> EventCategory value
int framework_category( const ncqe::Truth& t ) {
  switch ( t.cat ) {
    case 0: return kNCQE_FV;
    case 1: return kNCQE_NCOther;
    case 2: return kNCQE_CC;
    case 3: return kNCQE_NuBar;
    case 4: return t.nu ? kNCQE_OOFVNu : kNCQE_OOFVOther;
    case 5: return kNCQE_External;
    default: return kUnknown;   // true FV NC QE on 16O with an unexpected final-state lepton
  }
}

// Copy a weight vector into a fixed-length output vector, padding missing entries with NaN
// (UniverseMaker's safe_weight() maps NaN to 1, as the Python code does)
void fixed_length_copy( const std::vector<double>* in, size_t n, std::vector<double>& out ) {
  out.assign( n, kNaN );
  if ( !in ) return;
  if ( in->size() > n ) {
    throw std::runtime_error( "Weight vector longer than the expected number of universes" );
  }
  for ( size_t u = 0; u < in->size(); ++u ) out[u] = in->at( u );
}

int file_number_from_name( const std::string& path ) {
  std::smatch m;
  static const std::regex re( "MC_([0-9]+)\\.root$" );
  if ( std::regex_search( path, m, re ) ) return std::stoi( m[1] );
  return -1;
}

// Per-trigger result of the cluster selection
struct ClusterSel {
  bool selected = false;    // all cuts except the bunch-timing cuts
  bool in_win = false;      // +-1 sigma around a predicted nuNCQE bunch centre
  bool in26 = false;        // first 26 bunches
  double time = kNaN;       // cluster time (bunchTimes) [ns]
  double reco_e = kNaN;     // reconstructed energy [MeV]
  int n_selected = 0;       // number of selected clusters in the trigger (expected <= 1)
};

int write_placeholder( const std::string& out_name ) {
  TFile fo( out_name.c_str(), "recreate" );
  TTree to( "phaseIITriggerTree", "nuNCQE placeholder (beam-on / beam-off stand-in)" );
  Int_t category = kUnknown;
  Bool_t sel = false, sel_inwin = false;
  to.Branch( "category", &category, "category/I" );
  to.Branch( "sel_ncqe", &sel, "sel_ncqe/O" );
  to.Branch( "sel_ncqe_inwin", &sel_inwin, "sel_ncqe_inwin/O" );
  to.Fill();
  to.Write();
  fo.Close();
  std::cout << "wrote placeholder " << out_name << '\n';
  return 0;
}

} // namespace

int ncqe_prep( const std::string& in_name, const std::string& out_name, bool write_unisims,
  bool keep_all, int file_number )
{
  std::unique_ptr<TFile> f_in( TFile::Open( in_name.c_str(), "read" ) );
  if ( !f_in || f_in->IsZombie() ) {
    std::cerr << "SKIP " << in_name << ": cannot open\n";
    return 2;
  }
  TTree* S = nullptr;
  TTree* T = nullptr;
  f_in->GetObject( "phaseIITriggerTree", S );
  f_in->GetObject( "phaseIITankClusterTree", T );
  if ( !S || !T ) {
    std::cerr << "SKIP " << in_name << ": missing phaseIITriggerTree or phaseIITankClusterTree\n";
    return 2;
  }
  const Long64_t n_trig = S->GetEntries();
  const Long64_t n_clus = T->GetEntries();

  // ------------------------------------------------------------------ pass 1: trigger truth and veto
  S->SetBranchStatus( "*", false );
  for ( const char* b : { "trueNuIntxVtx_X", "trueNuIntxVtx_Y", "trueNuIntxVtx_Z", "trueNC",
    "trueQEL", "trueTargetZ", "trueFSLPdg", "MRDhitT", "FMVhitT" } ) S->SetBranchStatus( b, true );

  double vtx_x, vtx_y, vtx_z;
  int true_nc, true_qel, true_z, true_fsl;
  std::vector<double>* mrd_t = nullptr;
  std::vector<double>* fmv_t = nullptr;
  S->SetBranchAddress( "trueNuIntxVtx_X", &vtx_x );
  S->SetBranchAddress( "trueNuIntxVtx_Y", &vtx_y );
  S->SetBranchAddress( "trueNuIntxVtx_Z", &vtx_z );
  S->SetBranchAddress( "trueNC", &true_nc );
  S->SetBranchAddress( "trueQEL", &true_qel );
  S->SetBranchAddress( "trueTargetZ", &true_z );
  S->SetBranchAddress( "trueFSLPdg", &true_fsl );
  S->SetBranchAddress( "MRDhitT", &mrd_t );
  S->SetBranchAddress( "FMVhitT", &fmv_t );

  std::vector<ncqe::Truth> truth( n_trig );
  std::vector<char> veto_activity( n_trig, 0 );   // MRD or FMV hit in [0, 300] ns
  auto any_in_window = []( const std::vector<double>* v ) {
    if ( !v ) return false;
    for ( double t : *v ) if ( t >= 0. && t <= 300. ) return true;
    return false;
  };
  for ( Long64_t i = 0; i < n_trig; ++i ) {
    S->GetEntry( i );
    truth[i] = ncqe::truth_info( vtx_x, vtx_y, vtx_z, true_nc, true_qel, true_z, true_fsl );
    veto_activity[i] = any_in_window( mrd_t ) || any_in_window( fmv_t );
  }
  S->ResetBranchAddresses();

  // ------------------------------------------------------------------ cluster selection
  T->SetBranchStatus( "*", false );
  for ( const char* b : { "eventNumber", "clusterNumber", "clusterPE", "clusterChargeBalance",
    "clusterHits", "bunchTimes", "TankMRDCoinc", "NoVeto", "recoLeastSqVtxX", "recoLeastSqVtxY",
    "recoLeastSqVtxZ", "hitX", "hitY", "hitZ", "hitPE", "hitT" } ) T->SetBranchStatus( b, true );

  int ev_num, clus_num, tank_mrd, no_veto;
  unsigned int clus_hits;
  double clus_pe, clus_qb, bunch_t, rvx, rvy, rvz;
  std::vector<double>* hx = nullptr;
  std::vector<double>* hy = nullptr;
  std::vector<double>* hz = nullptr;
  std::vector<double>* hpe = nullptr;
  std::vector<double>* ht = nullptr;
  T->SetBranchAddress( "eventNumber", &ev_num );
  T->SetBranchAddress( "clusterNumber", &clus_num );
  T->SetBranchAddress( "clusterPE", &clus_pe );
  T->SetBranchAddress( "clusterChargeBalance", &clus_qb );
  T->SetBranchAddress( "clusterHits", &clus_hits );
  T->SetBranchAddress( "bunchTimes", &bunch_t );
  T->SetBranchAddress( "TankMRDCoinc", &tank_mrd );
  T->SetBranchAddress( "NoVeto", &no_veto );
  T->SetBranchAddress( "recoLeastSqVtxX", &rvx );
  T->SetBranchAddress( "recoLeastSqVtxY", &rvy );
  T->SetBranchAddress( "recoLeastSqVtxZ", &rvz );
  T->SetBranchAddress( "hitX", &hx );
  T->SetBranchAddress( "hitY", &hy );
  T->SetBranchAddress( "hitZ", &hz );
  T->SetBranchAddress( "hitPE", &hpe );
  T->SetBranchAddress( "hitT", &ht );
  TBranch* hit_branches[] = { T->GetBranch( "hitX" ), T->GetBranch( "hitY" ),
    T->GetBranch( "hitZ" ), T->GetBranch( "hitPE" ), T->GetBranch( "hitT" ) };
  std::vector<TBranch*> scalar_branches;
  for ( const char* b : { "eventNumber", "clusterNumber", "clusterPE", "clusterChargeBalance",
    "clusterHits", "bunchTimes", "TankMRDCoinc", "NoVeto", "recoLeastSqVtxX", "recoLeastSqVtxY",
    "recoLeastSqVtxZ" } ) scalar_branches.push_back( T->GetBranch( b ) );

  // Segments of the "single prompt cluster" requirement: a new segment starts at every cluster
  // with clusterNumber == 0, in tree order (np.cumsum(cn == 0) in the Python code)
  std::vector<long> seg( n_clus );
  std::vector<char> seg_has_other;   // indexed by segment id
  seg_has_other.push_back( 0 );
  long cur_seg = 0;
  for ( Long64_t k = 0; k < n_clus; ++k ) {
    for ( TBranch* b : scalar_branches ) b->GetEntry( k );
    if ( clus_num == 0 ) { ++cur_seg; seg_has_other.push_back( 0 ); }
    seg[k] = cur_seg;
    if ( clus_num != 0 && bunch_t < ncqe::kTimeToPromptEnd ) seg_has_other[cur_seg] = 1;
  }

  std::vector<ClusterSel> csel( n_trig );
  Long64_t n_candidates = 0;
  for ( Long64_t k = 0; k < n_clus; ++k ) {
    for ( TBranch* b : scalar_branches ) b->GetEntry( k );
    if ( ev_num < 0 || ev_num >= n_trig ) {
      throw std::runtime_error( "eventNumber outside the trigger tree in " + in_name );
    }
    // prompt cluster
    bool m = ( clus_num == 0 ) && !( bunch_t > ncqe::kTimeToPromptEnd ) && !( clus_hits < 10 );
    // no FMV / MRD activity
    m = m && !( tank_mrd == 1 || no_veto == 0 || veto_activity[ev_num] );
    // in spill
    m = m && !( bunch_t > ncqe::kBunchTimeEnd ) && !( bunch_t < ncqe::kBunchTimeStart );
    // single prompt cluster
    m = m && !seg_has_other[ seg[k] ];
    // q_b, reco FV, E_reco
    m = m && ncqe::core_cuts( clus_qb, rvx, rvy, rvz, clus_pe );
    if ( !m ) continue;
    ++n_candidates;
    // C_V direction cut (hit level)
    for ( TBranch* b : hit_branches ) b->GetEntry( k );
    if ( !ncqe::direction_ok( *hx, *hy, *hz, *hpe, *ht ) ) continue;

    auto& cs = csel[ev_num];
    ++cs.n_selected;
    cs.selected = true;
    cs.time = bunch_t;
    cs.reco_e = ncqe::reco_energy( clus_pe );
    cs.in_win = ncqe::in_bunch_window( bunch_t );
    cs.in26 = !( bunch_t > ncqe::kMCThirdSpillCut );
  }
  T->ResetBranchAddresses();

  // ------------------------------------------------------------------ pass 2: write kept triggers
  S->SetBranchStatus( "*", false );
  std::vector<std::string> weight_names = { "weight_TunedCentralValue_UBGenie", "weight_All0_UBGenie" };
  for ( const auto& k : kFluxKnobs ) weight_names.push_back( "weight_" + k );
  if ( write_unisims ) for ( const auto& s : kUnisimSlots ) weight_names.push_back( "weight_" + s + "_UBGenie" );
  for ( const auto& b : weight_names ) {
    if ( !S->GetBranch( b.c_str() ) ) {
      throw std::runtime_error( "Missing weight branch " + b + " in " + in_name );
    }
    S->SetBranchStatus( b.c_str(), true );
  }
  S->SetBranchStatus( "mcEntryNumber", true );

  std::map<std::string, std::vector<double>*> win;
  for ( const auto& b : weight_names ) {
    win[b] = nullptr;
    S->SetBranchAddress( b.c_str(), &win[b] );
  }
  int mc_entry;
  S->SetBranchAddress( "mcEntryNumber", &mc_entry );

  TFile fo( out_name.c_str(), "recreate" );
  TTree to( "phaseIITriggerTree", "ANNIE nuNCQE selection (annie_ncqe_prep)" );

  Int_t o_mc_entry, o_file_number, o_trigger_entry, o_category, o_py_cat;
  Bool_t o_sel_all, o_sel_inwin, o_sel, o_true_signal, o_true_fv, o_true_ext, o_true_nu;
  Double_t o_time, o_reco_e, o_dirt_mu, o_mrd_eff;
  std::vector<double> o_cv( 1 ), o_genie, o_flux, o_flux_tmp;
  std::vector<std::vector<double>> o_unisim( kUnisimSlots.size() );

  to.Branch( "mcEntryNumber", &o_mc_entry, "mcEntryNumber/I" );
  to.Branch( "mc_file_number", &o_file_number, "mc_file_number/I" );
  to.Branch( "trigger_entry", &o_trigger_entry, "trigger_entry/I" );
  to.Branch( "category", &o_category, "category/I" );
  to.Branch( "ncqe_py_category", &o_py_cat, "ncqe_py_category/I" );
  to.Branch( "sel_ncqe_nobunch", &o_sel_all, "sel_ncqe_nobunch/O" );
  to.Branch( "sel_ncqe_inwin", &o_sel_inwin, "sel_ncqe_inwin/O" );
  to.Branch( "sel_ncqe", &o_sel, "sel_ncqe/O" );
  to.Branch( "ncqe_true_signal", &o_true_signal, "ncqe_true_signal/O" );
  to.Branch( "true_fv", &o_true_fv, "true_fv/O" );
  to.Branch( "true_external", &o_true_ext, "true_external/O" );
  to.Branch( "true_nu_ncqe_o16", &o_true_nu, "true_nu_ncqe_o16/O" );
  to.Branch( "cluster_time", &o_time, "cluster_time/D" );
  to.Branch( "reco_energy", &o_reco_e, "reco_energy/D" );
  to.Branch( "DirtMu", &o_dirt_mu, "DirtMu/D" );
  to.Branch( "MRDEff", &o_mrd_eff, "MRDEff/D" );
  to.Branch( "weight_TunedCentralValue_UBGenie", &o_cv );
  to.Branch( "weight_All_UBGenie", &o_genie );
  to.Branch( "weight_flux_all", &o_flux );
  if ( write_unisims ) {
    for ( size_t s = 0; s < kUnisimSlots.size(); ++s ) {
      to.Branch( ( "weight_" + kUnisimSlots[s] + "_UBGenie" ).c_str(), &o_unisim[s] );
    }
  }

  Long64_t n_kept = 0, n_sel = 0, n_inwin = 0, n_cv = 0, n_signal = 0;
  Long64_t n_missing_cv = 0, n_cv_above_30 = 0, n_multi = 0, n_unassigned_sel = 0;
  for ( Long64_t i = 0; i < n_trig; ++i ) {
    const auto& tr = truth[i];
    const auto& cs = csel[i];
    bool true_signal = tr.fv && tr.nu;
    if ( !( cs.selected || true_signal || keep_all ) ) continue;
    S->GetEntry( i );

    o_mc_entry = mc_entry;
    o_file_number = file_number;
    o_trigger_entry = static_cast<Int_t>( i );
    o_category = framework_category( tr );
    o_py_cat = tr.cat;
    o_sel_all = cs.selected;
    o_sel_inwin = cs.selected && cs.in_win;
    o_sel = cs.selected && cs.in_win && cs.in26;
    o_true_signal = true_signal;
    o_true_fv = tr.fv;
    o_true_ext = tr.ext;
    o_true_nu = tr.nu;
    o_time = cs.time;
    o_reco_e = cs.reco_e;
    o_dirt_mu = tr.ext ? ncqe::kAlphaDirt : 1.;
    o_mrd_eff = 1.;

    // Tuned CV: missing -> 0 (the Python CV yields drop such events); +-inf kept, which
    // safe_weight() maps to 1 (the Python convention for the CV weight)
    const auto* cv_in = win.at( "weight_TunedCentralValue_UBGenie" );
    bool cv_missing = !cv_in || cv_in->empty();
    o_cv[0] = cv_missing ? 0. : cv_in->front();

    fixed_length_copy( win.at( "weight_All0_UBGenie" ), kNumGenieUniv, o_genie );

    // Product of the 13 flux knobs, universe by universe, in the Python order
    o_flux.assign( kNumFluxUniv, 1. );
    for ( size_t k = 0; k < kFluxKnobs.size(); ++k ) {
      fixed_length_copy( win.at( "weight_" + kFluxKnobs[k] ), kNumFluxUniv, o_flux_tmp );
      for ( size_t u = 0; u < kNumFluxUniv; ++u ) {
        o_flux[u] = ( k == 0 ) ? o_flux_tmp[u] : o_flux[u] * o_flux_tmp[u];
      }
    }

    if ( write_unisims ) {
      for ( size_t s = 0; s < kUnisimSlots.size(); ++s ) {
        fixed_length_copy( win.at( "weight_" + kUnisimSlots[s] + "_UBGenie" ), kNumUnisimUniv,
          o_unisim[s] );
      }
    }

    to.Fill();
    ++n_kept;
    if ( true_signal ) ++n_signal;
    if ( cs.selected ) {
      ++n_sel;
      if ( o_sel_inwin ) ++n_inwin;
      if ( o_sel ) ++n_cv;
      if ( cv_missing ) ++n_missing_cv;
      if ( !cv_missing && std::isfinite( o_cv[0] ) && o_cv[0] > 30. ) ++n_cv_above_30;
      if ( cs.n_selected > 1 ) ++n_multi;
      if ( o_category == kUnknown ) ++n_unassigned_sel;
    }
  }
  to.Write();
  fo.Close();

  // One summary line per file (parsed by the run script)
  std::cout << "PREP " << in_name << " triggers " << n_trig << " clusters " << n_clus
    << " candidates " << n_candidates << " selected " << n_sel << " inwin " << n_inwin
    << " cv " << n_cv << " true_signal " << n_signal << " kept " << n_kept
    << " missing_cv " << n_missing_cv << " cv_gt30 " << n_cv_above_30
    << " multi_selected " << n_multi << " unassigned " << n_unassigned_sel << '\n';
  if ( n_multi > 0 ) {
    std::cerr << "WARNING " << in_name << ": " << n_multi
      << " triggers with more than one selected cluster (only the last is stored)\n";
  }
  return 0;
}

int main( int argc, char* argv[] ) {
  std::vector<std::string> args( argv + 1, argv + argc );
  if ( args.size() == 2 && args[0] == "--placeholder" ) return write_placeholder( args[1] );

  if ( args.size() < 2 ) {
    std::cout << "Usage: " << argv[0]
      << " INPUT.root OUTPUT.root [--unisims] [--keep-all] [--file-number N]\n"
      << "       " << argv[0] << " --placeholder OUTPUT.root\n";
    return 1;
  }
  std::string in_name = args[0], out_name = args[1];
  bool unisims = false, keep_all = false;
  int file_number = file_number_from_name( in_name );
  for ( size_t a = 2; a < args.size(); ++a ) {
    if ( args[a] == "--unisims" ) unisims = true;
    else if ( args[a] == "--keep-all" ) keep_all = true;
    else if ( args[a] == "--file-number" && a + 1 < args.size() ) file_number = std::stoi( args[++a] );
    else { std::cerr << "Unknown option " << args[a] << '\n'; return 1; }
  }
  return ncqe_prep( in_name, out_name, unisims, keep_all, file_number );
}
