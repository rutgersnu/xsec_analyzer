#pragma once

// Systematics calculator for the ANNIE nuNCQE single-bin cross section.
//
// The observable evaluated in every universe is the extracted flux-averaged cross section itself,
// using the extraction of the PRD draft (T2K PRD 100, 112009, Eq. 5 with one beam mode):
//
//   N_obs = f (M_FV + M_OOFVnu) + M_bkg + D_CIT + D_sky                  (paper Eq. 25)
//   f     = (N_obs - D_CIT - D_sky - M_bkg) / (M_FV + M_OOFVnu)          (paper Eq. 26)
//   sigma = f <sigma_GENIE>,   <sigma_GENIE> = T / (POT_data Phi N_T)     (paper Eqs. 27-28)
//
// M_c are the selected MC yields per true category (true bins of the univmake configuration, matched
// by their "category == N" formulas), T is the number of true FV nuNCQE interactions (true bin of
// kNCQE_FV), and M_bkg = NC other + CC + nubar NCQE + OOFV other + External. N_obs, D_CIT and D_sky
// are data-driven constants read from a small configuration file.
//
// Conventions (NCQE_XS_analysis_v3, 08_systematics.py):
//   * External (dirt) is constrained by data, so its yield is held at the CV in every universe;
//   * GENIE universes recompute T (the signal model enters only through the efficiency);
//   * flux universes hold T, and therefore <sigma_GENIE> and the integrated flux, at the CV;
//   * covariance = mean over universes of (sigma_u - sigma_CV)^2 (RMS about the CV) when the
//     systematics configuration requests averaging.
// The option scale_oofv = 0 subtracts OOFV nuNCQE as background (the configuration used before
// 2026-10-09), for cross-checks only.

#include <cmath>
#include <fstream>
#include <map>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "EventCategory.hh"
#include "SystematicsCalculator.hh"

// Data-side inputs of the extraction (ncqe/ncqe_inputs.conf)
struct NCQEExtractionInputs {
  double n_obs = 0.;          // selected beam-on events
  double d_cit = 0.;          // beam-uncorrelated (CIT) background, scaled pre-beam count
  double d_cit_err = 0.;      // its statistical uncertainty
  double d_sky = 0.;          // skyshine background
  double d_sky_err = 0.;      // its fit uncertainty
  double flux_per_pot = 0.;   // FV-averaged integrated nu_mu + nu_e flux [nu / cm^2 / POT]
  double n_targets = 0.;      // 16O nuclei in the FV
  bool scale_oofv = true;     // scale OOFV nuNCQE with f (paper) or subtract it (earlier config)

  // Reads "key value" lines; '#' starts a comment
  static NCQEExtractionInputs from_file( const std::string& file_name ) {
    std::ifstream in( file_name );
    if ( !in ) throw std::runtime_error( "Cannot open NCQE inputs file " + file_name );
    std::map< std::string, double > kv;
    std::string line;
    while ( std::getline( in, line ) ) {
      auto hash = line.find( '#' );
      if ( hash != std::string::npos ) line.erase( hash );
      std::istringstream ss( line );
      std::string key;
      double value;
      if ( ss >> key >> value ) kv[ key ] = value;
    }
    auto get = [&]( const std::string& k ) {
      if ( !kv.count( k ) ) throw std::runtime_error( "Missing key " + k + " in " + file_name );
      return kv.at( k );
    };
    NCQEExtractionInputs r;
    r.n_obs = get( "n_obs" );
    r.d_cit = get( "d_cit" );
    r.d_cit_err = get( "d_cit_err" );
    r.d_sky = get( "d_sky" );
    r.d_sky_err = get( "d_sky_err" );
    r.flux_per_pot = get( "flux_per_pot" );
    r.n_targets = get( "n_targets" );
    r.scale_oofv = kv.count( "scale_oofv" ) ? ( kv.at( "scale_oofv" ) != 0. ) : true;
    return r;
  }
};

// Result of the extraction in one universe
struct NCQEExtraction {
  std::map< EventCategory, double > yields;   // selected yield per category (POT-scaled)
  double T = 0.;             // true FV nuNCQE (POT-scaled)
  double signal_like = 0.;   // denominator of Eq. 26
  double background = 0.;    // M_bkg
  double numerator = 0.;     // N_obs - D_CIT - D_sky - M_bkg
  double f = 0.;             // measured scale factor
  double sigma_genie = 0.;   // <sigma_GENIE> [1e-38 cm^2 / oxygen]
  double sigma = 0.;         // extracted <sigma> [1e-38 cm^2 / oxygen]
};

class NCQESystematicsCalculator : public SystematicsCalculator {

  public:

    NCQESystematicsCalculator( const std::string& input_respmat_file_name,
      const std::string& syst_cfg_file_name,
      const std::string& ncqe_inputs_file_name,
      const std::string& respmat_tdirectoryfile_name = "" );

    // Categories in the order used for printing
    static const std::vector< EventCategory >& categories() {
      static const std::vector< EventCategory > cats = { kNCQE_FV, kNCQE_OOFVNu,
        kNCQE_NCOther, kNCQE_CC, kNCQE_NuBar, kNCQE_OOFVOther, kNCQE_External };
      return cats;
    }

    // Extraction in a given universe. hold_truth_at_cv = true for flux universes.
    NCQEExtraction extract( const Universe& univ, bool hold_truth_at_cv = false ) const;

    // Observable: extracted cross section [1e-38 cm^2 / oxygen]
    virtual double evaluate_observable( const Universe& univ, int reco_bin,
      int flux_universe_index = -1 ) const override;

    // MC statistics: delta method (independent selected yields per category, plus the true FV
    // signal events that were not selected)
    virtual double evaluate_mc_stat_covariance( const Universe& univ,
      int reco_bin_a, int reco_bin_b ) const override;

    // use_ext = false: Poisson uncertainty on N_obs; use_ext = true: statistical uncertainty of the
    // CIT estimate (configure as "<name> EXTstat")
    virtual double evaluate_data_stat_covariance( int reco_bin_a,
      int reco_bin_b, bool use_ext ) const override;

    // Relative MC-statistics terms (squared terms add up to the total)
    std::map< std::string, double > mc_stat_terms( const Universe& univ ) const;

    // Relative uncertainty from the skyshine yield uncertainty (no universe; for reporting)
    double skyshine_rel_unc() const;

    inline const NCQEExtractionInputs& inputs() const { return inputs_; }

    // Zero-based true bin of a category
    inline int true_bin( EventCategory c ) const { return cat_to_true_bin_.at( c ); }

  protected:

    NCQEExtractionInputs inputs_;
    std::map< EventCategory, int > cat_to_true_bin_;
};

NCQESystematicsCalculator::NCQESystematicsCalculator(
  const std::string& input_respmat_file_name,
  const std::string& syst_cfg_file_name,
  const std::string& ncqe_inputs_file_name,
  const std::string& respmat_tdirectoryfile_name )
  : SystematicsCalculator( input_respmat_file_name, syst_cfg_file_name,
    respmat_tdirectoryfile_name ),
    inputs_( NCQEExtractionInputs::from_file( ncqe_inputs_file_name ) )
{
  // Match true bins to categories by their "category == N" formulas
  static const std::regex re( "^\\s*category\\s*==\\s*([0-9]+)\\s*$" );
  for ( size_t tb = 0u; tb < true_bins_.size(); ++tb ) {
    std::smatch m;
    const std::string& cuts = true_bins_.at( tb ).signal_cuts_;
    if ( std::regex_match( cuts, m, re ) ) {
      auto c = static_cast< EventCategory >( std::stoi( m[1] ) );
      cat_to_true_bin_[ c ] = static_cast< int >( tb );
    }
  }
  for ( EventCategory c : categories() ) {
    if ( !cat_to_true_bin_.count( c ) ) {
      throw std::runtime_error( "No true bin \"category == "
        + std::to_string( static_cast< int >( c ) ) + "\" in the univmake configuration" );
    }
  }
  if ( reco_bins_.size() != 1u ) {
    throw std::runtime_error( "The nuNCQE extraction expects exactly one reco bin" );
  }
}

NCQEExtraction NCQESystematicsCalculator::extract( const Universe& univ,
  bool hold_truth_at_cv ) const
{
  const Universe& cv = this->cv_universe();
  NCQEExtraction e;
  for ( EventCategory c : categories() ) {
    // External is constrained by data: always the CV yield
    const Universe& u = ( c == kNCQE_External ) ? cv : univ;
    e.yields[ c ] = u.hist_2d_->GetBinContent( this->true_bin( c ) + 1, 1 );
  }
  const Universe& u_truth = hold_truth_at_cv ? cv : univ;
  e.T = u_truth.hist_true_->GetBinContent( this->true_bin( kNCQE_FV ) + 1 );

  const auto& Y = e.yields;
  e.background = Y.at( kNCQE_NCOther ) + Y.at( kNCQE_CC ) + Y.at( kNCQE_NuBar )
    + Y.at( kNCQE_OOFVOther ) + Y.at( kNCQE_External );
  e.signal_like = Y.at( kNCQE_FV );
  if ( inputs_.scale_oofv ) e.signal_like += Y.at( kNCQE_OOFVNu );
  else e.background += Y.at( kNCQE_OOFVNu );

  e.numerator = inputs_.n_obs - e.background - inputs_.d_cit - inputs_.d_sky;
  e.f = e.numerator / e.signal_like;
  e.sigma_genie = e.T / ( total_bnb_data_pot_ * inputs_.flux_per_pot * inputs_.n_targets ) / 1e-38;
  e.sigma = e.f * e.sigma_genie;
  return e;
}

double NCQESystematicsCalculator::evaluate_observable( const Universe& univ,
  int /*reco_bin*/, int flux_universe_index ) const
{
  bool is_flux_universe = ( flux_universe_index >= 0 );
  return this->extract( univ, is_flux_universe ).sigma;
}

std::map< std::string, double > NCQESystematicsCalculator::mc_stat_terms(
  const Universe& univ ) const
{
  NCQEExtraction e = this->extract( univ );
  // Variance of a selected yield: summed squared weights (POT-scaled, incl. alpha for External)
  auto var = [&]( EventCategory c ) {
    double err = univ.hist_2d_->GetBinError( this->true_bin( c ) + 1, 1 );
    return err * err;
  };
  double var_bkg = var( kNCQE_NCOther ) + var( kNCQE_CC ) + var( kNCQE_NuBar )
    + var( kNCQE_OOFVOther ) + var( kNCQE_External );
  if ( !inputs_.scale_oofv ) var_bkg += var( kNCQE_OOFVNu );

  // T = S + U, with S the selected FV signal and U the true FV signal events not selected
  double err_T = univ.hist_true_->GetBinError( this->true_bin( kNCQE_FV ) + 1 );
  double var_U = err_T * err_T - var( kNCQE_FV );

  std::map< std::string, double > terms;
  terms[ "backgrounds" ] = std::sqrt( var_bkg ) / e.numerator;
  terms[ "FV_signal" ] = std::abs( 1. / e.T - 1. / e.signal_like ) * std::sqrt( var( kNCQE_FV ) );
  terms[ "OOFV_nuNCQE" ] = inputs_.scale_oofv
    ? std::sqrt( var( kNCQE_OOFVNu ) ) / e.signal_like : 0.;
  terms[ "truth_denominator" ] = std::sqrt( std::max( var_U, 0. ) ) / e.T;
  return terms;
}

double NCQESystematicsCalculator::evaluate_mc_stat_covariance( const Universe& univ,
  int /*reco_bin_a*/, int /*reco_bin_b*/ ) const
{
  double rel2 = 0.;
  for ( const auto& pair : this->mc_stat_terms( univ ) ) rel2 += pair.second * pair.second;
  double sigma = this->extract( univ ).sigma;
  return rel2 * sigma * sigma;
}

double NCQESystematicsCalculator::evaluate_data_stat_covariance( int /*reco_bin_a*/,
  int /*reco_bin_b*/, bool use_ext ) const
{
  // d sigma / d N_obs = - d sigma / d D_CIT = sigma / numerator
  NCQEExtraction e = this->extract( this->cv_universe() );
  double var_counts = use_ext ? inputs_.d_cit_err * inputs_.d_cit_err : inputs_.n_obs;
  double dsig = e.sigma / e.numerator;
  return var_counts * dsig * dsig;
}

double NCQESystematicsCalculator::skyshine_rel_unc() const {
  NCQEExtraction e = this->extract( this->cv_universe() );
  return inputs_.d_sky_err / e.numerator;
}
