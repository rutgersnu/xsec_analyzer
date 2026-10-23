/**
 * ncqe_xsec
 *
 * Single-bin nuNCQE cross-section extraction and uncertainty summary from a univmake output file.
 *
 * Usage:
 *   ncqe_xsec UNIVMAKE.root SYSTCALC.conf NCQE_INPUTS.conf FILE_PROPERTIES.txt [OUTPUT.json]
 *
 * Prints the CV yields per category, the extraction (f, <sigma_GENIE>, <sigma>), every covariance
 * defined in SYSTCALC.conf as a fractional uncertainty on <sigma>, the per-universe distribution of
 * each reweighting set (RMS about the CV, mean, standard deviation, 16th/84th percentiles; individual
 * shifts for two-universe unisims), the MC-statistics breakdown, and the CIT and skyshine terms.
 * No unfolding is used: with one reco bin the extraction is the scale factor of paper Eq. 26.
 */

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"

#include "NCQESystematicsCalculator.hh"

namespace {

bool ends_with( const std::string& s, const std::string& e ) {
  return s.size() >= e.size() && s.compare( s.size() - e.size(), e.size(), e ) == 0;
}

// Same names that UniverseMaker multiplies by the tuned CV weight
bool is_flux_weight( const std::string& name ) {
  return name == "weight_flux_all" || ends_with( name, "_FluxUnisim" )
    || name.find( "_PrimaryHadron" ) != std::string::npos;
}

// numpy.percentile with the default linear interpolation
double percentile( std::vector<double> v, double p ) {
  std::sort( v.begin(), v.end() );
  double idx = p / 100. * ( v.size() - 1 );
  size_t lo = static_cast<size_t>( std::floor( idx ) );
  size_t hi = std::min( lo + 1, v.size() - 1 );
  return v[lo] + ( v[hi] - v[lo] ) * ( idx - lo );
}

std::string json_str( const std::string& s ) {
  std::string out = "\"";
  for ( char c : s ) {
    if ( c == '"' || c == '\\' ) out += '\\';
    out += c;
  }
  return out + "\"";
}

std::string num( double x ) {
  if ( !std::isfinite( x ) ) return "null";
  std::ostringstream o;
  o << std::setprecision( 10 ) << x;
  return o.str();
}

} // namespace

int main( int argc, char* argv[] ) {
  if ( argc != 5 && argc != 6 ) {
    std::cout << "Usage: ncqe_xsec UNIVMAKE.root SYSTCALC.conf NCQE_INPUTS.conf"
      << " FILE_PROPERTIES.txt [OUTPUT.json]\n";
    return 1;
  }
  const std::string univ_file( argv[1] ), syst_cfg( argv[2] ), inputs_cfg( argv[3] ),
    fp_cfg( argv[4] );
  const std::string json_out = ( argc == 6 ) ? argv[5] : "";

  auto& fpm = FilePropertiesManager::Instance();
  fpm.load_file_properties( fp_cfg );

  NCQESystematicsCalculator calc( univ_file, syst_cfg, inputs_cfg );
  const auto& eci = EventCategoryInterpreter::Instance();
  const auto& in = calc.inputs();
  const Universe& cv = calc.cv_universe();
  NCQEExtraction e = calc.extract( cv );

  std::ostringstream js;
  js << "{\n";
  js << "  \"inputs\": {\"univmake\": " << json_str( univ_file ) << ", \"systcalc\": "
     << json_str( syst_cfg ) << ", \"ncqe_inputs\": " << json_str( inputs_cfg )
     << ", \"file_properties\": " << json_str( fp_cfg ) << ", \"n_obs\": " << num( in.n_obs )
     << ", \"d_cit\": " << num( in.d_cit ) << ", \"d_sky\": " << num( in.d_sky )
     << ", \"pot_data\": " << num( calc.total_bnb_data_pot_ ) << ", \"scale_oofv\": "
     << ( in.scale_oofv ? "true" : "false" ) << "},\n";

  // ------------------------------------------------------------------ CV extraction
  std::cout << std::fixed << std::setprecision( 4 );
  std::cout << "\n==================== nuNCQE extraction (CV universe) ====================\n";
  std::cout << "Data POT " << std::scientific << calc.total_bnb_data_pot_ << std::fixed
    << "   N_obs " << in.n_obs << "   D_CIT " << in.d_cit << "   D_sky " << in.d_sky
    << "   OOFV nuNCQE " << ( in.scale_oofv ? "scaled with f" : "subtracted as background" ) << '\n';
  js << "  \"cv\": {\"yields\": {";
  bool first = true;
  for ( EventCategory c : NCQESystematicsCalculator::categories() ) {
    std::cout << "  " << std::setw( 20 ) << std::left << eci.label( c ) << std::right
      << std::setw( 12 ) << e.yields.at( c ) << '\n';
    js << ( first ? "" : ", " ) << json_str( eci.label( c ) ) << ": " << num( e.yields.at( c ) );
    first = false;
  }
  std::cout << "  M_bkg " << e.background << "   signal-like " << e.signal_like
    << "   N_obs - D - M_bkg " << e.numerator << '\n';
  std::cout << "  T (true FV nuNCQE) " << e.T << "   <sigma_GENIE> " << e.sigma_genie << '\n';
  std::cout << "  f = " << e.f << "   <sigma> = " << e.sigma << " x 1e-38 cm^2 / oxygen\n";
  js << "}, \"M_bkg\": " << num( e.background ) << ", \"signal_like\": " << num( e.signal_like )
     << ", \"numerator\": " << num( e.numerator ) << ", \"T\": " << num( e.T )
     << ", \"sigma_genie\": " << num( e.sigma_genie ) << ", \"f\": " << num( e.f )
     << ", \"sigma\": " << num( e.sigma ) << "},\n";

  // Closure: every selected event should fall into one of the category true bins
  double sum_2d = 0.;
  for ( int tb = 0; tb < cv.hist_2d_->GetNbinsX(); ++tb ) sum_2d += cv.hist_2d_->GetBinContent( tb + 1, 1 );
  double reco = cv.hist_reco_->GetBinContent( 1 );
  std::cout << "  closure: selected (reco bin) " << reco << ", sum over true bins " << sum_2d
    << ", difference " << reco - sum_2d << '\n';
  js << "  \"closure\": {\"reco\": " << num( reco ) << ", \"sum_true_bins\": " << num( sum_2d ) << "},\n";

  // ------------------------------------------------------------------ covariances from the config
  auto cov_map = calc.get_covariances();
  // get_covariances() ends its progress output with a carriage return, not a newline
  std::cout << std::endl;
  std::cout << "\n==================== Fractional uncertainties on <sigma> (" << syst_cfg
    << ") ====================\n";
  js << "  \"covariances_percent\": {";
  first = true;
  for ( const auto& pair : *cov_map ) {
    double var = pair.second.cov_matrix_->GetBinContent( 1, 1 );
    double frac = 100. * std::sqrt( var ) / e.sigma;
    std::cout << "  " << std::setw( 24 ) << std::left << pair.first << std::right
      << std::setw( 10 ) << frac << " %\n";
    js << ( first ? "" : ", " ) << json_str( pair.first ) << ": " << num( frac );
    first = false;
  }
  js << "},\n";

  // ------------------------------------------------------------------ per-universe distributions
  std::cout << "\n==================== Universe distributions (sigma_u / sigma_CV - 1, %) ====================\n";
  std::cout << "  " << std::setw( 44 ) << std::left << "weight set" << std::right << std::setw( 6 )
    << "N" << std::setw( 10 ) << "RMS" << std::setw( 10 ) << "mean" << std::setw( 10 ) << "std"
    << std::setw( 10 ) << "p16" << std::setw( 10 ) << "p84" << '\n';
  js << "  \"universes\": {";
  first = true;
  for ( const auto& pair : calc.rw_universes_ ) {
    const std::string& name = pair.first;
    if ( name == calc.CV_UNIV_NAME || name == UNWEIGHTED_NAME ) continue;
    bool flux = is_flux_weight( name );
    std::vector<double> shifts;
    for ( const auto& u : pair.second ) {
      shifts.push_back( 100. * ( calc.extract( *u, flux ).sigma / e.sigma - 1. ) );
    }
    size_t n = shifts.size();
    double mean = 0., ms = 0.;
    for ( double s : shifts ) { mean += s; ms += s * s; }
    mean /= n;
    double rms = std::sqrt( ms / n );
    double sd = std::sqrt( std::max( ms / n - mean * mean, 0. ) );
    std::cout << "  " << std::setw( 44 ) << std::left << name << std::right << std::setw( 6 ) << n
      << std::setw( 10 ) << rms << std::setw( 10 ) << mean << std::setw( 10 ) << sd;
    if ( n > 2 ) {
      std::cout << std::setw( 10 ) << percentile( shifts, 16. ) << std::setw( 10 ) << percentile( shifts, 84. );
    }
    else {
      std::cout << "   shifts:";
      for ( double s : shifts ) std::cout << ' ' << s;
    }
    std::cout << ( flux ? "   [flux: T held at CV]" : "" ) << '\n';
    js << ( first ? "" : ", " ) << json_str( name ) << ": {\"n\": " << n << ", \"rms\": " << num( rms )
       << ", \"mean\": " << num( mean ) << ", \"std\": " << num( sd ) << ", \"flux\": "
       << ( flux ? "true" : "false" );
    if ( n > 2 ) js << ", \"p16\": " << num( percentile( shifts, 16. ) ) << ", \"p84\": " << num( percentile( shifts, 84. ) );
    js << ", \"shifts\": [";
    for ( size_t i = 0; i < n; ++i ) js << ( i ? ", " : "" ) << num( shifts[i] );
    js << "]}";
    first = false;
  }
  js << "},\n";

  // ------------------------------------------------------------------ MC statistics and data-driven terms
  std::cout << "\n==================== MC statistics (delta method, %) ====================\n";
  js << "  \"mc_stat_terms_percent\": {";
  double tot2 = 0.;
  first = true;
  for ( const auto& pair : calc.mc_stat_terms( cv ) ) {
    std::cout << "  " << std::setw( 24 ) << std::left << pair.first << std::right << std::setw( 10 )
      << 100. * pair.second << '\n';
    tot2 += pair.second * pair.second;
    js << ( first ? "" : ", " ) << json_str( pair.first ) << ": " << num( 100. * pair.second );
    first = false;
  }
  std::cout << "  " << std::setw( 24 ) << std::left << "total" << std::right << std::setw( 10 )
    << 100. * std::sqrt( tot2 ) << '\n';
  js << ", \"total\": " << num( 100. * std::sqrt( tot2 ) ) << "},\n";

  double stat = 100. * std::sqrt( in.n_obs ) / e.numerator;
  double cit = 100. * in.d_cit_err / e.numerator;
  double sky = 100. * calc.skyshine_rel_unc();
  std::cout << "\n  data statistics " << stat << " %   CIT statistics " << cit << " %   skyshine " << sky << " %\n";
  js << "  \"data_driven_percent\": {\"stat\": " << num( stat ) << ", \"cit_stat\": " << num( cit )
     << ", \"skyshine\": " << num( sky ) << "}\n}\n";

  if ( !json_out.empty() ) {
    std::ofstream out( json_out );
    out << js.str();
    std::cout << "\nwrote " << json_out << '\n';
  }
  return 0;
}
