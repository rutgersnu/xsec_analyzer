/**
 * Slice plotting for ANNIE.
 */

// univmake "/exp/annie/data/users/mastbaum/cc/ccinc_closure_univmake.root"
// systcalc "/exp/annie/app/users/mastbaum/cc/xsec_analyzer/ccinc/systcalc/systcalc_data.conf"
// slice "/exp/annie/app/users/mastbaum/cc/xsec_analyzer/ccinc/slice_config_ccinc.txt"

#include <algorithm>

#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "THStack.h"
#include "TLegend.h"
#include "TMatrixD.h"
#include "TPad.h"
#include "TStyle.h"

#include "FilePropertiesManager.hh"
#include "MCC9SystematicsCalculator.hh"
#include "SliceBinning.hh"
#include "SliceHistogram.hh"
#include "WienerSVDUnfolder.hh"
#include "AnnieGeometryTools.hh"


using NFT = NtupleFileType;

//#define USE_FAKE_DATA ""

TH1D* MatrixToTH1(const TMatrixD& m, const TString& name) {
  int n = m.GetNrows();
  TH1D* h1 = new TH1D(name, "", n, 0, n);

  for (int i=0, bin=1; i<n; i++, bin++) {
    h1->SetBinContent(bin, m(i, 0));
  }

  return h1;
}

TMatrixD* TH1ToMatrix(const TH1* h) {
  TMatrixD* m = new TMatrixD(h->GetNbinsX()-1, 1);
  for (int i=1; i<=h->GetNbinsX()-1; i++) {
    m->operator()(i-1, 0) = h->GetBinContent(i);
  }
  return m;
}


void annie_slice_plots(const char* univmake_file,
                       const char* systcalc_config,
                       const char* slice_config) {

  #ifdef USE_FAKE_DATA
    // Initialize the FilePropertiesManager and tell it to treat the NuWro
    // MC ntuples as if they were data
    auto& fpm = FilePropertiesManager::Instance();
    fpm.load_file_properties( "nuwro_file_properties.txt" );
  #endif

  auto* syst_ptr = new MCC9SystematicsCalculator(univmake_file, systcalc_config);
  auto& syst = *syst_ptr;

  // FIXME no hardcoded POT
  const double POTMC = 3.5e19; //12.744e20;
  const double POTBNB = 3.5e19; //11.682e20;

  // Get access to the relevant histograms owned by the SystematicsCalculator
  // object. These contain the reco bin counts that we need to populate the
  // slices below.
  TH1D* reco_bnb_hist = syst.data_hists_.at( NFT::kOnBNB ).get();
  TH1D* reco_ext_hist = syst.data_hists_.at( NFT::kExtBNB ).get();

  #ifdef USE_FAKE_DATA
    // Add the EXT to the "data" when working with fake data
    reco_bnb_hist->Add( reco_ext_hist );
  #endif

  TH2D* category_hist = syst.cv_universe().hist_categ_.get();

  // Total MC prediction in reco bin space
  TH1D* reco_mc_plus_ext_hist = dynamic_cast<TH1D*>(
    syst.cv_universe().hist_reco_.get()->Clone("reco_mc_plus_ext_hist"));
  reco_mc_plus_ext_hist->SetDirectory( nullptr );

  // Keys are covariance matrix types, values are CovMatrix objects that
  // represent the corresponding matrices
  auto* matrix_map_ptr = syst.get_covariances().release();
  auto& matrix_map = *matrix_map_ptr;

  auto* sb_ptr = new SliceBinning(slice_config);
  auto& sb = *sb_ptr;

  for ( size_t sl_idx = 0u; sl_idx < sb.slices_.size(); ++sl_idx ) {

    const auto& slice = sb.slices_.at( sl_idx );

    // We now have all of the reco bin space histograms that we need as input.
    // Use them to make new histograms in slice space.
    SliceHistogram* slice_bnb = SliceHistogram::make_slice_histogram(
      *reco_bnb_hist, slice, &matrix_map.at("BNBstats") );

    SliceHistogram* slice_mc_plus_ext = SliceHistogram::make_slice_histogram(
      *reco_mc_plus_ext_hist, slice, &matrix_map.at("total") );

    // Compute chi2 with the covariance matrix
    auto chi2_result = slice_bnb->get_chi2( *slice_mc_plus_ext );
    std::cout << "Slice " << sl_idx << ": \u03C7\u00b2 = "
      << chi2_result.chi2_ << '/' << chi2_result.num_bins_ << " bins,"
      << " p-value = " << chi2_result.p_value_ << '\n';

    // Build a stack of categorized central-value MC predictions plus the
    // extBNB contribution in slice space
    const auto& eci = EventCategoryInterpreter::Instance();

    THStack* slice_pred_stack = new THStack( "mc+ext", "" );

    const auto& cat_map = eci.label_map();

    // Canvas
    TCanvas* c1 = new TCanvas;
    TLegend* lg = new TLegend( 0.75, 0.6, 0.98, 0.98 );
    slice_bnb->hist_->SetLineColor( kBlack );
    slice_bnb->hist_->SetLineWidth( 1 );
    slice_bnb->hist_->SetMarkerStyle( kFullCircle );
    slice_bnb->hist_->SetMarkerSize( 0.8 );
    slice_bnb->hist_->SetStats( false );
    double ymax = std::max( slice_bnb->hist_->GetMaximum(),
      slice_mc_plus_ext->hist_->GetMaximum() ) * 1.2;
    slice_bnb->hist_->GetYaxis()->SetRangeUser( 0., ymax );

    // Go in reverse so that signal ends up on top. Note that this index is
    // one-based to match the ROOT histograms
    int cat_bin_index = cat_map.size();
    for ( auto iter = cat_map.crbegin(); iter != cat_map.crend(); ++iter )
    {
      EventCategory cat = iter->first;
      TH1D* temp_mc_hist = category_hist->ProjectionY( "temp_mc_hist",
        cat_bin_index, cat_bin_index );
      temp_mc_hist->SetDirectory( nullptr );

      SliceHistogram* temp_slice_mc = SliceHistogram::make_slice_histogram(
        *temp_mc_hist, slice  );

      eci.set_mc_histogram_style( cat, temp_slice_mc->hist_.get() );

      TH1D* h = dynamic_cast<TH1D*>(temp_slice_mc->hist_.get());
      slice_pred_stack->Add(h);

      lg->AddEntry(h, eci.label(cat).c_str(), "f");

      std::string cat_col_prefix = "MC" + std::to_string( cat );

      --cat_bin_index;
    }

    const auto &cm = &matrix_map.at("total");
    SliceHistogram* s4s = SliceHistogram::make_slice_histogram(
        *reco_mc_plus_ext_hist, slice, cm );

      // The SliceHistogram object already set the bin errors appropriately
      // based on the slice covariance matrix. Just change the bin contents
      // for the current histogram to be fractional uncertainties. Also set
      // the "uncertainties on the uncertainties" to zero.
      // TODO: revisit this last bit, possibly assign bin errors here
      for ( const auto& bin_pair : slice.bin_map_ ) {
        int global_bin_idx = bin_pair.first;
        double err = s4s->hist_->GetBinError( global_bin_idx );
        slice_bnb->hist_->SetBinError( global_bin_idx, err );
      }

    slice_bnb->hist_->SetTitle("");
    slice_bnb->hist_->SetXTitle("Reconstructed p_{#mu} (GeV/c)");
    slice_bnb->hist_->SetYTitle("Events");
    slice_bnb->hist_->Draw( "e" );
    slice_pred_stack->Draw( "hist same" );
    slice_mc_plus_ext->hist_->SetLineWidth( 1 );
    slice_mc_plus_ext->hist_->SetLineColor( kBlack );
    slice_mc_plus_ext->hist_->Draw( "same hist" );
    slice_bnb->hist_->Draw( "e same" );

    lg->AddEntry(slice_bnb->hist_.get(), "Data","ep");
    lg->AddEntry(slice_mc_plus_ext->hist_.get(), "MC","l");
    lg->Draw("same");

    TString nameo = Form("slice_hist_%02lu.pdf", sl_idx);
    c1->SaveAs(nameo);

    // Get the binning and axis labels for the current slice by cloning the
    // (empty) histogram owned by the Slice object
    TH1* slice_hist = dynamic_cast< TH1* >(
      slice.hist_->Clone("slice_hist") );

    slice_hist->SetDirectory( nullptr );

    // Keys are labels, values are fractional uncertainty histograms
    auto* fr_unc_hists = new std::map< std::string, TH1* >();
    auto& frac_uncertainty_hists = *fr_unc_hists;

    // Show fractional uncertainties computed using these covariance matrices
    // in the ROOT plot. All configured fractional uncertainties will be
    // included in the output pgfplots file regardless of whether they appear
    // in this vector.
    const std::vector< std::string > cov_mat_keys = {
      "total",
      "flux",
      "xsec_total",
      "MCstats",
      "BNBstats",
      "detVar_total",
      "POT",
      "numTargets"
    };
 
    // Loop over the various systematic uncertainties
    int color = 1;
    for ( const auto& pair : matrix_map ) {

      const auto& key = pair.first;
      const auto& cov_matrix = pair.second;

      SliceHistogram* slice_for_syst = SliceHistogram::make_slice_histogram(
        *reco_mc_plus_ext_hist, slice, &cov_matrix );

      // The SliceHistogram object already set the bin errors appropriately
      // based on the slice covariance matrix. Just change the bin contents
      // for the current histogram to be fractional uncertainties. Also set
      // the "uncertainties on the uncertainties" to zero.
      // TODO: revisit this last bit, possibly assign bin errors here
      for ( const auto& bin_pair : slice.bin_map_ ) {
        int global_bin_idx = bin_pair.first;
/*        if(key == "total"){
          slice_for_syst->hist_->SetBinContent( global_bin_idx, 0);//add together later
          slice_for_syst->hist_->SetBinError( global_bin_idx, 0. );
          continue;
        }*/
        double y = 0.;
        double err = 0.;
        if(key == "BNBstats"){
          y = slice_bnb->hist_->GetBinContent( global_bin_idx );//uses data CV
          err = slice_for_syst->hist_->GetBinError( global_bin_idx );//uses data CV
        } else if(key == "MCstats"){
          slice_for_syst->hist_->Scale(POTMC/POTBNB);
          y = slice_for_syst->hist_->GetBinContent( global_bin_idx );
          slice_for_syst->hist_->Scale(POTBNB/POTMC);
          err = slice_for_syst->hist_->GetBinError( global_bin_idx );
        } else {
          y = slice_for_syst->hist_->GetBinContent( global_bin_idx );
          err = slice_for_syst->hist_->GetBinError( global_bin_idx );
        }
        double frac = 0.;
        if ( y > 0. ) frac = err / y;
        slice_for_syst->hist_->SetBinContent( global_bin_idx, frac );
        slice_for_syst->hist_->SetBinError( global_bin_idx, 0. );
      }

      // Check whether the current covariance matrix name is present in
      // the vector defined above this loop. If it isn't, don't bother to
      // plot it, and just move on to the next one.
      auto cbegin = cov_mat_keys.cbegin();
      auto cend = cov_mat_keys.cend();
      auto iter = std::find( cbegin, cend, key );
      if ( iter == cend ) continue;

      frac_uncertainty_hists[ key ] = slice_for_syst->hist_.get();

      if ( color <= 9 ) ++color;
      if ( color == 5 ) ++color;
      if ( color >= 10 ) color += 10;

      slice_for_syst->hist_->SetLineColor( color );
      slice_for_syst->hist_->SetLineWidth( 3 );
    }

    TCanvas* c2 = new TCanvas;
    TLegend* lg2 = new TLegend( 0.25, 0.7, 0.75, 0.85 );

    auto* total_frac_err_hist = frac_uncertainty_hists.at("total");
    total_frac_err_hist->SetStats( false );
    total_frac_err_hist->GetYaxis()->SetRangeUser(
      0., total_frac_err_hist->GetMaximum() * 1.55 );
    total_frac_err_hist->SetLineColor( kBlack );
    total_frac_err_hist->SetLineWidth( 3 );
    total_frac_err_hist->SetTitle("");
    total_frac_err_hist->SetXTitle("Reconstructed p_{#mu} (GeV/c)");
    total_frac_err_hist->SetYTitle("Fractional uncertainty");
    total_frac_err_hist->Draw( "hist" );

    lg2->AddEntry( total_frac_err_hist, "total", "l" );

    for ( auto& pair : frac_uncertainty_hists ) {
      const auto& name = pair.first;
      TH1* hist = pair.second;
      // We already plotted the "total" one above
      if ( name == "total" ) continue;

      lg2->AddEntry( hist, name.c_str(), "l" );
      hist->Draw( "same hist" );
    }

    lg2->SetNColumns(2);
    lg2->Draw( "same" );

    TString name = Form("slice_syst_%02lu.pdf", sl_idx);
    c2->SaveAs(name);
  }

  // Analysis, unfolding, etc.
  // Smearing
  auto smearcept_ptr = syst.get_cv_smearceptance_matrix();

  TCanvas* c3 = new TCanvas();
  TH2D* hSmear = new TH2D(*smearcept_ptr);
  gStyle->SetOptStat(0);
  hSmear->Draw("hist colz text");
  hSmear->GetXaxis()->SetTitle("True p_{#mu} bin");
  hSmear->GetYaxis()->SetTitle("Reco p_{#mu} bin");
  c3->SaveAs("smearcept.pdf");

  // Unfolding
  std::unique_ptr<Unfolder> unfolder(new WienerSVDUnfolder(
    true, WienerSVDUnfolder::RegularizationMatrixType::kFirstDeriv));

  auto result = unfolder->unfold(syst);

  TCanvas* c4 = new TCanvas();
  TMatrixD* AC = result.add_smear_matrix_.get();
  TH2D* hAC = new TH2D(*AC);
  hAC->SetXTitle("True? p_{#mu} bin");
  hAC->SetYTitle("Unfolded p_{#mu} bin");
  hAC->Draw("colz text");
  c4->SaveAs("ac_matrix.pdf");

  TH2D* hUnfCov = new TH2D(*result.cov_matrix_);
  hUnfCov->SetXTitle("Unfolded p_{#mu} bin");
  hUnfCov->SetYTitle("Unfolded p_{#mu} bin");
  hUnfCov->Draw("colz text");
  c4->SaveAs("cov_matrix.pdf");

  TMatrixDDiag* cov_diag = new TMatrixDDiag(*result.cov_matrix_);

  // Scale by flux/N and bin widths
  double total_pot = syst.total_bnb_data_pot_;
  double integ_flux = integrated_numu_flux_in_FV(total_pot);
  double num_tgt = num_O_targets_in_FV();
  double conv_factor = (num_tgt * integ_flux) / 1e38;

  std::cout << "POT: " << total_pot << std::endl;
  std::cout << "integ_flux: " << integ_flux << std::endl;
  std::cout << "num_tgt: " << num_tgt << std::endl;
  std::cout << "conv_factor: " << conv_factor << std::endl;

  const auto& slice = sb.slices_.at(0);  // Slice 0
  int nbins = slice.hist_->GetNbinsX();
  TH1D* h_bin_widths_mev = \
    new TH1D("h_bin_widths_mev", "", nbins, 0, nbins);

  for (int bin=1; bin<=nbins+1; bin++) {
    h_bin_widths_mev->SetBinContent(bin, slice.hist_->GetBinWidth(bin));
  }

  TH1D* hSignal = MatrixToTH1(*result.unfolded_signal_, "asdf2");
  hSignal->Sumw2();
  hSignal->Scale(1.0 / conv_factor);
  hSignal->Divide(h_bin_widths_mev);
  hSignal->SetXTitle("p_{#mu} bin");
  hSignal->SetYTitle("d#sigma/dp_{#mu} (10^{-38} cm^{2}/(MeV/C)/^{16}O)");
  hSignal->Draw("e1");
  hSignal->SetLineColor(kBlack);
  hSignal->SetMarkerSize(0.05);
  hSignal->GetYaxis()->SetRangeUser(0, 0.005);

  TH1D* genie_cv_truth_t = syst.cv_universe().hist_true_.get();
  TMatrixD genie_cv_truth_m(
    *AC,
    TMatrixD::EMatrixCreatorsOp2::kMult,
    *TH1ToMatrix(genie_cv_truth_t)
  );
  TH1D* genie_cv_truth = MatrixToTH1(genie_cv_truth_m, "truth");
  genie_cv_truth->Sumw2();
  genie_cv_truth->Scale(1.0 / conv_factor);
  genie_cv_truth->Divide(h_bin_widths_mev);
  genie_cv_truth->SetLineColor(kRed);
  genie_cv_truth->SetLineWidth(2);
  genie_cv_truth->Draw("hist same");

  //TH1D* hTrueSignal = MatrixToTH1(*syst.get_cv_true_signal(), "asdf");
  //hTrueSignal->Sumw2();
  //hTrueSignal->Scale(1.0 / conv_factor);
  //hTrueSignal->Divide(h_bin_widths_mev);
  //hTrueSignal->SetLineColor(kBlue);
  //hTrueSignal->SetLineWidth(2);
  //hTrueSignal->Draw("hist same");

  hSignal->Draw("e1 same");

  c4->SetLeftMargin(0.14);
  c4->SaveAs("xs.pdf");
}


int main(int argc, char* argv[]) {
  if (argc < 4) {
    std::cout << "Usage: " << argv[0]
              << " UNIVMAKE.root SYSTCALC.conf SLICE_CFG.txt"
              << std::endl;
    return 0;
  }

  const char* univmake_file = argv[1];
  const char* systcalc_config = argv[2];
  const char* slice_config = argv[3];

  annie_slice_plots(univmake_file, systcalc_config, slice_config);

  return 0;
}


