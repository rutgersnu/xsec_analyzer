// Standard library includes
#include <algorithm>
#include <iomanip>

// ROOT includes
#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "THStack.h"
#include "TLegend.h"

// STV analysis includes
#include "FilePropertiesManager.hh"
#include "MCC9SystematicsCalculator.hh"
#include "includes/PlotUtils.hh"
#include "SliceBinning.hh"
#include "SliceHistogram.hh"

#include "WienerSVDUnfolder.hh"
//#include "TMatrixD.h"
//#include "TLatex.h"
//#include "RooStats/RooStatsUtils.h"
//#include "TDecompChol.h"
#include "TDecompSVD.h"
#include "includes/AnnieGeometryTools.hh"

using NFT = NtupleFileType;

//#define USE_FAKE_DATA ""
void visualizeMatrix(const TMatrixD & matrix, std::string title,  
 char *pdfname, std::string XaxisTitle, std::string YaxisTitle, bool Printnumbers  ) {
    // Get the number of rows and columns in the matrix
    Int_t numRows = matrix.GetNrows();
    Int_t numCols = matrix.GetNcols();
    char name[1024];
    sprintf(name, "%s" ,title.c_str()); 
    std::cout<<"inside::isualizeMatrix"<< "numRows = "<< numRows<< "numCols = "<< numCols<< std::endl;
    // Create a TH2D histogram to visualize the matrix
    TH2D* hist = new TH2D(uniq(), "Matrix Visualization", numCols, 0, numCols, numRows, 0, numRows);

    // Fill the histogram with the matrix elements
    for (Int_t i = 0; i < numRows; ++i) {
        for (Int_t j = 0; j < numCols; ++j) {
            hist->SetBinContent(j + 1, i + 1, matrix[i][j]); // SetBinContent takes bin indices starting from 1
        }
    }

    // Create a canvas and draw the histogram
    TCanvas* canvas = new TCanvas("canvas", "Matrix Visualization Canvas");
    
    gStyle->SetOptStat(0);
    //canvas.cd(); 
    hist->GetXaxis()->SetTitle(XaxisTitle.c_str()); 
    hist->GetYaxis()->SetTitle(YaxisTitle .c_str());
      
    hist->SetTitle(name);
    hist->Draw("COLZ"); // COLZ option for a colored 2D plot
    if( Printnumbers== true)
    {gStyle->SetPaintTextFormat("2.3f");
    hist->Draw("text same");}
    canvas -> Print(pdfname);
    // Run the event loop
    //canvas->Modified();
    //canvas->Update();
    //canvas->Draw(pdfname);

    // Keep the program running to interact with the plot
    canvas->WaitPrimitive();

    // Clean up
    delete hist;
    delete canvas;
}

TMatrixD covarianceToCorrelation(const TMatrixD& covarianceMatrix) {
    // Clone the covariance matrix
    TMatrixD clonedMatrix = TMatrixD(covarianceMatrix); // Using copy constructor to clone
    int n = clonedMatrix.GetNrows();
    TMatrixD correlationMatrix(n, n);

    for (int i = 0; i < n; ++i) {
        double var_i = clonedMatrix(i, i);
        for (int j = 0; j < n; ++j) {
            double cov_ij = clonedMatrix(i, j);
            double var_j = clonedMatrix(j, j);
            double correlation_ij = cov_ij / sqrt(var_i * var_j);
            correlationMatrix(i, j) = correlation_ij;
        }
    }

    return correlationMatrix;
}

void slice_plots() {

  #ifdef USE_FAKE_DATA
    // Initialize the FilePropertiesManager and tell it to treat the NuWro
    // MC ntuples as if they were data
    auto& fpm = FilePropertiesManager::Instance();
    fpm.load_file_properties( "nuwro_file_properties.txt" );
  #endif

  auto* syst_ptr = new MCC9SystematicsCalculator(
//    "/exp/annie/data/users/jminock/stv-analysis/stv-univmake-output-20k.root",
    "../stv_output/output_0pi_closure.root",
    "systcalc.conf" );
  auto& syst = *syst_ptr;

 std::string pdf_type = "Closure";
 std::string Pdf_name  = "CrossSection_ANNIE_" + pdf_type;
 std::string Name_DATATYPE = "Fake Data";
 char pdf_title[1024];
 char Plot_title[1024];
 char axisXtitle[1024];

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

  // Total MC+EXT prediction in reco bin space. Start by getting EXT.
  TH1D* reco_mc_plus_ext_hist = dynamic_cast< TH1D* >(
    reco_ext_hist->Clone("reco_mc_plus_ext_hist") );
  reco_mc_plus_ext_hist->SetDirectory( nullptr );

  // Add in the CV MC prediction
  reco_mc_plus_ext_hist->Add( syst.cv_universe().hist_reco_.get() );

  // Keys are covariance matrix types, values are CovMatrix objects that
  // represent the corresponding matrices
  auto* matrix_map_ptr = syst.get_covariances().release();
  auto& matrix_map = *matrix_map_ptr;

  auto* sb_ptr = new SliceBinning( "tutorial_slice_config.txt" );
  auto& sb = *sb_ptr;

//  auto* cov_mat = matrix_map.at( "total" ).cov_matrix_.get();

  std::cout<<"Now examining real data."<<std::endl;

  double total_pot = syst_ptr->total_bnb_data_pot_;
  std::cout << "total_pot: " << total_pot << std::endl;
  double integ_flux_annie = integrated_numu_flux_in_FV( total_pot );
  std::cout << "integ_flux: " << integ_flux_annie << std::endl;
  double num_O = num_O_targets_in_FV(); //5.25e28;
  std::cout << "num_O: " << num_O << std::endl;
  double conv_factor = (num_O * integ_flux_annie)/1e38;
  std::cout << "conv_factor: " << conv_factor << std::endl;

 std::unique_ptr< Unfolder > unfolder (new WienerSVDUnfolder( true, WienerSVDUnfolder::RegularizationMatrixType::kFirstDeriv ));
 
 UnfoldedMeasurement result = unfolder->unfold( syst );
 
 // Propagate all defined covariance matrices through the unfolding procedure
 const TMatrixD& err_prop = *result.err_prop_matrix_;
 TMatrixD err_prop_tr( TMatrixD::kTransposed, err_prop );
 
 std::map< std::string, std::unique_ptr<TMatrixD> > unfolded_cov_matrix_map;
 
  for ( size_t sl_idx = 0u; sl_idx < sb.slices_.size(); ++sl_idx ) {

    const auto& slice = sb.slices_.at( sl_idx );

    // We now have all of the reco bin space histograms that we need as input.
    // Use them to make new histograms in slice space.
//    SliceHistogram* slice_bnb = SliceHistogram::make_slice_histogram(
//      *reco_bnb_hist, slice, &matrix_map.at("BNBstats") );

//    SliceHistogram* slice_ext = SliceHistogram::make_slice_histogram(
//      *reco_ext_hist, slice, &matrix_map.at("EXTstats") );

    SliceHistogram* slice_mc_plus_ext = SliceHistogram::make_slice_histogram(
      *reco_mc_plus_ext_hist, slice, &matrix_map.at("total") );

    //auto chi2_result = slice_bnb->get_chi2( *slice_mc_plus_ext );
    //std::cout << "Slice " << sl_idx << ": \u03C7\u00b2 = "
    //  << chi2_result.chi2_ << '/' << chi2_result.num_bins_ << " bins,"
    //  << " p-value = " << chi2_result.p_value_ << '\n';

    // Build a stack of categorized central-value MC predictions plus the
    // extBNB contribution in slice space
    const auto& eci = EventCategoryInterpreter::Instance();
//    eci.set_ext_histogram_style( slice_ext->hist_.get() );

    THStack* slice_pred_stack = new THStack( "mc+ext", "" );
//    slice_pred_stack->Add( slice_ext->hist_.get() ); // extBNB

    const auto& cat_map = eci.label_map();

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

      slice_pred_stack->Add( temp_slice_mc->hist_.get() );

      std::string cat_col_prefix = "MC" + std::to_string( cat );

      --cat_bin_index;
    }

    TCanvas* c1 = new TCanvas;
/*    slice_bnb->hist_->SetLineColor( kBlack );
    slice_bnb->hist_->SetLineWidth( 3 );
    slice_bnb->hist_->SetMarkerStyle( kFullCircle );
    slice_bnb->hist_->SetMarkerSize( 0.8 );
    slice_bnb->hist_->SetStats( false );
    double ymax = std::max( slice_bnb->hist_->GetMaximum(),
      slice_mc_plus_ext->hist_->GetMaximum() ) * 1.07;
    slice_bnb->hist_->GetYaxis()->SetRangeUser( 0., ymax );

    slice_bnb->hist_->Draw( "e" );

    slice_pred_stack->Draw( "hist same" );
*/
    slice_mc_plus_ext->hist_->SetLineWidth( 3 );
    slice_mc_plus_ext->hist_->Draw( "same hist e" );

//    slice_bnb->hist_->Draw( "same e" );

    //std::string out_pdf_name = "plot_slice_";
    //if ( sl_idx < 10 ) out_pdf_name += "0";
    //out_pdf_name += std::to_string( sl_idx ) + ".pdf";
    //c1->SaveAs( out_pdf_name.c_str() );

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
    const std::vector< std::string > cov_mat_keys = { "total","flux","xsec_total","MCstats","detVar_total" };
 
/*    const std::vector< std::string > cov_mat_keys = { "total",
      "detVar_total", "flux", "reint", "xsec_total", "POT", "numTargets",
      "MCstats", "EXTstats", "BNBstats"
    };
*/
    // Loop over the various systematic uncertainties

 SliceHistogram* slice_unf = SliceHistogram::make_slice_histogram(*result.unfolded_signal_, slice, result.cov_matrix_.get() );
  
 TH1D *h_unfoledData_clone_Error = (TH1D*)slice_unf->hist_->Clone(uniq());
 
 
 int color=1;
 double Max = 0.4;

 for ( auto CovMatrixName:cov_mat_keys  ) {
   
     const auto& key = CovMatrixName;

 std::cout<<"Inside: matrix_map_error :: key::"<< key<< std::endl;

 auto CovMatrix_ = unfolded_cov_matrix_map[ CovMatrixName ].get(); 
 std::cout<<"Inside: matrix_map_error :: key::"<< key<< std::endl;
 


  TMatrixD CorrelationMatrix =  covarianceToCorrelation(*CovMatrix_);
  
  sprintf(Plot_title, "Correlation Matrix: %s",key.c_str());
  sprintf(pdf_title, "%s.pdf", Pdf_name.c_str());
  visualizeMatrix(CorrelationMatrix, Plot_title,  pdf_title, "Bin N", "Bin N" , 1);


 TH1D* h_Error = (TH1D*)h_unfoledData_clone_Error->Clone(uniq());

  
   // The SliceHistogram object already set the bin errors appropriately
   // based on the slice covariance matrix. Just change the bin contents
   // for the current histogram to be fractional uncertainties. Also set
   // the "uncertainties on the uncertainties" to zero.
   // TODO: revisit this last bit, possibly assign bin errors here
   for ( const auto& bin_pair : slice.bin_map_ ) {
     int global_bin_idx = bin_pair.first;
     
     const auto& bin_set = bin_pair.second;
       int Bin_matrix; 
         for (size_t element : bin_set) {
         // Use the element from the set
         std::cout << "global_bin_idx= "<< global_bin_idx <<  " element: " << element << std::endl;
         Bin_matrix = element;
     }
     
      double y = h_unfoledData_clone_Error->GetBinContent( global_bin_idx );
     double width = h_unfoledData_clone_Error->GetBinWidth( global_bin_idx );
//     width *= other_var_width;
     Double_t element = (*CovMatrix_)(Bin_matrix, Bin_matrix);
     double Var1 = (element);
     //double Var = (Var1 * 1e38 * 1e38) / (integ_flux*integ_flux * num_Ar*num_Ar * width*width);
    double Var = ( Var1 * 1e38 * 1e38) / (integ_flux_annie*integ_flux_annie * num_O*num_O * width*width);
    
      double frac = 0.;
     if ( y > 0. ) frac = sqrt(Var) / y;
     if(Max<frac && frac < 1.0){Max = frac;}
     //std::cout<<" y =  "<< y << " var1 = " << Var1 << " Var = "<< Var<< "Var sqrted = "<< (Var*Var) << " Frac = " << frac << std::endl;
     h_Error->SetBinContent( global_bin_idx, frac );
     h_Error->SetBinError( global_bin_idx, 0. );
   }

   frac_uncertainty_hists[ key ] = h_Error;
   // Check whether the current covariance matrix name is present in
   // the vector defined above this loop. If it isn't, don't bother to
   // plot it, and just move on to the next one.
   auto cbegin = cov_mat_keys.cbegin();
   auto cend = cov_mat_keys.cend();
   auto iter = std::find( cbegin, cend, key );
   if ( iter == cend ) continue;
  

   if ( color == 3 ) color = kGreen +2;
   if ( color == kGreen + 3) color = 4;
   if ( color == 5 ) color = 6;
   if ( color == 7 ) color = kYellow -2;
   if ( color == kYellow -1 ) color = kMagenta + 2;
   if (color == kMagenta + 3) color = kSpring - 2;
   if (color == kSpring - 1) color = kOrange + 2;
   if (color == kOrange + 3) color = kRed - 6;
   if (color == kRed - 5) color = kAzure + 8;
   if (color == kAzure + 8) color = kOrange - 3;

  if(key == "BNBstats" ||
  key == "EXTstats" ||
  key == "MCstats" ||
  key == "DataStats"){
  frac_uncertainty_hists[ key ]->SetLineStyle(2);}


   frac_uncertainty_hists[ key ]->SetLineColor( color );
   frac_uncertainty_hists[ key ]->SetLineWidth( 4 );
   color++;
 }
  

/*
 for ( const auto& matrix_pair : matrix_map ) {
  const std::string& matrix_key = matrix_pair.first;
  auto temp_cov_mat = matrix_pair.second.get_matrix();
  
  TMatrixD temp_mat( *temp_cov_mat, TMatrixD::EMatrixCreatorsOp2::kMult,
  err_prop_tr );
 
  unfolded_cov_matrix_map[ matrix_key ] = std::make_unique< TMatrixD >(
  err_prop, TMatrixD::EMatrixCreatorsOp2::kMult, temp_mat );
 
  TMatrixD CovMatrix_new(err_prop, TMatrixD::EMatrixCreatorsOp2::kMult,
   temp_mat );
 
  CovMatrix_new *= (1e38 * 1e38) / (integ_flux_annie*integ_flux_annie * num_O*num_O); 
 
  std::string input_string = matrix_key;
 //TMatrixD newMatrix(*CovMatrix_); 
// unfolded_CovMatrixMap.insert(std::pair(input_string, CovMatrix_new));
  std::cout << input_string << std::endl;
/////////////////////////////////////////////////////

      const auto& key = input_string;
      const auto& cov_matrix = CovMatrix_new;
  
      TH1D* h_Error = (TH1D*)h_unfoledData_clone_Error->Clone(uniq());


   // The SliceHistogram object already set the bin errors appropriately
   // based on the slice covariance matrix. Just change the bin contents
   // for the current histogram to be fractional uncertainties. Also set
   // the "uncertainties on the uncertainties" to zero.
   // TODO: revisit this last bit, possibly assign bin errors here
   
      for ( const auto& bin_pair : slice.bin_map_ ) {
        int global_bin_idx = bin_pair.first;
        const auto& bin_set = bin_pair.second;
        int Bin_matrix; 
        for (size_t element : bin_set) {
          // Use the element from the set
          std::cout << "global_bin_idx= "<< global_bin_idx <<  " element: " << element << std::endl;
          Bin_matrix = element;
        }
     
        std::cout << "global_bin_idx= "<< global_bin_idx << "  Bin_matrix ="<<Bin_matrix<< std::endl;
        double y = h_unfoledData_clone_Error->GetBinContent( global_bin_idx );
        double width = h_unfoledData_clone_Error->GetBinWidth( global_bin_idx );
        Double_t element = (*CovMatrix_)(Bin_matrix, Bin_matrix);
        double Var1 = (element);

        double factor1 = (integ_flux_annie*integ_flux_annie * num_O*num_O * width*width) ;
        std::cout<< "Var1 = "<<Var1  <<"  width ="<<width <<"   (integ_flux_annie*integ_flux_annie * num_O*num_O * width*width) = " << factor1<<std::endl;
        double Var = ( Var1 * 1e38 * 1e38) / (integ_flux_annie*integ_flux_annie * num_O*num_O * width*width);
    
        double frac = 0.;
        if ( y > 0. ) frac = sqrt(Var) / y;
        if(Max<frac && frac < 1.0){Max = frac;}

        //std::cout<<" y =  "<< y << " var1 = " << Var1 << " Var = "<< Var<< "Var sqrted = "<< (Var*Var) << " Frac = " << frac << std::endl;
        h_Error->SetBinContent( global_bin_idx, frac );
        h_Error->SetBinError( global_bin_idx, 0. );
      }

      frac_uncertainty_hists[ key ] = h_Error;
      // Check whether the current covariance matrix name is present in
      // the vector defined above this loop. If it isn't, don't bother to
      // plot it, and just move on to the next one.
      auto cbegin = cov_mat_keys.cbegin();
      auto cend = cov_mat_keys.cend();
      auto iter = std::find( cbegin, cend, key );
      if ( iter == cend ) continue;

      if ( color == 3 ) color = kGreen +2;
      if ( color == kGreen + 3) color = 4;
      if ( color == 5 ) color = 6;
      if ( color == 7 ) color = kYellow -2;
      if ( color == kYellow -1 ) color = kMagenta + 2;
      if (color == kMagenta + 3) color = kSpring - 2;
      if (color == kSpring - 1) color = kOrange + 2;
      if (color == kOrange + 3) color = kRed - 6;
      if (color == kRed - 5) color = kAzure + 8;
      if (color == kAzure + 8) color = kOrange - 3;

      if(key == "BNBstats" || key == "EXTstats" || key == "MCstats" || key == "DataStats"){
        frac_uncertainty_hists[ key ]->SetLineStyle(2);
      }

      frac_uncertainty_hists[ key ]->SetLineColor( color );
      if(key == "DataStats"){ frac_uncertainty_hists[ key ]->SetLineColor( 28 );}
      frac_uncertainty_hists[ key ]->SetLineWidth( 3 );
      color++;
    }

    TCanvas* c2 = new TCanvas;
    TLegend* lg2 = new TLegend( 0.7, 0.7, 0.9, 0.9 );

    auto* total_frac_err_hist = frac_uncertainty_hists.at( "total" );
    total_frac_err_hist->SetStats( false );
    total_frac_err_hist->GetYaxis()->SetRangeUser( 0.,
      total_frac_err_hist->GetMaximum() * 1.05 );
    total_frac_err_hist->SetLineColor( kBlack );
    total_frac_err_hist->SetLineWidth( 3 );
    total_frac_err_hist->Draw( "hist" );

    lg2->AddEntry( total_frac_err_hist, "total", "l" );

    for ( auto& pair : frac_uncertainty_hists ) {
      const auto& name = pair.first;
      TH1* hist = pair.second;
      // We already plotted the "total" one above
      if ( name == "total" ) continue;

      lg2->AddEntry( hist, name.c_str(), "l" );
      hist->Draw( "same hist" );

      std::cout << name << " frac err in bin #1 = "
        << hist->GetBinContent( 1 )*100. << "%\n";
    }

    lg2->Draw( "same" );

    std::cout << "Total frac error in bin #1 = "
      << total_frac_err_hist->GetBinContent( 1 )*100. << "%\n";
*/
  } // slices

}

int main() {
  slice_plots();
  return 0;
}
