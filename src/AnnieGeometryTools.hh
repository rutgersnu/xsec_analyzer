// Header of geometry tool 
//// Hist tools 
/// Christian Nguyen 
#pragma once

#include <algorithm>
#include <cmath>

// ROOT includes
#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "THStack.h"
#include "TLegend.h"
#include "TObjArray.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include <TGraph2D.h>

// STV analysis includes

#include "TLatex.h"
#include "TLine.h"
#include "TH2Poly.h"
#include "plotutils/UBTH2Poly.h" 
#include "plotutils/GridCanvas.hh"
 
#include "NamedCategory.hh"
#include "plotutils/HistUtils.h"
#include "TVector3.h"


TGraph* createBoundaryGraph(const std::string& plotType);

//ANNIE FV Boundaries
double FV_RAD = 100.0;
double A_FV_Y_MIN = -100.0;
double A_FV_Y_MAX =  100.0;
double A_Z_CTR = 168.1; //Z center of tank

template <typename Number> bool point_inside_AFV( Number x, Number y, Number z )
{
  if(y > A_FV_Y_MIN && y < A_FV_Y_MAX && FV_RAD > std::sqrt((z - A_Z_CTR)*(z - A_Z_CTR) + x*x) && (z - A_Z_CTR < FV_RAD)){
    return true;
  }
  else {
    return false;
  }
}

inline bool point_inside_AFV( const TVector3& pos ) {
  return point_inside_AFV( pos.X(), pos.Y(), pos.Z() );
}

// Returns the number of O nuclei inside the fiducial volume
inline double num_O_targets_in_FV() {
  double volume = FV_RAD * FV_RAD * M_PI * ( A_FV_Y_MAX - A_FV_Y_MIN ); // cm^3
  constexpr double m_mol_H2O = 18.015; // g/mol
  constexpr double N_Avogadro = 6.02214076e23; // mol^(-1)
  constexpr double mass_density_H2O = 1.; // g/cm^3

  double num_O = volume * mass_density_H2O * N_Avogadro / m_mol_H2O;
  return num_O;
}

// Returns the total BNB muon neutrino flux (numu / cm^2) in the fiducial
// volume as a function of a given beam exposure (measured in
// protons-on-target)
// NOTE: This is currently approximated using the flux in the *active volume*.
// TODO: Revisit this approximation
inline double integrated_numu_flux_in_FV( double pot ) {
  // Obtained using the histogram hEnumu_cv (the central-value numu flux
  // in the MicroBooNE active volume as a function of neutrino energy)
  // stored in /pnfs/uboone/persistent/uboonebeam/bnb_gsimple
  // /bnb_gsimple_fluxes_01.09.2019_463_hist/. The ROOT commmands executed
  // were
  // root [3] hEnumu_cv->Scale( 1/(4997.*5e8)/(256.35*233.) )
  // root [4] hEnumu_cv->Integral()
  // (double) 7.3762291e-10
  // See the README file in that same folder for details.
  constexpr double numu_per_cm2_per_POT_in_AV = 2.26256e-8;
  double flux = pot * numu_per_cm2_per_POT_in_AV; // numu / cm^2
  return flux;
}
