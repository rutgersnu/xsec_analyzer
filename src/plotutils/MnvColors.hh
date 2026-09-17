//#ifndef MNVCOLORS_HH
//#define MNVCOLORS_HH

#pragma once
//#ifndef __CINT__
#include <cassert>
#include <iostream>
#include <vector>

#include "TClass.h"
#include "TColor.h"
#include "TH1.h"
#include "TLatex.h"
#include "TList.h"
#include "TPad.h"
#include "TString.h"

namespace MnvColors {
enum EColorPalettes {
  kAlphabetPalette,
  kKellyPalette,
  k36Palette,
  kGlasbeyPalette,
  kBrewerSet1Palette,
  kBrewerDark2Palette,
  kTolBrightPalette,
  kTolMutedPalette,
  kTolLightPalette,
  kOkabeItoPalette,
  kOkabeItoLightPalette,
  kOkabeItoDarkPalette,
  kHeliumPalette,
  kNPalettes,
  kCCZeroPion
};

const std::vector<int>& GetColors(int palette = k36Palette);
void AutoColorHists(TPad* pad, int palette = k36Palette);
std::vector<TH1*> GetPadHists(TPad* pad);
const char* GetPaletteName(EColorPalettes e = k36Palette);
}  // namespace MnvColors

//#endif  // __CINT__
//#endif  // MNVCOLORS_H