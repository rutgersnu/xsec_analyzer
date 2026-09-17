//// Hist tools 
/// Christian Nguyen 
#pragma once

#include <algorithm>

// ROOT includes
#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "THStack.h"
#include "TLegend.h"
#include "TObjArray.h"
#include "TGraphErrors.h"
#include <TGraph2D.h>
#include <TGraph.h>
// STV analysis includes

//#include "SliceBinning.hh"
//#include "SliceHistogram.hh"
#include "TLatex.h"
#include "TLine.h"
#include "TH2Poly.h"
#include "UBTH2Poly.h" 
#include "GridCanvas.hh"
 
//#include "ConfigMakerUtils.hh"
#include "EventCategory.hh"
#include "NamedCategory.hh"
//#include "PlotUtils.hh"
#include "HistFolio_slim.hh"

////////////////////////////////////////////////////////////////////////////////
class Vertex_XYZ {
 public:
  double x;
  double y;
  double z;

  Vertex_XYZ(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z) {}

  bool operator==(const Vertex_XYZ& vector) const {
    return (x == vector.x && y == vector.y && z == vector.z);
  }

  Vertex_XYZ& operator=(const Vertex_XYZ& vector) {
    x = vector.x;
    y = vector.y;
    z = vector.z;
    return *this;
  }

  Vertex_XYZ operator+(const Vertex_XYZ& input) const {
    return Vertex_XYZ(x + input.x, y + input.y, z + input.z);
  }

  Vertex_XYZ operator-(const Vertex_XYZ& input) const {
    return Vertex_XYZ(x - input.x, y - input.y, z - input.z);
  }

  double DotProduct(const Vertex_XYZ& vector1, const Vertex_XYZ& vector2) {
    return vector1.x * vector2.x + vector1.y * vector2.y +
           vector1.z * vector2.z;
  }

  Vertex_XYZ CrossProduct(const Vertex_XYZ& vector1,
                          const Vertex_XYZ& vector2) {
    Vertex_XYZ inputVector;
    inputVector.x = vector1.y * vector2.z - vector1.z * vector2.y;
    inputVector.y = vector1.z * vector2.x - vector1.x * vector2.z;
    inputVector.z = vector1.x * vector2.y - vector1.y * vector2.x;
    return inputVector;
  }
};
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
TH1D* GetTH1DHist(TFile& fin, const char* name);
////////////////////////////////////////////////////////////////////////////////
TH1D* GetTH1DHist(TFile& fin, std::string name );
////////////////////////////////////////////////////////////////////////////////
TH2D* GetTH2DHist(TFile& fin, const char* name) ;
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
TH2D* GetTH2DHist(TFile& fin, std::string name) ;
////////////////////////////////////////////////////////////////////////////////

template <typename EnumType>
std::map<EnumType, TH1D*> MakeNamedCategoryMap_TH1D(TFile& fin, 
std::string histBase_name, const std::vector<NamedCategory<EnumType>> &categories )
{
std::map<EnumType, TH1D*> OutMap_result;
char Histtitle[1024];
   std::cout << "Processing categories:" << std::endl;
    for (const auto& namedCategory : categories) {
    std::string HistName= histBase_name + "_" + namedCategory.m_name;
    const char* name_char = namedCategory.m_name.c_str();
     const char* HISTname_char = HistName.c_str();
     TH1D* HistInput=(TH1D*)fin.Get(HISTname_char);
    HistInput->SetTitle(name_char);
    OutMap_result[namedCategory.m_value]=HistInput;
     std::cout << "Category: " << static_cast<int>(namedCategory.m_value) << ", Name: " <<  HistName << std::endl;
     
    }

return OutMap_result;

}//////////////////////////////// End of Function
TGraph  *Make_X_vs_Z_Tgraph_fromVector(std::vector<Vertex_XYZ> input);
////////////////////////////////////////////////////////////////////////////////
TGraph  *Make_Y_vs_Z_Tgraph_fromVector(std::vector<Vertex_XYZ>input);
////////////////////////////////////////////////////////////////////////////////
TGraph  *Make_R_vs_Z_Tgraph_fromVector(std::vector<Vertex_XYZ>input);
////////////////////////////////////////////////////////////////////////////////
TGraph  *Make_RR_vs_Z_Tgraph_fromVector(std::vector<Vertex_XYZ>input);
////////////////////////////////////////////////////////////////////////////////
TGraph  *Make_R_vs_Y_Tgraph_fromVector(std::vector<Vertex_XYZ>input);
////////////////////////////////////////////////////////////////////////////////
TGraph  *Make_RR_vs_Y_Tgraph_fromVector(std::vector<Vertex_XYZ>input);
////////////////////////////////////////////////////////////////////////////////
TGraph  *Make_RR_vs_X_Tgraph_fromVector(std::vector<Vertex_XYZ>input);
////////////////////////////////////////////////////////////////////////////////
TGraph  *Make_R_vs_X_Tgraph_fromVector(std::vector<Vertex_XYZ>input);
////////////////////////////////////////////////////////////////////////////////
TGraph  *Make_X_vs_Y_Tgraph_fromVector(std::vector<Vertex_XYZ> input);
