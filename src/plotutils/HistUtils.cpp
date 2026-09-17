// Christian Nguyen 
#include "HistUtils.h"


TH1D* GetTH1DHist(TFile& fin, const char* name) {

  TH1D* h=(TH1D*)fin.Get(name);
  if (h==0) {
      const char* fileName = fin.GetName();
  std::cout << "Could not get 1D TH1D hist :  " << name << "  From TFILE named : "<<fileName << "\n";
  
  }
  return h;
}//////////////////////////////// End of Function 
////////////////////////////////////////////////////////////////////////////////
TH1D* GetTH1DHist(TFile& fin, std::string name ){
  const char* name_char = name.c_str();
  TH1D* hist = GetTH1DHist(fin, name_char);
  return hist;
}//////////////////////////////// End of Function 
////////////////////////////////////////////////////////////////////////////////
TH2D* GetTH2DHist(TFile& fin, const char* name) {
  TH2D* h=(TH2D*)fin.Get(name);
  if (h==0) std::cout << "Could not get 2D TH2D hist " << name << "\n";
  return h;
}//////////////////////////////// End of Function 
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
TH2D* GetTH2DHist(TFile& fin, std::string name) {
  const char* name_char = name.c_str();
  TH2D* hist = GetTH2DHist(fin, name_char);
  return hist;
}//////////////////////////////// End of Function 

////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
TGraph  *Make_X_vs_Y_Tgraph_fromVector(std::vector<Vertex_XYZ> input)
{
  const int size = input.size();
  double YY[size];
  double XX[size];

  for(unsigned int j =0; j<size;j++){
  XX[j]=input.at(j).x;
  YY[j]=input.at(j).y;
  }

  TGraph *Tg_result = new TGraph(size,XX,YY);
  
   Tg_result->GetXaxis()->CenterTitle();
   Tg_result->GetYaxis()->CenterTitle();
   Tg_result->GetXaxis()->SetTitle("Z [cm]");
   Tg_result->GetYaxis()->SetTitle("Y [cm]");
   Tg_result->GetXaxis()->SetTitleSize(0.038);
   Tg_result->GetYaxis()->SetTitleSize(0.038);
   Tg_result->SetLineColor(2);
   Tg_result->SetLineWidth(4);
   Tg_result->SetMarkerColor(1);
   Tg_result->SetMarkerSize(1);
   Tg_result->SetMarkerStyle(20);
  

    return Tg_result;
}
////////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
TGraph  *Make_X_vs_Z_Tgraph_fromVector(std::vector<Vertex_XYZ> input)
{
  const int size = input.size();
  double ZZ[size];
  double XX[size];

  for(unsigned int j =0; j<size;j++){
  XX[j]=input.at(j).x;
  ZZ[j]=input.at(j).z;
  }

  TGraph *Tg_result = new TGraph(size,ZZ,XX);
   Tg_result->GetXaxis()->CenterTitle();
   Tg_result->GetYaxis()->CenterTitle();
   Tg_result->GetXaxis()->SetTitle("Z [cm]");
   Tg_result->GetYaxis()->SetTitle("X [cm]");
   Tg_result->GetXaxis()->SetTitleSize(0.038);
   Tg_result->GetYaxis()->SetTitleSize(0.038);
   Tg_result->SetLineColor(2);
   Tg_result->SetLineWidth(4);
   Tg_result->SetMarkerColor(1);
   Tg_result->SetMarkerSize(1);
   Tg_result->SetMarkerStyle(20);
   

    return Tg_result;
}
////////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
TGraph  *Make_Y_vs_Z_Tgraph_fromVector(std::vector<Vertex_XYZ>input)
{

  const int size = input.size();
  double YY[size];
  double ZZ[size];

  for(unsigned int j =0; j<size;j++){
  YY[j]=input.at(j).y;
  ZZ[j]=input.at(j).z;
  }

  TGraph *Tg_result = new TGraph(size,ZZ,YY);


  Tg_result->GetXaxis()->CenterTitle();
  Tg_result->GetYaxis()->CenterTitle();
  Tg_result->GetXaxis()->SetTitle("Z [mm]");
  Tg_result->GetYaxis()->SetTitle("Y [mm]");
  Tg_result->GetXaxis()->SetTitleSize(0.038);
  Tg_result->GetYaxis()->SetTitleSize(0.038);
  Tg_result->SetLineColor(2);
  Tg_result->SetLineWidth(4);
  Tg_result->SetMarkerColor(1);
  Tg_result->SetMarkerSize(1);
  Tg_result->SetMarkerStyle(20);


    return Tg_result;
}
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
TGraph  *Make_R_vs_Z_Tgraph_fromVector(std::vector<Vertex_XYZ>input)
{
  const int size = input.size();
  double ZZ[size];
  double RR[size];

  for(unsigned int j =0; j<size;j++){
  RR[j]=sqrt(pow(input.at(j).x,2)+pow(input.at(j).z,2));
  ZZ[j]=input.at(j).z;
  }

  TGraph *Tg_result = new TGraph(size,ZZ,RR);

  Tg_result->GetXaxis()->CenterTitle();
  Tg_result->GetYaxis()->CenterTitle();
  Tg_result->GetXaxis()->SetTitle("Z [mm]");
  Tg_result->GetYaxis()->SetTitle("R [mm]");
  Tg_result->GetXaxis()->SetTitleSize(0.038);
  Tg_result->GetYaxis()->SetTitleSize(0.038);
  Tg_result->SetLineColor(2);
  Tg_result->SetLineWidth(4);
  Tg_result->SetMarkerColor(1);
  Tg_result->SetMarkerSize(1);
  Tg_result->SetMarkerStyle(20);
    return Tg_result;

}
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////

TGraph  *Make_RR_vs_Z_Tgraph_fromVector(std::vector<Vertex_XYZ>input)
{

  const int size = input.size();
  double ZZ[size];
  double RR[size];

  for(unsigned int j =0; j<size;j++){
  RR[j]=pow(input.at(j).z,2) + pow(input.at(j).x,2);
  ZZ[j]=input.at(j).z;
  }

  TGraph *Tg_result = new TGraph(size,ZZ,RR);

  Tg_result->GetXaxis()->CenterTitle();
  Tg_result->GetYaxis()->CenterTitle();
  Tg_result->GetXaxis()->SetTitle("Z [mm]");
  Tg_result->GetYaxis()->SetTitle("R^{2} [mm]");
  Tg_result->GetXaxis()->SetTitleSize(0.038);
  Tg_result->GetYaxis()->SetTitleSize(0.038);
  Tg_result->SetLineColor(2);
  Tg_result->SetLineWidth(4);
  Tg_result->SetMarkerColor(1);
  Tg_result->SetMarkerSize(1);
  Tg_result->SetMarkerStyle(20);
    
    return Tg_result;
}
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////

TGraph  *Make_RR_vs_Y_Tgraph_fromVector(std::vector<Vertex_XYZ>input)
{

  const int size = input.size();
  double YY[size];
  double RR[size];

  for(unsigned int j =0; j<size;j++){
  RR[j]=pow(input.at(j).z,2) + pow(input.at(j).x,2);
  YY[j]=input.at(j).y;
  }

  TGraph *Tg_result = new TGraph(size,YY,RR);

  Tg_result->GetXaxis()->CenterTitle();
  Tg_result->GetYaxis()->CenterTitle();
  Tg_result->GetXaxis()->SetTitle("T [cm]");
  Tg_result->GetYaxis()->SetTitle("R^{2} [cm]");
  Tg_result->GetXaxis()->SetTitleSize(0.038);
  Tg_result->GetYaxis()->SetTitleSize(0.038);
  Tg_result->SetLineColor(2);
  Tg_result->SetLineWidth(4);
  Tg_result->SetMarkerColor(1);
  Tg_result->SetMarkerSize(1);
  Tg_result->SetMarkerStyle(20);


    return Tg_result;
}
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
TGraph  *Make_R_vs_Y_Tgraph_fromVector(std::vector<Vertex_XYZ>input)
{

  const int size = input.size();
  double YY[size];
  double RR[size];

  for(unsigned int j =0; j<size;j++){
  RR[j]=sqrt(pow(input.at(j).z,2) + pow(input.at(j).x,2));
  YY[j]=input.at(j).y;
  }

  TGraph *Tg_result = new TGraph(size,YY,RR);

  Tg_result->GetXaxis()->CenterTitle();
  Tg_result->GetYaxis()->CenterTitle();
  Tg_result->GetXaxis()->SetTitle("Y [mm]");
  Tg_result->GetYaxis()->SetTitle("R [mm]");
  Tg_result->GetXaxis()->SetTitleSize(0.038);
  Tg_result->GetYaxis()->SetTitleSize(0.038);
  Tg_result->SetLineColor(2);
  Tg_result->SetLineWidth(4);
  Tg_result->SetMarkerColor(1);
  Tg_result->SetMarkerSize(1);
  Tg_result->SetMarkerStyle(20);
  Tg_result->SetMinimum(0);

    return Tg_result;
}
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////

TGraph  *Make_RR_vs_X_Tgraph_fromVector(std::vector<Vertex_XYZ>input)
{

  const int size = input.size();
  double X[size];
  double RR[size];

  for(unsigned int j =0; j<size;j++){
  RR[j]=pow(input.at(j).z,2) + pow(input.at(j).x,2);
  X[j]=input.at(j).x;
  }

  TGraph *Tg_result = new TGraph(size,X,RR);

  Tg_result->GetXaxis()->CenterTitle();
  Tg_result->GetYaxis()->CenterTitle();
  Tg_result->GetXaxis()->SetTitle("X [mm]");
  Tg_result->GetYaxis()->SetTitle("R^{2} [mm]");
  Tg_result->GetXaxis()->SetTitleSize(0.038);
  Tg_result->GetYaxis()->SetTitleSize(0.038);
  Tg_result->SetLineColor(2);
  Tg_result->SetLineWidth(4);
  Tg_result->SetMarkerColor(1);
  Tg_result->SetMarkerSize(1);
  Tg_result->SetMarkerStyle(20);

    return Tg_result;
}
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////

TGraph  *Make_R_vs_X_Tgraph_fromVector(std::vector<Vertex_XYZ>input)
{

  const int size = input.size();
  double XX[size];
  double RR[size];

  for(unsigned int j =0; j<size;j++){
  RR[j]=sqrt(pow(input.at(j).z,2) + pow(input.at(j).x,2));
  XX[j]=input.at(j).x;
  }

  TGraph *Tg_result = new TGraph(size,XX,RR);

  Tg_result->GetXaxis()->CenterTitle();
  Tg_result->GetYaxis()->CenterTitle();
  Tg_result->GetXaxis()->SetTitle("X [cm]");
  Tg_result->GetYaxis()->SetTitle("R [cm]");
  Tg_result->GetXaxis()->SetTitleSize(0.038);
  Tg_result->GetYaxis()->SetTitleSize(0.038);
  Tg_result->SetLineColor(2);
  Tg_result->SetLineWidth(4);
  Tg_result->SetMarkerColor(1);
  Tg_result->SetMarkerSize(1);
  Tg_result->SetMarkerStyle(20);


    return Tg_result;
}
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////