/// Set of tool for plotting working working ANNIE's geometry
// Aurther :: christian nguyen 

#include "AnnieGeometryTools.hh"

TGraph* createBoundaryGraph(const std::string& plotType) {
    double radius = FV_RAD;  // cm
    double y_min = A_FV_Y_MIN;  // cm
    double y_max = A_FV_Y_MAX;   // cm
    double z_center = A_Z_CTR;  // cm

    std::vector<Vertex_XYZ> VertexPostion;

    int n_points = 10000;  // Number of points to plot

    // Generate points along the boundary
    for (int i = 0; i < n_points; ++i) {
        double angle = 2 * M_PI * i / n_points;
        double x = radius * std::cos(angle);
        double z = radius * std::sin(angle) + z_center;
        double y = y_min + (y_max - y_min) * i / n_points;

          if (y > y_min && y < y_max && radius > std::sqrt((z - z_center) * (z - z_center) + x * x) && (z - z_center < radius)) {
              Vertex_XYZ input{x,y,z};
              VertexPostion.push_back(input);
              std::cout<<"(x,y,z) = "<< "(" << x<< ", "<< y << ", "<< z<< " ) "<<std::endl;
              
            }
         }

        if (plotType == "X vs Y") {
            TGraph  *TGraph_output = Make_X_vs_Y_Tgraph_fromVector(VertexPostion);
            return TGraph_output;
        }
        else if (plotType == "X vs Z") {
            TGraph  *TGraph_output = Make_X_vs_Z_Tgraph_fromVector(VertexPostion);
            return TGraph_output;
        }
        else if (plotType == "Y vs Z") {
            TGraph  *TGraph_output = Make_Y_vs_Z_Tgraph_fromVector(VertexPostion);
            return TGraph_output;
        }
        else if (plotType == "R vs Y") {
            TGraph  *TGraph_output = Make_R_vs_Y_Tgraph_fromVector(VertexPostion);
            return TGraph_output;
        }
        else if (plotType == "R vs X") {
            TGraph  *TGraph_output = Make_R_vs_X_Tgraph_fromVector(VertexPostion);
            return TGraph_output;
        }
        else if (plotType == "R vs Z") {
            TGraph  *TGraph_output = Make_R_vs_Z_Tgraph_fromVector(VertexPostion);
            return TGraph_output;
        }
        else if (plotType == "RR vs Z") {
            TGraph  *TGraph_output = Make_RR_vs_Z_Tgraph_fromVector(VertexPostion);
            return TGraph_output;
        }
        else if (plotType == "RR vs X") {
            TGraph  *TGraph_output = Make_RR_vs_X_Tgraph_fromVector(VertexPostion);
            return TGraph_output;
        }
         else if (plotType == "RR vs Y") {
            TGraph  *TGraph_output = Make_RR_vs_Y_Tgraph_fromVector(VertexPostion);
            return TGraph_output;
        }
       else{TGraph  *TGraph_output =nullptr;\
       return TGraph_output;
       }
    
}
