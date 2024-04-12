/**
 * CC0pi reco corrections.
 *
 * ATM 2024/03/13
 * Based on Christian's Version 1/10/24
 * Panos's corrections, see MicroBooNE-doc-41275
 */

#ifndef __PMU_CORR__
#define __PMU_CORR__

class PMuCorrection {
public:

  static double correct(double p, bool contained, bool mu_quality) {
    if (contained) {
      return corr_contained(p, mu_quality);
    }
    else {
      return corr_uncontained(p);
    }
  }

protected:

  static double corr_contained(double p, bool mu_quality) {
    if (p <= 0.2 && !mu_quality) {
      return 0;  // FIXME
    }

    else if (p <= 0.2 && mu_quality) {
      return poly2(p, cR1_A, cR1_B, cR1_C);
    }
        
    else if (p > 0.2 && p < 1.5) {
      return poly2(p, cR2_A, cR2_B, cR2_C);
    }
                    
    else if (p >= 1.5) {
      return poly3(p, cR3_A, cR3_B, cR3_C, cR3_D);
    }

    else {
      std::cerr << "PmuCorrection::corr_contained: p out of range." << std::endl;
      assert(false);
    }
  }

  static double corr_uncontained(double p) {
    if (p <= 0.2) {
      return 0.0;  // FIXME
    }    

    else if (p > 0.2 && p < 1.2) {
      return poly2(p, uR1_A, uR1_B, uR1_C);
    }
                    
    else if (p >= 1.2 && p < 2.2 ) {
      return poly2(p, uR2_A, uR2_B, uR2_C);
    }

    else if ( 2.2 <= p ){
      return poly2(p, uR3_A, uR3_B, uR3_C);
    }

    else {
      std::cerr << "PmuCorrection::corr_uncontained: p out of range." << std::endl;
      assert(false);
    }
  }

private:

  static inline double poly2(double x, double A, double B, double C) {
   return A + B * x + C * x * x;
  }

  static inline double poly3(double x, double A, double B, double C, double D) {
    return A + B * x + C * x * x + D * x * x * x; 
  }

  // Fit parameters
  static constexpr double cR1_A =  0.59;
  static constexpr double cR1_B = -6.89;
  static constexpr double cR1_C = 20.63;
  static constexpr double cR2_A =  0.03;
  static constexpr double cR2_B = -0.010;
  static constexpr double cR2_C = -0.014;
  static constexpr double cR3_A =  1.96;
  static constexpr double cR3_B = -3.35;
  static constexpr double cR3_C =  1.97;
  static constexpr double cR3_D = -0.41;
  
  static constexpr double uR1_A =  0.18;
  static constexpr double uR1_B = -0.15;
  static constexpr double uR1_C = -0.02;
  static constexpr double uR2_A = -0.98;
  static constexpr double uR2_B =  1.63;
  static constexpr double uR2_C = -0.73;
  static constexpr double uR3_A =  1.9;
  static constexpr double uR3_B = -1.34;
  static constexpr double uR3_C =  0.03;
  static constexpr double uR4_A =  3.42119;
  static constexpr double uR4_B = -2.19198;
  static constexpr double uR4_C =  0.160376;
};

#endif  // __PMU_CORR__

