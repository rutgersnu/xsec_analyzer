#pragma once

// ANNIE nuNCQE event selection for the xsec_analyzer framework.
//
// Line-by-line port of the Python selection that produced the PRD-draft cross section
// (NCQE_XS_analysis_v3, branch adopted-extraction):
//   lib/selection.py       mc_select_file(), core_cuts(), in_bunch_window(), truth_masks(), _categories()
//   lib/reco_functions.py  reco_energy(), pairwise_relative_direction(), is_inside_rotated_parabola()
//   lib/constants.py       bunch centres, bunch sigma and offset, FV and tank dimensions
// Inequalities, cut order and floating-point operations follow the Python code so that the selected
// events are identical event by event. Two places need care to reproduce NumPy bit for bit:
//   * NumPy scalar x**2 calls libm pow(x, 2), so std::pow is used where the Python code squares a scalar
//     (arrays squared with ** use x*x, and so does this port);
//   * np.linalg.norm of a 3-vector is a BLAS dot product, sqrt(fma(z,z, fma(y,y, x*x))).
// Compile with -ffp-contract=off so that the compiler does not fuse other multiply-adds (NumPy does not).

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <vector>

namespace ncqe {

// ---------------------------------------------------------------------------------- geometry [m]
inline constexpr double kTankRadius = 1.524;
inline constexpr double kTankHalfHeight = 1.98;
inline constexpr double kFVRadius = 0.7;
inline constexpr double kFVHalfHeight = 1.0;

// Offsets applied to the true vertex (cm, detector coordinates) to centre it on the tank (m)
inline constexpr double kTrueYOffset = 0.1446;
inline constexpr double kTrueZOffset = 1.681;

// ---------------------------------------------------------------------------------- bunch timing [ns]
inline constexpr double kBunchSigma = 3.5;      // nuNCQE bunch sigma
inline constexpr double kBunchOffset = 0.785;   // simulated CCinc -> nuNCQE bunch offset

// CCinc bunch centres fitted in MC and in data (https://github.com/S81D/BNB_bunches)
inline constexpr std::array<double, 81> kMCCenters = {
  31.080241331783608, 49.86887398044972, 68.93296122138285, 87.99776529989828,
  106.92835749210317, 125.6744929249585, 144.55074229641292, 163.55533862061978,
  182.42054838303565, 201.3975667029726, 220.45443080446125, 239.35228866298255,
  258.3293431484236, 277.16438785242855, 296.09720802748456, 315.07114394995256,
  334.03324394912823, 352.9587576490855, 371.8466462675708, 390.79194334882425,
  409.63044291467685, 428.65567092584, 447.6549629320083, 466.57352377409313,
  485.51073635902037, 504.4792057236341, 523.4956574718212, 542.3691320377058,
  561.1099950563898, 580.2327318476815, 599.219010807661, 618.1019306658702,
  636.9854167698654, 655.951836439576, 674.9453451130016, 693.8124506487408,
  712.671642937589, 731.7183509540064, 750.6217703640303, 769.543817927331,
  788.6698605766949, 807.5171007720303, 826.3570279447015, 845.2861546004002,
  864.1547464834093, 883.1283566641179, 902.1248645210509, 920.9586744177927,
  939.8944893020667, 958.8604582494381, 977.9115116970311, 996.7815728120514,
  1015.6916816004485, 1034.593675449364, 1053.6218603405582, 1072.4974885887596,
  1091.536793059606, 1110.3787443226672, 1129.5039961808134, 1148.208660322583,
  1167.3275787742382, 1186.2991638425285, 1204.9871507926803, 1224.054624008908,
  1242.917134430439, 1261.9481194256011, 1280.8989737693325, 1299.7479318121127,
  1318.6949768605764, 1337.6816835239597, 1356.5741533701105, 1375.6189574036052,
  1394.2645061208234, 1413.405519458969, 1432.3530855833285, 1451.3251805522368,
  1470.3698038901168, 1489.1695033456674, 1508.0806190832418, 1526.995246745324,
  1546.0600858753075
};

inline constexpr std::array<double, 81> kDataCenters = {
  196.56224433440246, 215.66909099116478, 234.50541331974833, 252.77573081775793,
  271.82317272128086, 290.5525519751267, 309.47529727366754, 328.30485229514017,
  347.30963148293534, 366.2956121865753, 385.226842076395, 404.22901923484017,
  423.1284015795741, 442.1028957843268, 460.9563022567543, 479.9292098671415,
  498.67523613302967, 517.745556275794, 536.5634329571434, 555.5826172158363,
  574.5584125542986, 593.3412403321598, 612.4668820959361, 631.1017157126897,
  650.245272991258, 669.1190067717141, 688.0583663595144, 707.0581456469093,
  726.0269162795082, 744.773346737428, 763.7192011932474, 782.8556702712547,
  801.7117975261683, 820.6641961303449, 839.4297222624705, 858.4883851813762,
  877.452957835381, 896.3449195713947, 915.2971258549227, 934.355709628352,
  953.1996719890523, 972.2816369135228, 991.1165854922059, 1009.9643936927976,
  1028.951234108228, 1047.9972672500965, 1066.8793045738905, 1085.8351450523521,
  1104.763384381811, 1123.7468193612021, 1142.803255690862, 1161.5532573658047,
  1180.4810843925607, 1199.5057688888362, 1218.302045446269, 1237.1624520071487,
  1256.1779033910852, 1275.2306259397883, 1294.0361634248936, 1313.1200477268103,
  1331.9742215964689, 1350.8189134298868, 1369.786773197136, 1388.6371303935912,
  1407.6174152883736, 1426.6042835718315, 1445.7252579081544, 1464.4532145609044,
  1483.2454568431624, 1502.377311593686, 1521.258000819768, 1540.0844404855327,
  1559.0588621416, 1578.0394188515295, 1596.9215520806429, 1615.9925479998735,
  1635.0135581066502, 1653.8237509076696, 1672.8362535130436, 1691.772995793539,
  1710.7097762845738
};

// Derived times, evaluated in the same order as lib/constants.py
inline const double kBunchTimeStart = kMCCenters.front() - kBunchOffset - kBunchSigma;     // 26.795 ns
inline const double kBunchTimeEnd = kMCCenters.back() - kBunchOffset + kBunchSigma;        // 1548.775 ns
// MC time of the end of the 2 us prompt buffer, mimicking data
inline const double kTimeToPromptEnd = (2000. - kDataCenters.back()) + kMCCenters.back();  // 1835.35 ns
// Cut keeping the first 26 bunches (reduces skyshine)
inline const double kMCThirdSpillCut = kMCCenters.at(25) - kBunchOffset + kBunchSigma;    // 507.194 ns

// ---------------------------------------------------------------------------------- External (dirt) scaling
inline constexpr double kAlphaDirt = 0.27;      // from the CCinc sideband (paper Sec. V A)

// ---------------------------------------------------------------------------------- reconstruction helpers
inline double reco_energy( double pe ) { return ( pe + 2.90 ) / 5.90; }

// |(centre - offset) - t| < sigma for any bunch (selection.in_bunch_window)
inline bool in_bunch_window( double t ) {
  for ( double c : kMCCenters ) {
    double shifted = c - kBunchOffset;
    if ( std::abs( shifted - t ) < kBunchSigma ) return true;
  }
  return false;
}

// q_b < 0.6, vertex found, reco FV, 5 <= E_reco <= 12 MeV. Written as !(reject) so that NaN passes,
// as in the Python code. Vertex in metres.
inline bool core_cuts( double qb, double rvx, double rvy, double rvz, double pe ) {
  double r2 = rvx * rvx + rvz * rvz;
  double E = reco_energy( pe );
  return !( qb >= 0.6 ) && !( rvx < -50. ) && !( r2 > std::pow( kFVRadius, 2 ) )
    && !( std::abs( rvy ) > kFVHalfHeight ) && !( E < 5. ) && !( E > 12. );
}

// Charge-coherence (C_V) direction: PE-weighted sum of unit vectors between pairs of early hits
// (reco_functions.pairwise_relative_direction). Returns the normalised vector (0,0,0 if undefined).
inline std::array<double, 3> pairwise_relative_direction( const std::vector<double>& hx,
  const std::vector<double>& hy, const std::vector<double>& hz, const std::vector<double>& hpe,
  const std::vector<double>& ht )
{
  std::array<double, 3> zero = { 0., 0., 0. };
  size_t n_all = ht.size();
  if ( n_all == 0 ) return zero;

  // Stable sort by hit time (Python list.sort is stable)
  std::vector<size_t> order( n_all );
  std::iota( order.begin(), order.end(), 0 );
  std::stable_sort( order.begin(), order.end(),
    [&ht]( size_t a, size_t b ) { return ht[a] < ht[b]; } );

  double mini = *std::min_element( ht.begin(), ht.end() );
  std::vector<size_t> kept;
  for ( size_t k : order ) {
    if ( ( ht[k] - mini ) > 13. ) continue;
    kept.push_back( k );
  }
  size_t n = kept.size();
  if ( n < 2 ) return zero;

  std::array<double, 3> vec = { 0., 0., 0. };
  for ( size_t i = 0; i < n; ++i ) {
    for ( size_t j = i + 1; j < n; ++j ) {
      size_t a = kept[i], b = kept[j];
      double dx = hx[b] - hx[a];
      double dy = hy[b] - hy[a];
      double dz = hz[b] - hz[a];
      // NumPy scalars: dx**2 is pow(dx, 2)
      double r = std::sqrt( std::pow( dx, 2 ) + std::pow( dy, 2 ) + std::pow( dz, 2 ) );
      if ( r == 0. ) continue;
      double weight = hpe[a] * hpe[b];
      vec[0] += weight * ( dx / r );
      vec[1] += weight * ( dy / r );
      vec[2] += weight * ( dz / r );
    }
  }

  // np.linalg.norm of a 3-vector: BLAS ddot accumulates with fused multiply-adds
  double norm = std::sqrt( std::fma( vec[2], vec[2], std::fma( vec[1], vec[1], vec[0] * vec[0] ) ) );
  if ( norm == 0. ) return zero;
  return { vec[0] / norm, vec[1] / norm, vec[2] / norm };
}

// Rotated-parabola acceptance on the C_V direction (reco_functions.is_inside_rotated_parabola).
// The Python call passes (Z_val, Y_val) = (vec[2], vec[0]).
inline bool is_inside_rotated_parabola( double Z_val, double Y_val, double a = -0.45,
  double c = 0.4, double b = 0.0, double theta_deg = 275. )
{
  double theta = theta_deg * ( M_PI / 180. );   // np.deg2rad
  double u = Z_val * std::cos( theta ) + Y_val * std::sin( theta );
  double v = -Z_val * std::sin( theta ) + Y_val * std::cos( theta );
  return v < a * std::pow( u, 2 ) + b * u + c;
}

inline bool direction_ok( const std::vector<double>& hx, const std::vector<double>& hy,
  const std::vector<double>& hz, const std::vector<double>& hpe, const std::vector<double>& ht )
{
  auto cw = pairwise_relative_direction( hx, hy, hz, hpe, ht );
  return is_inside_rotated_parabola( cw[2], cw[0] );
}

// ---------------------------------------------------------------------------------- truth
struct Truth {
  bool ext = false;     // true vertex outside the tank (External / dirt)
  bool fv = false;      // true vertex in the fiducial volume
  bool nu = false;      // NC QE on 16O, final-state nu_mu or nu_e
  bool nubar = false;   // NC QE on 16O, final-state anti-nu_mu or anti-nu_e
  int cat = -1;         // Python categories: 0 nuNCQE(FV) 1 NCother 2 CC 3 nubarNCQE 4 OOFV 5 External, -1 none
};

// truth_masks() and _categories() for one trigger (true vertex in cm, detector coordinates)
inline Truth truth_info( double vtx_x_cm, double vtx_y_cm, double vtx_z_cm, int nc, int qel,
  int target_z, int fsl_pdg )
{
  Truth t;
  double vx = vtx_x_cm / 100.;
  double vy = vtx_y_cm / 100. + kTrueYOffset;
  double vz = vtx_z_cm / 100. - kTrueZOffset;
  double r2t = vx * vx + vz * vz;
  t.ext = ( r2t > std::pow( kTankRadius, 2 ) ) || ( std::abs( vy ) > kTankHalfHeight );
  t.fv = ( r2t < std::pow( kFVRadius, 2 ) ) && ( std::abs( vy ) < kFVHalfHeight );
  bool ncqe_o16 = ( nc == 1 ) && ( qel == 1 ) && ( target_z == 8 );
  t.nu = ncqe_o16 && ( fsl_pdg == 14 || fsl_pdg == 12 );
  t.nubar = ncqe_o16 && ( fsl_pdg == -14 || fsl_pdg == -12 );

  // Assignment order of _categories(); later assignments overwrite earlier ones
  int cat = -1;
  if ( t.ext ) cat = 5;
  bool f_ = !t.ext && t.fv;
  if ( f_ && nc != 1 ) cat = 2;
  if ( f_ && nc == 1 && qel != 1 ) cat = 1;
  if ( f_ && nc == 1 && qel == 1 && target_z != 8 ) cat = 1;
  bool fo = f_ && nc == 1 && qel == 1 && target_z == 8;
  if ( fo && ( fsl_pdg == 14 || fsl_pdg == 12 ) ) cat = 0;
  if ( fo && ( fsl_pdg == -14 || fsl_pdg == -12 ) ) cat = 3;
  if ( !t.ext && !t.fv ) cat = 4;
  t.cat = cat;
  return t;
}

} // namespace ncqe
