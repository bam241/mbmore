#include <gtest/gtest.h>

#include "StageConfig.h"

#include "agent_tests.h"
#include "context.h"
#include "facility_tests.h"

namespace mbmore {

// Benchmarked using a regression test (expected values calculated manually)
namespace stageconfig_test {
// Fixed for a cascade separating out U235 from U238 in UF6 gas
double M = 0.35202;   // kg/mol UF6
double dM = 0.003;  // kg/mol U238 - U235
double x = 1000;    // Pressure ratio (Glaser)

// General cascade assumptions
double flow_ratio = 2.0;
double eff = 1.0;
double cut = 0.5;

// Centrifgue parameters based on Glaser SGS 2009 paper
double v_a = 485;                                           // m/s
double height = 0.5;                                        // meters
double diameter = 0.15;                                     // meters
double feed_m = 15 * 60 * 60 / ((1e3) * 60 * 60 * 1000.0);  // kg/sec
double temp = 320.0;                                        // Kelvin

// Cascade params used in in calculating expected values
const double feed_assay = 0.0071;
const double prod_assay = 0.035;
const double waste_assay = 0.001;
const double feed_c = 739 / (30.4 * 24 * 60 * 60);    // kg/month -> kg/sec
const double product_c = 77 / (30.4 * 24 * 60 * 60);  // kg/month -> kg/sec
CentrifugeConfig centrifuge(v_a, height, diameter, feed_m, temp, eff, M, dM, x,
                            flow_ratio);

// del U=7.0323281e-08 alpha=1.16321
double delU = centrifuge.ComputeDeltaU(cut);

const double tol_assay = 1e-5;
const double tol_qty = 1e-6;
const double tol_num = 1e-2;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Find product assay from separation factor alpha
TEST(StageConfig_Test, TestAssays) {
  double cur_alpha = 1.4;
  double cur_f_assay = 0.007;

  StageConfig stage(cur_f_assay, feed_m, cut, delU, cur_alpha, 1e-16);
  stage.ProductAssay();

  // N_prime = alpha*R / ( 1+alpha*R)
  double target_prod_assay = 0.009773;
  double tol = 1e-6;

  EXPECT_NEAR(stage.product_assay(), target_prod_assay, tol);

  double n_stages = 5;
  double target_w_assay = 0.004227;
  stage.BetaByAlphaAndCut();
  stage.TailAssay();

  EXPECT_NEAR(stage.tail_assay(), target_w_assay, tol);
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Calculate ideal SWU params of single machinefeed_assay
// (separation potential delU and separation factor alpha)
TEST(StageConfig_Test, TestSWU) {
  double expected_U = 7.03232816847e-08;
  double tol = 1e-9;
  int k = 2;
  double P1_F[5] = {0, 0, 12.6e-6, 6.4e-6, 3.9e-6};
  CentrifugeConfig P1(320, 1.8, 0.1, P1_F[k], 320, 0.564, M, dM, x, k);
  StageConfig stg_P1;
  stg_P1.centrifuge = P1;
  //std::cout << "P1, du: " << P1.ComputeDeltaU(0.50)*3600*24*365.25/M*0.238 << " "; 
  //stg_P1.DU(P1.ComputeDeltaU(0.50)/M*0.238);
  std::cout << "P1, du: " << P1.ComputeDeltaU(0.50)*3600*24*365.25 << " "; 
  stg_P1.DU(P1.ComputeDeltaU(0.50));
  //stg_P1.DU(2.5/3600/24/365);
  stg_P1.cut(.5);
  stg_P1.feed_assay(0.007);
  stg_P1.AlphaByDU();
  stg_P1.BetaByAlphaAndCut();
  std::cout << stg_P1.feed_assay() << " " << stg_P1.alpha()*stg_P1.beta()<< std::endl;
  stg_P1.feed_assay(0.07);
  stg_P1.AlphaByDU();
  stg_P1.BetaByAlphaAndCut();
  std::cout << stg_P1.feed_assay() << " " << stg_P1.alpha()*stg_P1.beta()<< std::endl;

  double P2_F[5] = {0, 0, 15e-6, 7.7e-6, 4.6e-6};
  CentrifugeConfig P2(485, 1, 0.15, P2_F[k], 320, 0.465, M, dM, x, k);
  StageConfig stg_P2;
  stg_P2.centrifuge = P2;
  //std::cout << "P2, du: " << P2.ComputeDeltaU(0.50)*3600*24*365.25/M*0.238 << " "; 
  //stg_P2.DU(P2.ComputeDeltaU(0.50)/M*0.238);
  std::cout << "P2, du: " << P2.ComputeDeltaU(0.50)*3600*24*365.25 << " "; 
  stg_P2.DU(P2.ComputeDeltaU(0.50));
  //stg_P2.DU(6./3600/24/365);
  stg_P2.cut(.5);
  stg_P2.feed_assay(0.007);
  stg_P2.AlphaByDU();
  stg_P2.BetaByAlphaAndCut();
  std::cout << stg_P2.feed_assay() << " " << stg_P2.alpha()*stg_P2.beta()<< std::endl;
  stg_P2.feed_assay(0.07);
  stg_P2.AlphaByDU();
  stg_P2.BetaByAlphaAndCut();
  std::cout << stg_P2.feed_assay() << " " << stg_P2.alpha()*stg_P2.beta()<< std::endl;

  double A1_F[5] = {0, 0, 51.4e-6, 26.2e-6, 15.9e-6};
  CentrifugeConfig A1(600, 2, 0.2, A1_F[k], 320, 0.340, M, dM, x, k);
  StageConfig stg_A1;
  stg_A1.centrifuge = A1;
  //std::cout << "A1, du: " << A1.ComputeDeltaU(0.50)*3600*24*365.25/M*0.238 << " "; 
  //stg_A1.DU(A1.ComputeDeltaU(0.50)/M*0.238);
  std::cout << "A1, du: " << A1.ComputeDeltaU(0.50)*3600*24*365.25 << " "; 
  stg_A1.DU(A1.ComputeDeltaU(0.50));
  //stg_A1.DU(20.6/3600/24/365);
  stg_A1.cut(.5);
  stg_A1.feed_assay(0.007);
  stg_A1.AlphaByDU();
  stg_A1.BetaByAlphaAndCut();
  std::cout << stg_A1.feed_assay() << " " << stg_A1.alpha()*stg_A1.beta()<< std::endl;
  stg_A1.feed_assay(0.07);
  stg_A1.AlphaByDU();
  stg_A1.BetaByAlphaAndCut();
  std::cout << stg_A1.feed_assay() << " " << stg_A1.alpha()*stg_A1.beta()<< std::endl;

  double A2_F[5] = {0, 0, 214e-6, 109e-6, 66e-6};
  CentrifugeConfig A2(750, 5, 0.2, A2_F[k], 320, 0.263, M, dM, x, k);
  StageConfig stg_A2;
  stg_A2.centrifuge = A2;
  //std::cout << "A2, du: " << A2.ComputeDeltaU(0.50)*3600*24*365.25/M*0.238 << " "; 
  //stg_A2.DU(A2.ComputeDeltaU(0.50)/M*0.238);
  std::cout << "A2, du: " << A2.ComputeDeltaU(0.50)*3600*24*365.25 << " "; 
  stg_A2.DU(A2.ComputeDeltaU(0.50));
  //stg_A2.DU(97./3600/24/365);
  stg_A2.cut(.5);
  stg_A2.feed_assay(0.007);
  stg_A2.AlphaByDU();
  stg_A2.BetaByAlphaAndCut();
  std::cout << stg_A2.feed_assay() << " "<< stg_A2.alpha()*stg_A2.beta()<< std::endl;
  stg_A2.feed_assay(0.07);
  stg_A2.AlphaByDU();
  stg_A2.BetaByAlphaAndCut();
  std::cout << stg_A2.feed_assay() << " "<< stg_A2.alpha()*stg_A2.beta()<< std::endl;

  double A3_F[5] = {0, 0, 429e-6, 219e-6, 132e-6};
  CentrifugeConfig A3(750, 10, 0.6, A3_F[k], 320, 0.263, M, dM, x, k);
  StageConfig stg_A3;
  stg_A3.centrifuge = A3;
  //std::cout << "A3, du: " << A3.ComputeDeltaU(0.50)*3600*24*365.25/M*0.238 << " "; 
  //stg_A3.DU(A3.ComputeDeltaU(0.50)/M*0.238);
  std::cout << "A3, du: " << A3.ComputeDeltaU(0.50)*3600*24*365.25 << " "; 
  stg_A3.DU(A3.ComputeDeltaU(0.50));
  //stg_A3.DU(195./3600/24/365);
  stg_A3.cut(.5);
  stg_A3.feed_assay(0.007);
  stg_A3.AlphaByDU();
  stg_A3.BetaByAlphaAndCut();
  std::cout << stg_A3.feed_assay() << " " << stg_A3.alpha()*stg_A3.beta()<< std::endl;
  stg_A3.feed_assay(0.07);
  stg_A3.AlphaByDU();
  stg_A3.BetaByAlphaAndCut();
  std::cout << stg_A3.feed_assay() << " " << stg_A3.alpha()*stg_A3.beta()<< std::endl;
  
  
  StageConfig stage(feed_assay, feed_m, cut, delU, -1, 1e-16);

  double expected_alpha = 1.16321;
  double tol_alpha = 1e-2;
  EXPECT_NEAR(stage.alpha(), expected_alpha, tol_alpha);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Calculate the product assay for an ideal stage configuration.
TEST(StageConfig_Test, TestIdeal) {
  StageConfig stage_ideal(feed_assay, feed_m, cut, -1, -1, 1e-16);

  // Only setting precision for building ideal stage
  stage_ideal.precision(1e-3);
  stage_ideal.BuildIdealStg();
  stage_ideal.precision(1e-16);

  // All expected numbers were calculated using the methods used
  // and are trusted to be correct (regression test).
  double expected_alpha = 1.18181;
  double tol_alpha = 1e-2;

  double expected_cut = 0.4589269;
  double tol_cut = 1e-3;

  double expected_U = 7.4221362040947e-08;
  double tol_DU = 1e-9;

  EXPECT_NEAR(stage_ideal.alpha(), expected_alpha, tol_alpha);
  EXPECT_NEAR(stage_ideal.cut(), expected_cut, tol_cut);
  EXPECT_NEAR(stage_ideal.DU(), expected_U, tol_DU);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST(StageConfig_Test, ProductAssayByGamma) {
  double gamma = 1.3798316056650026;
  double target_product_assay = 0.00822;
  double theta_ = 0.46040372309;
  double feed_assay_ = 0.007;
  StageConfig stage(feed_assay_, feed_c, theta_, delU, -1, 1e-16);
  stage.ProductAssayByGamma(gamma);

  EXPECT_NEAR(target_product_assay,stage.product_assay(), 1e-5);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST(StageConfig_Test, AlphaByProductAssay) {
  double gamma = 1.3798316056650026;
  double target_product_assay = 0.00821;
  double theta_ = 0.46040372309;
  double feed_assay_ = 0.007;
  StageConfig stage(feed_assay_, feed_c, theta_, delU, -1, 1e-16);
  stage.feed_assay(0.1);
  stage.product_assay(0.3);
  double alpha_ = 0.3/(1-0.3)*(1-0.1)/0.1;
  stage.AlphaByProductAssay();

  EXPECT_NEAR(alpha_,stage.alpha(), 1e-5);
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Determine the output of the first enrich/strip stage of a cascade
// based on the design params for the cascade
TEST(StageConfig_Test, TestStages) {
  StageConfig stage(feed_assay, feed_c, cut, delU, -1, 1e-16);

  stage.ProductAssay();
  stage.MachinesNeededPerStage();

  //double expected_product_assay_s = 0.007782156959583;
  double expected_product_assay_s = 0.0082492;

  // Calculated using equations from 2009 Glaser paper
  int expected_n_mach_e = 19;

  EXPECT_NEAR(stage.n_machines(), expected_n_mach_e, tol_num);
  EXPECT_NEAR(stage.product_assay(), expected_product_assay_s, tol_assay);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Verify that, for a P1-type centrifuge, results are in agreement
// with those presented by Glaser paper.
TEST(StageConfig_Test, AlphaBeta) {
  // P1-type parameters
  double v_a = 320; // m/s
  double Z = 1.8; // m
  double d = 0.10; // m
  double eff = 0.564;
  double x = 1000.;
  double flow_ratio = 2.0;
  double temp = 320; // K
  double M = 0.352;   // kg/mol UF6
  double dM = 0.003;  // kg/mol U238 - U235
  double feed = 12.6 / 1000. / 1000.; // mg/s -> kg/s
  double cut = 0.5;
  double f_assay = 0.00720;

  CentrifugeConfig cent(v_a, Z, d, feed, temp, eff, M, dM, x, flow_ratio);

  double delU = cent.ComputeDeltaU(cut);

  double delU_e = 2.5 / (365.25 * 24 * 60 * 60);
  EXPECT_NEAR(delU, 1.0, 1e-8);

  StageConfig stg(f_assay, feed, cut, delU);

  double ab_e = 1.29;
  EXPECT_NEAR(stg.alpha() * stg.beta(), ab_e, 0.0001);

}

}  // namespace enrichfunctiontests
}  // namespace mbmore
