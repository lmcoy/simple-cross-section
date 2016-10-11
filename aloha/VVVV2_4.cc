// This File is Automatically generated by ALOHA
// The process calculated in this file is:
// Metric(1,4)*Metric(2,3) + Metric(1,3)*Metric(2,4) - 2*Metric(1,2)*Metric(3,4)
// 
#include "VVVV2_4.h"

void VVVV2_4(std::complex<double> V1[], std::complex<double> V2[], std::complex<double> V3[],
    std::complex<double> COUP, double M4, double W4, std::complex<double> V4[])
{
  std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> TMP17; 
  double OM4; 
  std::complex<double> denom; 
  std::complex<double> TMP26; 
  double P4[4]; 
  std::complex<double> TMP29; 
  std::complex<double> TMP25; 
  std::complex<double> TMP19; 
  std::complex<double> TMP13; 
  OM4 = 0.; 
  if (M4 != 0.)
    OM4 = 1./pow(M4, 2); 
  V4[0] = +V1[0] + V2[0] + V3[0]; 
  V4[1] = +V1[1] + V2[1] + V3[1]; 
  P4[0] = -V4[0].real(); 
  P4[1] = -V4[1].real(); 
  P4[2] = -V4[1].imag(); 
  P4[3] = -V4[0].imag(); 
  TMP25 = (V1[2] * P4[0] - V1[3] * P4[1] - V1[4] * P4[2] - V1[5] * P4[3]); 
  TMP26 = (V2[2] * P4[0] - V2[3] * P4[1] - V2[4] * P4[2] - V2[5] * P4[3]); 
  TMP29 = (V3[2] * P4[0] - V3[3] * P4[1] - V3[4] * P4[2] - V3[5] * P4[3]); 
  TMP19 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP17 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  TMP13 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  denom = COUP/(pow(P4[0], 2) - pow(P4[1], 2) - pow(P4[2], 2) - pow(P4[3], 2) -
      M4 * (M4 - cI * W4));
  V4[2] = denom * (OM4 * P4[0] * (-2. * cI * (TMP13 * TMP29) + cI * (TMP19 *
      TMP25 + TMP17 * TMP26)) + (-cI * (V1[2] * TMP19 + V2[2] * TMP17) + 2. *
      cI * (V3[2] * TMP13)));
  V4[3] = denom * (OM4 * P4[1] * (-2. * cI * (TMP13 * TMP29) + cI * (TMP19 *
      TMP25 + TMP17 * TMP26)) + (-cI * (V1[3] * TMP19 + V2[3] * TMP17) + 2. *
      cI * (V3[3] * TMP13)));
  V4[4] = denom * (OM4 * P4[2] * (-2. * cI * (TMP13 * TMP29) + cI * (TMP19 *
      TMP25 + TMP17 * TMP26)) + (-cI * (V1[4] * TMP19 + V2[4] * TMP17) + 2. *
      cI * (V3[4] * TMP13)));
  V4[5] = denom * (OM4 * P4[3] * (-2. * cI * (TMP13 * TMP29) + cI * (TMP19 *
      TMP25 + TMP17 * TMP26)) + (-cI * (V1[5] * TMP19 + V2[5] * TMP17) + 2. *
      cI * (V3[5] * TMP13)));
}

