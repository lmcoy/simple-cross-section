// This File is Automatically generated by ALOHA
// The process calculated in this file is:
// P(3,1)*Metric(1,2) - P(3,2)*Metric(1,2) - P(2,1)*Metric(1,3) +
// P(2,3)*Metric(1,3) + P(1,2)*Metric(2,3) - P(1,3)*Metric(2,3)
// 
#include "VVV1_0.h"

void VVV1_0(std::complex<double> V1[], std::complex<double> V2[], std::complex<double> V3[],
    std::complex<double> COUP, std::complex<double> & vertex)
{
  std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  std::complex<double> TMP10; 
  double P2[4]; 
  std::complex<double> TMP19; 
  std::complex<double> TMP17; 
  double P3[4]; 
  std::complex<double> TMP16; 
  std::complex<double> TMP15; 
  std::complex<double> TMP14; 
  std::complex<double> TMP9; 
  std::complex<double> TMP13; 
  std::complex<double> TMP18; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  TMP9 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  TMP19 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP18 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP15 = (P1[0] * V3[2] - P1[1] * V3[3] - P1[2] * V3[4] - P1[3] * V3[5]); 
  TMP14 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP17 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  TMP16 = (P2[0] * V3[2] - P2[1] * V3[3] - P2[2] * V3[4] - P2[3] * V3[5]); 
  TMP10 = (P3[0] * V1[2] - P3[1] * V1[3] - P3[2] * V1[4] - P3[3] * V1[5]); 
  TMP13 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  vertex = COUP * (TMP13 * (-cI * (TMP15) + cI * (TMP16)) + (TMP17 * (-cI *
      (TMP18) + cI * (TMP14)) + TMP19 * (-cI * (TMP9) + cI * (TMP10))));
}

