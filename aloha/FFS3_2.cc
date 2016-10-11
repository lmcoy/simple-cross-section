// This File is Automatically generated by ALOHA
// The process calculated in this file is:
// ProjP(2,1)
// 
#include "FFS3_2.h"

void FFS3_2(std::complex<double> F1[], std::complex<double> S3[], std::complex<double> COUP,
    double M2, double W2, std::complex<double> F2[])
{
  std::complex<double> cI = std::complex<double> (0., 1.); 
  double P2[4]; 
  std::complex<double> denom; 
  F2[0] = +F1[0] + S3[0]; 
  F2[1] = +F1[1] + S3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/(pow(P2[0], 2) - pow(P2[1], 2) - pow(P2[2], 2) - pow(P2[3], 2) -
      M2 * (M2 - cI * W2));
  F2[2] = denom * cI * S3[2] * (F1[4] * (P2[0] - P2[3]) + F1[5] * (+cI *
      (P2[2]) - P2[1]));
  F2[3] = denom * - cI * S3[2] * (F1[4] * (P2[1] + cI * (P2[2])) - F1[5] *
      (P2[0] + P2[3]));
  F2[4] = denom * cI * F1[4] * M2 * S3[2]; 
  F2[5] = denom * cI * F1[5] * M2 * S3[2]; 
}

