// This File is Automatically generated by ALOHA
// The process calculated in this file is:
// ProjP(2,1)
// 
#include "FFS3_1.h"

void FFS3_1(std::complex<double> F2[], std::complex<double> S3[], std::complex<double> COUP,
    double M1, double W1, std::complex<double> F1[])
{
  std::complex<double> cI = std::complex<double> (0., 1.); 
  double P1[4]; 
  std::complex<double> denom; 
  F1[0] = +F2[0] + S3[0]; 
  F1[1] = +F2[1] + S3[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  F1[2] = denom * - cI * S3[2] * (F2[4] * (P1[0] + P1[3]) + F2[5] * (P1[1] + cI
      * (P1[2])));
  F1[3] = denom * cI * S3[2] * (F2[4] * (+cI * (P1[2]) - P1[1]) + F2[5] *
      (P1[3] - P1[0]));
  F1[4] = denom * cI * F2[4] * M1 * S3[2]; 
  F1[5] = denom * cI * F2[5] * M1 * S3[2]; 
}

