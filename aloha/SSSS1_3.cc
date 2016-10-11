// This File is Automatically generated by ALOHA
// The process calculated in this file is:
// 1
// 
#include "SSSS1_3.h"

void SSSS1_3(std::complex<double> S1[], std::complex<double> S2[], std::complex<double> S4[],
    std::complex<double> COUP, double M3, double W3, std::complex<double> S3[])
{
  std::complex<double> cI = std::complex<double> (0., 1.); 
  double P3[4]; 
  std::complex<double> denom; 
  S3[0] = +S1[0] + S2[0] + S4[0]; 
  S3[1] = +S1[1] + S2[1] + S4[1]; 
  P3[0] = -S3[0].real(); 
  P3[1] = -S3[1].real(); 
  P3[2] = -S3[1].imag(); 
  P3[3] = -S3[0].imag(); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  S3[2] = denom * cI * S4[2] * S2[2] * S1[2]; 
}

