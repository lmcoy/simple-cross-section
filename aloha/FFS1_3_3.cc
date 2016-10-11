// This File is Automatically generated by ALOHA
// 
#include "FFS1_3_3.h"
#include "FFS1_3.h"
#include "FFS3_3.h"

void FFS1_3_3(std::complex<double> F1[], std::complex<double> F2[], std::complex<double>
    COUP1, std::complex<double> COUP2, double M3, double W3, std::complex<double> S3[])
{
  std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> denom; 
  std::complex<double> Stmp[3]; 
  double P3[4]; 
  int i; 
  FFS1_3(F1, F2, COUP1, M3, W3, S3); 
  FFS3_3(F1, F2, COUP2, M3, W3, Stmp); 
  i = 2; 
  while (i < 3)
  {
    S3[i] = S3[i] + Stmp[i]; 
    i++; 
  }
}
