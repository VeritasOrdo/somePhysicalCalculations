#include "../package/LorentzVector/LorentzVector.h"
#include<iostream>

int main() {
  LorentzVector v1(1, 2, 3, 4);
  LorentzVector v2(5, 6, 7, 8);
  double v3 = v1 * v2;
  std::cout << v3 << std::endl;
  return 0;
}