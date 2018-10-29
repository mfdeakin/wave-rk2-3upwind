
#include "mesh.hpp"

#include <stdio.h>

template <int _x, int _y>
void Mesh<_x, _y>::print_cva(const ControlVolumes &cva) noexcept {
  printf("  y \\ x |");
  for(int i = 0; i < cva.extent(0); i++) {
    printf("   % .3f,      ", median_x(i));
  }
  printf("\n--------+");
  for(int i = 0; i < cva.extent(0); i++) {
    printf("----------------");
  }
  printf("\n");
  for(int j = cva.extent(1) - 1; j >= 0; j--) {
    printf("% .3f  |", median_y(j));
    for(int i = 0; i < cva.extent(0); i++) {
      printf("   % .5e,", cva(i, j));
    }
    printf("\n");
  }
}

template <int _x, int _y>
void Mesh<_x, _y>::print() noexcept {
  print_cva(cva);
}

// Instantiate all the versions we need
template class Mesh<10, 10>;
template class Mesh<20, 20>;
template class Mesh<30, 30>;
template class Mesh<40, 40>;
template class Mesh<50, 50>;
template class Mesh<60, 60>;
template class Mesh<70, 70>;
template class Mesh<80, 80>;
template class Mesh<90, 90>;
template class Mesh<100, 100>;
template class Mesh<110, 110>;
template class Mesh<120, 120>;
template class Mesh<130, 130>;
template class Mesh<140, 140>;
template class Mesh<150, 150>;
template class Mesh<160, 160>;
template class Mesh<170, 170>;
template class Mesh<180, 180>;
template class Mesh<190, 190>;
template class Mesh<200, 200>;
template class Mesh<210, 210>;
template class Mesh<220, 220>;
template class Mesh<230, 230>;
template class Mesh<240, 240>;
template class Mesh<250, 250>;
template class Mesh<260, 260>;
template class Mesh<270, 270>;
template class Mesh<280, 280>;
template class Mesh<290, 290>;
template class Mesh<300, 300>;
template class Mesh<310, 310>;
template class Mesh<320, 320>;

