#include "force.h"

void newton(System *sys) {
  double phicutoff = 4.0/pow(2.5, 12) - 4.0/pow(2.5, 6);
  sys->potential = 0.0;
  for (int i = 0; i < 3 * sys->n_particles; i++)
    sys->force[i] = 0.0;
  for (int i = 0; i < sys->n_particles; i++) {
    for (int j = i + 1; j < sys->n_particles; j++) {
      double dr[3];
      for (int k = 0; k < 3; k++)
        dr[k] = sys->position[3*i+k] - sys->position[3*j+k];
      minimum_images(sys, dr);
      calculate_force(sys, i, j, dr);
    }
  }
}

__attribute__((always_inline,pure))
inline void minimum_images(System *sys, double *dr) {
  for (int k = 0; k < 3; k++) {
    if (dr[k] > 0.5 * sys->size)
      dr[k] -= sys->size;
    else if (dr[k] < -0.5 * sys->size)
      dr[k] += sys->size;
  }
}

__attribute__((always_inline,pure))
inline void calculate_force(System *sys, int i, int j, double *dr) {
  double distance = 0.0;
  for (int k = 0; k < 3; k++) {
    distance += dr[k] * dr[k];
  }
  if (distance < 2.5 * 2.5) {
    double rm2 = 1.0/distance;
    double rm6 = rm2 * rm2 * rm2;
    double rm12 = rm6 * rm6;
    double phi  = 4.0 * (rm12 - rm6);
    double dphi = 24.0*rm2*(2.0*rm12 - rm6);
    for (int k = 0; k < 3; k++) {
      sys->force[3*i+k] += dphi * dr[k];
      sys->force[3*j+k] -= dphi * dr[k];
    }
    sys->potential += phi - sys->phicut;
  }
}

void kinetic(System *sys) {
  double kin = 0.0;
  for (int i = 0; i < 3 * sys->n_particles; i++) {
    kin += sys->velocity[i] * sys->velocity[i];
  }
  sys->kinetic = kin/2;
}
