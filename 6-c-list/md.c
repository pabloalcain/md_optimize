#include "md.h"

float timedifference_msec(struct timeval t0, struct timeval t1)
{
    return (t1.tv_sec - t0.tv_sec) + (t1.tv_usec - t0.tv_usec) / 1000000.0f;
}

int main(int argc, char** argv) {
  FILE *file_time, *file_energy;
  struct timeval start, now;
  file_time = fopen("time.dat", "w");
  gettimeofday(&start, NULL);
  file_energy = fopen("energy.dat", "w");
  fprintf(file_energy, "#Potential, Kinetic, Total\n");
  System *sys = (System *) malloc(sizeof(System));
  sys->n_particles = 1000;
  sys->n_steps = 1000;
  sys->size = 13.0;
  sys->rcut = 2.5;
  sys->phicut = 4.0*(pow(2.5, -12) - pow (2.5, -6));
  init_system(sys);
  CellList *clist = (CellList *) malloc(sizeof(CellList));
  init_cells(clist, sys, 2.5);
  Integrator *integ = (Integrator *) malloc(sizeof(Integrator));
  integ->timestep = 0.0005;


  update_cells(clist, sys);
  newton(sys, clist);
  kinetic(sys);
  for (int i = 0; i < sys->n_steps; i++) {
    first_step(integ, sys);
    update_cells(clist, sys);
    newton(sys, clist);
    last_step(integ, sys);
    kinetic(sys);
    fprintf(file_energy, "%g, %g, %g\n", sys->potential, sys->kinetic, sys->potential+sys->kinetic);
    gettimeofday(&now, NULL);
    double elapsed = timedifference_msec(start, now);
    fprintf(file_time, "%i, %g\n", i, elapsed);
    printf("%i, %g\n", i, elapsed);
  }

  return 0;
}
