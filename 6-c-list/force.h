#ifndef FORCE_H
#define FORCE_H
#include "system.h"
#include "cell.h"

void newton(System *sys, CellList *clist);
void minimum_images(System *sys, double *dr);
void calculate_force(System *sys, int i, int j, double *dr);
void kinetic(System *sys);
#endif
