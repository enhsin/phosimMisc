#ifndef Par_H
#define Par_H

// general numerical accuracy
// int const SURFACE_POINTS = 8192;
// int const PERTURBATION_POINTS = 2048;
// int const SCREEN_SIZE = 1024;
// int const OPD_SCREEN_SIZE = 2048;

int const SURFACE_POINTS = 102400;
int const PERTURBATION_POINTS = 32;
int const SCREEN_SIZE = 1024;
int const OPD_SCREEN_SIZE = 255;

//for sersic
int const bin_number = 8000;
double const scale = 10.0;

//for raytrace
int const MAX_SOURCE = 150000;
int const MAX_SURF = 20;
int const MAX_LAYER = 100;
int const MAX_IMAGE = 100;
int const MAX_BOUNCE = 100;
int const NZERN = 21;
int const NCHEB = 21;
double const RAYTRACE_TOLERANCE = 1e-7;

// for trim
int const maxChip = 400;
int const maxCatalog = 20;

// for silicon
int const SILICON_STEPS = 1024;

#endif
