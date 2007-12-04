/*
 *  openfat -- A C multiaxial fatigue library
 *
 *  Copyright (C) 2007 Anthony Domi
 *  This file is part of openfat.
 *
 *  openfat is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  openfat is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with openfat; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/* $Id$ */


#ifndef UTILITIES_H
#define UTILITIES_H
#include "matrix.h"
#include <stdio.h>

extern int of_coordN(double p_crit, double tau_crit, double A, double B,
        double * pN, double * tN, double * ON);

extern double of_coefsecur(double p_crit, double tau_crit, double ON);
extern int of_cord(Matrix *tau_r, Matrix *tau_l, double *tau_r_m, double *tau_l_m);

#endif

