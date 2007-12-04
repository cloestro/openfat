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


#include <stdio.h>
#include "matrix.h"
#include <malloc.h>
#include <math.h>


int of_crossland(Matrix *Sigma, double sigma_d, double tau_d, double Re,
        double steptime, double R, double *Cs_i, double *Cs_ni);

int of_sines(Matrix *Sigma, double sigma_d, double tau_d, double Re,
        double steptime, double R, double *Cs_i, double *Cs_ni);

int main(int argc, char *argv[])
{
    unsigned int i,j;
    Matrix cro = of_new_matrix(361, 7);
    double csi, csni;

    for ( i=0; i < cro.rows; i++)
    {
        cro.tab[i][0] = (double)i;
        cro.tab[i][1] = 50.0 + 200.0 * sin( cro.tab[i][0] * (2.0*M_PI) / 360.0);
        cro.tab[i][2] = 25.0 + 100.0 * cos( cro.tab[i][0] * (2.0*M_PI) / 360.0);
        cro.tab[i][6] = -50.0 ;
    }
    
    if (!of_crossland(&cro, 210.0, 180.0, 500.0, 1, -1, &csi, &csni))
        printf("intr: %f, not intr: %f\n", csi, csni);
    if (!of_sines(&cro, 210.0, 180.0, 500.0, 1, 0.1, &csi, &csni))
        printf("intr: %f, not intr: %f\n", csi, csni);
     
    of_free_matrix(cro);
    return 0;
}

