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
#include <math.h>
#include "utilities.h"

extern int of_coordN(double p_crit, double tau_crit, double A, double B,
        double * pN, double * tN, double * ON)
{
    *pN = 0.0 ;
    *tN = 0.0 ;
    *ON = 0.0 ;

    if ( p_crit == 0.0 )
    {
        *tN = B;
        *ON = B;
    }
    else
    {
        if ( tau_crit <= - A * p_crit )
        {
            fprintf(stderr, "Error: safety coefficient not defined\n");
            return -1 ;
        }
        else
        {
            *pN = B / ( A + ( tau_crit / p_crit ) ) ;
            *tN = ( tau_crit / p_crit ) * ( *pN ) ;
            *ON = sqrt( pow(*pN, 2) + pow(*tN, 2 ) );
        }

    }
    return 0;
}
extern double of_coefsecur(double p_crit, double tau_crit, double ON)
{
    double OM = 0.0 ;
    double Cs = 0.0 ;

    OM = sqrt( pow(p_crit, 2) + pow(tau_crit, 2) );
    Cs = ON / OM;
    return Cs;


}






extern int of_cord(Matrix *tau_r, Matrix *tau_l, double *tau_r_moy, double *tau_l_moy)
{
    unsigned int i = 0, j = 0;
    unsigned int ic = 0, jc = 0;
    double norme = 0.0 ;
    double normemax = 0.0 ;

    if ( tau_r->rows == 0)
    {
        fprintf(stderr, "Error: number of rows of tau_r is zero\n");
        return -1 ;
    }

    if ( tau_l->rows == 0)
    {
        fprintf(stderr, "Error: number of rows of tau_l is zero\n");
        return -2 ;
    }

    if ( tau_r->cols != 1 )
    {
        fprintf(stderr, "Error: number of columns of tau_r is not one\n");
        return -3 ;
    }

    if ( tau_l->cols != 1)
    {
        fprintf(stderr, "Error: number of columns of tau_l is not one\n");
        return -4;
    }

    if ( tau_r->rows != tau_l->rows )
    {
        fprintf(stderr, "Error tau_r and tau_l have different sizes\n");
        return -5;
    }

    for ( i=0; i < tau_r->rows; i++)
    {
        for ( j = 0; j < tau_r->rows; j++)
        {
            norme = sqrt( pow( tau_r->tab[i][0] - tau_r->tab[j][0],2) + 
                    pow( tau_l->tab[i][0] - tau_l->tab[j][0] ,2) );
            
            if ( norme > normemax )
            {
                normemax = norme ;
                ic = i ;
                jc = j ;
            }
        }
    }

    *tau_r_moy = ( tau_r->tab[ic][0] + tau_r->tab[jc][0] ) / 2.0 ;
    *tau_l_moy = ( tau_l->tab[ic][0] + tau_l->tab[jc][0] ) / 2.0 ;
    return 0;
}

