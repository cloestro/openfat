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
#include <math.h>
#include "utilities.h"


//TODO: return the coordinates Sh_max, ... and
//compute intersecting point
//
int of_crossland(Matrix *Sigma, double sigma_d, double tau_d, double Re,
        double steptime, double R, double *Cs_i, double *Cs_ni)
{
    Matrix s_a_1, s_a_2, temp1;
    double sm11, sm22, sm33, sm12, sm23, sm31, T, Sh_max,
           tau_oct_a_max, A, B, Sh_N_intr, tau_N_intr, tau_oct_a_intr;
    unsigned int i, rows, cols, j;
    double tau_oct_a_el, sa11_el, sa22_el, sa33_el, sa12_el, sa23_el, sa31_el;
    double DJ2, DJ2max ;
    double ON_intr ;
    double Sh_N_non_intr, tau_N_non_intr, ON_non_intr;
    double Sha ;


    if ( Sigma->cols != 7)
    {
        fprintf(stderr, "Input table must have 7 columns (%i given)\n", 
                Sigma->cols);
        return -1;
    }
    //TODO: misc. checks on input

    rows = Sigma->rows;
    cols = Sigma->cols;
    
    sm11 = 0.0;
    sm22 = 0.0;
    sm33 = 0.0;
    sm23 = 0.0;
    sm31 = 0.0;
    sm12 = 0.0;
    T = Sigma->tab[rows-1][0];
    
    /* Trapeze method to obtain the mean Sigma */

    for (i=0; i < rows - 1; i++)
    {
        sm11 += ( Sigma->tab[i][1] + Sigma->tab[i+1][1] ) * 
            ( Sigma->tab[i+1][0] - Sigma->tab[i][0] ) / ( 2.0 * T ) ;

        sm22 += ( Sigma->tab[i][2] + Sigma->tab[i+1][2] ) * 
            ( Sigma->tab[i+1][0] - Sigma->tab[i][0] ) / ( 2.0 * T ) ;

        sm33 += ( Sigma->tab[i][3] + Sigma->tab[i+1][3] ) * 
            ( Sigma->tab[i+1][0] - Sigma->tab[i][0] ) / ( 2.0 * T ) ;

        sm23 += ( Sigma->tab[i][4] + Sigma->tab[i+1][4] ) * 
            ( Sigma->tab[i+1][0] - Sigma->tab[i][0] ) / ( 2.0 * T ) ;

        sm31 += ( Sigma->tab[i][5] + Sigma->tab[i+1][5] ) * 
            ( Sigma->tab[i+1][0] - Sigma->tab[i][0] ) / ( 2.0 * T ) ;
        
        sm12 += ( Sigma->tab[i][6] + Sigma->tab[i+1][6] ) * 
            ( Sigma->tab[i+1][0] - Sigma->tab[i][0] ) / ( 2.0 * T ) ;

 
    }


    s_a_1 = of_new_matrix(3, 3);
    s_a_2 = of_new_matrix(3, 3);
    temp1 = of_new_matrix(3, 3);
    DJ2max = 0.0 ;

    for ( i=0; i < rows; i++)
    {
        if ( i== 0)
        {
            Sh_max = (Sigma->tab[i][1] + Sigma->tab[i][2] + Sigma->tab[i][3] ) / 3.0;
            sa11_el = Sigma->tab[i][1] - sm11 ;
            sa22_el = Sigma->tab[i][2] - sm22 ;
            sa33_el = Sigma->tab[i][3] - sm33 ;
            sa23_el = Sigma->tab[i][4] - sm23 ;
            sa31_el = Sigma->tab[i][5] - sm31 ;
            sa12_el = Sigma->tab[i][6] - sm12 ;
            
            tau_oct_a_max = sqrt(pow((sa11_el - sa22_el),2) + 
                    pow(( sa22_el - sa33_el ), 2) + 
                    pow(( sa33_el - sa11_el), 2) +
                    6.0 * ( pow(sa12_el, 2) + pow(sa23_el, 2) + pow(sa31_el, 2) ) ) / 3.0;
        }
        else
        {
            if ( ( Sigma->tab[i][1] + Sigma->tab[i][2] + Sigma->tab[i][3] ) / 3.0 > Sh_max )
                Sh_max = ( Sigma->tab[i][1] + Sigma->tab[i][2] + Sigma->tab[i][3] ) / 3.0 ;

            sa11_el = Sigma->tab[i][1] - sm11 ;
            sa22_el = Sigma->tab[i][2] - sm22 ;
            sa33_el = Sigma->tab[i][3] - sm33 ;
            sa23_el = Sigma->tab[i][4] - sm23 ;
            sa31_el = Sigma->tab[i][5] - sm31 ;
            sa12_el = Sigma->tab[i][6] - sm12 ;
 
            tau_oct_a_el = sqrt(pow((sa11_el - sa22_el),2) + 
                    pow(( sa22_el - sa33_el ), 2) + 
                    pow(( sa33_el - sa11_el), 2) +
                    6.0 * ( pow(sa12_el, 2) + pow(sa23_el, 2) + pow(sa31_el, 2) ) ) / 3.0;
            if ( tau_oct_a_el > tau_oct_a_max )
                tau_oct_a_max = tau_oct_a_el ;

        }

        Sha =  ( sa11_el + sa22_el + sa33_el ) / 3.0 ;
        //( Sigma->tab[i][1] + Sigma->tab[i][2] + Sigma->tab[i][3] ) / 3.0 ;

        s_a_1.tab[0][0] = sa11_el - Sha ;
        s_a_1.tab[1][1] = sa22_el - Sha ; 
        s_a_1.tab[2][2] = sa33_el - Sha ; 
        s_a_1.tab[0][1] = sa12_el ;
        s_a_1.tab[1][0] = sa12_el ;
        s_a_1.tab[1][2] = sa23_el ;
        s_a_1.tab[2][1] = sa23_el ;
        s_a_1.tab[2][0] = sa31_el ;
        s_a_1.tab[0][2] = sa31_el ;


        for ( j=0; j < rows; j++)
        {
            sa11_el = Sigma->tab[j][1] - sm11 ;
            sa22_el = Sigma->tab[j][2] - sm22 ;
            sa33_el = Sigma->tab[j][3] - sm33 ;
            sa23_el = Sigma->tab[j][4] - sm23 ;
            sa31_el = Sigma->tab[j][5] - sm31 ;
            sa12_el = Sigma->tab[j][6] - sm12 ;
            
            Sha =  ( sa11_el + sa22_el + sa33_el ) / 3.0 ; 

            s_a_2.tab[0][0] = - sa11_el + Sha; 
            s_a_2.tab[1][1] = - sa22_el + Sha;
            s_a_2.tab[2][2] = - sa33_el + Sha; 
            s_a_2.tab[0][1] = - sa12_el ;
            s_a_2.tab[1][0] = - sa12_el ;
            s_a_2.tab[1][2] = - sa23_el ;
            s_a_2.tab[2][1] = - sa23_el ;
            s_a_2.tab[2][0] = - sa31_el ;
            s_a_2.tab[0][2] = - sa31_el ;

            of_sum_mat(&s_a_2, &s_a_1, &temp1);
            DJ2 = of_tracesquare( &temp1 ) / 2.0  ;

            if ( DJ2 > DJ2max )
                DJ2max = DJ2 ;
        }

    }

    tau_oct_a_intr = sqrt( 2.0 * DJ2max / 3.0) / 2.0 ;


    A = sqrt(3.0 / 2.0 ) * ( 1.0 - R ) * tau_d / sigma_d
        - sqrt(2.0) * ( 1.0 - R ) / 2.0 ;
    B = sqrt(2.0 / 3.0 ) * tau_d ;

    of_free_matrix(temp1);
    of_free_matrix(s_a_1);
    of_free_matrix(s_a_2);

    if ( A <= 0 )
    {
        fprintf(stderr, "Warning: the A coefficient is negative (%f)\n", A);
    }

   if (of_coordN(Sh_max, tau_oct_a_intr, A, B, &Sh_N_intr , 
               &tau_N_intr, &ON_intr) )
   {
       return -3 ;
   }
   if (of_coordN(Sh_max, tau_oct_a_max, A, B, &Sh_N_non_intr , 
               &tau_N_non_intr, &ON_non_intr) )
   {
       return -4 ;
   }

#ifdef DEBUG
   printf("A=%f, B=%f\n", A, B);
   printf("Shmax=%f, Tintr=%f, Tnintr=%f\n", Sh_max, tau_oct_a_intr,
           tau_oct_a_max);
   printf("ON_i=%f, ON_ni=%f\n", ON_intr, ON_non_intr);
#endif
    *Cs_i = of_coefsecur(Sh_max, tau_oct_a_intr, ON_intr);
    *Cs_ni = of_coefsecur(Sh_max, tau_oct_a_max, ON_non_intr);

    return 0;

}

