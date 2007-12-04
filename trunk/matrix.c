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

extern Matrix of_new_matrix(unsigned int rows, unsigned int cols)
{
    unsigned int i;
    Matrix new_mat ;
    new_mat.rows = rows ;
    new_mat.cols = cols ;
    new_mat.tab = (double **) malloc(rows * sizeof(double)) ;

    for(i=0; i < rows; i++)
    {
        (new_mat.tab)[i] = (double *) calloc(cols, sizeof(double));
    }

    return new_mat ;
}

extern void of_free_matrix(Matrix mat)
{
    unsigned int i;
    for(i = 0; i< mat.rows; i++)
        free( (mat.tab)[i]);
    free(mat.tab);
}
extern int of_set_m(Matrix *mat, unsigned int row, unsigned int col, double val)
{
    if ( row >= mat->rows)
    {
        fprintf(stderr, "set_m: row out of range (%i>=%i)\n", row, mat->rows);
        return -1;
    }
    if ( col >= mat->cols)
    {
        fprintf(stderr, "set_m: col out of range (%i>=%i)\n", col, mat->cols);
        return -2 ;
    }
    (mat->tab)[row][col] = val;
    return 0;
}
extern int of_get_m(Matrix *mat, unsigned int row, unsigned int col, double *val)
{
    if ( row >= mat->rows)
    {
        fprintf(stderr, "get_m: row out of range (%i>=%i)\n", row, mat->rows);
        return -1;
    }
    if ( col >= mat->cols)
    {
        fprintf(stderr, "get_m: col out of range (%i>=%i)\n", col, mat->cols);
        return -2 ;
    }
    *val = (mat->tab)[row][col];
    return 0;
}

extern void of_sum_mat(Matrix *mat1, Matrix *mat2, Matrix *mat)
{
    unsigned int rows, cols;
    unsigned int i,j;

    if ( mat1->rows > mat2->rows)
        rows = mat2->rows;
    else
        rows = mat1->rows;
    
    if ( mat1->cols > mat2->cols)
        cols = mat2->cols;
    else
        cols = mat1->cols;


    for(i=0; i< rows; i++)
    {
        for (j=0; j<cols; j++)
        {
            mat->tab[i][j] = mat1->tab[i][j] + mat2->tab[i][j];
        }
    }
}



extern void of_product_mat_el(Matrix *mat1, Matrix *mat2, Matrix *mat)
{
    unsigned int rows, cols;
    unsigned int i,j;

    if ( mat1->rows > mat2->rows)
        rows = mat2->rows;
    else
        rows = mat1->rows;
    
    if ( mat1->cols > mat2->cols)
        cols = mat2->cols;
    else
        cols = mat1->cols;


    for(i=0; i< rows; i++)
    {
        for (j=0; j<cols; j++)
        {
             mat->tab[i][j] = ( mat1->tab[i][j] ) * ( mat2->tab[i][j] );
        }
    }
}



extern int of_product_mat(Matrix *mat1, Matrix *mat2, Matrix * mat)
{
    unsigned int i,j, k;
    //TODO: check that mat1->cols == mat2->rows
    //mat = new_matrix(mat1->rows, mat2->cols);

    for(i=0; i < mat1->rows; i++)
    {
        for (j=0; j < mat2->cols; j++)
        {
            //set_m(mat, i, j, 0.0);
            mat->tab[i][j] = 0.0 ;
            for (k = 0; k < mat1->cols; k++)
               mat->tab[i][j] += ( ( mat1->tab[i][k] ) * ( mat2->tab[k][j] ) ); 
        }
    }
    return 0;
}


extern double of_trace(Matrix *mat)
{
    double t = 0.0 ;
    unsigned int n, i;
    if ( mat->rows > mat->cols)
        n = mat->cols;
    else
        n = mat->rows;
    
    for ( i=0; i < n; i++)
        t += mat->tab[i][i];

    return t;
}

extern double of_tracesquare(Matrix *mat)
{

    double t = 0.0 ;
    unsigned n, i, k;
    n = mat->rows;
    
    for ( i = 0; i < n; i++)
    {
        for ( k = 0; k < n; k++)
        {
            t += ( mat->tab[i][k] ) * ( mat->tab[k][i] ) ;
        }
    }
    return t ;
}

extern void of_print_mat(Matrix * mat)
{
    unsigned i, j;
    printf("[ ");
    for ( i=0; i< mat->rows; i++)
    {
        for ( j=0; j < mat->cols; j++)
        {
            printf("%f ", mat->tab[i][j]);
        }
        if ( i== mat->rows - 1 )
            printf("]");
        printf("\n");
    }
}

