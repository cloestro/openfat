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


#ifndef MATRIX_H
#define MATRIX_H
#include <stdarg.h>

typedef struct of_mat{ double ** tab; unsigned int rows; unsigned int cols; } Matrix;

extern Matrix of_new_matrix(unsigned int rows, unsigned int cols);

extern void of_free_matrix( Matrix my_matrix );

extern int of_set_m(Matrix * mat, unsigned int row, unsigned int col, double val);
extern int of_get_m(Matrix * mat, unsigned int row, unsigned int col, double * val);
extern void of_sum_mat(Matrix * mat1, Matrix * mat2, Matrix * mat);
extern void of_product_mat_el(Matrix * mat1, Matrix * mat2, Matrix *mat);
extern int of_product_mat(Matrix * mat1, Matrix * mat2, Matrix * mat);
extern double of_trace(Matrix * mat1);
extern double of_tracesquare(Matrix * mat);
extern void of_print_mat(Matrix * mat);
#endif

