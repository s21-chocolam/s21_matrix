#include "s21_matrix.h"

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int flag = OK;
  if (rows <= 0 || columns <= 0 || result == NULL) {
    return flag = INCORRECT_MATRIX;
  }

  result->rows = rows;
  result->columns = columns;
  result->matrix =
      (double **)calloc(rows, sizeof(double *) + columns * sizeof(double));

  if (result->matrix != NULL) {
    result->matrix[0] = (double *)(result->matrix + rows);
    for (int i = 1; i < rows; i++) {
      result->matrix[i] = result->matrix[0] + i * columns;
    }
  } else {
    flag = INCORRECT_MATRIX;
  }

  return flag;
}

void s21_remove_matrix(matrix_t *A) {
  if (A != NULL && A->matrix != NULL) {
    free(A->matrix);
    A->rows = 0;
    A->columns = 0;
    A->matrix = NULL;
  }
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int flag = SUCCESS;
  if (is_valid_matrix(A) || is_valid_matrix(B)) {
    return flag = FAILURE;
  }
  if (A->rows != B->rows || A->columns != B->columns) {
    flag = FAILURE;
  } else {
    for (int i = 0; i < A->rows && flag; i++) {
      for (int j = 0; j < A->columns && flag; j++) {
        if (fabs(A->matrix[i][j] - B->matrix[i][j]) > EPS) {
          flag = FAILURE;
        }
      }
    }
  }
  return flag;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int flag = OK;
  if (is_valid_matrix(A) || is_valid_matrix(B) || result == NULL) {
    return flag = INCORRECT_MATRIX;
  }
  if (A->rows == B->rows && A->columns == B->columns) {
    s21_create_matrix(A->rows, A->columns, result);
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
      }
    }
  } else {
    flag = CALCULATION_ERROR;
  }
  return flag;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int flag = OK;
  if (is_valid_matrix(A) || is_valid_matrix(B) || result == NULL) {
    return flag = INCORRECT_MATRIX;
  }
  if (A->rows == B->rows && A->columns == B->columns) {
    s21_create_matrix(A->rows, A->columns, result);
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
      }
    }
  } else {
    flag = CALCULATION_ERROR;
  }
  return flag;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  int flag = OK;
  if (is_valid_matrix(A) || result == NULL) {
    return flag = INCORRECT_MATRIX;
  }
  s21_create_matrix(A->rows, A->columns, result);
  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < A->columns; j++) {
      result->matrix[i][j] = A->matrix[i][j] * number;
    }
  }
  return flag;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int flag = OK;
  if (is_valid_matrix(A) || is_valid_matrix(B) || result == NULL) {
    return flag = INCORRECT_MATRIX;
  }
  if (A->rows == B->columns && A->columns == B->rows) {
    s21_create_matrix(A->rows, B->columns, result);
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < B->columns; j++) {
        result->matrix[i][j] = 0;
        for (int k = 0; k < A->columns; k++) {
          result->matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
        }
      }
    }
  } else {
    flag = CALCULATION_ERROR;
  }
  return flag;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  int flag = OK;
  if (is_valid_matrix(A) || result == NULL) {
    return flag = INCORRECT_MATRIX;
  }
  s21_create_matrix(A->columns, A->rows, result);
  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < A->columns; j++) {
      result->matrix[j][i] = A->matrix[i][j];
    }
  }
  return flag;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  int flag = OK;
  if (is_valid_matrix(A) || result == NULL) {
    return flag = INCORRECT_MATRIX;
  }
  if (A->columns == A->rows) {
    s21_create_matrix(A->columns, A->rows, result);
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        double determinant = 0;
        matrix_t reduced_matrix;
        s21_create_matrix(A->rows - 1, A->columns - 1, &reduced_matrix);
        create_reduced_matrix(i, j, A, &reduced_matrix);
        s21_determinant(&reduced_matrix, &determinant);
        result->matrix[i][j] = pow((-1), i + j) * determinant;
        s21_remove_matrix(&reduced_matrix);
      }
    }
  } else {
    flag = CALCULATION_ERROR;
  }
  return flag;
}

int s21_determinant(matrix_t *A, double *result) {
  int flag = OK;
  if (is_valid_matrix(A) || result == NULL) {
    return flag = INCORRECT_MATRIX;
  }
  if (A->rows != A->columns) {
    flag = CALCULATION_ERROR;
  } else {
    *result = get_determinant(A);
  }
  return flag;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  int flag = OK;
  double determinant = 0;
  if (is_valid_matrix(A) || result == NULL) {
    return flag = INCORRECT_MATRIX;
  } else {
    int check_determinant = s21_determinant(A, &determinant);
    if (check_determinant || determinant == 0) {
      flag = CALCULATION_ERROR;
    } else {
      matrix_t calc_compl;
      s21_calc_complements(A, &calc_compl);
      matrix_t transp_matrix;
      s21_transpose(&calc_compl, &transp_matrix);
      s21_mult_number(&transp_matrix, 1 / determinant, result);
      s21_remove_matrix(&calc_compl);
      s21_remove_matrix(&transp_matrix);
    }
  }
  return flag;
}

double get_determinant(matrix_t *A) {
  double determinant = 0;
  if (A->rows == 1) {
    determinant = A->matrix[0][0];
  } else {
    matrix_t reduced_matrix;
    s21_create_matrix(A->rows - 1, A->columns - 1, &reduced_matrix);
    for (int i = 0; i < A->rows; i++) {
      create_reduced_matrix(0, i, A, &reduced_matrix);
      if (i % 2 == 0) {
        determinant += A->matrix[0][i] * get_determinant(&reduced_matrix);
      } else {
        determinant -= A->matrix[0][i] * get_determinant(&reduced_matrix);
      }
    }
    s21_remove_matrix(&reduced_matrix);
  }
  return determinant;
}

void create_reduced_matrix(int rows, int columns, matrix_t *A,
                           matrix_t *result) {
  int result_rows = 0;
  for (int i = 0; i < A->rows; i++) {
    if (i == rows) {
      continue;
    }
    int result_columns = 0;
    for (int j = 0; j < A->columns; j++) {
      if (j == columns) {
        continue;
      }
      result->matrix[result_rows][result_columns] = A->matrix[i][j];
      result_columns += 1;
    }
    result_rows += 1;
  }
}

int is_valid_matrix(matrix_t *matrix) {
  int flag = 0;
  if (matrix == NULL || matrix->matrix == NULL || matrix->rows <= 0 ||
      matrix->columns <= 0) {
    flag = 1;
  } else {
    flag = 0;
  }
  return flag;
}