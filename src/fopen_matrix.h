#ifndef FOPEN_MATRIX_H
#define FOPEN_MATRIX_H

#include "matrix.h"

#include <cassert>
#include <cstdio>
#include <cstring>

namespace FopenMatrix {
  namespace internal {
#define MAX_LINE 100000
#define LABEL_LEN 64

    int fscanf2(FILE *file, short *x);

    int fscanf2(FILE *file, int *x);

    int fscanf2(FILE *file, float *x);

    int fscanf2(FILE *file, double *x);

    template<typename T>
    Matrix<T> load_matrix_from_file(FILE *fp, std::size_t reserved_count, bool ignore_first) {
      Matrix<T> matrix(reserved_count); // RVO

      char line[MAX_LINE];
      if (fgets(line, MAX_LINE, fp) == NULL)  throw - 1;

      const char delims[] = " \t\r\n";
      char *pch = strtok(line, delims);
      if (pch == NULL) throw - 1;
      if (ignore_first) pch = strtok(NULL, delims); // ignore the first value
      while (pch != NULL) {
        matrix.col_names.push_back(pch);
        pch = strtok(NULL, delims);
      }

      char value[LABEL_LEN];
      while (1 == fscanf(fp, "%s", value)) {
        matrix.row_names.push_back(value);
        std::vector<T> line_data(matrix.col_names.size());
        for (std::size_t i = 0; i < matrix.col_names.size(); i++) {
          fscanf2(fp, &line_data[i]);
        }
        matrix.data.emplace_back(line_data);
      }
      assert(matrix.row_names.size() == matrix.data.size());
      return matrix; // RVO
    }
  }

  template<typename T> Matrix<T> load_matrix(const char* file_name, std::size_t reserved_count, bool ignore_first = true) {
    FILE *fp = fopen(file_name, "r");
    if (NULL == fp) {
      printf("Failed to open '%s'", file_name);
      throw - 1;
    }
    Matrix<T> matrix = internal::load_matrix_from_file<T>(fp, reserved_count, ignore_first); // RVO
    fclose(fp);
    return matrix; // RVO
  }

  template<typename T> Matrix<T> load_matrix(const char* file_name, bool ignore_first = true) {
    return load_matrix<T>(file_name, 4096, ignore_first); // RVO
  }
}
#endif
