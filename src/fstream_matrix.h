#ifndef FSTREAM_MATRIX_H
#define FSTREAM_MATRIX_H

#include "matrix.h"

namespace FstreamMatrix {
  template<typename T> Matrix<T> load_matrix(const char* file_name, std::size_t reserved_count = 4096);

  namespace internal {
#include <fstream>
#include <cassert>

    template<typename T> Matrix<T> load_matrix_from_ifstream(std::ifstream &infile, std::size_t reserved_count) {
      Matrix<T> matrix(reserved_count); // RVO

      std::string line;
      if (!std::getline(infile, line)) throw - 1;

      std::istringstream iss(line);
      std::string value;
      iss >> value; // ignore this value

      matrix.col_names.reserve(reserved_count);
      matrix.row_names.reserve(reserved_count);
      matrix.data.reserve(reserved_count);

      while (iss >> value) {
        matrix.col_names.push_back(value);
      }

      while (std::getline(infile, line)) {
        std::istringstream iss(line);
        if (!(iss >> value)) {
          break;  // error
        }
        matrix.row_names.push_back(value);

        std::vector<float> line_data;
        line_data.reserve(matrix.col_names.size());
        float tmp;
        while (iss >> tmp) {
          line_data.push_back(tmp);
        }
        assert(matrix.col_names.size() == line_data.size());
        matrix.data.push_back(line_data);
      }
      assert(matrix.row_names.size() == matrix.data.size());
      return matrix; // RVO
    }
  }

  template<typename T>
  Matrix<T> FstreamMatrix::load_matrix(const char* file_name, std::size_t reserved_count) {
    std::ifstream infile(file_name);
    Matrix<T> matrix = internal::load_matrix_from_ifstream(infile, reserved_count); // RVO
    return matrix; // RVO
  }
}
#endif
