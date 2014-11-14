#ifndef _VECTOR_IO_H_
#define _VECTOR_IO_H_
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <utility>
#include <vector>
#include "macro_constant.h"
enum {
  kRow, kColumn
};
const char kSpace = '\t';

class TextFileWriter {
public:
  TextFileWriter(const char* file_name) : text_file_(file_name) {
    if (!text_file_.is_open()) {
      std::cerr << CURRENT_LINE << " => failed to open file " << file_name << std::endl;
      exit(0);
    }
  }

  template <class T>
  TextFileWriter& operator<<(T& v) {
    text_file_ << v;
    return *this;
  }

  ~TextFileWriter() {
    text_file_.close();
  }

private:
  std::ofstream text_file_;
};

class TextFileReader {
public:
  TextFileReader(const char* file_name) : text_file_(file_name)
  {
    if (!text_file_.is_open()) {
      std::cerr << CURRENT_LINE << " => failed to open file " << file_name << std::endl;
      exit(0);
    }
  }

  template <class T>
  TextFileReader& operator>>(T& v) {
    text_file_ >> v;
    return *this;
  }

  ~TextFileReader() {
    text_file_.close();
  }

private:
  std::ifstream text_file_;
};


class BinaryFileReader {
public:
  BinaryFileReader(const char* file_name) {
    std::ifstream in(file_name, std::ios::binary | std::ios::in);
    if (!in.is_open()) {
      std::cerr << "BinaryFileReader::BinaryFileReader() => failed to open input file " << file_name << std::endl;
      exit(0);
    }
    // Get file size
    in.seekg(0, in.end);
    int size = int(in.tellg());
    // Allocate buffer
    buf_.resize(size);
    ptr_ = &buf_[0];
    // Read file data to buffer
    in.seekg(0, in.beg);
    in.read(ptr_, size);
    in.close();
  }

  template <class T>
  inline void Read(T* data, int element_count) {
    if (!Valid()) {
      std::cerr << "BinaryFileReader::Read() => no more data to read from file" << std::endl;
      exit(0);
    }
    int size_in_bytes = sizeof(T) * element_count;
    memcpy(data, ptr_, size_in_bytes);
    ptr_ += size_in_bytes;
  }


  template <class T>
  inline void Read2DVector(std::vector<std::vector<T> >& vector) {
    int size1 = *((int*) ptr_);
    ptr_ += sizeof(int);
    vector.resize(size1);
    for (int i = 0; i < int(vector.size()); ++i) {
      int size2 = *((int*) ptr_);
      ptr_ += sizeof(int);
      vector[i].resize(size2);
      Read<T>(&vector[i][0], size2);
    }
  }

  template <class T>
  inline void Read1DVector(std::vector<T>& vector) {
    int size = *((int*) ptr_);
    ptr_ += sizeof(int);
    vector.resize(size);
    Read<T>(&vector[0], size);
  }

  bool Valid() {
    return (ptr_ - &buf_[0]) < (int) buf_.size();
  }

private:
  char* ptr_;
  std::vector<char> buf_;
};

template <class T>
inline void ReadBuffer(T* dest, char*& buf, int count) {
  memcpy(dest, buf, sizeof(T) * count);
  buf += sizeof(T) * count;
}//  fp = fopen(filename, "r+");
//  if (fp == NULL)	{
//    printf("ERROR: file %s not open.\n", filename);
//    return;
//  }
//  fscanf(fp, "%d %d %d %d\n", &vertex_num_, &temp_value, &temp_value, &bound);
//  P(bound, vertex_num_);
//  if (bound == 0) {
//    for (int i = 0; i < vertex_num_; i++) {
//      float temp_x0, temp_x1, temp_x2;
//      fscanf(fp, "%d %f %f %f\n", &temp_value, &temp_x0, &temp_x1, &temp_x2);
//      if (i == 0) one_indexed_ = (temp_value == 0) ? 0 : 1;
//      X[i * 3 + 0] = temp_x0;
//      X[i * 3 + 1] = temp_x1;
//      X[i * 3 + 2] = temp_x2;
//      //      P(dj::Vec3f(X + i * 3));
//      printf("en %d: %f, %f, %f\n", i, X[i * 3], X[i * 3 + 1], X[i * 3 + 2]);
//    }
//  } else {
//    for (int i = 0; i < vertex_num_; i++) {
//      float temp_x0, temp_x1, temp_x2;
//      fscanf(fp, "%d %f %f %f %d\n", &temp_value, &temp_x0, &temp_x1, &temp_x2, &temp_value);
//      if (i == 0) one_indexed_ = (temp_value == 0) ? 0 : 1;
//      X[i * 3 + 0] = temp_x0;
//      X[i * 3 + 1] = temp_x1;
//      X[i * 3 + 2] = temp_x2;
//      printf("en %d: %f, %f, %f\n", i, X[i * 3], X[i * 3 + 1], X[i * 3 + 2]);
//    }
//  }
//  fclose(fp);
//  exit(0);
  //for(int i=0; i<number; i++)
  //	printf("v %d: %f, %f, %f\n", i, X[i*3], X[i*3+1], X[i*3+2]);


template <class T>
inline void Write2DVector(std::ostream& out, std::vector<std::vector<T> >& vector) {
  int size1 = vector.size();
  out.write((char*) &size1, sizeof(int));
  for (int i = 0; i < int(vector.size()); ++i) {
    int size2 = (int) vector[i].size();
    out.write((char*) &size2, sizeof(int));
    if (size2 > 0) {
      out.write((char*) &vector[i][0], sizeof(T) * size2);
    }
  }
}

template <class T>
inline void Read2DVector(std::istream& in, std::vector<std::vector<T> >& vector) {
  int size1 = -1;
  in.read((char*) &size1, sizeof(int));
  vector.resize(size1);
  for (int i = 0; i < int(vector.size()); ++i) {
    int size2 = -1;
    in.read((char*) &size2, sizeof(int));
    vector[i].resize(size2);
    if (size2 > 0) {
      in.read((char*) &vector[i][0], sizeof(T) * size2);
    }
  }
}

template <class T>
inline void Write1DVector(std::ostream& out, std::vector<T>& vector) {
  int size = (int) vector.size();
  out.write((char*) &size, sizeof(int));
  out.write((char*) &vector[0], sizeof(T) * size);
}

template <class T>
inline void Read1DVector(std::istream& in, std::vector<T>& vector) {
  int size = -1;
  in.read((char*) &size, sizeof(T));
  vector.resize(size);
  in.read((char*) &vector[0], sizeof(T) * size);
}

template <class T>
inline bool WriteMatrix(int nx, int ny, T* matrix, std::ostream& out) {
  out << "[";
  for (int row = 0, offset = 0; row < nx; ++row) {
    for (int col = 0; col < ny; ++col, ++offset) {
      out << matrix[offset] << " ";
    }
    out << "; ";
  }
  out << "]" << std::endl;
  return true;
} //#WriteVector#

//------------------------------------------------------------------------------
// Write vector to a file
//------------------------------------------------------------------------------
// bool @WriteVector@(int size, T* vector, std::ostream& out, int format = kColumn)
template <class T>
bool WriteVector(int size, T* vector, std::ostream& out, int format = kColumn) {
  char space_char = (format == kColumn) ? '\n' : kSpace;
  out << size << std::endl;
  for (int i = 0; i < size; i++) {
    out << vector[i] << space_char;
  }
  if (space_char == kSpace) {
    out << std::endl;
  }
  return true;
} //#WriteVector#

// bool @WriteVector@(int size, T* vector, const char* file_name, int format = kColumn)
template <class T>
bool WriteVector(int size, T* vector, const char* file_name, int format = kColumn) {
  std::ofstream output_file(file_name);
  if (!output_file.is_open()) {
    std::cerr << "WriteVector() => failed to open file " << file_name << std::endl;
    return false;
  }
  WriteVector(size, vector, output_file, format);
  output_file.close();
  return true;
} // #WriteVector#

// bool @WriteVectorToMatlab@(int size, T* vector, std::ostream& out, int format = kColumn)
template <class T>
bool WriteVectorToMatlab(int size, T* vector, std::ostream& out, int format = kColumn) {
  char space_char = (format == kColumn) ? '\n' : kSpace;
  for (int i = 0; i < size; i++) {
    out << vector[i] << space_char;
  }
  if (space_char == kSpace) {
    out << std::endl;
  }
  return true;
} // #WriteVectorToMatlab#

// bool @WriteVectorToMatlab@(int size, T* vector, const char* file_name, int format = kColumn)
template <class T>
bool WriteVectorToMatlab(int size, T* vector, const char* file_name, int format = kColumn) {
  std::ofstream output_file(file_name);
  if (!output_file.is_open()) {
    std::cerr << "WriteVectorToMatlab() => failed to open file " << file_name << std::endl;
    return false;
  }
  WriteVectorToMatlab<T>(size, vector, output_file, format);
  output_file.close();
  return true;
} // #WriteVectorToMatlab#

//------------------------------------------------------------------------------
// Read vector from file
//------------------------------------------------------------------------------
// T* @ReadVector@(std::istream& in)
template <class T>
std::pair<int, T*> ReadVector(std::istream& in) {
  int vector_size;
  in >> vector_size;
  if (vector_size <= 0) {
    std::cerr << "ReadVector() => vector size on the first line is " << vector_size << std::endl;
    return std::pair<int, T*>(0, NULL);
  }
  T* vector = (T*) malloc(sizeof(T) * vector_size);
  for (int i = 0; i < vector_size; i++) {
    in >> vector[i];
  }
  return std::pair<int, T*>(vector_size, vector);
} // #ReadVector#

//  T* @ReadVector@(const char* file_name)
template <class T>
std::pair<int, T*> ReadVector(const char* file_name) {
  std::ifstream input_file(file_name);
  if (!input_file.is_open()) {
    std::cerr << "ReadVector() => failed to open input file " << file_name << std::endl;
    return std::pair<int, T*>(0, NULL);
  }
  return ReadVector<T>(input_file);
} // #ReadVector#

// std::pair<int, T*> @ReadVectorFromMatlab@(std::istream& in)
template <class T>
std::pair<int, T*> ReadVectorFromMatlab(std::istream& in) {
  std::vector<T> tmp_vector;
  while (!in.eof()) {
    T tmp;
    in >> tmp;
    if (in.eof()) break;
    tmp_vector.push_back(tmp);
  }
  T* vector = (T*) malloc(sizeof(T) * tmp_vector.size());
  memcpy(vector, &tmp_vector[0], sizeof(T) * tmp_vector.size());
  return std::pair<int, T*>(tmp_vector.size(), vector);
} // #ReadVectorFromMatlab#

// std::pair<int, T*> @ReadVectorFromMatlab@(const char* file_name)
template <class T>
std::pair<int, T*> ReadVectorFromMatlab(const char* file_name) {
  std::ifstream input_file(file_name);
  if (!input_file.is_open()) {
    std::cerr << "ReadVectorFromMatlab() => failed to open input file " << file_name << std::endl;
    return std::pair<int, T*>(0, NULL);
  }
  return ReadVectorFromMatlab<T>(input_file);
} // #ReadVectorFromMatlab#
#endif // _VECTOR_IO_H_
