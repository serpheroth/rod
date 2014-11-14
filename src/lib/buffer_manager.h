#ifndef BUFFER_MANAGER_H_
#define BUFFER_MANAGER_H_
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <new>
#include <iostream>
#include <vector>

template <class T>
class SmartArrayPointer
{
public:
  SmartArrayPointer(T* pointer) : pointer_(pointer) {}
  ~SmartArrayPointer(void)
  {
    delete[] pointer_;
  }
  inline T& operator[](int index)
  {
    return pointer_[index];
  }
private:
  T* pointer_;
};

class BufferManager
{
public:
  // @BufferManager@(void) {
  BufferManager(void)
  {
    buffer_size_.clear();
    pointers_.clear();
  } // #BufferManager#

  // @~BufferManager@(void) {
  ~BufferManager(void)
  {
    Clear();
  } // #~BufferManager#

  template <class T>
  void Free(T* pointer)
  {
    for (unsigned int i = 0; i < pointers_.size(); i++) {
      if (pointers_[i] == (void*) pointer) {
        free(pointer);
        buffer_size_.erase(buffer_size_.begin() + i);
        pointers_.erase(pointers_.begin() + i);
        break;
      }
    }
  }

  template <class T>
  void Push(T* pointer)
  {
    pointers_.push_back((void*) pointer);
    buffer_size_.push_back(0);
  }

  //  // T* @New@(unsigned int size) {
  //  template <class T>
  //  T* New(unsigned int size) {
  //    T* buffer = NULL:
  //    if (size == 1) {
  //     buffer = new (std::no_throw) T();
  //    } else {
  //     buffer = new (std::no_throw) T[size];
  //    }
  //    if (buffer == NULL) {
  //      std::cerr << "Failed to allocate memory in BufferManager::New()!" << std::endl;
  //      getchar();
  //    }
  //    return buffer;
  //  } // #New#
  //
  //  // T* @New@(unsigned int size) {
  //  template <class T>
  //  T* New(unsigned int size, const T& initial_value) {
  //    T* buffer = NULL:
  //    if (size == 1) {
  //     buffer = new (std::no_throw) T(initial_value);
  //    } else {
  //     buffer = new (std::no_throw) T[size];
  //    }
  //    if (buffer == NULL) {
  //      std::cerr << "Failed to allocate memory in BufferManager::New()!" << std::endl;
  //      getchar();
  //    } else if (size > 1) {
  //      for (int i = 0; i < size; i++) {
  //        buffer[i] = initial_value;
  //      }
  //    }
  //    return buffer;
  //  } // #New#

  // T* @Malloc2D@(int n)
  template <class T>
  T** Malloc2D(int n, T initial_value = T())
  {
    T** buffer_pointer = (T**) malloc(sizeof(T*) * n);
    T* buffer = (T*) malloc(sizeof(T) * n * n);
    for (int i = 0; i < n * n; i++) {
      *(buffer + i) = initial_value;
    }
    for (int i = 0; i < n; i++) {
      buffer_pointer[i] = buffer + n * i;
    }
    pointers_.push_back(buffer_pointer);
    pointers_.push_back(buffer);
    buffer_size_.push_back(sizeof(T*) * n);
    buffer_size_.push_back(sizeof(T) * n * n);
    return buffer_pointer;
  } // #Malloc2D#

  // T** @Malloc2D@(int x, int y)
  template <class T>
  T** Malloc2D(int x, int y, T initial_value = T())
  {
    T** buffer_pointer = (T**) malloc(sizeof(T*) * x);
    T* buffer = (T*) malloc(sizeof(T) * x * y);
    for (int i = 0; i < x * y; i++) {
      *(buffer + i) = initial_value;
    }
    for (int i = 0; i < x; i++) {
      buffer_pointer[i] = buffer + y * i;
    }
    pointers_.push_back(buffer_pointer);
    pointers_.push_back(buffer);
    buffer_size_.push_back(sizeof(T*) * x);
    buffer_size_.push_back(sizeof(T) * x * y);
    return buffer_pointer;
  } // #Malloc2D#


  // T* @Malloc@(unsigned int size) {
  template <class T>
  T* Malloc(unsigned int size)
  {
    void* buffer = malloc(sizeof(T) * size);
    if (buffer == NULL) {
      std::cerr << (size) << std::endl;
      std::cerr << "BufferManager::Malloc() => failed to allocate memory." << std::endl;
      getchar();
      exit(0);
    }
    memset(buffer, 0, sizeof(T) * size);
    pointers_.push_back(buffer);
    buffer_size_.push_back(sizeof(T) * size);
    return (T*) buffer;
  } // #Malloc#

  // T*  @Malloc@(unsigned int size, T initial_value) {
  template <class T>
  T*  Malloc(unsigned int size, const T initial_value)
  {
    void* buffer = malloc(sizeof(T) * size);
    if (buffer == NULL) {
      std::cerr << "BufferManager::Malloc() => failed to allocate memory." << std::endl;
      getchar();
      exit(0);
    }
    for (unsigned int i = 0; i < size; i++) {
      ((T*) buffer)[i] = initial_value;
    }
    pointers_.push_back(buffer);
    buffer_size_.push_back(sizeof(T) * size);
    return (T*) buffer;
  } // #Malloc#

  // T* @Malloc@(unsigned int size, T* initial_value) {
  template <class T>
  T*  Malloc(unsigned int size, T* initial_value)
  {
    void* buffer = malloc(sizeof(T) * size);
    if (buffer == NULL) {
      std::cerr << "BufferManager::Malloc() => failed to allocate memory." << std::endl;
      getchar();
      exit(0);
    }
    for (unsigned int i = 0; i < size; i++) {
      ((T*) buffer)[i] = *initial_value;
    }
    pointers_.push_back(buffer);
    buffer_size_.push_back(sizeof(T) * size);
    return (T*) buffer;
  } // #Malloc#

  int BufferSize() {
    int result = 0;
    for (unsigned i = 0; i < buffer_size_.size(); ++i) {
      result += buffer_size_[i];
    }
    return result;
  }

  // void @Clear@(void) {
  void Clear(void)
  {
    for (unsigned int i = 0; i < pointers_.size(); i++) {
      free(pointers_[i]);
    }
    buffer_size_.clear();
    pointers_.clear();
  } // #Clear#
private:
  BufferManager(const BufferManager&);
  BufferManager& operator=(const BufferManager&);
  std::vector<int> buffer_size_;
  std::vector<void*> pointers_;
};
#endif // BUFFER_MANAGER_H_
