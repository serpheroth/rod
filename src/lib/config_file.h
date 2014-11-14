#ifndef CONFIG_FILE_H
#define CONFIG_FILE_H
#include <map>
#include <vector>
#include <iostream>
#include <typeinfo>
#include <cstdlib>
#include "macro_constant.h"
#include "print_macro.h"

class ConfigFile;
template <class T> struct Adaptor {
  Adaptor(ConfigFile* conf) : conf_(conf) {}
  T operator()(std::string& str) {
    (void) str;
    return T();
  }

  std::map<std::string, std::vector<T> >* GetMap() {
    return NULL;
  }
private:
  ConfigFile* conf_;
};

template <> struct Adaptor<float>;
template <> struct Adaptor<double>;
template <> struct Adaptor<int>;
template <> struct Adaptor<std::string>;
//template <> struct Adaptor<bool>;

template <class T>
struct TypeTraits {
  typedef T value_type;
  typedef T return_type;
  static return_type Return(std::vector<T>& vec) {
    return vec[0];
  }
};

template <class T>
struct TypeTraits<T*> {
  typedef T value_type;
  typedef T* return_type;
  static return_type Return(std::vector<T>& vec) {
    return &vec[0];
  }
};

template <class T>
struct TypeTraits<std::vector<T>*> {
  typedef T value_type;
  typedef std::vector<T>* return_type;
  static return_type Return(std::vector<T>& vec) {
    return &vec;
  }
};


class ConfigFile {
private:
  template <class T>
  typename TypeTraits<T>::return_type
  Get(std::string& key, Adaptor<typename TypeTraits<T>::value_type>& conv) {
    typedef typename TypeTraits<T>::value_type value_type;
    typedef std::vector<value_type> vector;
    vector values;
    for (unsigned i = 0; i < raw_values_[key].size(); ++i) {
      int line_no = raw_values_[key][i];
      std::string& line = lines_[line_no];
      std::stringstream ss;
      value_type v;
      ss << line;
      ss >> v;
      values.push_back(v);
    }
    std::map<std::string, vector>* map = conv.GetMap();
    (*map)[key] = values;
    return TypeTraits<T>::Return((*map)[key]);
  }

public:
  ConfigFile(const char* config_file);
  void Save(const char* file_name);

  template <class T>
  typename TypeTraits<T>::return_type
  Get(const char* key) {
    typedef typename TypeTraits<T>::value_type value_type;
    typedef typename TypeTraits<T>::return_type return_type;
    std::string str_key(key);
    if (raw_values_.find(str_key) == raw_values_.end()) {
      std::cerr << CURRENT_LINE << " => property " << key << " does not exist in the conf file!!!" << std::endl;
      exit(0);
      return return_type(0);
    } else {
      Adaptor<value_type> adaptor(this);
      return Get<T>(str_key, adaptor);
    }
  }


  template <class T>
  void Put(const char *key, T value) {
    static_assert(std::is_pointer<T>::value == false, "Put only works for non-pointer type");
    std::string str_key(key);
    Adaptor<T> adaptor(this);
    T* v = Get<T*>(str_key, adaptor);
    v[0] = value;
  }

private:
  void Flush();
  template <class T>
  void Flush(std::map<std::string, std::vector<T> > &map,
             std::string (*Adaptor)(T));

  friend std::ostream& operator<<(std::ostream& out, ConfigFile& conf);


  void AddOption(std::string &key, std::vector<int> &values);
  friend  struct Adaptor<int>;
  //  friend  struct Adaptor<bool>;
  friend  struct Adaptor<float>;
  friend  struct Adaptor<double>;
  friend  struct Adaptor<std::string>;

  std::map<std::string, std::vector<int> > raw_values_;
  std::map<std::string, std::vector<float> > floats_;
  std::map<std::string, std::vector<double> > doubles_;
  std::map<std::string, std::vector<int> > ints_;
  std::map<std::string, std::vector<std::string> > strings_;
  //  std::map<std::string, std::vector<char> > bools_;
  std::vector<std::string> lines_;
};



template <> struct Adaptor<float> {
  Adaptor(ConfigFile* conf) : conf_(conf) {}
  float operator()(std::string& str) {
    std::stringstream ss;
    float result;
    ss << str;
    ss >> result;
    return result;
  }
  std::map<std::string, std::vector<float> >* GetMap() {
    return &conf_->floats_;
  }
private:
  ConfigFile* conf_;
};

template <> struct Adaptor<double> {
  Adaptor(ConfigFile* conf) : conf_(conf) {}
  double operator()(std::string& str) {
    return atof(str.c_str());
  }
  std::map<std::string, std::vector<double> >* GetMap() {
    return &conf_->doubles_;
  }
private:
  ConfigFile* conf_;
};

template <> struct Adaptor<int> {
  Adaptor(ConfigFile* conf) : conf_(conf) {}
  int operator()(std::string& str) {
    return atoi(str.c_str());
  }
  std::map<std::string, std::vector<int> >* GetMap() {
    return &conf_->ints_;
  }
private:
  ConfigFile* conf_;
};

template <> struct Adaptor<std::string> {
  Adaptor(ConfigFile* conf) : conf_(conf) {}
  std::string operator()(std::string& str) {
    return str;
  }
  std::map<std::string, std::vector<std::string> >* GetMap() {
    return &conf_->strings_;
  }
private:
  ConfigFile* conf_;
};

//template <> struct Adaptor<bool> {
//  Adaptor(ConfigFile* conf) : conf_(conf) {}
//  bool operator()(std::string& str)
//  {
//    return (str == "true") ? true : false;
//  }
//  std::map<std::string, std::vector<char> >* GetMap()
//  {
//    return &conf_->bools_;
//  }
//private:
//  ConfigFile* conf_;
//};

std::ostream& operator<<(std::ostream& out, ConfigFile& conf);

#endif // CONFIG_FILE_H
