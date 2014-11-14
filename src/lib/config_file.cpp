#include <iostream>
#include <iterator>
#include <fstream>
#include "config_file.h"
#include "macro_constant.h"
#include "print_macro.h"

std::ostream& operator<<(std::ostream& out, ConfigFile& conf)
{
  conf.Flush();
  std::ostream_iterator<std::string> out_iter(out, "\n");
  std::copy(conf.lines_.begin(), conf.lines_.end(), out_iter);
  return out;
}


// Truncate leading and trailing spaces
void TruncateSpace(std::string& str)
{
  int begin = 0;
  int end = str.length() - 1;
  for (; begin < int(str.length()) && isspace(str[begin]); ++begin) {
  }
  for (; end >= begin && isspace(str[end]); --end) {
  }
  str = str.substr(begin, end - begin + 1);
}

//std::string Bool2Str(bool value)
//{
//  return (value ? "true" : "false");
//}

//bool StartWith(std::string& str, char c)
//{
//  for (unsigned int i = 0; i < str.length(); ++i) {
//    if (str[i] != ' ' && str[i] != '\t') {
//      return str[i] == c;
//    }
//  }
//  return false;
//}

bool IsEmptyLine(std::string& line)
{
  for (int i = 0; i < int(line.size()); ++i) {
    if (!isspace(line[i])) {
      return false;
    }
  }
  return true;
}

char FirstNonSpaceChar(std::string& str)
{
  char first_char = 0;
  for (unsigned int i = 0; i < str.length(); ++i) {
    if (!isspace(str[i])) {
      first_char = str[i];
      break;
    }
  }
  //  P(str.length(), str, (int) first_char);
  return first_char;
}

ConfigFile::ConfigFile(const char* config_file)
{
  std::ifstream in(config_file);
  if (!in.is_open()) {
    std::cerr << CURRENT_LINE << "=> failed to open file " << config_file << std::endl;
    exit(0);
  }
  std::string key = "";
  std::vector<int> value_lines;
  while (!in.eof()) {
    std::string line;
    getline(in, line);
    lines_.push_back(line);
    TruncateSpace(line);
    char first_char = FirstNonSpaceChar(line);
    if (first_char == '#' || first_char == '{' || first_char == '}' || first_char == 0) { // comment or empty line
//      if (value_lines.size() > 0) {
//        AddOption(key, value_lines);
//      }
    } else if (first_char == '*') { // new key
      AddOption(key, value_lines);
      key = line.substr(1);
      TruncateSpace(key);
    } else { // value
      value_lines.push_back(lines_.size() - 1);
    }
  }
  AddOption(key, value_lines);
  in.close();
}

void ConfigFile::Save(const char *file_name)
{
  std::ofstream out(file_name);
  if (!out.is_open()) {
    std::cerr << CURRENT_LINE << "=> failed to open file " << file_name << std::endl;
    exit(0);
  }
  out << *this;
  out.close();
}

void ConfigFile::Flush()
{
  Flush<int>(ints_, NULL);
  Flush<float>(floats_, NULL);
  Flush<double>(doubles_, NULL);
  Flush<std::string>(strings_, NULL);
  //  Flush<bool>(bools_, Bool2Str);
}


void ConfigFile::AddOption(std::string& key, std::vector<int>& values)
{
  if (!IsEmptyLine(key)) {
    if (values.size() == 0) {
      std::cerr << CURRENT_LINE << " => " << key << " does not have corresponding values!!!" << std::endl;
      exit(0);
    } else if (raw_values_.find(key) != raw_values_.end()) {
      std::cerr << CURRENT_LINE << " => " << key << " has duplicated values!!!" << std::endl;
      exit(0);
    } else {
      raw_values_[key] = values;
      key = "";
      values.clear();
    }
  }
}


template <class T>
void ConfigFile::Flush(std::map<std::string, std::vector<T> > &map, std::string (*Converter)(T))
{
  typename std::map<std::string, std::vector<T> >::iterator iter = map.begin();
  for (; iter != map.end(); ++iter) {
    std::vector<int>& line_no = raw_values_[iter->first];
    ASSERT(iter->second.size() == line_no.size());
    for (unsigned i = 0; i < iter->second.size(); ++i) {
      if (Converter == NULL) {
        std::stringstream str;
        str << iter->second[i];
        lines_[line_no[i]] = str.str();
      } else {
        lines_[line_no[i]] = Converter(iter->second[i]);
      }
    }
  }
}
