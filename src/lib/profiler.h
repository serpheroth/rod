#ifndef PROFILER_H
#define PROFILER_H
#include <iostream>
#include <iomanip>
#include <map>
#include <vector>
#include <algorithm>
#include "macro_constant.h"

template <class Timer, class TimeType = float>
class Profiler
{
public:
  Profiler() {
    accumulated_time_.reserve(100);
    number_of_calls_.reserve(100);
  }

  ~Profiler() {
    PrintProfileInfo(std::cout);
  }

  int GetId(const char* name)
  {
    std::string str_name(name);
    std::map<std::string, int>::iterator pos = name2id_.find(str_name);
    if (pos == name2id_.end()) {
      name2id_[str_name] = profile_names_.size();
      profile_names_.push_back(str_name);
      last_start_time_.push_back(TimeType(0));
      accumulated_time_.push_back(TimeType(0));
      number_of_calls_.push_back(0);
      return profile_names_.size() - 1;
    } else {
      return name2id_[str_name];
    }
  }

  inline void StartTimer(int id)
  {
    last_start_time_[id] = timer_.GetTime();
  }

  inline void StartTimer(const char* prof_text) {
    int id = GetId(prof_text);
    StartTimer(id);
  }

  inline void EndTimer(const char* prof_text) {
    int id = GetId(prof_text);
    EndTimer(id);
  }

  inline TimeType EndTimer(int id)
  {
    TimeType elapsed_time = timer_.GetTime() - last_start_time_[id];
    accumulated_time_[id] += elapsed_time;
    number_of_calls_[id]++;
    return elapsed_time;
  }

  double GetNumOfCalls(int id)
  {
    return number_of_calls_[id];
  }

  double GetTimePerCall(int id)
  {
    if (number_of_calls_[id] != 0) {
      return double(accumulated_time_[id]) / number_of_calls_[id];
    } else {
      std::cerr << CURRENT_LINE << " => " << profile_names_[id] << " has never been profiled" << std::endl;
      exit(0);
    }
  }

  TimeType GetAccumulatedTime(int id)
  {
    return accumulated_time_[id];
  }

  void PrintProfileInfo(std::ostream& out = std::cout)
  {
    std::vector<std::pair<double, int> > all_info_;
    all_info_.reserve(accumulated_time_.size());
    double total_time = 0;
    for (unsigned i = 0; i < accumulated_time_.size(); ++i) {
      all_info_.push_back(std::make_pair(-double(accumulated_time_[i]), int(i)));
      total_time += accumulated_time_[i];
    }
    std::sort(all_info_.begin(), all_info_.end());
    out << std::setw(20) << std::left << "name"
        << std::setw(20) << std::left << "accumulated time"
        << std::setw(20) << std::left << "num. of calls"
        << std::setw(20) << std::left << "time per call"
        << std::setw(20) << std::left << "percentage"
        << std::endl;
    for (unsigned i = 0; i < all_info_.size(); ++i) {
      int idx = all_info_[i].second;
      out << std::setw(20) << std::left << profile_names_[idx]
          << std::setw(20) << std::left << accumulated_time_[idx]
          << std::setw(20) << std::left << number_of_calls_[idx]
          << std::setw(20) << std::left << double(accumulated_time_[idx]) / number_of_calls_[idx]
          << double(accumulated_time_[idx]) * 100.0 / total_time << "%"
          << std::endl;
    }
  }

private:
  Timer timer_;
  std::map<std::string, int> name2id_;
  std::vector<std::string> profile_names_;
  std::vector<TimeType> last_start_time_;
  std::vector<TimeType> accumulated_time_;
  std::vector<int> number_of_calls_;
};

template <class Timer, class TimeType>
inline std::ostream& operator<<(std::ostream& out, Profiler<Timer, TimeType>& profile) {
  profile.PrintProfileInfo(out);
  return out;
}

#endif // PROFILER_H
