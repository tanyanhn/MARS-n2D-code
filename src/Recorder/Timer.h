#pragma once
#include <chrono>
#include <iostream>
#include <map>
#include <string>

class Timer {
 public:
  // 构造函数，初始化计时器名称并记录开始时间
  explicit Timer(std::string name_)
      : name_(std::move(name_)), start_(std::chrono::steady_clock::now()) {
    if (!initialized_) return;
    auto iter = accumulatedTimes_.insert({this->name_, {0.0, 0}});
    if (iter.second) names_.push_back(name_);
  }

  // 析构函数，计算持续时间并累加到统计中
  ~Timer() {
    if (!initialized_) return;
    auto end = std::chrono::steady_clock::now();
    auto duration = end - start_;
    auto iter = accumulatedTimes_.find(this->name_);
    iter->second.first += (double)duration.count() / 1e9;
    iter->second.second++;
  }

  static void initial() { initialized_ = true; }
  // 静态方法输出所有计时器的累计时间
  static void printStatistics() {
    if (!initialized_) return;
    std::cout << "\n===== Timing Statistics =====\n";
    for (const auto& name : names_) {
      const auto& pair = accumulatedTimes_.at(name);
      const auto& seconds = pair.first;
      const auto& count = pair.second;
      std::cout << "[" << name << "]\n"
                << "  Total: " << seconds * 1000 << " ms\n"
                << "         " << count << " counts\n";
    }
    std::cout << "=============================\n";
    initialized_ = false;
  }

 private:
  std::string name_;
  std::chrono::steady_clock::time_point start_;

  // 使用C++17 inline静态成员存储累计时间
  inline static std::map<std::string, std::pair<double, int>> accumulatedTimes_;
  inline static std::vector<std::string> names_;
  inline static bool initialized_ = false;
};