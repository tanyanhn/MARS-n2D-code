#pragma once
#include <chrono>
#include <iostream>
#include <map>
#include <string>

class Timer {
 public:
  // 构造函数，初始化计时器名称并记录开始时间
  explicit Timer(std::string name_)
      : name(std::move(name_)), start(std::chrono::steady_clock::now()) {
    auto iter = accumulatedTimes.insert({this->name, {0.0, 0}});
    if (iter.second) names.push_back(name);
  }

  // 析构函数，计算持续时间并累加到统计中
  ~Timer() {
    auto end = std::chrono::steady_clock::now();
    auto duration = end - start;
    auto iter = accumulatedTimes.find(this->name);
    iter->second.first += (double)duration.count() / 1e9;
    iter->second.second++;
  }

  // 静态方法输出所有计时器的累计时间
  static void printStatistics() {
    std::cout << "\n===== Timing Statistics =====\n";
    for (const auto& name : names) {
      const auto& pair = accumulatedTimes.at(name);
      const auto& seconds = pair.first;
      const auto& count = pair.second;
      std::cout << "[" << name << "]\n"
                << "  Total: " << seconds * 1000 << " ms\n"
                << "         " << count << " counts\n";
    }
    std::cout << "=============================\n";
  }

 private:
  std::string name;
  std::chrono::steady_clock::time_point start;

  // 使用C++17 inline静态成员存储累计时间
  inline static std::map<std::string, std::pair<double, int>> accumulatedTimes;
  inline static std::vector<std::string> names;
};