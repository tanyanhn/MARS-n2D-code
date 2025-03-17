#pragma once

#include <any>
#include <map>
#include <nlohmann/json.hpp>

#include "InterfaceTracking/VectorFunction.h"

// 参数包装器
class ParamWrapper {
  std::any data_;

 public:
  template <typename... Args>
  ParamWrapper(Args&&... args)
      : data_(std::make_tuple(std::forward<Args>(args)...)) {}

  template <typename... Args>
  auto unwrap() const {
    return std::any_cast<std::tuple<Args...>>(data_);
  }
};

class Vector2DFactory {
  using CreatorFunc = std::function<VectorFunction<2>*(const ParamWrapper&)>;
  std::map<std::string, CreatorFunc> creators_;

 public:
  static Vector2DFactory& getInstance() noexcept {
    static Vector2DFactory instance;
    return instance;
  }

  template <typename T, typename... Args>
  bool registerClass() {
    creators_[T::staticClassName()] = [](const ParamWrapper& params) noexcept {
      return std::apply(
          [](Args... args) { return new T(std::forward<Args>(args)...); },
          params.unwrap<Args...>());
    };
    return true;
  }

  VectorFunction<2>* create(const std::string& name,
                            const ParamWrapper& params) const {
    auto it = creators_.find(name);
    return (it != creators_.end()) ? it->second(params) : nullptr;
  }
};

// JSON参数转换器
ParamWrapper parseParams(const nlohmann::json& j) {
  // 实际项目需要根据类型信息做更完善的解析
  if (j.is_array()) {
    if (j.size() == 2 && j[0].is_number_float() && j[1].is_number()) {
      return ParamWrapper(j[0].get<double>(), j[1].get<int>());
    }
    if (j.size() == 1 && j[0].is_number_float()) {
      return ParamWrapper(j[0].get<double>());
    }
  }
  return ParamWrapper();  // 默认空参数
}

// 自动注册宏
#define REGISTER_CLASS(CLASS, ...) \
  inline bool CLASS##_registered = \
      Vector2DFactory::getInstance().registerClass<CLASS, __VA_ARGS__>();