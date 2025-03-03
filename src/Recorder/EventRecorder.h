
#include <string>

#include "CPUTimeHelper.h"
#include "SpdlogWrapper.h"

#pragma once

class EventRecorder {
 private:
  std::string className_;
  std::string eventName_;
  DebugLevel level_;
  double duration_;
  CPUTimeHelper timer_;

 public:
  explicit EventRecorder(const std::string& eventName,
                         const std::string& className = "",
                         DebugLevel level = DebugLevel::INFO)
      : className_(className), eventName_(eventName), level_(level) {
    timer_.reset();
    std::string message;
    if (className == "") {
      message = "[Event] Event [" + eventName_ + "] started.";
    } else {
      message =
          "[Event] Event [" + className + "]::[" + eventName_ + "] started.";
    }
    Message(message, level_);
    increaseIndent();
  }

 public:
  void finish();
  void printSummary() const;
};
