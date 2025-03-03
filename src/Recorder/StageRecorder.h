

#include <stack>

#include "EventRecorder.h"
#include "Recorder/SpdlogWrapper.h"

class StageRecorder {
 private:
  std::string stageName_;
  double duration_;
  std::vector<EventRecorder> events_;
  std::stack<int> eventStack_;
  CPUTimeHelper timer_;
  DebugLevel level_;

 public:
  explicit StageRecorder(const std::string& stageName,
                         DebugLevel level = DebugLevel::INFO)
      : stageName_(stageName), duration_(0.0), level_(level) {
    std::string message =
        std::string("[Stage] Stage [") + stageName_ + "] started.";
    Message(message, level_);
    timer_.reset();
  }

 public:
  void pushLogEvent(const std::string& eventName,
                    const std::string& className = "",
                    DebugLevel level = DebugLevel::INFO);
  void popLogEvent();
  void finish();
  void printSummary() const;
};
