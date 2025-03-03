#include "StageRecorder.h"

#include <iostream>

#include "Recorder/SpdlogWrapper.h"

void StageRecorder::pushLogEvent(const std::string& eventName,
                                 const std::string& className,
                                 DebugLevel level) {
  increaseIndent();
  events_.emplace_back(eventName, className, level);
  eventStack_.push(events_.size() - 1);
}

void StageRecorder::popLogEvent() {
  if (eventStack_.size() == 0) {
    throw std::runtime_error(
        "StageRecorder::popLogEvent() called when event "
        "stack is empty");
  }
  int id = eventStack_.top();
  eventStack_.pop();
  events_.at(id).finish();
  decreaseIndent();
}

void StageRecorder::finish() {
  if (eventStack_.size() != 0) {
    throw std::runtime_error(
        "StageRecorder::finish() called when event "
        "stack is not empty");
  }
  duration_ = timer_();
  std::string message;
  message = "[Stage] Stage [" + stageName_ + "] finished, spent " +
            std::to_string(duration_) + " seconds";
  Message(message, level_);
}

void StageRecorder::printSummary() const {
  std::cout << "Stage [" << stageName_ << "] spent " << duration_
            << " seconds: " << std::endl;
  for (const auto& event : events_) {
    event.printSummary();
  }
}
