
#include "EventRecorder.h"

#include <iostream>

#include "Recorder/SpdlogWrapper.h"

void EventRecorder::finish() {
  decreaseIndent();
  duration_ = timer_();
  std::string message;
  if (className_ == "") {
    message = "[Event] Event [" + eventName_ + "] finished, spent " +
              std::to_string(duration_) + " seconds.";
  } else {
    message = "[Event] Event [" + className_ + "]::[" + eventName_ +
              "] finished, spent " + std::to_string(duration_) + " seconds.";
  }
  Message(message, level_);
}
void EventRecorder::printSummary() const {
  std::string message;
  if (className_ == "") {
    message =
        "Event [" + eventName_ + "] Duration: " + std::to_string(duration_);
  } else {
    message = "Event [" + className_ + "]::[" + eventName_ +
              "] Duration: " + std::to_string(duration_);
  }
  std::cout << "  " << message << std::endl;
}
