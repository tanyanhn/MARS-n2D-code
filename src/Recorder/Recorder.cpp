
#include "Recorder.h"

#include <iostream>
#include <memory>

#include "Recorder/SpdlogWrapper.h"

std::shared_ptr<Recorder> globalRecorder;
std::string exportDir;

void Recorder::pushLogStage(const std::string& stageName, DebugLevel level) {
  if (owned_) {
    ErrorMessage(" Recorder::pushLogStage(): Please pop the stage first.");
  }
  stages_.push_back(StageRecorder(stageName, level));
  owned_ = true;
}

void Recorder::popLogStage() {
  if (!owned_) {
    ErrorMessage(" Recorder::popLogStage(): Please push a stage first.");
  }
  stages_.back().finish();
  owned_ = false;
}

void Recorder::pushLogEvent(const std::string& eventName,
                            const std::string& className, DebugLevel level) {
  if (!owned_) {
    ErrorMessage(" Recorder::pushLogEvent(): No stage is pushed.");
  }
  stages_.back().pushLogEvent(eventName, className, level);
}

void Recorder::popLogEvent() {
  if (!owned_) {
    ErrorMessage(" Recorder::popLogEvent(): No stage is pushed.");
  }
  stages_.back().popLogEvent();
}

void Recorder::finish() { duration_ = timer_(); }

void Recorder::printSummary() const {
  std::string message = std::string("[Global] Total time: ") +
                        std::to_string(duration_) + " seconds";
  InfoMessage(message);
  InfoMessage("[Global] Printing Summary ...");
  for (const auto& stage : stages_) {
    stage.printSummary();
  }
}

void RecorderInitialize(const RecorderInfo info) {
  globalRecorder = std::make_shared<Recorder>();
  globalNumIndent = 0;
  setDebugLevel(info.level);
  setExportDir(info.exportDir);
  InfoMessage("[Global] Recorder initialized");
}

void RecorderFinalize() {
  InfoMessage("[Global] Recorder finalized");
  globalRecorder->finish();
  globalRecorder->printSummary();
}

void pushLogStage(const std::string& stageName, DebugLevel level) {
  globalRecorder->pushLogStage(stageName, level);
}

void popLogStage() { globalRecorder->popLogStage(); }

void pushLogEvent(const std::string& eventName, const std::string& className,
                  DebugLevel level) {
  globalRecorder->pushLogEvent(eventName, className, level);
}

void popLogEvent() { globalRecorder->popLogEvent(); }

void setExportDir(const std::string& dir) {
  if (dir == "") {
    return;
  }
  exportDir = dir;
  InfoMessage("Global", "ExportDir set to " + exportDir);
}

std::string getExportDir() { return exportDir; }
