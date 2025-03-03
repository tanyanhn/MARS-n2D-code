/*
MIT License

Copyright (c) 2024 Yixiao Qian

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <list>

#include "Recorder/CPUTimeHelper.h"
#include "Recorder/SpdlogWrapper.h"
#include "StageRecorder.h"

#pragma once

class Recorder {
 private:
  std::list<StageRecorder> stages_;
  bool owned_ = false;
  CPUTimeHelper timer_;
  double duration_;

 public:
  Recorder() {
    timer_.reset();
    owned_ = false;
  }

 public:
  void pushLogStage(const std::string& stageName,
                    DebugLevel level = DebugLevel::INFO);
  void popLogStage();
  void pushLogEvent(const std::string& eventName,
                    const std::string& className = "",
                    DebugLevel level = DebugLevel::INFO);
  void popLogEvent();
  void finish();
  void printSummary() const;
};

extern std::shared_ptr<Recorder> globalRecorder;
extern std::string exportDir;

struct RecorderInfo {
  DebugLevel level = DebugLevel::INFO;
  std::string exportDir = "";
};

void RecorderInitialize(const RecorderInfo info = RecorderInfo());
void RecorderFinalize();
void pushLogStage(const std::string& stageName,
                  DebugLevel level = DebugLevel::INFO);
void popLogStage();
void pushLogEvent(const std::string& eventName,
                  const std::string& className = "",
                  DebugLevel level = DebugLevel::INFO);
void popLogEvent();
void setExportDir(const std::string& dir);
std::string getExportDir();
