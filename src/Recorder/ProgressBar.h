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

#pragma once
#include "indicators.h"

class ProgressBar {
 private:
  int totalIter_;
  int currentIter_;
  int currentProgress_;
  indicators::ProgressBar bar_;

 public:
  ProgressBar(int total, std::string name)
      : totalIter_(total), currentIter_(0), currentProgress_(0) {
    using namespace indicators;
    bar_.set_option(option::ShowElapsedTime(true));
    bar_.set_option(option::ShowRemainingTime(true));
    bar_.set_option(option::PostfixText(name));
  };

  void update() {
    ++currentIter_;
    while ((int)(100.0 * currentIter_ / totalIter_) > currentProgress_) {
      ++currentProgress_;
      bar_.tick();
    }
  }
};
