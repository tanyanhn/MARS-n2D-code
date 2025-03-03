
#include <spdlog/spdlog.h>

#include <string>

#pragma once

enum DebugLevel {
  TRACE = 0,
  DEBUG = 1,
  INFO = 2,
  WARNING = 3,
  ERROR = 4,
  CRITICAL = 5,
  OFF = 6
};

extern DebugLevel globalDebugLevel;
extern int globalNumIndent;

inline std::string getIndentedString(const std::string& str) {
  std::string indent(globalNumIndent * 2, ' ');
  return indent + str;
}

inline void increaseIndent() { globalNumIndent += 1; }
inline void decreaseIndent() { globalNumIndent -= 1; }

inline void TraceMessage(const std::string& trace) {
  spdlog::trace(getIndentedString(trace));
}
inline void DebugMessage(const std::string& debug) {
  spdlog::debug(getIndentedString(debug));
}

inline void DebugMessage(const std::string& tag, const std::string& debug) {
  std::string message = std::string("[") + tag + "] " + debug;
  spdlog::debug(getIndentedString(message));
}
inline void InfoMessage(const std::string& info) {
  spdlog::info(getIndentedString(info));
}

inline void InfoMessage(const std::string& tag, const std::string& info) {
  std::string message = std::string("[") + tag + "] " + info;
  spdlog::info(getIndentedString(message));
}
inline void WarningMessage(const std::string& warning) {
  spdlog::warn(getIndentedString(warning));
}
inline void ErrorMessage(const std::string& error) {
  spdlog::error(getIndentedString(error));
}

inline void CriticalMessage(const std::string& critical) {
  spdlog::critical(getIndentedString(critical));
}

inline void setDebugLevel(DebugLevel level) {
  InfoMessage("[Global] Debug level set to " + std::to_string(level));
  spdlog::set_level(static_cast<spdlog::level::level_enum>(level));
  globalDebugLevel = level;
}

inline void Message(const std::string& message, DebugLevel level) {
  auto logger = spdlog::default_logger();
  logger->log(static_cast<spdlog::level::level_enum>(level),
              getIndentedString(message));
}

inline void Message(const std::string& message, int level) {
  auto logger = spdlog::default_logger();
  logger->log(static_cast<spdlog::level::level_enum>(level),
              getIndentedString(message));
}

inline bool checkDebugLevel(DebugLevel level) {
  return level >= globalDebugLevel;
}