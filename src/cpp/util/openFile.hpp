#pragma once

#include <cassert>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>

std::pair<std::string, std::ofstream> openFile(std::string fileName) {
  std::string curPath = std::filesystem::current_path().string();
  assert(curPath.substr(curPath.size() - 8, 8) == "/src/cpp");
  std::string path = curPath.substr(0, curPath.size() - 8) + "/" + fileName;
  std::ofstream file(path);
  if (!file.is_open()) {
    std::cerr << "Error: file not found\n";
    exit(1);
  }
  return {path, std::move(file)};
}
