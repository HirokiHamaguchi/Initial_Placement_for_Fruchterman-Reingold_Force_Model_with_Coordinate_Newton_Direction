#pragma once

#include <cassert>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>

std::pair<std::string, std::ofstream> openFile(std::string path) {
  std::ofstream file(path);
  if (!file.is_open()) {
    std::cerr << "Error: file not found on \"openFile\" (path: " << path << ")"
              << std::endl;
    exit(1);
  }
  return {path, std::move(file)};
}
