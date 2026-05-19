#pragma once

#include <cassert>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>

std::pair<std::string, std::ofstream> openFile(std::string path)
{
  assert('A' <= path[0] && path[0] <= 'z'); // relative path
  std::string curPath = std::filesystem::current_path().string();
  size_t folderPathSize = curPath.find("Initial_Placement_for_FR");
  assert(folderPathSize != std::string::npos);
  curPath = curPath.substr(
      0, folderPathSize + std::string("Initial_Placement_for_FR").size());
  std::ofstream file(curPath + "/" + path);
  if (!file.is_open())
  {
    std::cerr << "Error: file not found on \"openFile\" (path: " << path << ")"
              << std::endl;
    exit(1);
  }
  return {path, std::move(file)};
}
