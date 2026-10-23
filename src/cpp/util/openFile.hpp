#pragma once

#include <cassert>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>

std::filesystem::path findProjectRoot()
{
  // Identify the checkout by its directories so renaming it preserves file paths.
  auto root = std::filesystem::current_path();
  while (!std::filesystem::is_directory(root / "data") ||
         !std::filesystem::is_directory(root / "src"))
  {
    if (root == root.parent_path())
      throw std::runtime_error("Project root with data/ and src/ not found");
    root = root.parent_path();
  }
  return root;
}

std::pair<std::string, std::ofstream> openFile(std::string path)
{
  assert('A' <= path[0] && path[0] <= 'z'); // relative path
  std::ofstream file(findProjectRoot() / path);
  if (!file.is_open())
  {
    std::cerr << "Error: file not found on \"openFile\" (path: " << path << ")"
              << std::endl;
    exit(1);
  }
  return {path, std::move(file)};
}
