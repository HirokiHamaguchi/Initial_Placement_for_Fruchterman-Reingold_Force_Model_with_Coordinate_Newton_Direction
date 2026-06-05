#include <algorithm>
#include <cassert>
#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "../../../src/cpp/solve/solve.cpp"

struct StanfordGraph
{
  std::string stem;
  std::string displayName;
  size_t n;
  size_t m;
  std::vector<size_t> row;
  std::vector<size_t> col;
  std::vector<double> data;
};

static std::string trim(std::string s)
{
  auto isSpace = [](unsigned char c)
  { return std::isspace(c) != 0; };
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), [&](unsigned char c)
                                  { return !isSpace(c); }));
  s.erase(std::find_if(s.rbegin(), s.rend(), [&](unsigned char c)
                       { return !isSpace(c); })
              .base(),
          s.end());
  return s;
}

static StanfordGraph readStanfordGraph(const std::string &relativePath)
{
  std::ifstream file(relativePath);
  if (!file.is_open())
  {
    std::cerr << "Error: file not found: " << relativePath << std::endl;
    std::exit(1);
  }

  std::string line;
  std::getline(file, line);
  std::getline(file, line);
  std::string displayName = trim(line);
  if (!displayName.empty() && displayName[0] == '#')
    displayName = trim(displayName.substr(1));

  std::getline(file, line);
  std::istringstream meta(line);
  std::string hashTag, nodesLabel, edgesLabel;
  size_t n = 0, m = 0;
  meta >> hashTag >> nodesLabel >> n >> edgesLabel >> m;
  assert(hashTag == "#");
  assert(nodesLabel == "Nodes:");
  assert(edgesLabel == "Edges:");

  std::getline(file, line);

  std::unordered_map<size_t, size_t> idToIndex;
  idToIndex.reserve(n * 2);
  std::vector<std::pair<size_t, size_t>> edges;
  edges.reserve(m);
  while (std::getline(file, line))
  {
    line = trim(line);
    if (line.empty() || line[0] == '#')
      continue;
    std::istringstream iss(line);
    size_t u = 0, v = 0;
    iss >> u >> v;
    if (u == 0 || v == 0)
      continue;
    auto [itU, insertedU] = idToIndex.emplace(u, idToIndex.size());
    auto [itV, insertedV] = idToIndex.emplace(v, idToIndex.size());
    if (u == v)
      continue;
    edges.emplace_back(std::minmax(itU->second, itV->second));
  }

  std::sort(edges.begin(), edges.end());
  edges.erase(std::unique(edges.begin(), edges.end()), edges.end());

  StanfordGraph graph;
  graph.stem = std::filesystem::path(relativePath).stem().string();
  graph.displayName = displayName;
  graph.n = idToIndex.size();
  graph.m = edges.size();
  graph.row.reserve(graph.m);
  graph.col.reserve(graph.m);
  graph.data.reserve(graph.m);
  for (const auto &[u, v] : edges)
  {
    graph.row.push_back(u);
    graph.col.push_back(v);
    graph.data.push_back(1.0);
  }

  return graph;
}

static Eigen::VectorXd makeInitialPosition(size_t n, int seed)
{
  Grid grid(static_cast<int>(n), seed);
  return grid.toPosition();
}

int main()
{
  std::filesystem::create_directories("doc/main/large/out");

  const int seed = 1;
  const std::vector<std::string> inputFiles = {
      "data/StanfordNetworkAnalysisProject/com-amazon.ungraph.txt",
      "data/StanfordNetworkAnalysisProject/com-dblp.ungraph.txt",
      "data/StanfordNetworkAnalysisProject/roadNet-PA.txt",
  };

  for (const auto &inputFile : inputFiles)
  {
    StanfordGraph graph = readStanfordGraph(inputFile);

    std::vector<size_t> row = graph.row;
    std::vector<size_t> col = graph.col;
    std::vector<double> data = graph.data;
    Problem problem;
    problem.n = graph.n;
    problem.m = graph.m;
    problem.k = 1.0 / std::sqrt(graph.n);
    problem.row = std::move(row);
    problem.col = std::move(col);
    problem.data = std::move(data);
    problem.matrixName = graph.stem;
    problem.modelType = ForceModelType::FR;
    problem.makeAdj();
    problem.createModel();

    Eigen::VectorXd initPos = makeInitialPosition(problem.n, seed);
    std::vector<Eigen::VectorXd> finalPositions;
    solve_init(problem, true, seed, finalPositions);
    assert(!finalPositions.empty());

    std::vector<Eigen::VectorXd> positions = {initPos, finalPositions.back()};
    const std::string outPath = "doc/main/large/out/" + graph.stem + "_FR.out";
    problem.printOutput(positions, outPath);

    std::cout << graph.stem << "\n";
    std::cout << "  displayName: " << graph.displayName << "\n";
    std::cout << "  nodes: " << graph.n << "\n";
    std::cout << "  edges: " << graph.m << "\n";
    std::cout << "  output: " << outPath << "\n";
  }

  return 0;
}
