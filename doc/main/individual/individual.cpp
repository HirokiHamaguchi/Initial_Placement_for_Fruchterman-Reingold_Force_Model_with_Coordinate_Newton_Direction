#include "../../../src/cpp/solve/solve.cpp"
#include "../../../src/cpp/util/problem.hpp"

#include <iostream>
#include <string>
#include <vector>

namespace
{
  const int NUM_SEEDS = 3;
  const int MAX_ITER = 200;   // Default max iterations
  const int MAX_ITER_2 = 500; // Compute the optimal solution with the seed 1
  const int OPT_SEED = 1;

  std::string getModelStr(ForceModelType model)
  {
    switch (model)
    {
    case ForceModelType::FR:
      return "FR";
    case ForceModelType::HC:
      return "HC";
    case ForceModelType::Eades:
      return "Eades";
    default:
      return "Unknown";
    }
  }

  std::vector<Method> getMethods(ForceModelType model)
  {
    if (model == ForceModelType::FR)
    {
      return {
          Method::FR,
          Method::CN_FR,
          Method::L_BFGS,
          Method::CN_L_BFGS,
      };
    }

    return {
        Method::L_BFGS,
        Method::CN_L_BFGS,
    };
  }

  struct ExperimentConfig
  {
    std::vector<std::string> matrixNames;
    ForceModelType model;
    std::string suffix;   // output file suffix
    std::string histName; // history file name
  };

  void appendHistory(
      std::string &histStr,
      const std::vector<std::pair<double, double>> &hist)
  {
    histStr += std::to_string(hist.size()) + "\n";

    for (auto [score, _] : hist)
      histStr += std::to_string(score) + " ";
    histStr += "\n";

    for (auto [_, time] : hist)
      histStr += std::to_string(time) + " ";
    histStr += "\n";

    auto [score, time] = hist.back();

    histStr += "Elapsed_time: " + std::to_string(time) + "\n";
    histStr += "Score: " + std::to_string(score) + "\n";

    std::cout << "Elapsed time: " << time << " seconds" << std::endl;
    std::cout << "Score: " << score << std::endl;
  }

  void runExperiment(const ExperimentConfig &config)
  {
    const auto modelStr = getModelStr(config.model);
    const auto methods = getMethods(config.model);

    std::cout << "Model: " << modelStr << std::endl;

    std::string histStr =
        std::to_string(config.matrixNames.size()) + " " +
        std::to_string(methods.size()) + "\n";

    for (const auto &matrixName : config.matrixNames)
    {
      Problem problem(matrixName, config.model, modelStr == "HC" || modelStr == "Eades");

      histStr += matrixName + "\n";
      std::cout << matrixName << std::endl;

      for (const auto &method : methods)
      {
        histStr += MethodStr[method] + "\n";
        std::cout << MethodStr[method] << std::endl;

        for (int seed = 1; seed <= NUM_SEEDS; seed++)
        {
          auto [hist, positions] = solve(method, problem, true, seed, MAX_ITER);

          appendHistory(histStr, hist);

          if (seed == OPT_SEED)
          {
            problem.printOutput(
                positions,
                "doc/main/individual/out/" +
                    problem.matrixName + "_" +
                    MethodStr[method] +
                    config.suffix +
                    ".out");
          }
        }

        if (method != Method::CN_L_BFGS)
          continue;

        auto [_, optPositions] = solve(method, problem, true, OPT_SEED, MAX_ITER_2);

        problem.printOutput(
            optPositions,
            "doc/main/individual/out/opt_" +
                problem.matrixName +
                config.suffix +
                ".out");
      }
    }

    const std::string fileName = "doc/main/individual/" + config.histName;
    auto [histPath, fileForHist] = openFile(fileName);

    fileForHist << histStr;
    fileForHist.close();

    std::cout << "Hist path: " << histPath << std::endl;
  }
} // namespace

int main()
{
  const std::vector<std::string> matrixNames1 = {
      "cycle300",
      "jagmesh1",
      "btree9",
      "1138_bus",
      "dwt_1005",
      "dwt_2680",
      "3elt",
  };

  const std::vector<std::string> matrixNames2 = {
      "jagmesh8",    // 1141
      "bcsstk14",    // 1806 1766
      "bcsstk15",    // 3948 3942
      "bcsstk16",    // 4884 4810
      "USpowerGrid", // 4941
      "bcspwr10",    // 5300
      "wiki-Vote",   // 8297 7066
  };

  const std::vector<ExperimentConfig> experiments = {
      {
          matrixNames1,
          ForceModelType::FR,
          "_FR",
          "hist_FR.txt",
      },
      {
          matrixNames1,
          ForceModelType::Eades,
          "_Eades",
          "hist_Eades.txt",
      },
      {
          matrixNames1,
          ForceModelType::HC,
          "_HC",
          "hist_HC.txt",
      },
      {
          matrixNames2,
          ForceModelType::FR,
          "_2",
          "hist_2.txt",
      },
  };

  for (const auto &experiment : experiments)
    runExperiment(experiment);

  return 0;
}
