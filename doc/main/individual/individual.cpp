#include "../../../src/cpp/solve/solve.cpp"
#include "../../../src/cpp/util/problem.hpp"

#include <iostream>
#include <string>
#include <vector>

namespace
{
  int getIterations(ForceModelType model)
  {
    switch (model)
    {
    case ForceModelType::HC:
      return 300;
    case ForceModelType::Eades:
      return 300;
    default:
      return 150;
    }
  }

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
    int num_seeds;
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
    const auto MAX_ITER = getIterations(config.model);
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

        for (int seed = 1; seed <= config.num_seeds; seed++)
        {
          auto [hist, positions] = solve(method, problem, true, seed, MAX_ITER);

          appendHistory(histStr, hist);

          if (seed == 1)
          {
            std::string outputFile = "doc/main/individual/out/" + problem.matrixName + "_" + MethodStr[method] + config.suffix + ".out";
            problem.printOutput(positions, outputFile);
          }
        }
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
      // "dwt_2680",
      // "3elt",
  };

  const std::vector<std::string> matrixNames2 = {
      "cylinder_30_30", //   900
      "gr_30_30",       //   900
      "sierpinski_06",  //  1095
      "jagmesh8",       //  1141
      "USpowerGrid",    //  4941
      "wiki-Vote",      //  7066 (original: 8297)
      "crack",          // 10240
  };

  const std::vector<ExperimentConfig> experiments = {
      // {
      //     matrixNames1,
      //     ForceModelType::FR,
      //     "_FR",
      //     "hist_FR.txt",
      //     5,
      // },
      {
          matrixNames1,
          ForceModelType::Eades,
          "_Eades",
          "hist_Eades.txt",
          1,
      },
      {
          matrixNames1,
          ForceModelType::HC,
          "_HC",
          "hist_HC.txt",
          1,
      },
      // {
      //     matrixNames2,
      //     ForceModelType::FR,
      //     "_2",
      //     "hist_2.txt",
      //     5,
      // },
  };

  for (const auto &experiment : experiments)
    runExperiment(experiment);

  return 0;
}
