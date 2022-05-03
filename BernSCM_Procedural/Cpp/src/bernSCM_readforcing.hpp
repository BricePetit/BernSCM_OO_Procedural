//
// Created by Brice Petit on 25-02-22.
//

#ifndef CPP_BERNSCM_READFORCING_HPP
#define CPP_BERNSCM_READFORCING_HPP

#include "Globals.hpp"

#include <filesystem>
#include <fstream>
#include <string>

using namespace std;
namespace fs = std::filesystem;

string getPath(const string &folder_name);

void readForcing();

#endif //CPP_BERNSCM_READFORCING_HPP
