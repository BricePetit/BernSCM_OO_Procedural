/**
 * Brice Petit 000439675
 * Master Thesis
 * Transformation of BernSCM using UML and OOP
 * CPP version
 */

#ifndef CPP_UTILSFILE_HPP
#define CPP_UTILSFILE_HPP

#include "Earth.hpp"
#include "Globals.hpp"

#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

using namespace std;
namespace fs = std::filesystem;

class UtilsFile {
public:
    explicit UtilsFile(string ID);

    ~UtilsFile();

    static string getPath(const string &folder_name);

    string getFileName(const double &t2x, const bool &t_dep, const bool &co2_dep, const string &scenario,
                       const double &dt);

    double** readForcing(double **forcing, const string &scenario);

    void output(Earth &earth);

private:
    string _id;
    const int _do_interpol = 1;
};


#endif //CPP_UTILSFILE_HPP
