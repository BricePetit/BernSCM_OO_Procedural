//
// Created by Brice Petit on 25-02-22.
//

#include "bernSCM_readforcing.hpp"

/**
 * Get the complete path and go to the folder corresponding to the folder name.
 *
 * @param folder_name   The name of the folder where we want to go.
 * @return              Return the path.
 */
string getPath(const string &folder_name) {
    string separator, path;
    path = fs::current_path().string();
    for (int i = 0; i < 3; i++) {
        path.pop_back();
    }
    path += folder_name;
#ifdef _WIN32
    separator = "\\";
#else
    separator = "/";
#endif
    return path + separator;
}

/**
 * Read forcing file.
 */
void readForcing() {
    string line, tmp;
    double frecord[NFORC];
    int i = 0, j, k, l;
    ifstream my_file;
    my_file.open(getPath("forcing") + "forcing_" + scenario + ".dat");
    // We count the number of lines in the file.
    while (getline(my_file, line)) {
        if (line[0] != '#' && line[1] != '#') {
            FORCING_ROWS++;
        }
    }
    my_file.clear();
    my_file.seekg(0);
    // Define the number of rows needed.
    forcing = new double *[FORCING_ROWS];
    // We create our forcing matrix.
    while (getline(my_file, line)) {
        if (line[0] != '#' && line[1] != '#') {
            l = 0;
            // We search each value in the line that we read.
            for (k = 0; k <= line.size(); k++) {
                if (line[k] != ' ' && k < line.size()) {
                    tmp += line[k];
                } else if ((!tmp.empty() && line[k] == ' ') || (k == line.size())) {
                    frecord[l] = stod(tmp);
                    l++;
                    tmp = "";
                }
            }
            forcing[i] = new double[NFORC];
            for (j = 0; j < NFORC; j++) {
                if (abs(frecord[j] - NA) < 1e-3) {
                    forcing[i][j] = NA;
                } else if (j == JACO2) {
                    forcing[i][j] = frecord[j] * PPMTOGT;
                } else {
                    forcing[i][j] = frecord[j];
                }
            }
            i++;
        }
    }
    my_file.close();
}
