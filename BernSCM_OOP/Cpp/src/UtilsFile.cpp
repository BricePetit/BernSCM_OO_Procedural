/**
 * Brice Petit 000439675
 * Master Thesis
 * Transformation of BernSCM using UML and OOP
 * CPP version
 */

#include "UtilsFile.hpp"

/**
 * Default constructor of the util class.
 *
 * @param ID    Additional simulation identifier.
 */
UtilsFile::UtilsFile(string ID) : _id(move(ID)) {}

/**
 * Destructor of the class UtilsFile. We use the default constructor because there is no
 * complex attribute to destroy.
 */
UtilsFile::~UtilsFile() = default;

/**
 * Get the complete path and go to the folder corresponding to the folder name.
 *
 * @param folder_name   The name of the folder where we want to go.
 * @return              Return the path.
 */
string UtilsFile::getPath(const string &folder_name) {
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
 * Compute the name of the file.
 *
 * @param t2x       Climate sensitivity.
 * @param t_dep     Temperature dependence.
 * @param co2_dep   Co2 dependence.
 * @param scenario  The name of the scenario file without forcing_ (at the beginning) and .dat (at the end).
 * @param dt        Time step.
 * @return          Return the name of the file where we want to write the solution.
 */
string UtilsFile::getFileName(const double &t2x, const bool &t_dep, const bool &co2_dep, const string &scenario,
                              const double &dt) {
    string x_imp, x_lin, x_sens;
    char xdt[16], temp[16];
    int vork, nachk;
    if (IMPLICIT_O) {
        if (IMPLICIT_L) {
            x_imp = "I";    // Implicit step
        } else {
            x_imp = "Io";   // Ocean-only implicit step
        }
    }
    if (LINEAR) {
        x_lin = "Q";
    }
    vork = max(1, int(log10(1e-5 + dt)) + 1);
    nachk = int(-floor(log10(1e-5 + dt - int(dt))));
    if (nachk > 4) {
        sprintf(xdt, "%*d", vork, int(dt));
    } else {
        sprintf(xdt, "%*.*f", vork + nachk + 1, nachk, dt);
    }

    if (t_dep) {
        x_sens = "_t";
    } else {
        x_sens = "_t0";
    }

    if (co2_dep) {
        x_sens += "_f";
    } else {
        x_sens += "_f0";
    }

    vork = int(log10(1e-1 + t2x * 10e0)) + 1;

    sprintf(temp, "%*d", vork, int(t2x * 10 + 0.5e0));
    x_sens += "_CS";
    x_sens += temp;
    if (!this->_id.empty()) {
        this->_id = "_" + this->_id + "_";
    } else {
        this->_id = "_";
    }
    return scenario + "_D" + xdt + x_lin + x_imp + this->_id + "BernSCM" + x_sens + ".dat";
}

/**
 * Read forcing file.
 *
 * @param forcing   The matrix of all forcing values.
 * @param scenario  The name of the scenario file without forcing_ (at the beginning) and .dat (at the end).
 * @return          Return all forcing values.
 */
double **UtilsFile::readForcing(double **forcing, const string &scenario) {
    string line, tmp;
    double frecord[NFORC];
    int i = 0, j, k, l;
    ifstream my_file;
    my_file.open(this->getPath("forcing") + "forcing_" + scenario + ".dat");
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
    return forcing;
}

/**
 * Function to write the output of our simulation.
 *
 * @param earth Object Earth containing the atmosphere, land and ocean model.
 */
void UtilsFile::output(Earth &earth) {
    int i, tmp_index;
    char tmp_str[32];
    double *fa{}, *frh{}, *npp_out{}, *land_c_uptake{}, *ocean_c_uptake{}, *fossil_emissions{}, *f_deep{}, *midtime{};
    string desc_file, output_path, file_name;
    ofstream my_file;
    // Create the time to the format d/m/Y H:M:S
    time_t now = time(0);
    tm *ltm = localtime(&now);
    string time = to_string(ltm->tm_mday) + "/" + to_string(ltm->tm_mon) + "/" + to_string(ltm->tm_year) + " " +
                  to_string(ltm->tm_hour) + ":" + to_string(ltm->tm_min) + ":" + to_string(ltm->tm_sec);
    for (i = 1; i < earth.getNtime(); i++) {
        if (earth.getRfElem(i) == NA || earth.getRfncElem(i) == NA || earth.getRfcElem(i) == NA) {
            earth.setRfbElem(i, NA);
        } else {
            earth.setRfbElem(i, earth.getRfElem(i) - earth.getRfncElem(i) - earth.getRfcElem(i));
        }
    }
    fa = new double[earth.getNtime()]{0};
    frh = new double[earth.getNtime()]{earth.getNpp0()};
    npp_out = new double[earth.getNtime()]{earth.getNpp0()};
    land_c_uptake = new double[earth.getNtime()]{0};
    ocean_c_uptake = new double[earth.getNtime()]{0};
    fossil_emissions = new double[earth.getNtime()]{0};
    f_deep = new double[earth.getNtime()]{0};
    midtime = new double[earth.getNtime()]{0};
    // C budget as in code
    // The following fluxes correspond to midtime (box-centered).
    for (i = 1; i < earth.getNtime(); i++) {
        fa[i] = (earth.getMaElem(i) - earth.getMaElem(i - 1)) / earth.getDt();
        fossil_emissions[i] = (earth.getECo2Elem(i) + earth.getECo2Elem(i - 1)) / 2e0;
        if (LINEAR) {
            ocean_c_uptake[i] = (earth.getFoElem(i) + earth.getFoElem(i - 1)) / 2e0;
            npp_out[i] = (earth.getFnppElem(i) + earth.getFnppElem(i - 1)) / 2e0;
        } else {
            ocean_c_uptake[i] = earth.getFoElem(i);
        }
        land_c_uptake[i] = (earth.getMlElem(i) - earth.getMlElem(i - 1)) / earth.getDt();
        frh[i] = earth.getFnppElem(i) - land_c_uptake[i];
        earth.setFbElem(i, fossil_emissions[i] - fa[i] - ocean_c_uptake[i] - land_c_uptake[i]);
    }
    if (!LINEAR) {
        npp_out = earth.getFnpp();
    }
    if (this->_do_interpol) {
        for (i = 0; i < earth.getNtime() - 1; i++) {
            fa[i] = (fa[i] + fa[i + 1]) / 2;
            ocean_c_uptake[i] = (ocean_c_uptake[i] + ocean_c_uptake[i + 1]) / 2;
            land_c_uptake[i] = (land_c_uptake[i] + land_c_uptake[i + 1]) / 2;
            fossil_emissions[i] = (fossil_emissions[i] + fossil_emissions[i + 1]) / 2;
            earth.setFbElem(i, fossil_emissions[i] - fa[i] - land_c_uptake[i] - ocean_c_uptake[i]);
            frh[i] = land_c_uptake[i] - npp_out[i];
        }
        tmp_index = earth.getNtime() - 1;
        earth.setFbElem(tmp_index, fossil_emissions[tmp_index] - fa[tmp_index] - land_c_uptake[tmp_index] -
                                   ocean_c_uptake[tmp_index]);
        frh[tmp_index] = land_c_uptake[tmp_index] - npp_out[tmp_index];
    } else {
        for (i = 0; i < earth.getNtime(); i++) {
            midtime[i] = earth.getTimeElem(i) - earth.getDt() / 2;
        }
    }
    desc_file = "# Creation date: " + time \
                    + "\n" \
                      "#\n" \
                      "# Bern Simple Climate Model (BernSCM) version 1.0\n" \
                      "# = = = = = = = = = = = = = = = = = = = = = = = =\n" \
                      "#\n" \
                      "# Authors:\n" \
                      "# Kuno Strassmann, Fortunat Joos\n" \
                      "# Climate and Environmental Physics\n" \
                      "# Sidlerstr 5\n" \
                      "# CH-3012 Bern\n" \
                      "# kuno.strassmann@alumni.ethz.ch, joos@climate.unibe.ch\n" \
                      "# Tel.: 0041-31-631 44 61\n" \
                      "#\n" \
                      "# Description:\n" \
                      "#\n" \
                      "# The model is based on impulse response representation (Joos et al., 1996) of the ocean (surface-to-deep\n" \
                      "# transport) and the terrestrial biosphere (accumulation and decay of NPP carbon). Nonlinear ocean surface carbon\n" \
                      "# chemistry and NPP changes are treated by separate parametrizations (see description of C cycle components below).\n" \
                      "#\n" \
                      "# References:\n" \
                      "#\n" \
                      "# - K. Strassmann and F. Joos, 2017: The Bern Simple Climate Model: an extensible and fully\n" \
                      "# documented open source reimplementation of the Bern reduced form model for global carbon\n" \
                      "# cycle-climate simulations. Submitted to Geophysical Model Development.\n" \
                      "# - F. Joos, M. Bruno, R. Fink, T. F. Stocker, U. Siegenthaler, C. Le Quere, and\n" \
                      "# J. L. Sarmiento, 1996: An efficient and accurate representation of complex oceanic\n" \
                      "# and biospheric models of anthropogenic carbon uptake. Tellus, 48B:397-417\n" \
                      "# - Joos, F. and M. Bruno, 1996: Pulse response functions are cost-efficient tools to model the link\n" \
                      "# between carbon emissions, atmospheric CO2 and global warming. Physics and Chemistry of the Earth\n" \
                      "# 21:471-476.\n" \
                      "#\n" \
                      "#\n" \
                      "#\n" \
                      "# Numerical solution\n" \
                      "# = = = = = = = = =\n" \
                      "# Time step: " + to_string(earth.getDt()) + "yr\n";
    if (IMPLICIT_L || IMPLICIT_O) {
        desc_file += "# Implicite step: \n";
        if (IMPLICIT_L) {
            desc_file += "# - Land C exchange\n";
        }
        if (IMPLICIT_O) {
            desc_file += "# - Ocean C and heat exchange\n";
        }
    }
    if (LINEAR) {
        desc_file += "# Discretization: piecewise linear\n";
    } else {
        desc_file += "# Discretization: piecewise constant\n";
    }
    desc_file += "#\n" \
                 "#\n" \
                 "#\n" \
                 "# Simulation setup\n" \
                 "# = = = = = = = = =\n" \
                 "#\n" \
                 "# Forcing scenario: " \
                 + earth.getScenario() \
                 + "\n" \
                   "#\n" \
                   "# Climate Sensitivity: " \
                 + to_string(earth.getT2x()) \
                 + " degrees C per doubling of atm. CO2 \n" \
                   "#\n" \
                   "# Carbon Cycle:\n" \
                   "# Process sensitivity to atmospheric CO2:\n" \
                   "# - Ocean CO2 uptake (see ocean component)\n";
    if (earth.isCo2Dep()) {
        desc_file += "# - Land C exchange (see land component)\n";
    }
    desc_file += "# Process sensitivity to global mean SAT:\n";
    if (earth.isTDep()) {
        desc_file += "# - Ocean CO2 uptake (see ocean component)\n";
        if (LAND_T_DEP) {
            desc_file += "# - Land C exchange (see land component)\n";
        }
        desc_file += "#\n#\n";
    } else {
        desc_file += "# none\n#\n#\n";
    }
    desc_file += earth.getOceanDoc();
    desc_file += earth.getLandDoc();
    desc_file += "#time_(yr) glob_temp_dev_(°C) RF_tot_(W/m²) RF_CO2_(W/m²) RF_nonCO2_(W/m²) RF_budget_(W/m²) " \
                 "ocean_heat_uptake_(PW) co2_atm_(ppm) co2_seasurf_(ppm) atm_CO2_change_(GtC/yr) " \
                 "fossil_CO2_em_(GtC/yr) budget_C_uptake_(GtC/yr) ocean_C_uptake_(GtC/yr) land_C_uptake_(GtC/yr) " \
                 "NPP_(GtC/yr) RH_(GtC/yr) LandC_(GtC) dDIC_(µmol/kg) fdeep_(GtC/yr)";
    if (!this->_do_interpol) {
        desc_file += " midtime_(yr)";
    }
    desc_file += "\n#time glob_temp_dev RF_tot RF_CO2 RF_nonCO2 RF_budget ocean_heat_uptake co2_atm co2_seasurf " \
                 "atm_CO2_change fossil_CO2_em budget_C_uptake ocean_C_uptake land_C_uptake NPP RH LandC dDIC fdeep";
    if (!this->_do_interpol) {
        desc_file += " midtime";
    }
    desc_file += "\n";
    for (i = 0; i < earth.getNtime(); i++) {
        sprintf(tmp_str, "%20.10f", earth.getTimeElem(i));
        desc_file += tmp_str;
        memset(tmp_str, 0, sizeof tmp_str);
        sprintf(tmp_str, "%20.10f", earth.getTempElem(i));
        desc_file += tmp_str;
        memset(tmp_str, 0, sizeof tmp_str);
        sprintf(tmp_str, "%20.10f", earth.getRfElem(i));
        desc_file += tmp_str;
        memset(tmp_str, 0, sizeof tmp_str);
        sprintf(tmp_str, "%20.10f", earth.getRfcElem(i));
        desc_file += tmp_str;
        memset(tmp_str, 0, sizeof tmp_str);
        sprintf(tmp_str, "%20.10f", earth.getRfncElem(i));
        desc_file += tmp_str;
        memset(tmp_str, 0, sizeof tmp_str);
        sprintf(tmp_str, "%20.10f", earth.getRfbElem(i));
        desc_file += tmp_str;
        memset(tmp_str, 0, sizeof tmp_str);
        sprintf(tmp_str, "%20.10f", earth.getFhElem(i));
        desc_file += tmp_str;
        memset(tmp_str, 0, sizeof tmp_str);
        sprintf(tmp_str, "%20.10f", earth.getMaElem(i) / PPMTOGT);
        desc_file += tmp_str;
        memset(tmp_str, 0, sizeof tmp_str);
        sprintf(tmp_str, "%20.10f", earth.getCo2Atm0() + earth.getDpcsElem(i));
        desc_file += tmp_str;
        memset(tmp_str, 0, sizeof tmp_str);
        sprintf(tmp_str, "%20.10f", fa[i]);
        desc_file += tmp_str;
        memset(tmp_str, 0, sizeof tmp_str);
        sprintf(tmp_str, "%20.10f", fossil_emissions[i]);
        desc_file += tmp_str;
        memset(tmp_str, 0, sizeof tmp_str);
        sprintf(tmp_str, "%20.10f", earth.getFbElem(i));
        desc_file += tmp_str;
        memset(tmp_str, 0, sizeof tmp_str);
        sprintf(tmp_str, "%20.10f", ocean_c_uptake[i]);
        desc_file += tmp_str;
        memset(tmp_str, 0, sizeof tmp_str);
        sprintf(tmp_str, "%20.10f", land_c_uptake[i]);
        desc_file += tmp_str;
        memset(tmp_str, 0, sizeof tmp_str);
        sprintf(tmp_str, "%20.10f", npp_out[i]);
        desc_file += tmp_str;
        memset(tmp_str, 0, sizeof tmp_str);
        sprintf(tmp_str, "%20.10f", frh[i]);
        desc_file += tmp_str;
        memset(tmp_str, 0, sizeof tmp_str);
        sprintf(tmp_str, "%20.10f", earth.getMlElem(i));
        desc_file += tmp_str;
        memset(tmp_str, 0, sizeof tmp_str);
        sprintf(tmp_str, "%20.10f", earth.getMsElem(i));
        desc_file += tmp_str;
        memset(tmp_str, 0, sizeof tmp_str);
        sprintf(tmp_str, "%20.10f", f_deep[i]);
        desc_file += tmp_str;
        memset(tmp_str, 0, sizeof tmp_str);
        if (!this->_do_interpol) {
            sprintf(tmp_str, "%20.10f", midtime[i]);
            desc_file += tmp_str;
            memset(tmp_str, 0, sizeof tmp_str);
        }
        desc_file += "\n";
    }
    output_path = this->getPath("output");
    file_name = this->getFileName(earth.getT2x(), earth.isTDep(), earth.isCo2Dep(), earth.getScenario(), earth.getDt());
    output_path.append(file_name);
    my_file.open(output_path, fstream::out);
    my_file << desc_file;
    my_file.close();

    delete[] fa;
    delete[] frh;
    delete[] npp_out;
    delete[] land_c_uptake;
    delete[] ocean_c_uptake;
    delete[] fossil_emissions;
    delete[] f_deep;
    delete[] midtime;
}