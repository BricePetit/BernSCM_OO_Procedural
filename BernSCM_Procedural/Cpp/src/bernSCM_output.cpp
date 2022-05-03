//
// Created by Brice Petit on 25-02-22.
//

#include "bernSCM_output.hpp"

/**
 * Compute the name of the file.
 *
 * @return  Return the name of the file where we want to write the solution.
 */
string getFileName() {
    string x_imp, x_lin, x_sens;
    char xdt[16], temp_char[16];
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
    vork = max(1, int(log10(1e-5 + DELTA_T)) + 1);
    nachk = int(-floor(log10(1e-5 + DELTA_T - int(DELTA_T))));
    if (nachk > 4) {
        sprintf(xdt, "%*d", vork, int(DELTA_T));
    } else {
        sprintf(xdt, "%*.*f", vork + nachk + 1, nachk, DELTA_T);
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

    sprintf(temp_char, "%*d", vork, int(t2x * 10 + 0.5e0));
    x_sens += "_CS";
    x_sens += temp_char;
    if (!ID.empty()) {
        ID = "_" + ID + "_";
    } else {
        ID = "_";
    }
    return scenario + "_D" + xdt + x_lin + x_imp + ID + "BernSCM" + x_sens + ".dat";
}

/**
 * Function to write the output of our simulation.
 */
void output() {
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
    for (i = 1; i < ntime; i++) {
        if (rf[i] == NA || rfnc[i] == NA || rfc[i] == NA) {
            rfb[i] = NA;
        } else {
            rfb[i] = rf[i] - rfnc[i] - rfc[i];
        }
    }
    fa = new double[ntime]{0};
    frh = new double[ntime]{NPP0};
    npp_out = new double[ntime]{NPP0};
    land_c_uptake = new double[ntime]{0};
    ocean_c_uptake = new double[ntime]{0};
    fossil_emissions = new double[ntime]{0};
    f_deep = new double[ntime]{0};
    midtime = new double[ntime]{0};
    // C budget as in code
    // The following fluxes correspond to midtime (box-centered).
    for (i = 1; i < ntime; i++) {
        fa[i] = (ma[i] - ma[i - 1]) / DELTA_T;
        fossil_emissions[i] = (e_co2[i] + e_co2[i - 1]) / 2e0;
        if (LINEAR) {
            ocean_c_uptake[i] = (fo[i] + fo[i - 1]) / 2e0;
            npp_out[i] = (fnpp[i] + fnpp[i - 1]) / 2e0;
        } else {
            ocean_c_uptake[i] = fo[i];
        }
        land_c_uptake[i] = (ml[i] - ml[i - 1]) / DELTA_T;
        frh[i] = fnpp[i] - land_c_uptake[i];
        fb[i] =  fossil_emissions[i] - fa[i] - ocean_c_uptake[i] - land_c_uptake[i];
    }
    if (!LINEAR) {
        npp_out = fnpp;
    }
    if (DO_INTERPOL) {
        for (i = 0; i < ntime - 1; i++) {
            fa[i] = (fa[i] + fa[i + 1]) / 2;
            ocean_c_uptake[i] = (ocean_c_uptake[i] + ocean_c_uptake[i + 1]) / 2;
            land_c_uptake[i] = (land_c_uptake[i] + land_c_uptake[i + 1]) / 2;
            fossil_emissions[i] = (fossil_emissions[i] + fossil_emissions[i + 1]) / 2;
            fb[i] =  fossil_emissions[i] - fa[i] - land_c_uptake[i] - ocean_c_uptake[i];
            frh[i] = land_c_uptake[i] - npp_out[i];
        }
        tmp_index = ntime - 1;
        fb[tmp_index] =  fossil_emissions[tmp_index] - fa[tmp_index] - land_c_uptake[tmp_index] -
                            ocean_c_uptake[tmp_index];
        frh[tmp_index] = land_c_uptake[tmp_index] - npp_out[tmp_index];
    } else {
        for (i = 0; i < ntime; i++) {
            midtime[i] = _time[i] - DELTA_T / 2;
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
                      "# Time step: " + to_string(DELTA_T) + "yr\n";
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
                 + scenario \
                 + "\n" \
                   "#\n" \
                   "# Climate Sensitivity: " \
                 + to_string(t2x) \
                 + " degrees C per doubling of atm. CO2 \n" \
                   "#\n" \
                   "# Carbon Cycle:\n" \
                   "# Process sensitivity to atmospheric CO2:\n" \
                   "# - Ocean CO2 uptake (see ocean component)\n";
    if (co2_dep) {
        desc_file += "# - Land C exchange (see land component)\n";
    }
    desc_file += "# Process sensitivity to global mean SAT:\n";
    if (t_dep) {
        desc_file += "# - Ocean CO2 uptake (see ocean component)\n";
        if (LAND_T_DEP) {
            desc_file += "# - Land C exchange (see land component)\n";
        }
        desc_file += "#\n#\n";
    } else {
        desc_file += "# none\n#\n#\n";
    }
    desc_file += OCEAN_DOC;
    desc_file += LAND_DOC;
    desc_file += "#time_(yr) glob_temp_dev_(°C) RF_tot_(W/m²) RF_CO2_(W/m²) RF_nonCO2_(W/m²) RF_budget_(W/m²) " \
                 "ocean_heat_uptake_(PW) co2_atm_(ppm) co2_seasurf_(ppm) atm_CO2_change_(GtC/yr) " \
                 "fossil_CO2_em_(GtC/yr) budget_C_uptake_(GtC/yr) ocean_C_uptake_(GtC/yr) land_C_uptake_(GtC/yr) " \
                 "NPP_(GtC/yr) RH_(GtC/yr) LandC_(GtC) dDIC_(µmol/kg) fdeep_(GtC/yr)";
    if (!DO_INTERPOL) {
        desc_file += " midtime_(yr)";
    }
    desc_file += "\n#time glob_temp_dev RF_tot RF_CO2 RF_nonCO2 RF_budget ocean_heat_uptake co2_atm co2_seasurf " \
                 "atm_CO2_change fossil_CO2_em budget_C_uptake ocean_C_uptake land_C_uptake NPP RH LandC dDIC fdeep";
    if (!DO_INTERPOL) {
        desc_file += " midtime";
    }
    desc_file += "\n";
    for (i = 0; i < ntime; i++) {
        sprintf(tmp_str, "%20.10f", _time[i]);
        desc_file += tmp_str;
        memset(tmp_str, 0, sizeof tmp_str);
        sprintf(tmp_str, "%20.10f", temp[i]);
        desc_file += tmp_str;
        memset(tmp_str, 0, sizeof tmp_str);
        sprintf(tmp_str, "%20.10f", rf[i]);
        desc_file += tmp_str;
        memset(tmp_str, 0, sizeof tmp_str);
        sprintf(tmp_str, "%20.10f", rfc[i]);
        desc_file += tmp_str;
        memset(tmp_str, 0, sizeof tmp_str);
        sprintf(tmp_str, "%20.10f", rfnc[i]);
        desc_file += tmp_str;
        memset(tmp_str, 0, sizeof tmp_str);
        sprintf(tmp_str, "%20.10f", rfb[i]);
        desc_file += tmp_str;
        memset(tmp_str, 0, sizeof tmp_str);
        sprintf(tmp_str, "%20.10f", fh[i]);
        desc_file += tmp_str;
        memset(tmp_str, 0, sizeof tmp_str);
        sprintf(tmp_str, "%20.10f", ma[i] / PPMTOGT);
        desc_file += tmp_str;
        memset(tmp_str, 0, sizeof tmp_str);
        sprintf(tmp_str, "%20.10f", co2_atm0 + dpcs[i]);
        desc_file += tmp_str;
        memset(tmp_str, 0, sizeof tmp_str);
        sprintf(tmp_str, "%20.10f", fa[i]);
        desc_file += tmp_str;
        memset(tmp_str, 0, sizeof tmp_str);
        sprintf(tmp_str, "%20.10f", fossil_emissions[i]);
        desc_file += tmp_str;
        memset(tmp_str, 0, sizeof tmp_str);
        sprintf(tmp_str, "%20.10f", fb[i]);
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
        sprintf(tmp_str, "%20.10f", ml[i]);
        desc_file += tmp_str;
        memset(tmp_str, 0, sizeof tmp_str);
        sprintf(tmp_str, "%20.10f", ms[i]);
        desc_file += tmp_str;
        memset(tmp_str, 0, sizeof tmp_str);
        sprintf(tmp_str, "%20.10f", f_deep[i]);
        desc_file += tmp_str;
        memset(tmp_str, 0, sizeof tmp_str);
        if (!DO_INTERPOL) {
            sprintf(tmp_str, "%20.10f", midtime[i]);
            desc_file += tmp_str;
            memset(tmp_str, 0, sizeof tmp_str);
        }
        desc_file += "\n";
    }
    output_path = getPath("output");
    file_name = getFileName();
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
