#ifdef __cplusplus
#include "xsreader.h"
#endif

/**
 * @var xs_directory
 * @brief A character array with the cross-section library directory.
 * @details By default this is in the xs-lib folder
 */
static std::string xs_directory = "pinspec/xs-lib/";


/**
 * @brief Sets the directory for the cross-section library.
 * @param directory character array for the cross-section library directory
 */
void setXSLibDirectory(const char* directory) {
    xs_directory = directory;
    log_printf(INFO, "Set the cross-section library directory "
                        "to: %s", xs_directory.c_str());
    return;
}


/**
 * @brief Returns the directory for the cross-section library directory.
 * @return a character array with the cross-section library directory
 */
const char* getXSLibDirectory() {
    return xs_directory.c_str();
}


/**
 * @brief Restores the cross-section library to its original, installed state.
 * @details This function copies the original cross-section data files installed
 *          with PINSPEC from the pinspec/xs-lib/BackupXS folder into the main
 *          cross-section library directory pinspec/xs-lib. This function is
 *          primarily when a user is using PINSPEC's SLBW module to create
 *          artifical cross-sections from Reich-Moore date. This function is
 *          called anytime the python pinspec module is imported. 
 * @return returns an int representing whether or not the copy was successful
 */
int restoreXSLibrary() {
    std::string cmd = std::string("cp ") + xs_directory + 
                                                "/BackupXS/* " + xs_directory;
    int ret = system(cmd.c_str());
    return ret;
}


/**
 * @brief Parses an input file of cross-section data
 * @details Loads the energy values (as eV) and the cross-section values 
 *          (barns) into two float arrays. It takes in the energies in MeV
            since that is the most typical unit for ENDF nuclear data, but 
            converts the values into eV during parsing.
 * @param file the filename for the data
 * @param energies a pointer to a float array for the energies (MeV)
 * @param xs_values a pointer to a float array fo the xs values (barns)
 * @return the number of data points
 */
int parseCrossSections(const char* file, float* energies, float* xs_values) {

    /* Instantiate I/O variables */
    std::ifstream input_file(file, std::ios::in);
    std::string buff;
    int count = 0;

    /* get the first line of the input file */
    getline(input_file, buff);

    /* Parse over each line in the file */
    while(getline(input_file, buff)) {
        /* Load data into the arrays for this line */
        sscanf(buff.c_str(), "%f,%f", &energies[count], &xs_values[count]);
        count++;
    }

    /* Close the file and return the number of data points */
    input_file.close();

    return count;
}


/**
 * @brief Counts the number of lines in an cross-section input file
 * @param filename the file of interest
 * @return the number of lines in the file
 */
int getNumCrossSectionDataPoints(const char* filename) {

    /* Instantiate I/O variables */
    std::ifstream input_file(filename, std::ios::in);
    std::string line;
    int num_xs_values = -1;

    /* Loop over each line and update counter */
    while(getline(input_file, line))
        num_xs_values++;

    /* Close the file and return the number of lines */
    input_file.close();
    return num_xs_values;
}
