/**
 * @file main.cpp
 * @brief Main file for T Operator Product Generation.
 * 
 * This program generates all permutations of input patterns, considers modifications
 * by toggling the rows, includes the transposed version of the patterns, and processes
 * binary patterns from files. It utilizes OpenMP for parallel processing and handles
 * specific cosets based on current T count and free multiply depth. The main function
 * orchestrates reading input files, generating permutations, and executing the main
 * loop for T count iteration and pattern processing.
 * 
 * @author Connor Mooney
 * @author Michael Jarret
 * @author Andrew Glaudell
 * @author Jacob Weston
 * @author Mingzhen Tian
 * @version 6/1/21
 */

#include <chrono>
#include <set>
#include <fstream>
#include <vector>
#include <omp.h>
#include <thread>
#include <dirent.h> // Directory Entry
#include "Globals.hpp"
#include "utils.hpp"

using namespace std;

/**
 * @brief Generates all permutations of the input pattern, its transposed version, and
 *        modifications by toggling the rows, inserting them into a set.
 * @param pat The pattern to permute and modify.
 * @return A set containing all unique permutations and modifications of the pattern.
 */

static std::set<pattern> permutation_set(pattern & patternToPermute)
{   
    std::set<pattern> perms;
    // Initialize a set to store the initial and transposed pattern
    std::set<pattern> inits;
    inits.insert(patternToPermute); // Insert the original pattern
    inits.insert(patternToPermute.transpose()); // Insert the transposed pattern

    for(pattern p : inits) {
        perms.insert(p); 
        // Initialize an array to represent the rows of the pattern
        int row[6] = {0, 1, 2, 3, 4, 5};

        // Generate all permutations of the pattern
        while (std::next_permutation(row, row + 6))
        {
            // Create a new pattern based on the current permutation
            pattern perm_of_orig;
            for (int c = 0; c < 6; c++)
            {
                for (int r = 0; r < 6; r++)
                    perm_of_orig.arr[c][r] = p.arr[c][row[r]];
            }
            // Order the new pattern lexicographically and insert it into the set
            perm_of_orig.lexicographic_order();
            perms.insert(perm_of_orig);
            // Iterate over all possible combinations of row modifications
            for(unsigned int counter = 0; counter < (1 << 6); counter++) {
                pattern mod_of_perm = perm_of_orig; // Start with a copy of the original permutation
                for(int j = 0; j < 6; j++) {
                    if((counter >> j) & 1) {
                        // Modify the j-th row if the j-th bit of 'counter' is set
                        mod_of_perm.mod_row(j);
                    }
                }
                // Order the modified pattern lexicographically and insert it into the set
                mod_of_perm.lexicographic_order();
                perms.insert(mod_of_perm);
            }
        }
    }
    return perms;
}

/**
 * @brief Inserts all permutations of a given pattern into the global pattern set.
 * @param patternToInsert The pattern to permute and insert.
 */
static void insert_all_permutations(pattern & patternToInsert)
{   
    std::set<pattern> perms_of_p = permutation_set(patternToInsert); 
    pattern_set.insert(perms_of_p.begin(),perms_of_p.end());
}
/**
 * @brief Erases all permutations of a given pattern from the global pattern set.
 * @param p The pattern whose permutations are to be erased.
 */

static void erase_all_permutations(pattern &p)
{
    std::set<pattern> perms_of_p = permutation_set(p);
    for(const pattern &erase : perms_of_p) pattern_set.erase(erase);
}

/**
 * @brief Reads binary patterns from a file, processes each line as a pattern, orders
 *        it lexicographically, and inserts all its permutations into a set.
 * @param pattern_file_path Path to the pattern file to be read.
 */
static void read_pattern_file(std::string pattern_file_path)
{
    if(pattern_file_path.empty()) return;

    std::cout << "[Read] Reading patterns from " << pattern_file << std::endl;
    
    std::ifstream patternFile(pattern_file_path); // More descriptive variable name

    if (!patternFile.is_open())
    {
        std::cerr << "Failed to open pattern file: " << pattern_file_path << std::endl;
        return;
    }

    bool isCSV = pattern_file_path.substr(pattern_file_path.find_last_of(".") + 1) == "csv";

    

     // Renamed for clarity
    std::string line;
    while (std::getline(patternFile, line))
    {
        std::string binaryString = isCSV ? utils::convert_csv_line_to_binary(line) : line;
        pattern currentPattern(binaryString); // More descriptive name
        currentPattern.lexicographic_order();
        // std::cout << currentPattern;
        insert_all_permutations(currentPattern);
    }

    patternFile.close(); // Close file after processing
    // Handle special case of the identity pattern
    pattern identityPattern = pattern::identity();
    pattern_set.erase(identityPattern);
    pattern_set.erase(identityPattern.pattern_mod());
    std::cout << "[Finished] Loaded " << pattern_set.size() << " non-identity patterns." << std::endl;
}

/**
 * @brief Checks if the given pattern is valid based on specific criteria.
 * @param pat The pattern to check.
 * @return True if the pattern is valid, false otherwise.
 */

static bool is_valid_pattern(pattern pat) {
    SO6 S, St;
    for(int col = 0; col < 6; col++) {
        for(int row = 0; row < 6; row++) {
            S[col][row][0]=pat.arr[col][row].first;
            S[col][row][1]=pat.arr[col][row].second;
        }
    }
    for(int col = 0; col < 6; col++) {
        for(int row = 0; row < 6; row++) {
            St[col][row][0]=pat.arr[row][col].first;
            St[col][row][1]=pat.arr[row][col].second;
        }
    }

    S = St*S;
    for(int col = 0; col < 6; col++) {
        for(int row = 0; row < 6; row++) {
            if(!(St[col][row] == 0)) return false;
        }
    }
    return true;
}

/**
 * @brief Reads binary patterns from a file and processes them for case analysis.
 * @param string path The path to the file containing case patterns.
*/
static void read_case_file(std::string path)
{
    // Attempt to open the file at the given path for reading.
    std::ifstream file(path); 

    // Check if the file was successfully opened.
    if (!file.is_open())
    {
        // Print an error message to the standard error stream and return early from the function.
        std::cerr << "Failed to open template file: " << path << std::endl;
        return;
    }

    std::string line;// Holds each line read from the file. 

    // Read the file line by line.
    while (std::getline(file, line))
    {
        // Directly use the line as a binary string representation of the pattern.
        std::string binaryString = line;

        // Check if the line has the expected length for a valid template
        if(line.length()!=36) {
            // If not, print an error message and exit the program with a failure status.
            std::cout << "Invalid template!\n";
            std::exit(EXIT_FAILURE);
        }
        // Convert the binary string to a pattern object.
        pattern currentPattern(binaryString); 
         // Order the current pattern according to a specific case ordering.
        currentPattern.case_order();

         // Add the ordered pattern to a global collection for later use.
        cases.push_back(currentPattern);
        // Optionally print information about the current pattern for debugging or logging.
        std::cout << currentPattern.case_string();
        std::cout << "Case number " << currentPattern.case_number() << "\n";
    }
    file.close();  // Close the file after all lines have been processed.
}

/** 
 * @brief Gets the current high-resolution time point.
 * @param time_point representing the current time.
 */
static std::chrono::_V2::high_resolution_clock::time_point now()
{
    return chrono::_V2::high_resolution_clock::now();
}

/**
 * @brief Erases the pattern of an SO6 from the global pattern set if it exists.
 * @param so6Pattern The SO6 pattern to be erased.
 * @return True if the pattern was found and erased, false otherwise.
 */
static bool erase_pattern(SO6 & so6Pattern) {
    // Convert the SO6 object to a pattern. 
    pattern pat = so6Pattern.to_pattern();
    //Track if the pattern was successfully erased.
    bool ret = false;

    //Check if the pattern exists in the global pattern set
    if (pattern_set.find(pat) != pattern_set.end()) {
        // Lock the global pattern set to ensure thread safety
        omp_set_lock(&lock);
        // Double check after grabbing the lock
        if (pattern_set.find(pat) != pattern_set.end()) {
            erase_all_permutations(pat);
            // Indicate that the pattern was successfully erased.
            ret = true;
        } 
         // Release the lock after modifying the set to allow other threads to access the set.
        omp_unset_lock(&lock);
    }
    // Return true if the pattern was erased, false otherwise.
    return ret;
}

/**
 * @brief Records the circuit string of an SO6 to the output file.
 * @param s The SO6 whose circuit string is to be recorded.
 * @param outputFileStream Output file stream where the circuit string is written.
 */
static void record_pattern(SO6 & s, std::ofstream& outputFileStream) {
    // Lock a shared mutex before accessing the shared output file stream.
    omp_set_lock(&lock);
    // Write the string representation of the SO6 pattern's circuit to the output file.
    outputFileStream << s.circuit_string() << std::endl;
    // Unlock the mutex after the write operation is complete, allowing other threads to access the output file stream.
    omp_unset_lock(&lock);
}

/**
 * @brief Erases the pattern of an SO6 object from a global set and records it to a file.
 * @param SO6 &s The SO6 object to process.
 * @param ofstream& of The file stream to write to.
 */
static void erase_and_record_pattern(SO6 &s, std::ofstream& of) {
    if(erase_pattern(s)) record_pattern(s,of);
}

/**
* @brief Reads a ".dat" file containing strings of gates circuits and processes them.
* @param file_name The name of the ".dat" file to read.
*/
static void read_dat(std::string file_name) {
    std::string line; // To store each line read from the file
    std::ifstream file(file_name); // Open the file for reading
    std::ofstream null_stream("/dev/null");// Stream for discarding output 

    if (file.is_open()) { // Check if the file was successfully opened.
        while (getline(file, line)) { // Read the file line by line.
            SO6 s = SO6::reconstruct_from_circuit_string(line); // Reconstruct SO6 object from each line's circuit string.
            std::cout << "current size: " << pattern_set.size() << "\n";
            erase_pattern(s); // Attempt to erase the pattern corresponding to the SO6 object from the global pattern set.
            std::cout << s.circuit_string() << "\n";  // Output the circuit string of the SO6 object.
        }
    }
}

/**
* @brief Calculates the time elapsed since a given time point.
* @param time_point &s The start time point.
*/
static std::string time_since(std::chrono::_V2::high_resolution_clock::time_point &s)
{
    chrono::duration<double> duration = now() - s;
    int64_t time = chrono::duration_cast<chrono::milliseconds>(duration).count();
    if (time < 1000)
        return std::to_string(time).append("ms");
    if (time < 60000)
        return std::to_string((float)time / 1000).substr(0, 5).append("s");
    if (time < 3600000)
        return std::to_string((float)time / 60000).substr(0, 5).append("min");
    if (time < 86400000)
        return std::to_string((float)time / 3600000).substr(0, 5).append("hr");
    return std::to_string((float)time / 86400000).substr(0, 5).append("days");
}

/**
 * @brief Method to report the start of a T count iteration
 * @param T the number of the current T count
 */
static void report_begin_T_count(const int T) {
    std::cout << " ||\t[Start] Beginning T=" << T << std::endl;
    tcount_init_time = now();
}

/**
 * @brief Method to report percentage completed
 * @param c current integer countr
 * @param s total size
 */
inline static void report_percent_complete(const uint64_t &c, const uint64_t s)
{
    if ((c & 0x7F) == 0) // if c is divisible by 128
    {
        std::cout << "\033[A\033[A\r ||\t↪ [Progress] Processing .....    "
                  << (100*c/s) << "\%" << "\n ||\t↪ [Patterns] "
                  << pattern_set.size() << " patterns remain." << std::endl;
    }
}

/**
 * @brief Function to report completion
 * @param matrices_found how many matrices were found
 */
static void finish_io(const uint &matrices_found, const bool b, std::ofstream &of) {
    std::cout << "\033[A\033[A\r ||\t↪ [Progress] Processing .....    100\%" << std::endl; 
    std::cout << " ||\t↪ [Patterns] " << pattern_set.size() << " patterns remain." << std::endl;
    if(b) std::cout << " ||\t↪ [Finished] Found " << matrices_found << " new matrices in " << time_since(tcount_init_time) << "\n ||" << std::endl;
    else std::cout << " ||\t↪ [Finished] Completed in " << time_since(tcount_init_time) << "\n ||" << std::endl;
    of.close();
}

/**
 * @brief Prepares the output file stream for recording results based on the current T count iteration.
 * 
 * This function opens an output file specific to the current T count iteration and prepares it for writing.
 * It reports the beginning of the T count iteration, including information about the depth of free multiplication
 * and other relevant details. The function ensures that the output file is opened successfully and is ready for
 * data recording. If the file cannot be opened, the program exits.
 *
 * The function also handles logging the initial state of the process, including the total number of patterns remaining,
 * and prints status messages to the console to inform the user about the current state of the operation.
 *
 * @param t The current T count iteration. Used to name the output file and control the logging output.
 * @param stored_depth_max The maximum depth stored from previous iterations. Used for reporting.
 * @param target_T_count The target T count the algorithm aims to achieve. Used for calculating the number of generating sets and for reporting.
 * 
 * @return A standard output file stream (std::ofstream) opened and ready for writing the results of the current T count iteration.
 */
static std::ofstream prepare_T_count_io(const int t, uint8_t &stored_depth_max, uint8_t &target_T_count);

static std::ofstream prepare_T_count_io(const int t, uint8_t &stored_depth_max, uint8_t &target_T_count) {
    // Calculate the depth of multiplication that can be performed freely without exceeding the target T count.
    int free_multiply_depth = utils::free_multiply_depth(target_T_count,stored_depth_max);
   // Special logging for the start of T=1 processing.
    if (t == 1) {
        std::cout << "\n\n[Begin] Generating T=1 through T=" << static_cast<int>(stored_depth_max);
        std::cout << " iteratively, but will only save ";

        if (free_multiply_depth >= 2) {
            std::cout << (free_multiply_depth == 2 ? "T=2" : "2≤T≤" + std::to_string(free_multiply_depth));
            std::cout << " and ";
        }
        std::cout << "T=" << static_cast<int>(stored_depth_max) << " in memory\n ||\n";
    }

    report_begin_T_count(t);
     // Prepare the filename based on the current T count and open the file.
    std::string file_string = "./data/" + to_string(t) + ".dat";
    std::ofstream of;
    of.open(file_string, ios::out | ios::trunc);
    // Exit if the file couldn't be opened.
    if(!of.is_open()) std::exit(0);
    // Notify that the file is opened for saving results.
    std::cout << " ||\t↪ [Save] Opening file " << file_string << "\n";
    // Additional logging for T counts immediately after stored depth and beyond.
    if(t == stored_depth_max+1) 
        std::cout << " ||\t↪ [Rep] Left multiplying everything by T₀\n";
    if(t > stored_depth_max+1) 
        std::cout << " ||\t↪ [Rep] Using generating_set[" << t-stored_depth_max-1 << "]\n";
    // Initial progress and patterns remaining.
    std::cout << " ||\t↪ [Progress] Processing .....    0\%\n"
                << " ||\t↪ [Patterns] " << pattern_set.size() << " patterns remain." << std::endl;
    return of; // Return the file stream for further use.
}

/**
 * @brief Store specific cosets T_0{curr} based on the current T count and free multiply depth.
 * This method saves a subset of the SO6 objects to the generating set, which are used in later iterations.
 * 
 * @param curr_T_count The current T count in the main computation loop.
 * @param free_multiply_depth The depth until which free multiplication is performed.
 * @param num_generating_sets The total number of generating sets.
 * @param current The current set of SO6 objects.
 * @param generating_set Reference to an array of vectors of SO6 objects to store the generated sets.
 */
void storeCosets(int curr_T_count, 
                 std::set<SO6>& current, std::vector<SO6> &generating_set)
{
    // Calculate the number of generating sets needed
    int ngs = utils::num_generating_sets(target_T_count,stored_depth_max);
    // Check if the current T count is within the range of required generating sets.
    if (curr_T_count < ngs)
    {
         // Log the action of saving the current set
        std::cout << "\033[A\r ||\t↪ [Save] Saving coset T₀{T=" << curr_T_count + 1 << "} as generating_set[" << curr_T_count << "]\n ||" << std::endl;
         // Copy the current set of SO6 objects into the generating set vector.
        generating_set = std::vector<SO6>(current.begin(),current.end());
        // Erase any SO6 objects from the generating set that do not meet a specific criterion
        generating_set.erase(std::remove_if(generating_set.begin(), generating_set.end(),
                                [](SO6& S) {
                                    return (S.circuit_string().back() == '0');
                                }),
                    generating_set.end());

 // Left multiply each SO6 object in the generating set by T₀
        for(SO6 &S : generating_set) {
                S = S.left_multiply_by_T(0);
        }
    }
}

/**
 * @brief The main function of the program. Initializes parameters, reads input files,
 *        and executes the main processing loop.
 * @param argc Number of command-line arguments.
 * @param argv Array of command-line arguments.
 * @return Exit status of the program.
 */
int main(int argc, char **argv)
{
    auto program_init_time = now();          // Begin timekeeping
    Globals::setParameters(argc, argv);      // Initialize parameters to command line argument
    Globals::configure();                    // Configure the globals to remove inconsistencies
    read_pattern_file(pattern_file);
    // Initialize sets for current and previous patterns.
    set<SO6> prior, current = std::set<SO6>({root});

    // Determine the number of generating sets needed.
    int ngs = utils::num_generating_sets(target_T_count, stored_depth_max);

    std::vector<SO6> generating_set[ngs];    
     // Main loop for processing T counts up to the stored depth maximum.
    for (int curr_T_count = 0; curr_T_count < stored_depth_max; ++curr_T_count)
    {
        set<SO6> next;
        // Prepare file I/O for the current T count.
        std::ofstream of = prepare_T_count_io(curr_T_count+1,stored_depth_max,target_T_count);

        uint64_t count = 0, interval_size = 15*current.size();

        for (int T = 0; T < 15; T++)
        {                
            for(const SO6& S : current)
            {
                SO6 toInsert = S.left_multiply_by_T(T);
                if (next.insert(toInsert).second)  // Want to insert toInsert no matter what, if we haven't seen this
                { 
                    if(T==0) erase_and_record_pattern(toInsert, of); // Only erase the pattern if the last T operation was done by 0, since we'll do things up to permutation?
                }
                report_percent_complete(++count, interval_size);
            }
        }

        utils::setDifference(next,prior);
         // Prepare for the next iteration.
        utils::rotate_and_clear(prior, current, next); 

        finish_io(current.size(), true, of);
        storeCosets(curr_T_count, current, generating_set[curr_T_count]);
    }
    
    std::set<SO6>().swap(prior); // Swap to clear
    std::cout << " ||\n[End] Stored T=" << (int)stored_depth_max << " as current to generate T=" << stored_depth_max + 1 << " through T=" << (int)target_T_count << "\n" << std::endl;

    std::vector<SO6> to_compute = utils::convert_to_vector_and_clear(current);

    std::cout << "[Report] Current patterns: " << pattern_set.size() << std::endl;

    std::cout << "[Begin] Beginning brute force multiply.\n ||" << std::endl;
    uint64_t set_size = to_compute.size();
    uint64_t interval_size = std::ceil(set_size / THREADS); // Equally divide among threads, not sure how to balance but each should take about the same time

    for (int curr_T_count = stored_depth_max; curr_T_count < target_T_count; ++curr_T_count)
    {    
        std::ofstream of = prepare_T_count_io(curr_T_count+1,stored_depth_max, target_T_count);

        std::vector<std::ofstream> file_stream(THREADS);
        // for(int i=0; i < THREADS; i++) {
        //     std::string file_name = "./data/reductions/thread" + std::to_string(i) + ".dat";
        //     file_stream[i].open(file_name, std::ios::out | std::ios::trunc);
        // }

        omp_init_lock(&lock);
        #pragma omp parallel for schedule(static, interval_size) num_threads(THREADS)
        for (uint64_t i = 0; i < set_size; i++)
        {
            int current_thread = omp_get_thread_num();
            const SO6 &S = to_compute.at(i); 
            if (omp_get_thread_num() == 0)
                report_percent_complete(i % interval_size, interval_size);

            if (curr_T_count == stored_depth_max)
            {
                SO6 N = S.left_multiply_by_T(0);
                if(!cases_flag) {
                    erase_and_record_pattern(N, of);
                    continue;
                }

                // for(pattern P : cases) {
                //     SO6 post = P*N;
                //     if(post.getLDE() == -1) {
                //         file_stream[current_thread] << N.circuit_string() << std::endl;
                //     }
                // }
            }
 
            for (const SO6 &G : generating_set[curr_T_count-stored_depth_max - 1])
            {
                SO6 N = G*S; 
                if(!cases_flag) {
                    erase_and_record_pattern(N, of);
                    continue;
                }
            }
        }
        omp_destroy_lock(&lock);
        finish_io(0, false, of);
        for(auto &stream : file_stream) stream.close();
    }
    std::cout << " ||\n[Finished] Free multiply complete.\n\n[Time] Total time elapsed: " << time_since(program_init_time) << std::endl;
    return 0;
}
