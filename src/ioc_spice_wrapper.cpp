#include "ioccultcalc/ioc_spice_wrapper.h"
#include "spice_omp_wrapper.h"
#include <sstream>
#include <cstring>

namespace ioccultcalc {

IOCSpiceWrapper::IOCSpiceWrapper(int num_threads)
    : num_threads_(num_threads), initialized_(false) {
    
    int status = spice_omp_init(num_threads);
    if (status != 0) {
        char error_msg[1024];
        spice_omp_get_last_error(error_msg, sizeof(error_msg));
        throw SpiceException(std::string("Failed to initialize IOC_SPICE: ") + error_msg);
    }
    
    initialized_ = true;
}

IOCSpiceWrapper::~IOCSpiceWrapper() {
    if (initialized_) {
        spice_omp_unload_all();
    }
}

void IOCSpiceWrapper::loadKernel(const std::string& kernel_path) {
    int status = spice_omp_furnsh(kernel_path.c_str());
    checkError(status, "Loading kernel: " + kernel_path);
}

void IOCSpiceWrapper::loadKernels(const std::vector<std::string>& kernel_paths) {
    for (const auto& path : kernel_paths) {
        loadKernel(path);
    }
}

StateVector IOCSpiceWrapper::getState(
    const std::string& target,
    double et,
    const std::string& observer,
    const std::string& frame,
    const std::string& abcorr
) {
    SpiceDouble positions[1][6];
    SpiceDouble et_array[1] = { et };
    
    int status = spice_omp_multi_state(
        target.c_str(),
        observer.c_str(),
        frame.c_str(),
        abcorr.c_str(),
        et_array,
        1,
        positions
    );
    
    checkError(status, "Getting state for " + target);
    
    return StateVector(
        positions[0][0], positions[0][1], positions[0][2],
        positions[0][3], positions[0][4], positions[0][5]
    );
}

std::vector<StateVector> IOCSpiceWrapper::getMultiState(
    const std::string& target,
    const std::vector<double>& et_times,
    const std::string& observer,
    const std::string& frame,
    const std::string& abcorr
) {
    if (et_times.empty()) {
        return {};
    }
    
    int n_times = static_cast<int>(et_times.size());
    std::vector<SpiceDouble> et_array(et_times.begin(), et_times.end());
    
    // Allocate output array
    std::vector<SpiceDouble[6]> positions(n_times);
    
    int status = spice_omp_multi_state(
        target.c_str(),
        observer.c_str(),
        frame.c_str(),
        abcorr.c_str(),
        et_array.data(),
        n_times,
        positions.data()
    );
    
    checkError(status, "Getting multiple states for " + target);
    
    // Convert to StateVector
    std::vector<StateVector> result;
    result.reserve(n_times);
    for (int i = 0; i < n_times; ++i) {
        result.emplace_back(
            positions[i][0], positions[i][1], positions[i][2],
            positions[i][3], positions[i][4], positions[i][5]
        );
    }
    
    return result;
}

std::vector<StateVector> IOCSpiceWrapper::getMultiBodyState(
    const std::vector<std::string>& targets,
    double et,
    const std::string& observer,
    const std::string& frame,
    const std::string& abcorr
) {
    if (targets.empty()) {
        return {};
    }
    
    int n_targets = static_cast<int>(targets.size());
    
    // Convert target names to C strings
    std::vector<const char*> target_ptrs;
    target_ptrs.reserve(n_targets);
    for (const auto& target : targets) {
        target_ptrs.push_back(target.c_str());
    }
    
    // Allocate output array
    std::vector<SpiceDouble[6]> positions(n_targets);
    
    int status = spice_omp_multi_body_state(
        target_ptrs.data(),
        n_targets,
        et,
        observer.c_str(),
        frame.c_str(),
        abcorr.c_str(),
        positions.data()
    );
    
    checkError(status, "Getting multi-body states");
    
    // Convert to StateVector
    std::vector<StateVector> result;
    result.reserve(n_targets);
    for (int i = 0; i < n_targets; ++i) {
        result.emplace_back(
            positions[i][0], positions[i][1], positions[i][2],
            positions[i][3], positions[i][4], positions[i][5]
        );
    }
    
    return result;
}

double IOCSpiceWrapper::utcToET(const std::string& utc_str) {
    SpiceDouble et;
    int status = spice_omp_str2et(utc_str.c_str(), &et);
    checkError(status, "Converting UTC to ET: " + utc_str);
    return et;
}

std::string IOCSpiceWrapper::etToUTC(double et) {
    // Note: IOC_SPICE doesn't have et2utc yet, use CSPICE directly
    // This is a placeholder - you may need to add et2utc to IOC_SPICE
    char utc_str[128];
    et2utc_c(et, "ISOC", 3, sizeof(utc_str), utc_str);
    return std::string(utc_str);
}

void IOCSpiceWrapper::checkError(int status, const std::string& operation) {
    if (status != 0) {
        char error_msg[1024];
        spice_omp_get_last_error(error_msg, sizeof(error_msg));
        
        std::ostringstream oss;
        oss << "IOC_SPICE error during " << operation << ": " << error_msg;
        throw SpiceException(oss.str());
    }
}

} // namespace ioccultcalc
