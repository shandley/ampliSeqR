#' Detect and profile system hardware for optimal configuration
#'
#' @description
#' Detects available system hardware resources including CPU cores, memory,
#' GPU availability, and architecture-specific features. Results can be used
#' to configure parallel backends for optimized performance.
#'
#' @param verbose Logical indicating whether to print hardware detection results.
#' Default is TRUE.
#'
#' @return A list containing the detected hardware information:
#' \itemize{
#'   \item{\code{os}: Operating system details}
#'   \item{\code{cpu_cores}: Number of physical CPU cores}
#'   \item{\code{logical_cores}: Number of logical CPU cores (with hyperthreading)}
#'   \item{\code{architecture}: CPU architecture (e.g., "x86_64", "arm64")}
#'   \item{\code{total_memory}: Total available RAM in GB}
#'   \item{\code{free_memory}: Currently free RAM in GB}
#'   \item{\code{gpu_available}: Whether a compatible GPU is available}
#'   \item{\code{gpu_info}: GPU details if available}
#'   \item{\code{optimal_threads}: Recommended number of threads for parallelization}
#'   \item{\code{recommended_backend}: Recommended parallel backend type}
#' }
#'
#' @examples
#' \dontrun{
#' # Get hardware profile
#' hw_profile <- detectHardware()
#'
#' # Use it to configure parallel backend
#' backend <- configureParallelBackend(hw_profile)
#' }
#'
#' @export
detectHardware <- function(verbose = TRUE) {
  # Get basic system information
  os_type <- Sys.info()["sysname"]
  os_release <- Sys.info()["release"]
  
  # Detect CPU cores
  cpu_cores <- ._detectCPUCores()
  logical_cores <- parallel::detectCores(logical = TRUE)
  physical_cores <- ._detectPhysicalCores()
  
  # Detect CPU architecture
  architecture <- ._detectArchitecture()
  
  # Detect available memory
  memory_info <- getAvailableMemory()
  
  # Detect GPU (if available)
  gpu_info <- ._detectGPU()
  
  # Determine optimal thread configuration
  optimal_threads <- ._calculateOptimalThreads(
    physical_cores = physical_cores,
    logical_cores = logical_cores,
    memory_info = memory_info,
    gpu_available = gpu_info$available
  )
  
  # Determine recommended backend
  recommended_backend <- ._recommendBackend(
    os_type = os_type,
    gpu_available = gpu_info$available,
    architecture = architecture
  )
  
  # Create hardware profile
  hardware_profile <- list(
    os = list(
      type = os_type,
      release = os_release
    ),
    cpu_cores = physical_cores,
    logical_cores = logical_cores,
    architecture = architecture,
    total_memory = memory_info$total,
    free_memory = memory_info$available,
    gpu_available = gpu_info$available,
    gpu_info = if (gpu_info$available) gpu_info$details else NULL,
    optimal_threads = optimal_threads,
    recommended_backend = recommended_backend
  )
  
  # Print information if verbose mode is on
  if (verbose) {
    ._printHardwareProfile(hardware_profile)
  }
  
  return(hardware_profile)
}

#' Configure parallel backend based on hardware profile
#'
#' @description
#' Sets up an appropriate parallel backend based on hardware detection results.
#' Configures BiocParallel or future backends with optimized settings for
#' the detected hardware.
#'
#' @param hardware_profile List containing hardware information as returned by
#'   \code{detectHardware()}.
#' @param strategy Character string specifying the parallelization strategy.
#'   One of "auto" (automatic selection), "multicore" (forking), "multisession" 
#'   (separate R sessions), or "gpu" (GPU acceleration if available).
#'   Default is "auto".
#' @param memory_limit Numeric value specifying the maximum memory (in GB) to use.
#'   Default is NULL, which means the function will determine a safe limit
#'   (usually 80% of available memory).
#'
#' @return A configured parallel backend object (either BiocParallel or future)
#'   ready to use in downstream functions.
#'
#' @examples
#' \dontrun{
#' # Detect hardware and configure backend
#' hw_profile <- detectHardware(verbose = FALSE)
#' backend <- configureParallelBackend(hw_profile)
#'
#' # Use backend in functions
#' # (Future enhancement: integrate with dada2 pipeline)
#' }
#'
#' @importFrom magrittr %>%
#' @export
configureParallelBackend <- function(hardware_profile, 
                                     strategy = c("auto", "multicore", "multisession", "gpu"),
                                     memory_limit = NULL) {
  strategy <- match.arg(strategy)
  
  if (is.null(hardware_profile)) {
    stop("Hardware profile is required. Run detectHardware() first.")
  }
  
  # Choose strategy if "auto" is selected
  if (strategy == "auto") {
    strategy <- hardware_profile$recommended_backend
  }
  
  # Set optimal thread count
  thread_count <- hardware_profile$optimal_threads
  
  # Set memory limits
  if (is.null(memory_limit)) {
    # Use 80% of available memory as a safe default
    memory_limit <- hardware_profile$free_memory * 0.8
  }
  
  # Set up the requested backend
  if (strategy == "gpu" && hardware_profile$gpu_available) {
    backend <- ._configureGPUBackend(
      hardware_profile, 
      memory_limit
    )
  } else if (strategy == "multicore" && hardware_profile$os$type != "Windows") {
    # BiocParallel:
    if (requireNamespace("BiocParallel", quietly = TRUE)) {
      workers <- min(thread_count, BiocParallel::multicoreWorkers())
      backend <- BiocParallel::MulticoreParam(
        workers = workers,
        # Set task size and other parameters based on memory availability
        tasks = ._estimateSafeTasks(workers, memory_limit, "default")
      )
    } 
    # future:
    else if (requireNamespace("future", quietly = TRUE)) {
      future::plan(future::multicore, workers = thread_count)
      backend <- list(
        type = "future_multicore",
        workers = thread_count,
        memory_limit = memory_limit
      )
    } 
    else {
      # Fallback if neither package is available
      backend <- list(
        type = "base_parallel",
        workers = thread_count,
        cluster = parallel::makeCluster(thread_count)
      )
      warning("Neither BiocParallel nor future found, using base parallel package.")
    }
  } else {
    # Default to multisession (works on all platforms)
    # BiocParallel:
    if (requireNamespace("BiocParallel", quietly = TRUE)) {
      backend <- BiocParallel::SnowParam(
        workers = thread_count,
        type = "SOCK"
      )
    } 
    # future:
    else if (requireNamespace("future", quietly = TRUE)) {
      future::plan(future::multisession, workers = thread_count)
      backend <- list(
        type = "future_multisession",
        workers = thread_count,
        memory_limit = memory_limit
      )
    } 
    else {
      # Fallback if neither package is available
      backend <- list(
        type = "base_parallel",
        workers = thread_count,
        cluster = parallel::makeCluster(thread_count)
      )
      warning("Neither BiocParallel nor future found, using base parallel package.")
    }
  }
  
  # Set memory safeguards for BiocParallel backends
  if (inherits(backend, "BiocParallelParam")) {
    BiocParallel::bpworkers(backend) <- thread_count
    # Set memory constraints if possible and needed
    if (memory_limit > 0 && inherits(backend, "SnowParam")) {
      # Divide memory limit by number of workers
      memory_per_worker <- memory_limit / thread_count
      # Convert GB to bytes
      memory_per_worker_bytes <- memory_per_worker * 1024^3
      
      # Set memory limit (this is a soft limit, might not be honored by all OS)
      BiocParallel::bpbackend(backend) <- BiocParallel::SerialParam(
        progressbar = TRUE
      )
    }
  }
  
  return(backend)
}

#' Get available system memory
#'
#' @description
#' Detects total and available memory on the system using platform-specific methods.
#'
#' @return A list with total and available memory in GB:
#' \itemize{
#'   \item{\code{total}: Total system memory in GB}
#'   \item{\code{available}: Available memory in GB}
#' }
#'
#' @examples
#' \dontrun{
#' # Get memory information
#' mem_info <- getAvailableMemory()
#' print(paste("Total memory:", mem_info$total, "GB"))
#' print(paste("Available memory:", mem_info$available, "GB"))
#' }
#'
#' @export
getAvailableMemory <- function() {
  os_type <- Sys.info()["sysname"]
  
  # Default values in case detection fails
  total_memory <- NA
  available_memory <- NA
  
  # macOS memory detection
  if (os_type == "Darwin") {
    tryCatch({
      # Get total memory
      hw_mem_cmd <- system("sysctl hw.memsize", intern = TRUE)
      total_memory <- as.numeric(sub(".*: ", "", hw_mem_cmd)) / (1024^3)
      
      # Get available memory using vm_stat
      vm_stat <- system("vm_stat", intern = TRUE)
      
      # Parse page size
      page_size_line <- vm_stat[1]
      page_size <- as.numeric(gsub("[^0-9]", "", page_size_line))
      
      # Parse free pages
      free_pages_line <- grep("Pages free:", vm_stat, value = TRUE)
      free_pages <- as.numeric(gsub("[^0-9]", "", free_pages_line))
      
      # Parse inactive pages
      inactive_pages_line <- grep("Pages inactive:", vm_stat, value = TRUE)
      inactive_pages <- as.numeric(gsub("[^0-9]", "", inactive_pages_line))
      
      # Calculate available memory
      available_pages <- free_pages + inactive_pages
      available_memory <- (available_pages * page_size) / (1024^3)
      
    }, error = function(e) {
      warning("Failed to detect memory on macOS: ", e$message)
      # Default to using gc information as fallback
      mem <- gc(full = TRUE, reset = TRUE)
      total_memory <- sum(mem[, 1] + mem[, 2]) / (1024^2)
      available_memory <- total_memory * 0.8  # Assume 80% is available
    })
  }
  # Linux memory detection 
  else if (os_type == "Linux") {
    tryCatch({
      mem_info <- system("cat /proc/meminfo", intern = TRUE)
      
      # Extract total memory
      total_line <- grep("MemTotal:", mem_info, value = TRUE)
      total_kb <- as.numeric(gsub("[^0-9]", "", total_line))
      total_memory <- total_kb / (1024^2)
      
      # Extract available memory
      available_line <- grep("MemAvailable:", mem_info, value = TRUE)
      if (length(available_line) > 0) {
        available_kb <- as.numeric(gsub("[^0-9]", "", available_line))
      } else {
        # Fallback for older kernels
        free_line <- grep("MemFree:", mem_info, value = TRUE)
        cached_line <- grep("Cached:", mem_info, value = TRUE)[1]
        buffers_line <- grep("Buffers:", mem_info, value = TRUE)
        
        free_kb <- as.numeric(gsub("[^0-9]", "", free_line))
        cached_kb <- as.numeric(gsub("[^0-9]", "", cached_line))
        buffers_kb <- as.numeric(gsub("[^0-9]", "", buffers_line))
        
        available_kb <- free_kb + cached_kb + buffers_kb
      }
      
      available_memory <- available_kb / (1024^2)
      
    }, error = function(e) {
      warning("Failed to detect memory on Linux: ", e$message)
      # Fallback
      mem <- gc(full = TRUE, reset = TRUE)
      total_memory <- sum(mem[, 1] + mem[, 2]) / (1024^2)
      available_memory <- total_memory * 0.8  # Assume 80% is available
    })
  }
  # Windows memory detection
  else if (os_type == "Windows") {
    tryCatch({
      if (requireNamespace("utils", quietly = TRUE)) {
        mem_info <- utils::memory.size(max = TRUE)
        total_memory <- mem_info / 1024
        
        mem_info_current <- utils::memory.size(max = FALSE)
        available_memory <- (mem_info - mem_info_current) / 1024
      } else {
        # Fallback if utils package is not available
        win_mem <- system("wmic OS get FreePhysicalMemory,TotalVisibleMemorySize /Value", intern = TRUE)
        
        total_line <- grep("TotalVisibleMemorySize", win_mem, value = TRUE)
        total_kb <- as.numeric(sub(".*=", "", total_line))
        total_memory <- total_kb / (1024^2)
        
        free_line <- grep("FreePhysicalMemory", win_mem, value = TRUE)
        free_kb <- as.numeric(sub(".*=", "", free_line))
        available_memory <- free_kb / (1024^2)
      }
    }, error = function(e) {
      warning("Failed to detect memory on Windows: ", e$message)
      # Fallback
      mem <- gc(full = TRUE, reset = TRUE)
      total_memory <- sum(mem[, 1] + mem[, 2]) / (1024^2)
      available_memory <- total_memory * 0.8  # Assume 80% is available
    })
  }
  # Other OS types
  else {
    warning("Unsupported OS type for detailed memory detection: ", os_type)
    # Use gc as a fallback
    mem <- gc(full = TRUE, reset = TRUE)
    total_memory <- sum(mem[, 1] + mem[, 2]) / (1024^2)
    available_memory <- total_memory * 0.8  # Assume 80% is available
  }
  
  # Make sure values are reasonable
  if (is.na(total_memory) || total_memory <= 0) {
    total_memory <- 4  # Assume 4GB as default
  }
  if (is.na(available_memory) || available_memory <= 0) {
    available_memory <- total_memory * 0.8  # Assume 80% is available
  }
  
  return(list(
    total = total_memory,
    available = available_memory
  ))
}

#' Estimate memory requirements per task
#'
#' @description
#' Estimates the amount of memory needed for various amplicon sequencing operations.
#' This helps functions allocate appropriate resources and avoid memory oversubscription.
#'
#' @param task_type Character string specifying the task type. Options are:
#'   "default", "filtering", "error_learning", "denoising", "merging", or "chimera_removal".
#' @param read_size Numeric value specifying the average size of read data in bytes.
#'   Default is 1e6 (1 million bytes).
#'
#' @return A numeric value representing the estimated memory requirement in GB.
#'
#' @examples
#' \dontrun{
#' # Estimate memory for filtering 10M reads
#' mem_needed <- estimateMemoryPerTask("filtering", read_size = 1e7)
#' print(paste("Memory needed:", mem_needed, "GB"))
#' }
#'
#' @export
estimateMemoryPerTask <- function(task_type = "default", read_size = 1e6) {
  # Convert bytes to GB
  size_gb <- read_size / (1024^3)
  
  # Memory estimation factors based on task type
  # These are approximate multipliers based on empirical observations
  memory_factors <- list(
    default = 2.0,            # Base case
    filtering = 1.5,          # Quality filtering
    error_learning = 4.0,     # DADA2 error rate learning
    denoising = 3.0,          # ASV inference
    merging = 2.0,            # Merging paired reads
    chimera_removal = 2.5     # Chimera detection and removal
  )
  
  # Get factor for the requested task type
  factor <- memory_factors[[task_type]]
  if (is.null(factor)) {
    warning("Unknown task type: ", task_type, ". Using default factor.")
    factor <- memory_factors[["default"]]
  }
  
  # Calculate estimated memory need
  estimated_memory <- size_gb * factor
  
  # Set reasonable minimum
  if (estimated_memory < 0.1) {
    estimated_memory <- 0.1  # Minimum 100MB
  }
  
  return(estimated_memory)
}

# Internal helper functions

#' Detect CPU cores using platform-specific methods
#'
#' @return Integer number of logical CPU cores
#' @keywords internal
._detectCPUCores <- function() {
  tryCatch({
    # Try to use parallel package first
    if (requireNamespace("parallel", quietly = TRUE)) {
      return(parallel::detectCores(logical = TRUE))
    }
    
    # Fallback methods by OS
    os_type <- Sys.info()["sysname"]
    
    if (os_type == "Darwin") {
      cores <- as.integer(system("sysctl -n hw.ncpu", intern = TRUE))
    } else if (os_type == "Linux") {
      cores <- as.integer(system("nproc", intern = TRUE))
    } else if (os_type == "Windows") {
      cores <- as.integer(Sys.getenv("NUMBER_OF_PROCESSORS"))
    } else {
      warning("Unknown OS type for core detection: ", os_type)
      cores <- 2  # Conservative fallback
    }
    
    if (is.na(cores) || cores <= 0) {
      cores <- 2  # Default if detection fails
    }
    
    return(cores)
  }, error = function(e) {
    warning("Failed to detect CPU cores: ", e$message)
    return(2)  # Conservative fallback
  })
}

#' Detect physical CPU cores
#'
#' @return Integer number of physical CPU cores
#' @keywords internal
._detectPhysicalCores <- function() {
  os_type <- Sys.info()["sysname"]
  
  tryCatch({
    if (os_type == "Darwin") {
      # macOS
      physical_cores <- as.integer(system("sysctl -n hw.physicalcpu", intern = TRUE))
    } else if (os_type == "Linux") {
      # Linux - count unique core IDs
      cores_cmd <- system("cat /proc/cpuinfo | grep 'core id' | sort -u | wc -l", intern = TRUE)
      physical_cores <- as.integer(cores_cmd)
      
      # Fallback if the above doesn't work
      if (is.na(physical_cores) || physical_cores == 0) {
        sockets_cmd <- system("cat /proc/cpuinfo | grep 'physical id' | sort -u | wc -l", intern = TRUE)
        cores_per_socket_cmd <- system("cat /proc/cpuinfo | grep 'cpu cores' | head -1", intern = TRUE)
        sockets <- as.integer(sockets_cmd)
        cores_per_socket <- as.integer(sub(".*: ", "", cores_per_socket_cmd))
        
        if (!is.na(sockets) && !is.na(cores_per_socket) && sockets > 0 && cores_per_socket > 0) {
          physical_cores <- sockets * cores_per_socket
        } else {
          # Last resort - use half of logical cores as estimate
          physical_cores <- max(1, floor(parallel::detectCores(logical = TRUE) / 2))
        }
      }
    } else if (os_type == "Windows") {
      # Windows - this is approximate, as Windows doesn't have a direct way to get physical cores
      # Use WMI or try to estimate from logical cores
      logical_cores <- as.integer(Sys.getenv("NUMBER_OF_PROCESSORS"))
      physical_cores <- max(1, floor(logical_cores / 2))  # Rough estimate assuming hyperthreading
    } else {
      warning("Unknown OS type for physical core detection: ", os_type)
      logical_cores <- parallel::detectCores(logical = TRUE)
      physical_cores <- max(1, floor(logical_cores / 2))  # Rough estimate
    }
    
    # Ensure we return a reasonable number
    if (is.na(physical_cores) || physical_cores <= 0) {
      physical_cores <- 1
    }
    
    return(physical_cores)
  }, error = function(e) {
    warning("Failed to detect physical CPU cores: ", e$message)
    # Fallback to a conservative estimate
    logical_cores <- parallel::detectCores(logical = TRUE)
    return(max(1, floor(logical_cores / 2)))
  })
}

#' Detect CPU architecture
#'
#' @return Character string with CPU architecture
#' @keywords internal
._detectArchitecture <- function() {
  os_type <- Sys.info()["sysname"]
  
  tryCatch({
    if (os_type == "Darwin") {
      # macOS
      arch_cmd <- system("uname -m", intern = TRUE)
      arch <- arch_cmd[1]
      
      # Check if running on Apple Silicon via Rosetta
      if (arch == "x86_64") {
        # Check if running under Rosetta on Apple Silicon
        sysctl_cmd <- system("sysctl -n sysctl.proc_translated 2>/dev/null || echo 0", intern = TRUE)
        is_translated <- as.integer(sysctl_cmd[1])
        
        if (is_translated == 1) {
          arch <- "arm64_rosetta"  # Running ARM code via Rosetta
        }
      }
    } else if (os_type == "Linux") {
      # Linux
      arch_cmd <- system("uname -m", intern = TRUE)
      arch <- arch_cmd[1]
    } else if (os_type == "Windows") {
      # Windows
      arch <- Sys.getenv("PROCESSOR_ARCHITECTURE")
      if (arch == "") {
        arch_cmd <- system("wmic OS get OSArchitecture", intern = TRUE)
        arch <- arch_cmd[2]  # Skip header line
        
        # Map Windows architecture names to standard ones
        if (grepl("64", arch)) {
          arch <- "x86_64"
        } else if (grepl("32", arch)) {
          arch <- "x86"
        } else if (grepl("ARM", arch, ignore.case = TRUE)) {
          arch <- "arm64"
        }
      }
    } else {
      # Other
      arch <- R.version$arch
    }
    
    # Clean up the architecture string
    arch <- tolower(trimws(arch))
    
    # Map to standardized architecture names
    if (arch %in% c("x86_64", "amd64", "x64")) {
      return("x86_64")
    } else if (arch %in% c("i386", "i686", "x86")) {
      return("x86")
    } else if (arch %in% c("arm64", "aarch64", "arm64_rosetta")) {
      return(arch)
    } else {
      return(arch)  # Return as-is if unknown
    }
  }, error = function(e) {
    warning("Failed to detect architecture: ", e$message)
    return(R.version$arch)  # Fallback to R's architecture
  })
}

#' Detect GPU availability and information
#'
#' @return A list with GPU information:
#' \itemize{
#'   \item{\code{available}: Logical indicating whether a supported GPU is available}
#'   \item{\code{details}: List with GPU details if available}
#' }
#' @keywords internal
._detectGPU <- function() {
  # Default return value
  result <- list(
    available = FALSE,
    details = NULL
  )
  
  # Check for common GPU compute packages
  gpu_available <- FALSE
  
  # Check for TensorFlow/Keras with GPU
  if (requireNamespace("tensorflow", quietly = TRUE)) {
    tryCatch({
      tf_gpu <- tensorflow::tf$config$list_physical_devices("GPU")
      if (length(tf_gpu) > 0) {
        gpu_available <- TRUE
        result$details <- list(
          type = "tensorflow",
          device_count = length(tf_gpu),
          devices = sapply(tf_gpu, function(x) x$name)
        )
      }
    }, error = function(e) {
      # TensorFlow not working or no GPU
    })
  }
  
  # Check for CUDA availability
  if (!gpu_available && requireNamespace("RcppCNPy", quietly = TRUE)) {
    tryCatch({
      cuda_available <- system("nvidia-smi -L", intern = TRUE, ignore.stdout = FALSE, ignore.stderr = TRUE)
      if (length(cuda_available) > 0) {
        gpu_available <- TRUE
        result$details <- list(
          type = "cuda",
          device_count = length(cuda_available),
          devices = cuda_available
        )
      }
    }, error = function(e) {
      # NVIDIA tools not available or no GPU
    })
  }
  
  # Check for ROCm (AMD) availability
  if (!gpu_available) {
    tryCatch({
      rocm_available <- system("rocm-smi --showproductname", intern = TRUE, ignore.stdout = FALSE, ignore.stderr = TRUE)
      if (length(rocm_available) > 0) {
        gpu_available <- TRUE
        result$details <- list(
          type = "rocm",
          device_count = length(rocm_available),
          devices = rocm_available
        )
      }
    }, error = function(e) {
      # ROCm tools not available or no GPU
    })
  }
  
  # Check for Mac GPU/Metal
  if (!gpu_available && Sys.info()["sysname"] == "Darwin") {
    arch <- ._detectArchitecture()
    if (arch == "arm64") {
      # Apple Silicon likely has GPU
      gpu_available <- TRUE
      result$details <- list(
        type = "metal",
        device_count = 1,
        devices = "Apple Silicon Integrated GPU"
      )
    }
  }
  
  result$available <- gpu_available
  return(result)
}

#' Calculate optimal thread count based on hardware profile
#'
#' @param physical_cores Integer indicating the number of physical CPU cores
#' @param logical_cores Integer indicating the number of logical CPU cores
#' @param memory_info List with memory information
#' @param gpu_available Logical indicating whether a GPU is available
#'
#' @return Integer representing the optimal number of threads to use
#' @keywords internal
._calculateOptimalThreads <- function(physical_cores, logical_cores, memory_info, gpu_available) {
  # Start with all logical cores
  optimal_threads <- logical_cores
  
  # Adjust based on available memory
  # Rule of thumb: 1-2GB per thread for demanding bioinformatics tasks
  memory_based_threads <- floor(memory_info$available / 2)
  
  # For IO-bound operations where CPU might be waiting for disk, 
  # using more threads than cores can be beneficial
  io_optimal_threads <- min(logical_cores * 1.5, logical_cores + 4)
  
  # For CPU-bound operations, staying at or below logical core count is better
  cpu_optimal_threads <- logical_cores
  
  # If GPU is available, we might want to reserve some CPU for GPU coordination
  if (gpu_available) {
    gpu_optimal_threads <- max(1, logical_cores - 2)  # Reserve 2 cores for GPU coordination
  } else {
    gpu_optimal_threads <- logical_cores
  }
  
  # Consider task type (generalized for now, but could be more specific)
  # For high-memory DADA2 operations, be more conservative
  dada2_optimal_threads <- min(physical_cores, memory_based_threads, logical_cores)
  
  # Take the most conservative approach for now
  optimal_threads <- min(
    optimal_threads,
    memory_based_threads,
    dada2_optimal_threads,
    gpu_optimal_threads
  )
  
  # Ensure we use at least 1 thread
  optimal_threads <- max(1, optimal_threads)
  
  # Round to integer
  optimal_threads <- as.integer(optimal_threads)
  
  return(optimal_threads)
}

#' Recommend the appropriate parallel backend type
#'
#' @param os_type Character string indicating the operating system type
#' @param gpu_available Logical indicating whether a GPU is available
#' @param architecture Character string indicating the CPU architecture
#'
#' @return Character string with the recommended backend type
#' @keywords internal
._recommendBackend <- function(os_type, gpu_available, architecture) {
  # Start with default recommendation
  recommendation <- "multisession"
  
  # If GPU is available and we have proper GPU support
  if (gpu_available && requireNamespace("tensorflow", quietly = TRUE)) {
    recommendation <- "gpu"
  }
  # On Unix-like systems, multicore is often more efficient
  else if (os_type %in% c("Darwin", "Linux")) {
    recommendation <- "multicore"
  }
  
  return(recommendation)
}

#' Estimate safe number of tasks for a given number of workers and memory limit
#'
#' @param workers Integer indicating the number of worker processes
#' @param memory_limit Numeric indicating the memory limit in GB
#' @param task_type Character string indicating the task type
#'
#' @return Integer representing the safe number of tasks
#' @keywords internal
._estimateSafeTasks <- function(workers, memory_limit, task_type) {
  # Get memory per task estimate
  mem_per_task <- estimateMemoryPerTask(task_type)
  
  # Calculate how many tasks can run with the given memory
  max_concurrent_tasks <- floor(memory_limit / mem_per_task)
  
  # Ensure we don't schedule more concurrent tasks than we have workers
  safe_tasks <- min(max_concurrent_tasks, workers * 2)
  
  # Make sure we have at least one task
  safe_tasks <- max(1, safe_tasks)
  
  return(as.integer(safe_tasks))
}

#' Configure a GPU-specific backend
#'
#' @param hardware_profile List containing hardware information
#' @param memory_limit Numeric value indicating the memory limit in GB
#'
#' @return A configured GPU backend object
#' @keywords internal
._configureGPUBackend <- function(hardware_profile, memory_limit) {
  gpu_type <- hardware_profile$gpu_info$type
  
  if (gpu_type == "tensorflow" && requireNamespace("tensorflow", quietly = TRUE)) {
    # Set up TensorFlow backend
    # Note: This is a placeholder for future integration
    # For now, just return a basic configuration
    backend <- list(
      type = "tensorflow_gpu",
      device_count = hardware_profile$gpu_info$device_count,
      memory_limit = memory_limit
    )
  } else if (gpu_type == "cuda" && requireNamespace("RcppCNPy", quietly = TRUE)) {
    # Set up CUDA backend
    backend <- list(
      type = "cuda_gpu",
      device_count = hardware_profile$gpu_info$device_count,
      memory_limit = memory_limit
    )
  } else if (gpu_type == "metal" && Sys.info()["sysname"] == "Darwin") {
    # Set up Apple Metal backend
    backend <- list(
      type = "metal_gpu",
      device_count = 1,
      memory_limit = memory_limit
    )
  } else {
    warning("GPU detected but no supported GPU framework found. Falling back to CPU-based backend.")
    # Fall back to multisession
    if (requireNamespace("BiocParallel", quietly = TRUE)) {
      backend <- BiocParallel::SnowParam(
        workers = hardware_profile$optimal_threads,
        type = "SOCK"
      )
    } else if (requireNamespace("future", quietly = TRUE)) {
      future::plan(future::multisession, workers = hardware_profile$optimal_threads)
      backend <- list(
        type = "future_multisession",
        workers = hardware_profile$optimal_threads,
        memory_limit = memory_limit
      )
    } else {
      backend <- list(
        type = "base_parallel",
        workers = hardware_profile$optimal_threads,
        cluster = parallel::makeCluster(hardware_profile$optimal_threads)
      )
    }
  }
  
  return(backend)
}

#' Print hardware profile information in a readable format
#'
#' @param hardware_profile List containing hardware information
#'
#' @return None (invisible NULL)
#' @keywords internal
._printHardwareProfile <- function(hardware_profile) {
  cat("\n=== Hardware Profile ===\n")
  cat(sprintf("OS: %s %s\n", hardware_profile$os$type, hardware_profile$os$release))
  cat(sprintf("CPU: %d physical cores, %d logical cores\n", 
              hardware_profile$cpu_cores, hardware_profile$logical_cores))
  cat(sprintf("Architecture: %s\n", hardware_profile$architecture))
  cat(sprintf("Memory: %.1f GB total, %.1f GB available\n", 
              hardware_profile$total_memory, hardware_profile$free_memory))
  
  if (hardware_profile$gpu_available) {
    cat("GPU: Available\n")
    if (!is.null(hardware_profile$gpu_info)) {
      cat(sprintf("  - Type: %s\n", hardware_profile$gpu_info$type))
      cat(sprintf("  - Devices: %d\n", hardware_profile$gpu_info$device_count))
    }
  } else {
    cat("GPU: Not available\n")
  }
  
  cat(sprintf("Optimal thread count: %d\n", hardware_profile$optimal_threads))
  cat(sprintf("Recommended backend: %s\n", hardware_profile$recommended_backend))
  cat("========================\n")
  
  return(invisible(NULL))
}