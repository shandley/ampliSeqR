test_that("detectHardware returns correctly structured hardware profile", {
  # Skip on CRAN to avoid environment-dependent test failures
  skip_on_cran()
  
  # Run hardware detection
  hw_profile <- detectHardware(verbose = FALSE)
  
  # Test structure
  expect_type(hw_profile, "list")
  expect_named(hw_profile, c("os", "cpu_cores", "logical_cores", "architecture", 
                             "total_memory", "free_memory", "gpu_available", 
                             "gpu_info", "optimal_threads", "recommended_backend"))
  
  # Test OS info
  expect_type(hw_profile$os, "list")
  expect_named(hw_profile$os, c("type", "release"))
  
  # Test CPU info
  expect_type(hw_profile$cpu_cores, "integer")
  expect_type(hw_profile$logical_cores, "integer")
  expect_gte(hw_profile$cpu_cores, 1)
  expect_gte(hw_profile$logical_cores, hw_profile$cpu_cores)
  
  # Test memory info
  expect_type(hw_profile$total_memory, "double")
  expect_type(hw_profile$free_memory, "double")
  expect_gte(hw_profile$total_memory, 0)
  expect_gte(hw_profile$free_memory, 0)
  
  # Test architecture
  expect_type(hw_profile$architecture, "character")
  
  # Test GPU info
  expect_type(hw_profile$gpu_available, "logical")
  
  # Test optimal threads
  expect_type(hw_profile$optimal_threads, "integer")
  expect_gte(hw_profile$optimal_threads, 1)
  
  # Test recommended backend
  expect_type(hw_profile$recommended_backend, "character")
  expect_true(hw_profile$recommended_backend %in% c("multicore", "multisession", "gpu"))
})

test_that("getAvailableMemory returns memory information", {
  # Skip on CRAN to avoid environment-dependent test failures
  skip_on_cran()
  
  mem_info <- getAvailableMemory()
  
  # Test structure
  expect_type(mem_info, "list")
  expect_named(mem_info, c("total", "available"))
  
  # Test values
  expect_type(mem_info$total, "double")
  expect_type(mem_info$available, "double")
  expect_gte(mem_info$total, 0)
  expect_gte(mem_info$available, 0)
  expect_lte(mem_info$available, mem_info$total * 1.1)  # Allow 10% margin for fluctuations
})

test_that("estimateMemoryPerTask returns reasonable estimates", {
  # Test default case
  default_mem <- estimateMemoryPerTask()
  expect_type(default_mem, "double")
  expect_gte(default_mem, 0)
  
  # Test with different task types
  filtering_mem <- estimateMemoryPerTask("filtering", read_size = 1e7)
  expect_type(filtering_mem, "double")
  expect_gte(filtering_mem, 0)
  
  # Test with large read size
  large_mem <- estimateMemoryPerTask("default", read_size = 1e9)
  expect_gt(large_mem, default_mem)
  
  # Test with unknown task type (should use default)
  unknown_mem <- estimateMemoryPerTask("unknown_task")
  expect_type(unknown_mem, "double")
  expect_gte(unknown_mem, 0)
})

test_that("configureParallelBackend returns a backend object", {
  # Skip on CRAN to avoid environment-dependent test failures
  skip_on_cran()
  
  # Mock hardware profile
  hw_profile <- list(
    os = list(type = Sys.info()["sysname"], release = Sys.info()["release"]),
    cpu_cores = 2,
    logical_cores = 4,
    architecture = R.version$arch,
    total_memory = 8,
    free_memory = 4,
    gpu_available = FALSE,
    gpu_info = NULL,
    optimal_threads = 2,
    recommended_backend = "multisession"
  )
  
  # Configure backend with mock profile
  backend <- configureParallelBackend(hw_profile, strategy = "multisession")
  
  # Test structure (type depends on available packages)
  expect_true(
    inherits(backend, "BiocParallelParam") || 
    inherits(backend, "list")
  )
  
  # Test with memory limit
  backend_with_limit <- configureParallelBackend(hw_profile, memory_limit = 2)
  expect_true(
    inherits(backend_with_limit, "BiocParallelParam") || 
    inherits(backend_with_limit, "list")
  )
  
  # Test with NULL hardware profile (should error)
  expect_error(configureParallelBackend(NULL))
})

# Test that internal helper functions are working
test_that("internal helper functions work correctly", {
  # Skip on CRAN to avoid environment-dependent test failures
  skip_on_cran()
  
  # These are internal functions, so we need the triple colon
  physical_cores <- ampliSeqR:::._detectPhysicalCores()
  expect_type(physical_cores, "integer")
  expect_gte(physical_cores, 1)
  
  architecture <- ampliSeqR:::._detectArchitecture()
  expect_type(architecture, "character")
  
  gpu_info <- ampliSeqR:::._detectGPU()
  expect_type(gpu_info, "list")
  expect_named(gpu_info, c("available", "details"))
  expect_type(gpu_info$available, "logical")
})