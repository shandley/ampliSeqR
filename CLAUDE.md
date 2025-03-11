# ampliSeqR Development Guide

## Commands
- Build package: `devtools::build()` or `R CMD build .`
- Check package: `devtools::check()` or `R CMD check package.tar.gz`
- Install package: `devtools::install()` or `R CMD INSTALL package.tar.gz`
- Document package: `devtools::document()`
- Run all tests: `devtools::test()` 
- Run single test: `devtools::test_file("tests/testthat/test-file.R")`
- Style code: `styler::style_pkg()`
- Lint code: `lintr::lint_package()`

## Style Guidelines
- **Imports**: Use `@importFrom pkg function` in roxygen2, not `@import`
- **Formatting**: 2-space indent, 80 char line limit, `<-` for assignment
- **Naming**: 
  - Functions: lowerCamelCase verbs (`calculateDistance()`)
  - Variables: snake_case nouns (`sequence_data`)
  - Internal functions: prefix with dot (`._helper()`)
- **Error handling**: Use `stop()` for errors, `warning()`, `message()` appropriately
- **Documentation**: Document all exported functions with roxygen2, include examples

## Project-Specific Notes
- Optimize performance-critical code paths
- Consider GPU acceleration for appropriate algorithms
- Include benchmarking in tests where relevant