using Pkg
using Coverage

# Run your tests with coverage enabled
Pkg.test("ControlDashboard"; coverage = true)

# Process the coverage data
coverage_results = process_folder("src") # Or the specific folder containing your source code

# Print a summary
LCOV.writefile("lcov.info", coverage_results)
