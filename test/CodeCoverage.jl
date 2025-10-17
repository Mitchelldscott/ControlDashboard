using Pkg;

Pkg.add("Coverage");
Pkg.add("Dates");

Pkg.instantiate();

using Coverage, Dates

# Run your tests with coverage enabled
Pkg.test("ControlDashboard"; coverage = true)

# Process the coverage data
coverage_results = process_folder("src") # Or the specific folder containing your source code

# Write the coverage results to the new file
LCOV.writefile("lcov.info", coverage_results)
