# ControlDashboard.jl

[![Build Status](https://github.com/MitchellDScott/ControlDashboard/actions/workflows/run_tests_with_coverage.yml/badge.svg)](https://github.com/MitchellDScott/ControlDashboard/actions/workflows/run_tests_with_coverage.yml)
[![Coverage](https://codecov.io/gh/mitchelldscott/ControlDashboard/branch/master/graph/badge.svg?token=)](https://app.codecov.io/gh/mitchelldscott/ControlDashboard)
<!-- [![Docs](https://img.shields.io/badge/docs-html-blue.svg)](https://MitchellDScott.github.io/ControlDashboard/) -->

A visualization toolkit for control system simulations and interactive design in Julia.

![quadcopter design tool](./examples/stabilization_demo.png)

## Getting Started

### 1. **Install Julia**

* Download and install from [julialang.org/downloads](https://julialang.org/downloads/).
* Recommended: Julia ≥ 1.11.

### 2. **Set up the environment**

```bash
git clone https://github.com/MitchellDScott/ControlDashboard.jl.git
cd ControlDashboard
julia --project=.
```

### 3. **Activate and instantiate dependencies**

```julia
import Pkg
Pkg.activate(".")
Pkg.instantiate()
```

### 4. **Run tests**

In a julia script/console:

```julia
include("test/runtests.jl")
```

or from the Julia REPL

```julia
julia> ]
(@v1.11) pkg> activate .
(ControlDashboard) pkg> test
```

### 5. **Run examples**

```bash
julia --project=. examples/sinusoid.jl
```

```bash
julia --project=. examples/quadcopter_attitude_stabilization.jl
```
