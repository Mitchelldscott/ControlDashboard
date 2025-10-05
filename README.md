# ControlDashboard

A Dash + Plotly package to setup interactive simulations

![quadcopter design tool](./examples/stabilization_demo.png)

Here’s a concise **“Getting Started (Developing)”** section for your `ControlDashboard` README:

---

## 🧭 Getting Started (Developing)

### 1. **Install Julia**

* Download and install from [julialang.org/downloads](https://julialang.org/downloads/).
* Recommended: Julia ≥ 1.10.

### 2. **Set up the environment**

```bash
git clone https://github.com/<your-username>/ControlDashboard.jl.git
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

```julia
include("test/runtests.jl")
```

```julia
$> julia --project=.
julia> ]
(ControlDashboard) pkg> test
```

### 5. **Run examples**

```bash
julia --project=. examples/quadcopter_attitude_stabilization.jl
```
