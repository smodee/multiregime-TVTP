# multiregime-TVTP

**Multi-Regime Time-Varying Transition Probability Models**

A comprehensive R implementation of regime-switching models with time-varying transition probabilities. This repository extends the two-regime framework from Sendstad, Chronopoulos, & Li (2025) to support arbitrary K-regime models with multiple transition dynamics.

## Overview

This repository implements several regime-switching models where the transition probabilities between states can vary over time based on different dynamics:

- **Constant Model**: Fixed transition probabilities (baseline model)
- **TVP Model**: Autoregressive dynamics where transitions depend on lagged observations
- **Exogenous Model**: Transitions driven by external variables
- **GAS Model**: Score-driven dynamics based on the predictive likelihood (Generalized Autoregressive Score)

The key innovation is the extension from 2-regime to K-regime models while maintaining numerical stability and computational efficiency.

## Repository Structure

```
multiregime-TVTP/
│
├── models/                    # Core model implementations
│   ├── model_constant.R       # Constant transition probability model
│   ├── model_TVP.R            # Time-varying transition probability (autoregressive)
│   ├── model_exogenous.R      # Exogenous-driven transition model
│   └── model_GAS.R            # Generalized Autoregressive Score model
│
├── helpers/                   # Utility and helper functions
│   ├── utility_functions.R    # Basic utilities (logit, logistic, etc.)
│   ├── transition_helpers.R   # Transition matrix operations
│   ├── parameter_transforms.R # Parameter extraction and transformation
│   └── score_functions.R      # GAS score calculations and scaling
│
├── simulation/                # Simulation and testing scripts
│   ├── sim_GAS.R             # GAS model simulations
│   ├── sim_exogenous.R       # Exogenous model simulations
│   └── sim_comparisons.R     # Model comparison studies
│
├── comparisons/              # Model comparison and analysis
│   └── compare_tvtp_models.R # Comprehensive model comparisons
│
├── tests/                    # Unit tests (to be implemented)
│
└── results/                  # Output directory for results (gitignored)
```

## Installation

### Prerequisites

Required R packages:
```r
# Core dependencies
install.packages(c("MASS", "optimx", "numDeriv"))

# For GAS model (Gauss-Hermite quadrature)
install.packages("statmod")

# For parallel processing (optional but recommended)
install.packages(c("future", "future.apply"))

# For plotting and analysis
install.packages(c("ggplot2", "reshape2"))
```

### Setup

1. Clone the repository:
```bash
git clone https://github.com/yourusername/multiregime-TVTP.git
cd multiregime-TVTP
```

2. Source the main functions:
```r
# Load all models
source("models/model_constant.R")
source("models/model_TVP.R")
source("models/model_exogenous.R")
source("models/model_GAS.R")

# Load comparison tools
source("comparisons/compare_tvtp_models.R")
```

## Quick Start

### Basic Usage Example

```r
# Generate sample data from a 3-regime TVP model
K <- 3  # Number of regimes
N <- 1000  # Length of time series
M <- 1  # Number of simulations

# Set true parameters
mu_true <- c(-2, 0, 2)           # Regime means
sigma2_true <- c(0.5, 1, 1.5)    # Regime variances
init_trans_true <- rep(0.2, 6)   # Initial transition probabilities (K*(K-1))
A_true <- rep(0.1, 6)             # Autoregressive coefficients

# Generate data
data <- dataTVPCD(M, N, mu_true, sigma2_true, init_trans_true, A_true)
y <- data[1,]  # Use first simulation

# Estimate model parameters
estimate <- estimate_tvp_model(
  y = y,
  K = 3,
  verbose = TRUE
)

# View results
print(estimate$parameters)
print(estimate$diagnostics)
```

### Model Comparison

```r
# Compare different models on the same data
comparison <- compare_tvtp_models(
  y = y,
  K = 3,
  models = c("Constant", "TVP", "GAS"),
  verbose = TRUE
)

# View comparison results
print(comparison$summary)
print(comparison$rankings)
```

## Model Details

### 1. Constant Model
Fixed transition probabilities that don't change over time.

**Parameters:** `mu` (K), `sigma2` (K), `trans_prob` (K*(K-1))

**Functions:**
- `dataCD()` - Generate data
- `Rfiltering_CD()` - Filter/likelihood calculation
- `estimate_constant_model()` - Parameter estimation

### 2. TVP Model (Autoregressive)
Transition probabilities depend on lagged observations:
```
f[t+1] = omega + A * y[t]
p[t+1] = logistic(f[t+1])
```

**Parameters:** `mu` (K), `sigma2` (K), `init_trans` (K*(K-1)), `A` (K*(K-1))

**Functions:**
- `dataTVPCD()` - Generate data
- `Rfiltering_TVP()` - Filter/likelihood calculation
- `estimate_tvp_model()` - Parameter estimation

### 3. Exogenous Model
Transition probabilities driven by external variables:
```
f[t+1] = omega + A * X_exo[t]
p[t+1] = logistic(f[t+1])
```

**Parameters:** `mu` (K), `sigma2` (K), `init_trans` (K*(K-1)), `A` (K*(K-1))

**Functions:**
- `dataTVPXExoCD()` - Generate data
- `Rfiltering_TVPXExo()` - Filter/likelihood calculation
- `estimate_exo_model()` - Parameter estimation

### 4. GAS Model (Score-Driven)
Transition probabilities evolve based on the score of the predictive likelihood:
```
f[t+1] = omega + A * s[t] + B * (f[t] - omega)
p[t+1] = logistic(f[t+1])
```

**Parameters:** `mu` (K), `sigma2` (K), `init_trans` (K*(K-1)), `A` (K*(K-1)), `B` (K*(K-1))

**Functions:**
- `dataGASCD()` - Generate data
- `Rfiltering_GAS()` - Filter/likelihood calculation
- `estimate_gas_model()` - Parameter estimation

## Advanced Features

### Multi-Start Optimization
All estimation functions support multiple starting points for robust optimization:

```r
estimate <- estimate_tvp_model(
  y = y,
  K = 3,
  n_starts = 10,  # Use 10 random starting points
  parallel = TRUE, # Enable parallel processing
  cores = 4        # Use 4 CPU cores
)
```

### Custom Parameter Bounds
Specify bounds for parameter estimation:

```r
estimate <- estimate_gas_model(
  y = y,
  K = 3,
  lower_bounds = list(mu = -10, sigma2 = 0.01, A = 0, B = 0),
  upper_bounds = list(mu = 10, sigma2 = 10, A = 1, B = 1)
)
```

### Simulation Studies
Run comprehensive simulation studies:

```r
# GAS model simulation study
results <- run_gas_simulation_study(
  num_repetitions = 100,
  sample_sizes = c(500, 1000, 2000),
  K = 3,
  seed = 123,
  output_dir = "results/gas_simulation"
)

# Plot results
plot_gas_simulation_results(results, output_file = "results/gas_plots.pdf")
```

## Helper Functions

### Transition Matrix Operations
- `transition_matrix()` - Convert vector to transition matrix
- `convert_to_valid_probs()` - Ensure valid stochastic matrix
- `validate_transition_matrix()` - Check matrix validity

### Parameter Transformations
- `transform_parameters()` - Transform to unconstrained space
- `untransform_parameters()` - Transform back to natural space
- `mean_from_par()`, `sigma2_from_par()`, etc. - Extract components

### Score Functions (GAS)
- `calculate_gas_score()` - Compute scaled score
- `setup_gauss_hermite_quadrature()` - Setup numerical integration
- `apply_moore_penrose_scaling()` - Apply score scaling

## Performance Considerations

### Numerical Stability
- All models use log-likelihood calculations to avoid underflow
- Transition probabilities are constrained to valid ranges
- Variance parameters are log-transformed during optimization

### Computational Efficiency
- Parallel processing for multi-start optimization
- Vectorized operations where possible
- Efficient matrix operations for transition calculations

### Scaling to Large K
- Memory usage scales as O(K²) for transition matrices
- Computation time scales approximately as O(K³N) for filtering
- Recommended maximum K ≈ 10 for practical applications

## Testing and Validation

### Running Tests
```r
# Test individual models
source("tests/test_models.R")
run_model_tests()

# Test parameter recovery
source("tests/test_parameter_recovery.R")
test_parameter_recovery(K = 3, N = 1000)

# Test numerical stability
source("tests/test_stability.R")
test_numerical_stability()
```

### Validation Against Original Implementation
The repository includes `simulation.R` which contains the original 2-regime implementation. You can validate the extended models by comparing results for K=2:

```r
# Compare with original for K=2
source("tests/validate_against_original.R")
validate_two_regime_case()
```

## Known Issues and Limitations

### Current Limitations
1. **GAS Model Scaling**: The Moore-Penrose pseudo-inverse scaling for GAS models can be numerically unstable. The simple scaling method is recommended for production use.

2. **Large K Performance**: For K > 10, optimization becomes challenging due to the high-dimensional parameter space.

3. **Identification Issues**: Some parameter combinations may lead to label switching or weak identification.

### Work in Progress
- [ ] Implement comprehensive unit tests
- [ ] Add diagnostic tools for model selection
- [ ] Implement regime prediction functionality
- [ ] Add support for multivariate observations
- [ ] Develop interactive visualization tools

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Submit a pull request with a clear description

## Citation

If you use this code in your research, please cite:

**Original methodology (2-regime):**
```bibtex
@article{sendstad2025value,
  title={The Value of Turning-Point Detection for Optimal Investment},
  author={Sendstad, L.H. and Chronopoulos, M. and Li, X.},
  journal={Journal Name},
  year={2025}
}
```

**GAS methodology:**
```bibtex
@article{bazzi2017time,
  title={Time-varying transition probabilities for Markov regime switching models},
  author={Bazzi, M. and Blasques, F. and Koopman, S.J. and Lucas, A.},
  journal={Journal of Time Series Analysis},
  volume={38},
  number={3},
  pages={458--478},
  year={2017}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

Copyright (c) 2025 Samuel Modée

## Contact

Samuel Modée - Repository maintainer

For questions, issues, or contributions, please open an issue on GitHub.

## Acknowledgments

This implementation extends the work of Sendstad, Chronopoulos, & Li (2025) and incorporates ideas from Bazzi et al. (2017) for score-driven dynamics.
