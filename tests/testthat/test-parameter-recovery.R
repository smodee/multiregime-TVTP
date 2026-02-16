# Parameter Recovery Tests for multiregimeTVTP
# ==============================================
#
# These tests are OFF by default. To run them:
#   Sys.setenv(MULTIREGIME_PARAM_RECOVERY = "true")
#   devtools::test(filter = "parameter-recovery")
#
# Or from command line (PowerShell):
#   $env:MULTIREGIME_PARAM_RECOVERY="true"
#   Rscript -e "devtools::test(filter='parameter-recovery')"
#
# Assertion strategy (based on exploratory results, 5 reps x 16 scenarios):
#   - mu:         ASSERT at 15% relative tolerance
#   - sigma2:     ASSERT at 20% relative tolerance (robust to seed variability)
#   - A, B:       LOG ONLY — A recovery is fundamentally hard (0-60% pass rate)
#   - trans_prob: LOG ONLY — informational
#   - GAS:        LOG ONLY for all components — known sigma2/A identifiability issue (#29)
#
# Single replication per scenario with fixed seed, parallel=TRUE for speed.

skip_if(
  !identical(Sys.getenv("MULTIREGIME_PARAM_RECOVERY"), "true"),
  message = "Parameter recovery tests disabled. Set MULTIREGIME_PARAM_RECOVERY=true to run."
)
skip_on_cran()

# ---- Level 1: K=2, diag_probs=TRUE, equal_variances=TRUE --------------------

describe("Level 1: K=2 diagonal equal-variance parameter recovery", {

  it("recovers constant model parameters", {
    scenario <- recovery_build_scenario(1L, "constant")
    result   <- recovery_run_single(scenario, seed = 101L)

    expect_true(result$accuracy$mu$pass,
                info = paste("mu:", result$accuracy$mu$detail))
    expect_true(result$accuracy$sigma2$pass,
                info = paste("sigma2:", result$accuracy$sigma2$detail))

    message("L1 constant trans_prob: ", result$accuracy$trans_prob$detail)
  })

  it("recovers TVP model parameters", {
    scenario <- recovery_build_scenario(1L, "tvp")
    result   <- recovery_run_single(scenario, seed = 102L)

    expect_true(result$accuracy$mu$pass,
                info = paste("mu:", result$accuracy$mu$detail))
    expect_true(result$accuracy$sigma2$pass,
                info = paste("sigma2:", result$accuracy$sigma2$detail))

    message("L1 TVP trans_prob: ", result$accuracy$trans_prob$detail)
    message("L1 TVP A: ", result$accuracy$A$detail)
  })

  it("recovers exo model parameters", {
    scenario <- recovery_build_scenario(1L, "exo")
    result   <- recovery_run_single(scenario, seed = 103L)

    expect_true(result$accuracy$mu$pass,
                info = paste("mu:", result$accuracy$mu$detail))
    expect_true(result$accuracy$sigma2$pass,
                info = paste("sigma2:", result$accuracy$sigma2$detail))

    message("L1 exo trans_prob: ", result$accuracy$trans_prob$detail)
    message("L1 exo A: ", result$accuracy$A$detail)
  })

  it("GAS model estimation completes (known issues, see #29)", {
    scenario <- recovery_build_scenario(1L, "gas")
    result   <- recovery_run_single(scenario, seed = 104L)

    # Assert estimation completed without error
    expect_true(TRUE, info = "GAS estimation completed successfully")

    # Log only — GAS has known sigma2/A identifiability issues (issue #29)
    message("GAS L1 mu: ", result$accuracy$mu$detail)
    message("GAS L1 sigma2: ", result$accuracy$sigma2$detail)
    message("GAS L1 A: ", result$accuracy$A$detail)
    message("GAS L1 B: ", result$accuracy$B$detail)
  })
})

# ---- Level 2: K=2, diag_probs=FALSE, equal_variances=FALSE ------------------

describe("Level 2: K=2 off-diagonal distinct-variance parameter recovery", {

  it("recovers constant model parameters", {
    scenario <- recovery_build_scenario(2L, "constant")
    result   <- recovery_run_single(scenario, seed = 201L)

    expect_true(result$accuracy$mu$pass,
                info = paste("mu:", result$accuracy$mu$detail))
    expect_true(result$accuracy$sigma2$pass,
                info = paste("sigma2:", result$accuracy$sigma2$detail))

    message("L2 constant trans_prob: ", result$accuracy$trans_prob$detail)
  })

  it("recovers TVP model parameters", {
    scenario <- recovery_build_scenario(2L, "tvp")
    result   <- recovery_run_single(scenario, seed = 202L)

    expect_true(result$accuracy$mu$pass,
                info = paste("mu:", result$accuracy$mu$detail))
    expect_true(result$accuracy$sigma2$pass,
                info = paste("sigma2:", result$accuracy$sigma2$detail))

    message("L2 TVP trans_prob: ", result$accuracy$trans_prob$detail)
    message("L2 TVP A: ", result$accuracy$A$detail)
  })

  it("recovers exo model parameters", {
    scenario <- recovery_build_scenario(2L, "exo")
    result   <- recovery_run_single(scenario, seed = 203L)

    expect_true(result$accuracy$mu$pass,
                info = paste("mu:", result$accuracy$mu$detail))
    expect_true(result$accuracy$sigma2$pass,
                info = paste("sigma2:", result$accuracy$sigma2$detail))

    message("L2 exo trans_prob: ", result$accuracy$trans_prob$detail)
    message("L2 exo A: ", result$accuracy$A$detail)
  })

  it("GAS model estimation completes (known issues, see #29)", {
    scenario <- recovery_build_scenario(2L, "gas")
    result   <- recovery_run_single(scenario, seed = 204L)

    # Assert estimation completed without error
    expect_true(TRUE, info = "GAS estimation completed successfully")

    # Log only — GAS has known sigma2/A identifiability issues (issue #29)
    message("GAS L2 mu: ", result$accuracy$mu$detail)
    message("GAS L2 sigma2: ", result$accuracy$sigma2$detail)
    message("GAS L2 A: ", result$accuracy$A$detail)
    message("GAS L2 B: ", result$accuracy$B$detail)
  })
})

# ---- Level 3: K=3, diag_probs=FALSE, equal_variances=TRUE -------------------

describe("Level 3: K=3 off-diagonal equal-variance parameter recovery", {

  it("recovers constant model parameters", {
    scenario <- recovery_build_scenario(3L, "constant")
    result   <- recovery_run_single(scenario, seed = 301L)

    expect_true(result$accuracy$mu$pass,
                info = paste("mu:", result$accuracy$mu$detail))
    expect_true(result$accuracy$sigma2$pass,
                info = paste("sigma2:", result$accuracy$sigma2$detail))

    message("L3 constant trans_prob: ", result$accuracy$trans_prob$detail)
  })

  it("recovers TVP model parameters", {
    scenario <- recovery_build_scenario(3L, "tvp")
    result   <- recovery_run_single(scenario, seed = 302L)

    expect_true(result$accuracy$mu$pass,
                info = paste("mu:", result$accuracy$mu$detail))
    expect_true(result$accuracy$sigma2$pass,
                info = paste("sigma2:", result$accuracy$sigma2$detail))

    message("L3 TVP trans_prob: ", result$accuracy$trans_prob$detail)
    message("L3 TVP A: ", result$accuracy$A$detail)
  })

  it("recovers exo model parameters", {
    scenario <- recovery_build_scenario(3L, "exo")
    result   <- recovery_run_single(scenario, seed = 303L)

    expect_true(result$accuracy$mu$pass,
                info = paste("mu:", result$accuracy$mu$detail))
    expect_true(result$accuracy$sigma2$pass,
                info = paste("sigma2:", result$accuracy$sigma2$detail))

    message("L3 exo trans_prob: ", result$accuracy$trans_prob$detail)
    message("L3 exo A: ", result$accuracy$A$detail)
  })

  it("GAS model estimation completes (known issues, see #29)", {
    scenario <- recovery_build_scenario(3L, "gas")
    result   <- recovery_run_single(scenario, seed = 304L)

    # Assert estimation completed without error
    expect_true(TRUE, info = "GAS estimation completed successfully")

    # Log only — GAS has known sigma2/A identifiability issues (issue #29)
    message("GAS L3 mu: ", result$accuracy$mu$detail)
    message("GAS L3 sigma2: ", result$accuracy$sigma2$detail)
    message("GAS L3 A: ", result$accuracy$A$detail)
    message("GAS L3 B: ", result$accuracy$B$detail)
  })
})

# ---- Level 4: K=3, diag_probs=FALSE, equal_variances=FALSE ------------------

describe("Level 4: K=3 off-diagonal distinct-variance parameter recovery", {

  it("recovers constant model parameters", {
    scenario <- recovery_build_scenario(4L, "constant")
    result   <- recovery_run_single(scenario, seed = 401L)

    expect_true(result$accuracy$mu$pass,
                info = paste("mu:", result$accuracy$mu$detail))
    expect_true(result$accuracy$sigma2$pass,
                info = paste("sigma2:", result$accuracy$sigma2$detail))

    message("L4 constant trans_prob: ", result$accuracy$trans_prob$detail)
  })

  it("recovers TVP model parameters", {
    scenario <- recovery_build_scenario(4L, "tvp")
    result   <- recovery_run_single(scenario, seed = 402L)

    expect_true(result$accuracy$mu$pass,
                info = paste("mu:", result$accuracy$mu$detail))
    expect_true(result$accuracy$sigma2$pass,
                info = paste("sigma2:", result$accuracy$sigma2$detail))

    message("L4 TVP trans_prob: ", result$accuracy$trans_prob$detail)
    message("L4 TVP A: ", result$accuracy$A$detail)
  })

  it("recovers exo model parameters", {
    scenario <- recovery_build_scenario(4L, "exo")
    result   <- recovery_run_single(scenario, seed = 403L)

    expect_true(result$accuracy$mu$pass,
                info = paste("mu:", result$accuracy$mu$detail))
    expect_true(result$accuracy$sigma2$pass,
                info = paste("sigma2:", result$accuracy$sigma2$detail))

    message("L4 exo trans_prob: ", result$accuracy$trans_prob$detail)
    message("L4 exo A: ", result$accuracy$A$detail)
  })

  it("GAS model estimation completes (known issues, see #29)", {
    scenario <- recovery_build_scenario(4L, "gas")
    result   <- recovery_run_single(scenario, seed = 404L)

    # Assert estimation completed without error
    expect_true(TRUE, info = "GAS estimation completed successfully")

    # Log only — GAS has known sigma2/A identifiability issues (issue #29)
    message("GAS L4 mu: ", result$accuracy$mu$detail)
    message("GAS L4 sigma2: ", result$accuracy$sigma2$detail)
    message("GAS L4 A: ", result$accuracy$A$detail)
    message("GAS L4 B: ", result$accuracy$B$detail)
  })
})
