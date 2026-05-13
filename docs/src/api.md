# API Reference

## Data types

```@docs
TransientSample
TransientDataset
ModalBasis
FiniteDifference
EnvelopeODE
RecoveryRateCurve
PolynomialSurface
BifurcationPoint
BifurcationDiagram
```

!!! note "Recovery-rate method dispatch"
    `FiniteDifference` and `EnvelopeODE` are concrete subtypes of the abstract
    type `BifurcationForecasting.RecoveryRateMethod`.  Pass an instance to
    [`estimate_recovery_rate`](@ref) or the `method` keyword of
    [`forecast`](@ref) to select the estimation strategy.

## High-level pipeline

```@docs
forecast
```

## ERA and modal filtering

```@docs
era_modal_basis
modal_project
bandpass_filter
```

## Signal processing

```@docs
find_peaks
discard_initial_transient
```

## Recovery-rate estimation

```@docs
estimate_recovery_rate
```

## Polynomial surface

```@docs
fit_polynomial_surface
evaluate
extract_a2
```

## Bifurcation diagram

```@docs
extract_diagram
trace_flutter_boundary
```

## State velocity method

```@docs
forecast_state_velocity
```

## Sampling utilities

```@docs
lhs_samples
amplitude_grid
```

## Diagnostics

```@docs
convergence_check
check_condition
```
