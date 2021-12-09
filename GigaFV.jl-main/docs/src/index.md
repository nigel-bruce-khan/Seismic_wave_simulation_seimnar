# GigaFV.jl Documentation

```@contents
```

# Main
```@docs
GigaFV.evaluate_rhs
GigaFV.main
```
# I/O
## Configuration
```@docs
GigaFV.Configuration
```

## VTK/Paraview output
```@docs
GigaFV.VTKPlotter
GigaFV.plot
GigaFV.save
```

## Error Writer
```@docs
GigaFV.evaluate_error
```

# Equations
```@docs
GigaFV.Equation
GigaFV.make_equation
GigaFV.Scenario
GigaFV.make_scenario
GigaFV.interpolate_initial_dofs
GigaFV.get_nvars
GigaFV.get_nparams
GigaFV.get_variable_name
GigaFV.is_periodic_boundary
GigaFV.evaluate_boundary
GigaFV.get_initial_values
GigaFV.is_analytical_solution
GigaFV.evaluate_flux
GigaFV.max_eigenval
GigaFV.@declare_dofs
```

# Grid
```@docs
GigaFV.FaceType
GigaFV.Cell
GigaFV.Grid
GigaFV.Face
GigaFV.get_neighbor
GigaFV.make_mesh
GigaFV.make_grid
GigaFV.globalposition
GigaFV.localposition
GigaFV.volume
GigaFV.area
GigaFV.inverse_jacobian
```

## Basis
```@docs
GigaFV.get_quadpoints
```

# Kernels
Many kernels take both global matrices and buffers.
Try to use them to avoid costly re-computations or memory
allocations.

## Global Matrices
```@docs
GigaFV.GlobalMatrices
```

## Time
```@docs
GigaFV.TimeIntegrator
GigaFV.make_timeintegrator
GigaFV.step
GigaFV.ExplicitEuler
GigaFV.SSPRK2
GigaFV.SSPRK3
```
## Surface
```@docs
GigaFV.BuffersFaceIntegral
GigaFV.local_lax_friedrichs
GigaFV.evaluate_face_integral
```

## Index
```@index
```