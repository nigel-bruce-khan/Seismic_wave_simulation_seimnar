module GigaFV
using WriteVTK
using Printf
using Logging
using LinearAlgebra
import YAML

include("configuration.jl")
include("basis.jl")
include("equations.jl")
include("grid.jl")
include("kernels/surface.jl")
include("kernels/time.jl")
include("plotters.jl")
include("error_writer.jl")
include("global_matrices.jl")

"""
    evaluate_rhs(eq, scenario, filter, globals, du, dofs, grid)

Evalutes the right-hand-side of the equation `eq` for 
scenario `scenario`, with filter `filter`, 
collection of global matrices `globals`, update
`du`, degrees of freedom `dofs` and grid `grid`.

Updates `du` in place.
"""
function evaluate_rhs(eq, scenario, globals, du, dofs, grid)
    buffers_face = Array{BuffersFaceIntegral,1}(undef, Threads.nthreads())
    for i=1:Threads.nthreads()
        buffers_face[i] = BuffersFaceIntegral(get_ndofs(eq))
    end

    nvar = get_ndofs(eq)
    du .= 0.0

    # compute max eigenvalue per thread, reduce later
    maxeigenvals = Array{Float64,1}(undef, Threads.nthreads())
    maxeigenvals .= -Inf

    # Evaluate all fluxes, compute eigenvalues
    Threads.@threads for i in eachindex(grid.cells)
        @views cell = grid.cells[i]
        @views data = dofs[:, cell.dataidx]
        @views flux = grid.flux[:, :,cell.dataidx]

        cureigenval = max(
            max_eigenval(eq, dofs, 1),
            max_eigenval(eq, dofs, 1)
        )
        maxeigenvals[Threads.threadid()] = max(maxeigenvals[Threads.threadid()], cureigenval)

        evaluate_flux(eq, data, flux)
    end
    Threads.@threads for i in eachindex(grid.cells)
        @views buffer_face = buffers_face[Threads.threadid()]

        @views cell = grid.cells[i]
        @views data = dofs[:, cell.dataidx]
        @views flux = grid.flux[:,:,cell.dataidx]

        # Here we also need to compute the maximum eigenvalue of each cell
        # and store it for each cell (needed for timestep restriction later!)
        faces = (left, top, right, bottom)
        # faces = instances(Face)
        for (i, neigh) in enumerate(cell.neighbors)
            @views dofsneigh = dofs[:,neigh.dataidx]
            @views fluxneigh = grid.flux[:,:,neigh.dataidx]

            isboundary = cell.facetypes[i] == boundary
            if isboundary
                normalidx = globals.normalidxs[faces[i]]

                # For boundary cells, we operate directly on the neighbour data
                evaluate_boundary(eq, scenario, faces[i], normalidx, data, buffer_face.boundary_state)
                evaluate_flux(eq, buffer_face.boundary_state, buffer_face.boundary_flux)
            end
            @views evaluate_face_integral(eq, 
                globals,
                data,
                if isboundary buffer_face.boundary_state else dofsneigh end,
                flux,
                if isboundary buffer_face.boundary_flux else fluxneigh end,
                buffer_face,
                cell,
                faces[i],
                du[:,cell.dataidx])
            
        end

        @views du[:,cell.dataidx] .= du[:,cell.dataidx] / volume(cell)
    end
    
    grid.maxeigenval = maximum(maxeigenvals)
    @show grid.maxeigenval
end

"""
    main(configfile::String)

Runs a FV-simulation with configuration from `configfile`.
"""
function main(configfile::String)
    config = Configuration(configfile)
    eq = make_equation(config)
    scenario = make_scenario(config)
    grid = make_grid(config, eq, scenario)
    integrator = make_timeintegrator(config, grid)
    globals = GlobalMatrices(2)

    filename = "output/plot"

    # Init everything
    for cell in grid.cells
        @views interpolate_initial_dofs(eq, scenario, grid.dofs[:,cell.dataidx],cell)
    end

    plotter = VTKPlotter(eq, scenario, grid, filename)

    grid.time = 0
    timestep = 0
    next_plotted = config.plot_start
     
    while grid.time < config.end_time
        if timestep > 0
            time_start = time()
            dt = config.cellsize[1] * config.courant * 1/grid.maxeigenval
            # Only step up to either end or next plotting
            dt = min(dt, next_plotted-grid.time, config.end_time - grid.time)
            @assert dt > 0

            @info "Running timestep" timestep dt grid.time
            step(integrator, grid, dt) do du, dofs, time
                evaluate_rhs(eq, scenario, globals, du, dofs, grid)
            end

            grid.time += dt
            time_end = time()
            time_elapsed = time_end - time_start
            @info "Timestep took" time_elapsed
        else
            # Compute initial eigenvalue (needed for dt)
            grid.maxeigenval = -1
            for cell in grid.cells
                @views celldata = grid.dofs[:,cell.dataidx]
                for normalidx=1:2
                    cureigenval = max_eigenval(eq, celldata, normalidx)
                    grid.maxeigenval = max(grid.maxeigenval, cureigenval)
                end
            end
        end

        if abs(grid.time - next_plotted) < 1e-10
            @info "Writing output" grid.time
            plot(plotter)
            next_plotted = grid.time + config.plot_step
        end
        timestep += 1
    end
    save(plotter)
    if is_analytical_solution(eq, scenario)
        evaluate_error(eq, scenario, grid, grid.time)
    end
end

end
