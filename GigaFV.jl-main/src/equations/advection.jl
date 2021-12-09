struct Advection <: Equation end
@declare_dofs Advection [:rho1, :rho2, :rho3]

struct PlanarWaves <: Scenario
end

function is_periodic_boundary(equation::Advection, scenario::PlanarWaves)
    true
end

function get_initial_values(eq::Advection, scenario::PlanarWaves, global_position; t=0.0)
    pxg, pyg = global_position


    # Solution is ρ(x,t) = ρ(x - at, 0)
    # So we just adjust the global_position

    # Takes periodic boundary conditions into account
    # Not elegant, but works ;)
    function adj_coord(x)
        while (x < 0)
            x += 1.0
        end
        while (x > 1.0)
            x -= 1.0
        end
        x
    end

    pxg = adj_coord(pxg - t)
    pyg = adj_coord(pyg - t)

    # And then use the initial condition.
    k = (2 * pi) # wave number
    val1 = sin(k * (pxg + pyg))
    val2 = sin(k * (pyg))
    val3 = 1
    [val1, val2, val3]
end

function is_analytical_solution(equation::Advection, scenario::PlanarWaves)
    true
end


function evaluate_flux(eq::Advection, celldofs, cellflux)
    s = AdvectionShortcuts()
    cellflux .= 0.0
    for pointidx=1:size(celldofs, 1)
        for dimension=1:2
            if dimension == 1
                cellflux[1,1] = 1 * celldofs[s.rho1]
                cellflux[1,2] = 1 * celldofs[s.rho2]
                cellflux[1,3] = 1 * celldofs[s.rho3]
            else
                @assert(dimension == 2)
                cellflux[2,1] = 1 * celldofs[s.rho1]
                cellflux[2,2] = 1 * celldofs[s.rho2]
                cellflux[2,3] = 1 * celldofs[s.rho3]
            end
        end
    end
end

function max_eigenval(eq::Advection, celldata, normalidx)
    1.0
end
