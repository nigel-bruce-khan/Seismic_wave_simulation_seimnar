struct Telegraphers <: Equation end
@declare_dofs Telegraphers [:I, :V, :L, :C] 2

function evaluate_flux(eq::Telegraphers, celldofs, cellflux)
    s = TelegraphersShortcuts()
    cellflux .= 0.0
    L = celldofs[s.L]
    C = celldofs[s.C]
    for dimension=1:2
        if dimension == 1
            cellflux[1,s.I] = 1/L * celldofs[s.V]
            cellflux[1,s.V] = 1/C * celldofs[s.I]
        else
            # Only defined in x-dir!
            cellflux[2,1] = 0.0
            cellflux[2,2] = 0.0
        end
    end
end

function max_eigenval(eq::Telegraphers, celldofs, normalidx)
    s = TelegraphersShortcuts()
    L = celldofs[s.L]
    C = celldofs[s.C]
    return sqrt(1/(C*L))    
end


function get_initial_values(eq::Telegraphers, scenario::PlanarWaves, global_position; t=0.0)
    pxg, pyg = global_position
    f1(x) = sin(2*pi * x)
    f2(x) = 0.1 * cos(2*pi * x)

    C = 2.0
    L = 1.5
    c = sqrt(1/(L*C)) # wavespeed
    Z_0 = sqrt(L/C) # impedance

    V = f1(pxg - c*t) + f2(pxg + c*t)
    I = f1(pxg - c*t)/Z_0 - f2(pxg + c*t)/Z_0

    [I, V, L, C]
end

function is_periodic_boundary(eq::Telegraphers, scenario::PlanarWaves)
    true
end

function is_analytical_solution(eq::Telegraphers, scenario::PlanarWaves)
    true
end

