struct ElasticWave <: Equation end
@declare_dofs ElasticWave [:ox, :oy, :oxy, :u, :v, :lam, :mu, :rho] 3

function evaluate_flux(eq::ElasticWave, celldofs, cellflux)
    
    s = ElasticWaveShortcuts()
    cellflux .= 0.0
   
    for dimension=1:2
        if dimension == 1
            cellflux[1,1] = -1 * celldofs[s.u] * (celldofs[s.lam] + (2* celldofs[s.mu]))
            cellflux[1,2] = -1 * celldofs[s.u] * celldofs[s.lam] 
            cellflux[1,3] = -1 * celldofs[s.v] * celldofs[s.mu]
            cellflux[1,4] = -1 * celldofs[s.ox] / celldofs[s.rho] 
            cellflux[1,5] = -1 * celldofs[s.oxy] / celldofs[s.rho]
        else
            @assert(dimension == 2)
            cellflux[2,1] = -1 * celldofs[s.v] * celldofs[s.lam]
            cellflux[2,2] = -1 * celldofs[s.v] * (celldofs[s.lam] + (2* celldofs[s.mu])) 
            cellflux[2,3] = -1 * celldofs[s.u] * celldofs[s.mu]
            cellflux[2,4] = -1 * celldofs[s.oxy] / celldofs[s.rho] 
            cellflux[2,5] = -1 * celldofs[s.oy] / celldofs[s.rho]
        end
    end
end

function max_eigenval(eq::ElasticWave, celldofs, normalidx)
    s = ElasticWaveShortcuts()
    lam = celldofs[s.lam]
    mu = celldofs[s.mu]
    rho = celldofs[s.rho]
    return sqrt((lam + 2*mu)/rho)    
end


function get_initial_values(eq::ElasticWave, scenario::PlanarWaves, global_position; t=0.0)
    pxg, pyg = global_position
    
    f1(x) = 2 * pi * sin(2 * pi * x)
   
    
    lam = 2.2
    mu = 1.3
    rho = 1.2
    cp = sqrt((lam + 2 * mu) / rho)   
    cs = sqrt(mu / rho) 
    
    """
    # only x-side
    ox = (lam+2*mu) * (-f1(pxg-cp*t))
    oy = 0
    oxy = 0
    u = cp * f1(pxg-cp*t)
    v = 0
    """
    
    """
    # only y-side
	ox = 0
    oy = (lam + 2 * mu) * (-f1(pyg-cp*t))
    oxy = 0
    u = 0
    v = cp * f1(pyg-cp*t)
    """
    
    ox = (lam + 2 * mu) * (- f1(pxg + pyg - cp * t) + f1(pxg + pyg - cs * t)) + lam * (- f1(pxg + pyg - cp * t) - f1(pxg + pyg - cs * t))
    oy = lam * (- f1(pxg + pyg - cp * t) + f1(pxg + pyg - cs * t)) + (lam + 2 * mu) * (- f1(pxg + pyg - cp * t) - f1(pxg + pyg - cs * t))
    oxy = mu * (- 2 * f1(pxg + pyg - cp * t))
    u = cp * f1(pxg + pyg - cp * t) - cs * f1(pxg + pyg - cs * t)
    v = cp * f1(pxg + pyg - cp * t) + cs * f1(pxg + pyg - cs * t)
    
    [ox, oy, oxy, u, v, lam, mu, rho]
end

function is_periodic_boundary(eq::ElasticWave, scenario::PlanarWaves)
    true
end

function is_analytical_solution(eq::ElasticWave, scenario::PlanarWaves)
    true
end

