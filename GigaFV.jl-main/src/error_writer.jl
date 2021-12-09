"""
    evaluate_error(eq::Equation, scenario::Scenario, grid::Grid, t::Float64)

Prints out the L1/L2/L∞ integral norms of the error for
equation `eq`, scenario (with analytical solution) `scenario`, grid
`grid` at time `t`.
Error is defined as difference between analytical solution and numerical
approximation.
"""
function evaluate_error(eq::Equation, scenario::Scenario, grid::Grid, t::Float64)
    @assert is_analytical_solution(eq, scenario)
    l1_error = similar(grid.dofs, get_ndofs(eq))
    l2_error = similar(l1_error)
    linf_error = similar(l1_error)
    l1_error .= 0
    l2_error .= 0
    linf_error .= 0


    # TODO "Project solution to basis" = comp. average
    # Compare cell average with exact average
    
    #dg_order = size(grid.basis,1)
    dg_order = 1
    quad_order = 3
    quadpoints, quadweights = get_quadpoints(quad_order)
    #quad_order = max(dg_order + 2, 6)
    #quad_basis = Basis(quad_order, grid.basis.dimensions)
    number_of_quadpoints = length(quadpoints)

    for cell in grid.cells
        # Constant cell averages -> only one value
        for var=1:get_ndofs(eq)
            @views value = grid.dofs[var, cell.dataidx]
            for i=1:number_of_quadpoints
                for j=1:number_of_quadpoints
                    x = quadpoints[i]
                    y = quadpoints[j]
                    quadweight = quadweights[i] * quadweights[j]
                    weight = volume(cell) * quadweight
                    pxg, pyg = globalposition(cell, (x,y))
                    analytical = get_initial_values(eq, scenario, (pxg,pyg), t=t)
                    error = abs(value - analytical[var])
                    l1_error[var] += weight * error
                    l2_error[var] += weight * error^2
                    linf_error[var] = max(linf_error[var], error)
                end
            end
        end
    end

    l2_error .= sqrt.(l2_error)
    println("Errors for each variable.")
    println("Var\tL1\tL2\tL∞")
    for i=1:get_ndofs(eq)
        var_name = string(get_variable_name(eq,i))
        @printf "%s\t%15.6e\t%15.6e\t%15.6e\n" var_name l1_error[i] l2_error[i] linf_error[i]
    end

end