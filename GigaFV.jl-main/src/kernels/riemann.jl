function telegrapher(eq, dofs, dofsneigh, flux, fluxneigh, dx, normalidx, normalsign, numericalflux)
    if normalidx != 1
        numericalflux .= 0.0
        return
    end

    maxeigenval_left = max_eigenval(eq, dofs, normalidx)
    maxeigenval_right = max_eigenval(eq, dofsneigh, normalidx)
    maxeigenval = max(maxeigenval_left, maxeigenval_right)

    s = TelegraphersShortcuts()
    C = dofs[s.C]
    L = dofs[s.L]
    Il = dofs[s.I]
    Ir = dofsneigh[s.I]
    Vl = dofs[s.V]
    Vr = dofsneigh[s.V]
    Istar = 0.5 .* (Il .+ Ir) + 0.5 * sqrt(C/L) * normalsign * (Vl - Vr)
    Vstar = 0.5 .* (Vl .+ Vr) + 0.5 * sqrt(L/C) * normalsign * (Il - Ir)

    numericalflux[s.I] = 1/L * Vstar
    numericalflux[s.V] = 1/C * Istar
end