"""
    struct BuffersFaceIntegral

Stores all temporary buffers needed during face integrals.
Avoids costly re-allocation.
"""
struct BuffersFaceIntegral
    boundary_state::Array{Float64,1}
    boundary_flux::Array{Float64,2}
    numericalflux::Array{Float64,1}

    function BuffersFaceIntegral(ndofs)
        boundary_state = zeros(ndofs)
        boundary_flux = zeros(2,ndofs)
        numericalflux = zeros(ndofs)
        new(boundary_state, boundary_flux, numericalflux)
    end
end

"""
    local_lax_friedrichs(eq, dofs, dofsneigh, flux, fluxneigh, dx, normalidx, normalsign, numericalflux)

Computes the local Lax-Friedrichs (or Rusanov) numerical flux for degrees of freedom
`dofs`, degrees of freedoms of the neighbor `dofsneigh`, flux `flux`, flux of neighbor
`fluxneigh`, cellsize `dx`.
All quantities are assumed to be represented by a basis on the reference line.
The face is parametrized by a index `normalidx`, where 1 stands for a face in x-direction
and 2 for a face in y-direction.
The sign of the outer normal of the face is given by `normalsign`.

The numerical flux is stored in `numericalflux`.
Method also returns the maximal eigenvalue.
"""
function local_lax_friedrichs(eq, dofs, dofsneigh, flux, fluxneigh, dx, normalidx, normalsign, numericalflux)
    maxeigenval_left = max_eigenval(eq, dofs, normalidx)
    maxeigenval_right = max_eigenval(eq, dofsneigh, normalidx)
    maxeigenval = max(maxeigenval_left, maxeigenval_right)
    
    s = ElasticWaveShortcuts()
    oxl = dofsneigh[s.ox]
    oxr = dofs[s.ox]
    oyl = dofsneigh[s.oy]
    oyr = dofs[s.oy]
    oxyl = dofsneigh[s.oxy]
    oxyr = dofs[s.oxy]
    ul = dofsneigh[s.u]
    ur = dofs[s.u]
    vl = dofsneigh[s.v]
    vr = dofs[s.v]
    lam = 2.2
    mu = 1.3
    rho = 1.2
    cp = sqrt((lam + 2*mu)/rho)   
    cs = sqrt(mu/rho)
    xterm3 = oyl - lam * (oxl/(lam + 2*mu))
	yterm3 = oxl - lam * (- (vl/cp) + (oyl+vl/cp*(lam+2*mu))/(lam+2*mu))
	
	
	xoxstar = (lam+2mu)/2  * ((oxl/(lam+2*mu) + ul/cp) + ((oxr)/(lam+2*mu) - ur/cp)) 
	xoystar = lam/2 * ((oxl/(lam+2*mu) + ul/cp) + ((oxr)/(lam+2*mu) - ur/cp)) + xterm3
	xoxystar = mu/2 * ((oxyr/mu - vr/cs) + (oxyl/mu + vl/cs))
	xustar = cp/2 * ((oxl/(lam+2*mu) + ul/cp) - ((oxr)/(lam+2*mu) - ur/cp)) 
	xvstar = cs/2 * (-(oxyr/mu - vr/cs) + (oxyl/mu + vl/cs))
	
	yoxstar = lam * ((oyl + (vl/cp * (lam + 2*mu)))/(2*(lam + 2*mu)) + (-vr/cp + (oyr + vr*(lam+2*mu)/cp)/(2*(lam+2*mu)))) + yterm3
	yoystar = (lam+2mu) * ((oyl + (vl/cp * (lam + 2*mu)))/(2*(lam + 2*mu)) + (-vr/cp + (oyr + vr*(lam+2*mu)/cp)/(2*(lam+2*mu))))
	yoxystar = mu * ((oxyl + mu*(ul-cs*oxyl/mu)/(mu+cs))/(mu) - (ur - cs*oxyr/mu)/(mu+cs))
	yustar = cs * ((oxyl + mu*(ul-cs*oxyl/mu)/(mu+cs))/(mu) + (ur - cs*oxyr/mu)/(mu+cs))
	yvstar = cp * ((oyl + (vl/cp * (lam + 2*mu)))/(2*(lam + 2*mu)) - (-vr/cp + (oyr + vr*(lam+2*mu)/cp)/(2*(lam+2*mu))))
	
	if normalidx == 1
	
			"""
			numericalflux[s.ox] = (lam+2*mu) * ustar
            numericalflux[s.oy] = lam * ustar
            numericalflux[s.oxy] = mu * vstar
            numericalflux[s.u] = 1/rho * oxstar
            numericalflux[s.v] = 1/rho * oxystar
            
			
            numericalflux[s.ox] = 0
            numericalflux[s.oy] = 0
            numericalflux[s.oxy] = 0
            numericalflux[s.u] = - cp .* ((normalsign .* (oxr .- oxl) .* cp ./ (lam .+ 2 .* mu)) .+ (normalsign .* (oxyr .- oxyl)))
            numericalflux[s.v] = - cs .* ((normalsign .* (oxyr .- oxyl) .* cs ./ mu) .+ (normalsign .* (vr .- vl)))
      		"""
      		
      		numericalflux[s.ox] = (lam+2*mu) * xustar
            numericalflux[s.oy] = lam * xustar
            numericalflux[s.oxy] = mu * xvstar
            numericalflux[s.u] = 1/rho * xoxstar
            numericalflux[s.v] = 1/rho * xoxystar
            
	elseif normalidx == 2
			  """
			numericalflux[s.ox] = 0
            numericalflux[s.oy] = 0
            numericalflux[s.oxy] = 0
            numericalflux[s.u] = 0
            numericalflux[s.v] = 0
            """
			numericalflux[s.ox] = lam * yvstar
            numericalflux[s.oy] = (lam+2*mu) * yvstar
            numericalflux[s.oxy] = mu * yustar
            numericalflux[s.u] = 1/rho * yoxystar
            numericalflux[s.v] = 1/rho * yoystar
            
          
	   		"""
            numericalflux[s.ox] = - lam .* normalsign .* (vr .- vl) 
            numericalflux[s.oy] = - (lam .+ 2 .* mu) .* normalsign .* (vr .- vl) 
            numericalflux[s.oxy] = - mu .* normalsign .* (ur .- ul) 
            numericalflux[s.u] = - ((cs .* cs) ./ mu) .* normalsign .* (oxyr .- oxyl) 
            numericalflux[s.v] = - ((cp .* cp) ./ (lam .+ 2 .* mu)) .* normalsign .* (oyr .- oyl) 
            """
	end
	
	
	
end

"""
    evaluate_face_integral(eq,globals, dofs, dofsneigh, flux, fluxneigh, buffers, cell, face, celldu)


Computes face integrals for equation `eq`, cell `cell` on face `face`.
Buffers are passed in `buffers`.
Result is stored in `celldu`.
"""
function evaluate_face_integral(eq, globals, dofs, dofsneigh, flux, fluxneigh, buffers, cell, face, celldu)
    normalidx = globals.normalidxs[face]
    normalsign = globals.normalsigns[face]

    # Compute Riemann solver (in normal)
    local_lax_friedrichs(eq, dofs, dofsneigh, flux, fluxneigh, cell.size[1], normalidx, normalsign, buffers.numericalflux)

    # Hack: Set numerical flux for materials to zero
    @views buffers.numericalflux[get_nvars(eq) + 1 : end] .= 0

    celldu .-= normalsign * area(cell) * buffers.numericalflux
end

