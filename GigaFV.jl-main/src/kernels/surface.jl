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
    oxl = dofs[s.ox]
    oxr = dofsneigh[s.ox]
    oyl = dofs[s.oy]
    oyr = dofsneigh[s.oy]
    oxyl = dofs[s.oxy]
    oxyr = dofsneigh[s.oxy]
    ul = dofs[s.u]
    ur = dofsneigh[s.u]
    vl = dofs[s.v]
    vr = dofsneigh[s.v]
    lam = dofs[s.lam]
    mu = dofs[s.mu]
    rho = dofs[s.rho]
    cp = sqrt((lam + 2*mu)/rho)   
    cs = sqrt(mu/rho)
    xterm3 = oyr - lam * (oxr/(lam + 2*mu))
	yterm3 = oxr - lam * (oyr/(lam + 2*mu))
	
	xoxstar = (lam+2*mu)/2 * ((oxl+oxr)/(lam+2*mu) + normalsign*(ul - ur)/cp)
	xoystar = lam/2 * ((oxl+oxr)/(lam+2*mu) + normalsign*(ul - ur)/cp) + xterm3
	xoxystar = mu/2 * ((oxyr+oxyl)/mu + normalsign*(vl-vr)/cs)
	xustar = cp/2 * (normalsign*(oxl-oxr)/(lam+2*mu) + (ul+ur)/cp)
	xvstar = cs/2 * (normalsign*(oxyl-oxyr)/mu + (vr+vl)/cs)
	
	yoxstar = lam * ((oyl + oyr)/(2*(lam + 2*mu)) + normalsign*(vl-vr)/(2*cp)) + yterm3
	yoystar = (lam+2*mu) * ((oyl + oyr)/(2*(lam + 2*mu)) + normalsign*(vl-vr)/(2*cp))
	yoxystar = mu * (oxyl/mu + normalsign*(ul-ur)/(2*cs) + normalsign*(oxyr-oxyl)/(mu*2))
	yustar = cs/2 * (normalsign*(oxyl-oxyr)/mu + (ul+ur)/cs)
	yvstar = cp * (normalsign*(oyl - oyr)/(2*(lam + 2*mu)) + (vl+vr)/(2*cp))
	

	if normalidx == 1
			"""
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
      		
      	"""
      		numericalflux[s.ox] = (lam+2*mu) * xustar
            numericalflux[s.oy] = 0
            numericalflux[s.oxy] = 0
            numericalflux[s.u] = 1/rho * xoxstar
            numericalflux[s.v] = 0
           """ 
            
         """
          	numericalflux[s.ox] = 0
            numericalflux[s.oy] = 0
            numericalflux[s.oxy] = 0
            numericalflux[s.u] = 0
            numericalflux[s.v] = 0
            """
            
	elseif normalidx == 2
			 """ 
			numericalflux[s.ox] = 0
            numericalflux[s.oy] = 0
            numericalflux[s.oxy] = 0
            numericalflux[s.u] = 0
            numericalflux[s.v] = 0
            """
              """
            numericalflux[s.ox] = 0
            numericalflux[s.oy] = (lam+2*mu) * yvstar
            numericalflux[s.oxy] = 0
            numericalflux[s.u] = 0
            numericalflux[s.v] = 1/rho * yoystar
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

