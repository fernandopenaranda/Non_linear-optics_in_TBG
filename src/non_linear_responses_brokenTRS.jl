

function injection_current_noTRS(a, b, c, p, list; valley_resolved = false, kws...) 
    println("Mass method: ", return_kwargs(kws, :mass_term))  
    return nonlinear_terms_noTRS(a, b, c, p, list, :INJECTION, valley_resolved; kws...)    
end

function nonlinear_terms_noTRS(a, b, c, p, list, method::Symbol, valley_resolved::Bool;
    save = true, η = 10^-3, evals = 2, kws...)
   println("Number of iterations: ", evals)
   vals = zeros(Float64, length(list))
   errs = similar(vals)
   vals, errs = integral_nonlinear_substitution_noTRS(list, a, b, c, ParamsBM(p, ν = p.ν), η, evals, method; kws...)
   vals2, errs2 = integral_nonlinear_substitution_noTRS(list, a, b, c, ParamsBM(p, ν = -p.ν), η, evals, method; kws...)
       saveconductivity(p, list, vals, vals2, errs, η, evals, a, b, c, :real; kws...)
       saveconductivity(p, list, vals2, vals, errs, η, evals, a, b, c, :real; kws...)
       saveconductivity(p, list, vals+vals2, errs, η, evals, a, b, c, :real; kws...)
    return vals, vals2
end

function integral_nonlinear_substitution_noTRS(ωlist::Array, a, b, c, p, η, evals, method; kws...)
    check_method(method)
    M, xmin, xmax = int_boundaries(p)
    #                                    integral 
    # _________________________________________________________________________________________
    #ifelse(method == :SHIFT, shift_current_noTRS_ω, ifelse(method == :INJECTION, injection_current_noTRS_ω, cd_ω))
    integrand(q) = real(injection_current_noTRS_ω(ωlist, a, b, c, q, p, η; kws...))
    
    val = zeros(eltype(ωlist))
    err = zeros(eltype(ωlist))
    quad_func = (x,v) -> v[:] = integrand(x)
    val, err = hcubature(length(ωlist), quad_func, [-xmax[1]/4, xmin[2]], [xmax[1]/2 - xmax[1]/4, xmax[2]]; 
        reltol = 1e-8, abstol=0, maxevals=evals) 
        #                                    constants 
        # _________________________________________________________________________________________
    angstroms_to_nm = 1/10   
    spin_dof = 2
    if method == :SHIFT || method == :INJECTION                                  
        cnst =  spin_dof * C * ħ_ev_s * angstroms_to_nm * (1/(2*pi*p.a0))^2
        # same constant for shift and injection methods
        #UNITS: μA /(s*V^2) nm (eV s) / eV = nm μA/V^2 
    else 
        cnst = C_cd * spin_dof * angstroms_to_nm * (1/(2*pi*p.a0))^2 
    end
        # the units after this are (nm μA/V s) remember that this is multiplied by  
        #[E_j q_k]  = V/nm [h ω c] =  V/nm [2π ϵ c] = V * Jules/ s #CHECK
    println("val: ", val * cnst)
    return val * cnst, err # the prefactor is the jacobian     
end

"""
    `injection_current_noTRS_ω(ωlist, a, b, c, q, p, η; jdos = false, kws...)`
returns the JDOS (if `jdos = true`) or the injection current due to a real and symmetric 
component of the η tensor. Note that this contribution (summed over the two valleys) is only 
finite when TRS is broken  with chosen method (else) at all ω in 
`ωlist` for a given momentum `q` which will be integrated over the BZ.
    See `σab_inter_linear(a, b, p, ωlist; kws...)`.
    returns Imag(η_abc)
    where η_abc = ∂tJ
"""
function injection_current_noTRS_ω(list, a, b, c, q, p, η; jdos = false, T = 0, vsfreq = true, kws...)  
        if vsfreq
            h = bistritzer_hamiltonian(p, q; kws...)
            ϵs, ψs = eigen(Matrix(h))
            if jdos == true
                return [sum_nondiag(f(ϵs, 0.0, T) .* -lorentz(ϵs, ω, η)) for ω in list] 
            else 
                mat = injection_current_integrand_noTRS(ϵs, ψs, p, q, a, b, c, T, 0.0; kws...)  
                return [ħ_ev_s * π/2 * sum_lowerdiag((mat) .* (-lorentz(ϵs, ω, η) .- 0*lorentz(ϵs, -ω, η))) 
                    for ω in list]
            end
        else # path to variation with μ
            integrand_arr = spzeros(ComplexF64,length(list))
            ω = 0.005 # 
            for (i, μ) in enumerate(list)
                h = bistritzer_hamiltonian(ParamsBM(p, mu = μ), q; kws...)
                ϵs, ψs = eigen(Matrix(h))
                mat = injection_current_integrand_noTRS(ϵs, ψs, ParamsBM(p, mu = μ), q, a, b, c, T, 0.0; kws...) 
                integrand_arr[i] = sum_lowerdiag((mat) .* (-lorentz(ϵs, ω, η) .- lorentz(ϵs, -ω, η)))
            end
                return (ħ_ev_s * π/2) .* integrand_arr
        end
        # This is the ∂tJ = ηijk thus 1/γ does not appear as a prefactor. Units: Å^3 /eV
        # Both Lorentzian should have the same size and the sum must go over n>m.
        # Note however that this could give rise to finite conductivity at ω = 0 due to the broadening in the Lorentzians. 
end


injection_current_integrand_noTRS(ϵs, ψs, p, q, a, b, c, T, doping; kws...) =  
    f(ϵs, doping, T) .* Imat_injection_noTRS(a, b, c, ϵs, ψs, p, q; kws...)

""" returns the symmetrized product sum """
function Imat_injection_noTRS(a, b, c, ϵs, ψs, p, q; method = :SHIFT, kws...)
    unbounded = (:x, :y)
    bounded_dirs = findall(x -> x ==:z, (a,b,c))
    if  a == :z
        throw(ArgumentError("We only seek in-plane currents since the z 
          direction is bounded"))
    elseif a ∈ unbounded && b ∈ unbounded && c ∈ unbounded
        return Imat_injection_uuu_noTRS(a, b, c, ϵs, ψs, p, q, method; kws...)
    else
        throw(ArgumentError(" Tensor component not currently implemented")) 
    end
end

"""
 a, b, c correspond to unbounded directions
 """
function Imat_injection_uuu_noTRS(a, b, c, ϵs, ψs, p, q, method; kws...)
    dha = d_bistritzer_hamiltonian(p, a)
    Δa = Δ(ψs, dha)
    _, ra, Δa = omega_r_Δ(ϵs, ψs, d_bistritzer_hamiltonian(p, a))
    if a == b == c 
        # ra = r(ϵs, ψs, dha)
        # return  real.(ra .* t(ra))
        return Δa .* t(real.(ra .* t(ra))) # Comment this this is the real syymmetric part that gives opposite contribution in different valleys whose sum is only nonzero if TRS is broken
    else
        throw(ArgumentError(" Tensor component not currently implemented"))
    end
end
