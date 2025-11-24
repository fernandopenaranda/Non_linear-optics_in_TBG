
"""
    `injection_current_TRS_ω(ωlist, a, b, c, q, p, η; jdos = false, kws...)`
returns the JDOS (if `jdos = true`) or the BPGE with chosen method (else) at all ω in 
`ωlist` for a given momentum `q` which will be integrated over the BZ.
    See `σab_inter_linear(a, b, p, ωlist; kws...)`.
    returns Imag(η_abc)
    where η_abc = ∂tJ
"""
function injection_current_TRS_ω(ωlist, a, b, c, q, p, η; jdos = false, vsfreq = true, T = 0, kws...)  
    h = bistritzer_hamiltonian(p, q; kws...)
    ϵs, ψs = eigen(Matrix(h))
    if  vsfreq
        if jdos == true
            return [sum_nondiag(f(ϵs, 0.0, T) .* lorentz(ϵs, ω, η)) for ω in ωlist] 
        else 
         mat = injection_current_integrand(ϵs, ψs, p, q, a, b, c, T, 0.0; kws...) 
         return [ħ_ev_s * π/2 * sum_lowerdiag(imag(mat) .* (lorentz(ϵs, ω, η) .- lorentz(ϵs, -ω, η)))  #this is the ∂tJ = ηijk thus 1/γ does not appear as a prefactor
            for ω in ωlist]  # Units: Å^3 /eV
         # Both Lorentzian should have the same size and the sum must go over n>m.
         # Note however that this could give rise to finite conductivity at ω = 0 due to the
         # broadening in the Lorentzians. 
        end
    else # route to variation with μ
        integrand_arr = spzeros(ComplexF64,length(list))
        ω = 0.005 # CHANGE BY HAND
        for (i, μ) in enumerate(list)
            h = bistritzer_hamiltonian(ParamsBM(p, mu = μ), q; kws...)
            ϵs, ψs = eigen(Matrix(h))
            mat = injection_current_integrand_noTRS(ϵs, ψs, ParamsBM(p, mu = μ), q, a, b, c, T, 0.0; kws...) 
            integrand_arr[i] = sum_lowerdiag((mat) .* (-lorentz(ϵs, ω, η) .- lorentz(ϵs, -ω, η)))
        end
            return (ħ_ev_s * π/2) .* integrand_arr
    end
end
    
injection_current_integrand(ϵs, ψs, p, q, a, b, c, T, doping; kws...) =  
    f(ϵs, doping, T) .* Imat_injection(a, b, c, ϵs, ψs, p, q; kws...)

""" returns the symmetrized product sum """
function Imat_injection(a, b, c, ϵs, ψs, p, q; method = :SHIFT, kws...)
    unbounded = (:x, :y)
    bounded_dirs = findall(x -> x ==:z, (a,b,c))
    if  a == :z
        throw(ArgumentError("We only seek in-plane currents since the z 
          direction is bounded"))
    elseif a ∈ unbounded && b ∈ unbounded && c ∈ unbounded
        return Imat_injection_uuu(a, b, c, ϵs, ψs, p, q, method; kws...)
    elseif length(bounded_dirs) == 1
        return Imat_injection_1b(a, b, c, ϵs, ψs, p, method, bounded_dirs[1]; kws...)
    else throw(ArgumentError(" xzz and yzz tensor components are not implemented")) end
end


"""
 a, b, c correspond to unbounded directions
 """
function Imat_injection_uuu(a, b, c, ϵs, ψs, p, q, method; kws...)
    dha = d_bistritzer_hamiltonian(p, a)
    Δa = Δ(ψs, dha)
    if a == b == c 
        ra = r(ϵs, ψs, dha)
        return Δa .* t(1im .*imag.(ra .* t(ra)))   
    elseif b == c
        rb = r(ϵs, ψs, d_bistritzer_hamiltonian(p, b; kws...)) 
        return Δa .* t(1im .*imag.(rb .* t(rb)))
    else
        rb = r(ϵs, ψs, d_bistritzer_hamiltonian(p, b; kws...))
        rc = r(ϵs, ψs, d_bistritzer_hamiltonian(p, c; kws...))
        return Δa .* t(1im .*imag.(rc .* t(rb)))
    end
end

function Imat_injection_1b(a, b, c, ϵs, ψs, p, method, index; kws...)
    if index == 1
        return throw(ArgumentError("Jz is not allowed")) #Imat_injection_buu(a, b, c, ϵs, ψs, p, method)
    elseif index == 2
        return Imat_injection_ubu(a, b, c, ϵs, ψs, p, method; kws...)
    else 
        return Imat_injection_uub(a, b, c, ϵs, ψs, p, method; kws...)
    end
end

""" 
b is the bounded direction (i.e. :z)
sigma_uzu = -C ∫_k ∑_nm fnm * Im[zmn rnm^c]
"""
function Imat_injection_ubu(a, b, c, ϵs, ψs, p, method; kws...)
    dha = d_bistritzer_hamiltonian(p, a; kws...)
    Δa = Δ(ψs, dha)
    rb = rz(ψs, length(ϵs))
    if a == c
        return 1im .* Δa .* t(imag.(r(ϵs, ψs, dha) .* t(rb)))
    else
        dhc = d_bistritzer_hamiltonian(p, c; kws...)
        return 1im .* Δa .* t(imag.(r(ϵs, ψs, dhc) .* t(rb)))
    end
end

"""
 c is the bounded direction (i.e. :z)
"""
function Imat_injection_uub(a, b, c, ϵs, ψs, p, method; kws...)
    dha = d_bistritzer_hamiltonian(p, a; kws...)
    Δa = Δ(ψs, dha)
    rc = rz(ψs, length(ϵs))
    if a == b
        return 1im .* Δa .* t(imag.(rc .* t(r(ϵs, ψs, dha))))
    else
        dhb = d_bistritzer_hamiltonian(p, b; kws...)
        return 1im .* Δa .* t(imag.(rc .* t(r(ϵs, ψs, dhb))))
    end
end