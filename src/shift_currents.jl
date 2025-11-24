"""
    `shift_current_TRS_ω(ωlist, a, b, c, q, p, η; jdos = false, kws...)`
returns the JDOS (if `jdos = true`) or the BPGE with chosen method (else) at all ω in 
`ωlist` for a given momentum `q` which will be integrated over the BZ.
    See `σab_inter_linear(a, b, p, ωlist; kws...)`.
"""
function shift_current_TRS_ω(ωlist, a, b, c, q, p, η; vsfreq = true, omega = 0.005, jdos = false, T = 0, kws...) 
    if return_kwargs(kws, :ph) == true
        charge_neutrality = 0.
    else
        charge_neutrality = 0.00164448
    end
    if vsfreq
        h = bistritzer_hamiltonian(p, q; kws...) # 1
        # h = bistritzer_hamiltonian(ParamsBM(p, mu = charge_neutrality ), q; kws...) # 2
        ϵs, ψs = eigen(Matrix(h))
        μeff =  p.mu-charge_neutrality
        if jdos == true
            return [sum_nondiag( -f(ϵs, 0, T) .* lorentz(ϵs, ω, η) ) for ω in ωlist]      # 1 Correct
            # return [sum_nondiag( -f(ϵs, μeff, T) .* lorentz(ϵs, ω, η) ) for ω in ωlist] # 2 Correct also
        else  
            mat = shift_current_integrand(ϵs, ψs, p, q, a, b, c, T, 0; kws...) # 1
            # mat = shift_current_integrand(ϵs, ψs, p, q, a, b, c, T, μeff; kws...) # 2
            return [- π/2 * sum_lowerdiag(mat .* (-lorentz(ϵs, ω, η) .- lorentz(ϵs, -ω, η))) 
                for ω in ωlist]  # Units: Å^3 /eV
        end 

    else # route to variation with μ
        integrand_arr = spzeros(ComplexF64,length(ωlist))
        for (i, μ) in enumerate(ωlist)
            h = bistritzer_hamiltonian(ParamsBM(p, mu = μ), q; kws...)
            ϵs, ψs = eigen(Matrix(h))
            mat = injection_current_integrand_noTRS(ϵs, ψs, ParamsBM(p, mu = μ), q, a, b, c, T, 0.0; kws...) 
            integrand_arr[i] = sum_lowerdiag((mat) .* (-lorentz(ϵs, omega, η) .- lorentz(ϵs, -omega, η)))
        end
            return (ħ_ev_s * π/2) .* integrand_arr
    end
end

shift_current_integrand(ϵs, ψs, p, q, a, b, c, T, doping; kws...) =  
    f(ϵs, doping, T) .* Imat_shift(a, b, c, ϵs, ψs, p, q; kws...)

""" returns the symmetrized product sum """
function Imat_shift(a, b, c, ϵs, ψs, p, q; method = :SHIFT, kws...)
    unbounded = (:x, :y)
    bounded_dirs = findall(x -> x ==:z, (a,b,c))
    if  a == :z
        throw(ArgumentError("We only seek in-plane currents since the z 
          direction is bounded"))
    elseif a ∈ unbounded && b ∈ unbounded && c ∈ unbounded
        return Imat_shift_uuu(a, b, c, ϵs, ψs, p, q, method; kws...)
    elseif length(bounded_dirs) == 1
        return Imat_shift_1b(a, b, c, ϵs, ψs, p, method, bounded_dirs[1])
    else throw(ArgumentError(" xzz and yzz tensor components are not implemented")) end
end

"""
 a, b, c correspond to unbounded directions
 """
function Imat_shift_uuu(a, b, c, ϵs, ψs, p, q, method; kws...)
    omega, ra, Δa = omega_r_Δ(ϵs, ψs, d_bistritzer_hamiltonian(p, a))
    r_cov = similar(ra)
    # println(ra)
    if a == b == c #FAST
        r_cov = r_covariant(omega, ra, Δa, ra, Δa)
        return whichBPGE(method, r_cov, ra, r_cov ,ra)
    else 
        _, rb, Δb = omega_r_Δ(ϵs, ψs, d_bistritzer_hamiltonian(p, b; kws...))
        _, rc, Δc = omega_r_Δ(ϵs, ψs, d_bistritzer_hamiltonian(p, c; kws...))
        r_cov_ca = r_covariant(omega, rc, Δc, ra, Δa)
        r_cov_ba = r_covariant(omega, rb, Δb, ra, Δa)
        return whichBPGE(method, r_cov_ca, rb, r_cov_ba, rc)
        
    end
end

function Imat_shift_1b(a, b, c, ϵs, ψs, p, method, index; kws...)
    if index == 1
        return Imat_shift_buu(a, b, c, ϵs, ψs, p, method; kws...)
    elseif index == 2
        return Imat_shift_ubu(a, b, c, ϵs, ψs, p, method; kws...)
    else 
        return Imat_shift_uub(a, b, c, ϵs, ψs, p, method; kws...)
    end
end

"""
 a is the bounded direction (i.e. :z)
 """
function Imat_shift_buu(a, b, c, ϵs, ψs, p, method; kws...)
    zmat = rz(ψs, length(ϵs)) 
    omega, ra, Δa = omega_r_Δ(ϵs, ψs, d_bistritzer_hamiltonian(p, a; kws...))
    _, rb, Δb = omega_r_Δ(ϵs, ψs, d_bistritzer_hamiltonian(p, b; kws...))
    r_cov_ba = r_covariant(omega, rb, Δb, ra, Δa)
    z_cov_ca = z_covariant(zmat, ra)     
    return whichBPGE(method, z_cov_ca, rb, r_cov_ba, zmat)
end

""" 
b is the bounded direction (i.e. :z)
sigma_uzu = -C ∫_k ∑_nm fnm * Im[znm;a rnm^c + r^znm;a znm] *δ(ω - ω0)
"""
function Imat_shift_ubu(a, b, c, ϵs, ψs, p, method; kws...)
    zmat = rz(ψs, length(ϵs)) 
    omega, ra, Δa = omega_r_Δ(ϵs, ψs, d_bistritzer_hamiltonian(p, a; kws...))
    _, rc, Δc = omega_r_Δ(ϵs, ψs, d_bistritzer_hamiltonian(p, c; kws...))
    r_cov_ca = r_covariant(omega, rc, Δc, ra, Δa)
    z_cov_a = z_covariant(zmat, ra)     
    return whichBPGE(method, zmat, r_cov_ca, rc, z_cov_a)
end
"""
 c is the bounded direction (i.e. :z)
 """
function Imat_shift_uub(a, b, c, ϵs, ψs, p, method; kws...)
    zmat = rz(ψs, length(ϵs)) 
    omega, ra, Δa = omega_r_Δ(ϵs, ψs, d_bistritzer_hamiltonian(p, a; kws...))
    _, rb, Δb = omega_r_Δ(ϵs, ψs, d_bistritzer_hamiltonian(p, b; kws...))
    r_cov_ba = r_covariant(omega, rb, Δb, ra, Δa)
    z_cov_a = z_covariant(zmat, ra)     
    return whichBPGE(method, z_cov_a, rb, r_cov_ba, zmat)
end