using Combinatorics
"""
oc_band is the first occupied level
unoc_band is the first unoccupied level
num_band_subset is the max number of bands considered 
        so transition between num_band_subset/2 bands
        can happen.
"""
struct Transitions
    oc_band::Int64
    unoc_band::Int64
end

transitions(trans_array::Array) = Transitions(trans_array[1], trans_array[2])
transitions(trans_array::Tuple) = Transitions(trans_array[1], trans_array[2])

struct Transitions_Array
    involved_bands::Int64
end

function index_convertor(transition::Transitions, ϵs)
    idx = findfirst(x -> x ≥ 0, ϵs)
    return idx -transition.oc_band, (transition.unoc_band-1) + idx
end

function index_convertor(transition::Tuple, ϵs)
    idx = findfirst(x -> x ≥ 0, ϵs)
    return idx -transition[1], (transition[2]-1) + idx
end

function shift_current_TRS(trans::Union{Transitions, Transitions_Array}, a, b, c, p, ωlist; 
        valley_resolved = false, kws...) 
    # println("Mass method: ", return_kwargs(kws, :mass_term))  
    return nonlinear_terms_TRS(trans, a, b, c, p, ωlist, :SHIFT, valley_resolved; kws...)    
end

function nonlinear_terms_TRS(trans::Transitions_Array, a, b, c, p, ωlist, method::Symbol, valley_resolved::Bool;
    save = true, η = 10^-3, evals = 2, kws...)
    part = ifelse(method == :SHIFT, :real, :imag) # injection current not yet implemented
    println("Number of iterations: ", evals)
 
    # Integral function
    integralfunc(valley; kws...) = integral_nonlinear_substitution_fast(trans, ωlist, a, b, c, ParamsBM(p, ν = valley),
                η, evals, method; kws...)
    # Compute
    trans_mat_pos, trans_indices =  integralfunc(p.ν; kws...)
    trans_mat_neg, _ = integralfunc(-p.ν; kws...)
    # println("Integral values: ", "freq: ", ωlist, " values pos:" , trans_mat_pos, "values neg", trans_mat_neg)
    # Store
    saveconductivity(trans, p, ωlist, trans_mat_pos, trans_mat_neg, trans_indices, η, evals, a, b, c, part, mass_term = :none; kws...)
end

function nonlinear_terms_TRS(trans::Transitions, a, b, c, p, ωlist, method::Symbol, valley_resolved::Bool;
    save = true, η = 10^-3, evals = 2, kws...)
    println("Number of iterations: ", evals)
    vals = zeros(Float64, length(ωlist)); errs = similar(vals)
    # Integral function
    integralfunc(valley; kws...) = integral_nonlinear_substitution_fast(trans, ωlist, a, b, c, ParamsBM(p, ν = valley),
                η, evals, method; kws...)

    vals, errs =  integralfunc(p.ν; kws...)
    while isnan(vals[1]) == true
        print("NaN encountered, retrying...")
        vals, errs = integralfunc(p.ν; kws...)
    end
    println("success!")

    vals2 = integralfunc(-p.ν; kws...)[1]
    while !isa(vals2[1] , Number) == true
        print("NaN encountered, retrying...")
        vals2 = integralfunc(-p.ν; kws...)[1]
    end
    println("success!")
    println("Integral values: ", "freq: ", ωlist, " values pos:" , vals, "values neg", vals2)

    saveconductivity(p, ωlist, vals, vals2, errs, η, evals, a, b, c, 
        ifelse(method == :SHIFT, :real, :imag); kws...)
    saveconductivity(p, ωlist, vals+vals2, errs, η, evals, a, b, c, 
        ifelse(method == :SHIFT, :real, :imag); kws...)
end

function integral_nonlinear_substitution_fast(trans::Transitions_Array, ωlist::Array, a, b, c, p, η, evals, method; 
    focus = :none, kws...)
    trans_indices = vec(collect(Base.Iterators.product(Base.Iterators.repeated(1:trans.involved_bands, 2)...)))
    trans_mat =  Matrix{Float64}(undef, length(trans_indices), length(ωlist))
    integralfunc(transition; kws...) = integral_nonlinear_substitution_fast(transition, 
        ωlist, a, b, c, p, η, evals, method; focus = :none, kws...)
    println("Number of transitions computed is:", length(trans_indices))
    for k in 1:length(trans_indices)
        println("Computed ", k, " of ", length(trans_indices), " transitions")
        trans_mat[k,:] = integralfunc(transitions(trans_indices[k]); kws...)[1]
        while isnan(trans_mat[k, 1]) == true # change
            print("NaN encountered, retrying...")
            trans_mat[k,:] = integralfunc(transitions(trans_indices[k]); kws...)[1]
        end
    end
    println("success!")
    return trans_mat, trans_indices
end

function integral_nonlinear_substitution_fast(trans::Transitions, ωlist::Array, a, b, c, p, η, evals, method; 
    focus = :none, kws...)
    M, xmin, xmax = int_boundaries(p)
    if focus == :none
        xminp = xmin
        xmaxp = [xmax[1]/2, xmax[2]]
    elseif focus == :Γ 
        xminp = [-xmax[1]/2, xmax[1]/2]
        xmaxp = [xmin[2]/2, xmax[2]/2]
    else nothing end 
    check_method(method)
    integral_method = ifelse(method == :SHIFT, shift_current_TRS_ω, injection_current_TRS_ω)
    part = real
    integrand(q) = part(integral_method(trans, ωlist, a, b, c, q, p, η; kws...))      
    val = zeros(eltype(ωlist))
    err = zeros(eltype(ωlist))
    quad_func = (x,v) -> v[:] = integrand(x)
    val, err = hcubature(length(ωlist), quad_func, xminp, xmaxp; reltol = 1e-8, abstol=0, maxevals=evals) 
    angstroms_to_nm = 1/10   
    spin_dof = 2                                  
    cnst =  spin_dof * C * ħ_ev_s * angstroms_to_nm * (1/(2*pi*p.a0))^2
    return val * cnst, err
end

"""
    `shift_current_TRS_ω(ωlist, a, b, c, q, p, η; jdos = false, kws...)`
returns the JDOS (if `jdos = true`) or the BPGE with chosen method (else) at all ω in 
`ωlist` for a given momentum `q` which will be integrated over the BZ.
    See `σab_inter_linear(a, b, p, ωlist; kws...)`.
"""
function shift_current_TRS_ω(transition::Transitions, ωlist, a, b, c, q, p, η, T = 0; jdos = false, kws...) 
    h = bistritzer_hamiltonian(p, q; kws...)
    ϵs, ψs = eigen(Matrix(h))
    i, j = index_convertor(transition, ϵs)
    if jdos == true
        return [( -f(ϵs, 0, T) .* lorentz(ϵs, ω, η) )[j,i] for ω in ωlist] 
        #return [sum_nondiag( -f(ϵs, 0, T) .* lorentz(ϵs, ω, η) ) for ω in ωlist] 
    else 
        mat = shift_current_integrand(ϵs, ψs, p, q, a, b, c, T, 0; kws...)
        return [-π/2 * (mat .* (-lorentz(ϵs, ω, η) .- lorentz(ϵs, -ω, η)))[j, i] 
            for ω in ωlist]  # Units: Å^3 /eV
    end 
end

function injection_current_TRS_ω(transition::Transitions, ωlist, a, b, c, q, p, η, T = 0; jdos = false, kws...) 
    h = bistritzer_hamiltonian(p, q; kws...)
    ϵs, ψs = eigen(Matrix(h))
    i, j = index_convertor(transition, ϵs)
    if jdos == true
        return [( -f(ϵs, 0, T) .* lorentz(ϵs, ω, η) )[j,i] for ω in ωlist] 
        #return [sum_nondiag( -f(ϵs, 0, T) .* lorentz(ϵs, ω, η) ) for ω in ωlist] 
    else 
        mat = injection_current_integrand(ϵs, ψs, p, q, a, b, c, T, 0; kws...)
        return [-π/2 * (mat .* (-lorentz(ϵs, ω, η) .- lorentz(ϵs, -ω, η)))[j, i] 
            for ω in ωlist]  # Units: Å^3 /eV
    end 
end

# We could have implemented a ShiftInvert(ArnoldiMethod(mat::SparseArrays)) to get the central eigenvalues
# but seamed overkill. An alternative is to get the sqrt(eigen), eiges = eigen (h^2) more allocations and
# generally slower.


# # More intelligent approach

# function integral_nonlinear_substitution_fast(trans::Transitions_Array, ωlist::Array, a, b, c, p, η, evals, method; 
#     focus = :none, kws...)
#     trans_indices = vec(collect(Base.Iterators.product(Base.Iterators.repeated(1:trans.involved_bands, 2)...)))
#     trans_mat =  Matrix{Float64}(undef, length(trans_indices), length(ωlist))
#     integralfunc(transition; kws...) = integral_nonlinear_substitution_fast(transition, 
#         ωlist, a, b, c, p, η, evals, method; focus = :none, kws...)
#     println("Number of transitions computed is:", length(trans_indices))
#     trans_mat = integralfunc(transitions(trans_indices); kws...)[1] #
#     while isnan(trans_mat[k, 1]) == true # change
#         print("NaN encountered, retrying...")
#         trans_mat = integralfunc(transitions(trans_indices); kws...)[1] #
#     end
#     println("success!")
#     return trans_mat, trans_indices
# end



# function integral_nonlinear_substitution_fast(trans::Vector, ωlist::Array, a, b, c, p, η, evals, method; 
#     focus = :none, kws...)
#     M, xmin, xmax = int_boundaries(p)
#     if focus == :none
#         xminp = xmin
#         xmaxp = [xmax[1]/2, xmax[2]]
#     elseif focus == :Γ 
#         xminp = [-xmax[1]/2, xmax[1]/2]
#         xmaxp = [xmin[2]/2, xmax[2]/2]
#     else nothing end 
#     num_trans = length(trans)
#     check_method(method)
#     integral_method = ifelse(method == :SHIFT, shift_current_TRS_ω, injection_current_TRS_ω)
#     part = real
#     integrand(q) = part(integral_method(trans, ωlist, a, b, c, q, p, η; kws...))      
#     val = zeros(eltype(ωlist))
#     err = zeros(eltype(ωlist))
#     quad_func = (x,v) -> v[:] = integrand(x)
#     val, err = hcubature(length(ωlist), quad_func, xminp, xmaxp; reltol = 1e-8, abstol=0, maxevals=evals) 
#     angstroms_to_nm = 1/10   
#     spin_dof = 2                                  
#     cnst =  spin_dof * C * ħ_ev_s * angstroms_to_nm * (1/(2*pi*p.a0))^2
#     return val * cnst, err
# end

# """
#     `shift_current_TRS_ω(ωlist, a, b, c, q, p, η; jdos = false, kws...)`
# returns the JDOS (if `jdos = true`) or the BPGE with chosen method (else) at all ω in 
# `ωlist` for a given momentum `q` which will be integrated over the BZ.
#     See `σab_inter_linear(a, b, p, ωlist; kws...)`.
# """
# function shift_current_TRS_ω(trans_inds::Vector, ωlist, a, b, c, q, p, η, T = 0; jdos = false, kws...) 
#     h = bistritzer_hamiltonian(p, q; kws...)
#     ϵs, ψs = eigen(Matrix(h))
#     for (pos, (i_ind,j_ind) ) in enumerate(trans_indices)
#         i, j = index_convertor(Transitions(i_ind,j_ind), ϵs)
#         push!(mat_positions, CartesianIndex(j,i))
#     end
#     if jdos == true
#         return [sum_nondiag( -f(ϵs, 0, T) .* lorentz(ϵs, ω, η) ) for ω in ωlist] 
#     else 
#         v =  [mat_positions[i] for i in 1:length(mat_positions)]
#         mat = shift_current_integrand(ϵs, ψs, p, q, a, b, c, T, 0; kws...)
#         return [-π/2 * (mat .* (-lorentz(ϵs, ω, η) .- lorentz(ϵs, -ω, η)))[v]
#             for ω in ωlist]  # Units: Å^3 /eV
#     end 
# end
