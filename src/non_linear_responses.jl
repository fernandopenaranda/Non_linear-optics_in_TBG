include("bistritzermodel.jl")
include("shift_currents.jl")
include("injection_currents.jl")
include("atlas_multithread_calculations.jl")
include("transision_resolved_nl_photocurrents.jl")

#module shiftcurrent
using Arpack
using Cubature
using ProgressMeter
using Base.Threads
using Distributed
using Dierckx
using PhysicalConstants
using PhysicalConstants.CODATA2018
using Unitful

const k_B = (PhysicalConstants.CODATA2018.k_B |> u"eV/mK").val
const ħ = PhysicalConstants.CODATA2018.ħ
const e = PhysicalConstants.CODATA2018.e
const C = ((e^3 / ħ^2) |> u"μA/V^2/s").val
const C_cd = ((e^2/ħ) |> u"μA/V").val
const ħ_ev_s = (ħ |> u"eV*s").val


# ------------------------------------------------------------------------------------------
# σabc_inter_nonlinear
# ------------------------------------------------------------------------------------------
"""
    `shift_current_TRS(a, b, c, p, ωlist; save = false, η = 10^-3, evals = 6400)`
returns the shift current in a the BM model. The TRS symmetry ensures that all absortive 
contribution to the shift current comes from the LPGE term in sigma_abc 
(i.e. that which is real and symmetric (b <-> c))
Ref: PHYSICAL REVIEW RESEARCH 2, 012017(R) (2020) generalized to include bounded dimensions.

Presets: 
- `ωlist` list of frequencies `ω` 
-`p::ParamsBM` is a struct that contains the system information.  
 - kws... 
    - `mass_term = :none` dos Santos model with no mass term
                 = `:normal` sigma_z mass in both layers
                 = `:layered`sigma_z mass in a single layer
                 = `:kang_vafek_mass` term 
    - `η` defines the energy broadening of the Lorentzian peaks.
    - `evals` sets the mesh size in the BZ integration
    - `save` establishes whether the gen data should be stored or not automatically 
"""
function shift_current_TRS(a, b, c, p, ωlist; valley_resolved = false, kws...) 
    println("PH: ", return_kwargs(kws, :ph))
    println("Mass method: ", return_kwargs(kws, :mass_term))  
    return nonlinear_terms_TRS(a, b, c, p, ωlist, :SHIFT, valley_resolved; kws...)    
end

"""
    `injection_current_TRS(a, b, c, p, ωlist; save = false, η = 10^-3, evals = 6400)`
returns the dc injection current in a the BM model. 
The TRS symmetry ensures that all absortive
contribution to the injection current  comes from the CPGE term in η_abc 
(i.e. that which is real and symmetric (b <-> c)). 
The contribution to the optical susceptibility is obtained by dividing the result by
 1/(-i ω_Σ), which yields a divergence as ω_Σ = ω - ω = 0, 
 I define γ  = 0.1 meV to be the relaxation rate as in
(DOI: 10.1103/PhysRevResearch.4.013209).

Main Ref: PHYSICAL REVIEW RESEARCH 2, 012017(R) (2020) generalized to include bounded 
dimensions.

Presets: 
- `ωlist` list of frequencies `ω` 
-`p::ParamsBM` is a struct that contains the system information.  
 - kws... 
    - `mass_term = :none` dos Santos model with no mass term
                 = `:normal` sigma_z mass in both layers
                 = `:layered`sigma_z mass in a single layer
                 = `:kang_vafek_mass` term 
    - `η` defines the energy broadening of the Lorentzian peaks.
    - `evals` sets the mesh size in the BZ integration
    - `save` establishes whether the gen data should be stored or not automatically 
"""
function injection_current_TRS(a, b, c, p, ωlist; valley_resolved = false, kws...) 
    println("Mass method: ", return_kwargs(kws, :mass_term))  
    return nonlinear_terms_TRS(a, b, c, p, ωlist, :INJECTION, valley_resolved; kws...)    
end

"""
circular dichroinism Jᵢ = αᵢⱼₖ(ω) Eⱼqₖ. Units of α are  [nm A /V s]
"""
function circular_dichroinism(a, b, c, p, ωlist; valley_resolved = false, kws...) 
    println("Mass method: ", return_kwargs(kws, :mass_term))  
    return nonlinear_terms_TRS(a, b, c, p, ωlist, :CD, valley_resolved; kws...)    
end

function nonlinear_terms_TRS(a, b, c, p, ωlist, method::Symbol, valley_resolved::Bool;
     save = true, η = 10^-3, evals = 2, kws...)
    println("Number of iterations: ", evals)
    vals = zeros(Float64, length(ωlist))
    errs = similar(vals)
    vals, errs = integral_nonlinear_substitution(ωlist, a, b, c, ParamsBM(p, ν = p.ν), η, evals, method; kws...)
    println(vals[1])
    while isnan(vals[1]) == true
        print("NaN encountered, retrying...")
        vals, errs = integral_nonlinear_substitution(ωlist, a, b, c, ParamsBM(p, ν = p.ν), η, evals, method; kws...)
    end
    println("success!")
    if valley_resolved == false
        vals2 = integral_nonlinear_substitution(ωlist, a, b, c, ParamsBM(p, ν = -p.ν),
             η, evals, method; kws...)[1]
        while !isa(vals2[1] , Number) == true
            print("NaN encountered, retrying...")
            vals2 = integral_nonlinear_substitution(ωlist, a, b, c, ParamsBM(p, ν = -p.ν),
            η, evals, method; kws...)[1]
        end
        println("success!")
        saveconductivity(p, ωlist, vals, vals2, errs, η, evals, a, b, c, 
            ifelse(method == :SHIFT, :real, ifelse(method == :INJECTION, :imag, :cd)); kws...)
        saveconductivity(p, ωlist, vals+vals2, errs, η, evals, a, b, c, 
            ifelse(method == :SHIFT, :real, ifelse(method == :INJECTION, :imag, :cd)); kws...)
    else
        saveconductivity(p, ωlist, vals, errs, η, evals, a, b, c, 
            ifelse(method == :SHIFT, :real, ifelse(method == :INJECTION, :imag, :cd)); kws...)
    end
    # vals2, errs2 = integral_nonlinear_substitution(ωlist, a, b, c, ParamsBM(p, ν = -1 ), η, evals, method; kws...)
 end

 """ 
    `integral_nonlinear_substitution(ωlist, a, b, p, η, evals, method)`
computes the integral of the nonlinear optical conductivity
over twice the Brilloin zone (thus we must multiply the result by an 1/2)
It performs an adaptive integral in k space doubling the k mesh till the prescribed
tolerance is achieved using hcubature techniques. Since we are computing the integral 

  -> shift_current_TRS(:x, :x, :x,  ParamsBM(p, θ = 1.12, nmax = 3, ν  = -1), ωlist, η = 10^-3, 
        evals = 1000, jdos = false, mass_term = :layered, valley_resolved = true)

If focus != :none (default) it computes the integral around a high symmetry point
(not at the whole Brilluoin Zone)

"""
function integral_nonlinear_substitution(ωlist::Array, a, b, c, p, η, evals, method, evaluation = :fast; kws...)
    if evaluation == :fast
        return integral_nonlinear_substitution_fast(ωlist, a, b, c, p, η, evals, method; kws...)     # CHANGE TO FAST 
    elseif evaluation == :slow 
        return integral_nonlinear_substitution_slow(ωlist, a, b, c, p, η, evals, method; kws...)
    else throw(ArgumentError("evaluation can only be :slow or :fast"))
    end
end

"""
The integration strategy is h_adaptive cubature from Cubature.jl.
Literally quoting h_cubature documentation [see https://juliahub.com/ui/Packages/General/Cubature]:
"h-adaptive integration works by recursively subdividing the integration domain into smaller and smaller regions, 
applying the same fixed-order (fixed number of points) integration rule within each sub-region and subdividing a region 
if its error estimate is too large. (Technically, we use a Gauss-Kronrod rule in 1d and a Genz-Malik rule in higher dimensions.) 
This is well-suited for functions that have localized sharp features [[as our deltas]] in a portion of the domain, 
because it will adaptively add more points in this region while using a coarser set of points elsewhere." 
The algorithms used in the Cubature.jl package are:
 -> "Vectorization of one dimensional quadrature codes," pp. 230--238 in Numerical Integration. Recent Developments, 
    Software and Applications, G. Fairweather and P. M. Keast, eds., NATO ASI Series C203, Dordrecht (1987),
 -> "Parallel Globally Adaptive Algorithms for Multi-dimensional Integration," J. M. Bull and T. L. Freeman, 
    http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.42.6638 (1994). 
"""
function integral_nonlinear_substitution_fast(ωlist::Array, a, b, c, p, η, evals, method; focus = :none, kws...)
    M, xmin, xmax = int_boundaries(p)
    if focus == :none
        xminp = xmin
        xmaxp = [xmax[1]/2, xmax[2]]
    elseif focus == :Γ 
        xminp = [-xmax[1]/2, xmax[1]/2]
        xmaxp = [xmin[2]/2, xmax[2]/2]
    else nothing end 
    check_method(method)
    integral_method = ifelse(method == :SHIFT, shift_current_TRS_ω, ifelse(method == :INJECTION, injection_current_TRS_ω, cd_ω))
    part = real
    integrand(q) = part(integral_method(ωlist, a, b, c, q, p, η; kws...))       # par with respect to pm k cause same valley contribution
    
    val = zeros(eltype(ωlist))
    err = zeros(eltype(ωlist))
    quad_func = (x,v) -> v[:] = integrand(x)
    # [-xmax[1]/4, xmin[2]], [xmax[1]/2 - xmax[1]/4, xmax[2]]
    val, err = hcubature(length(ωlist), quad_func, [-xmax[1]/4, xmin[2]], [xmax[1]/2 - xmax[1]/4, xmax[2]]; reltol = 1e-8, abstol=0, maxevals=evals)   # non-vectorized method

    # val, err = hcubature(length(ωlist), quad_func, xminp, xmaxp; reltol = 1e-8, abstol=0, maxevals=evals)   # non-vectorized method
    # quad_func_v(x::Array{Float64,2}, v) = for i = 1:size(v,2) ### SLOW DISREGARDED
    #     for j = 1:size(v,1) #omegalist
    #         v[j,i] = integrand(x[:,i])[j]
    #     end
    # end
    # val, err = hcubature_v(length(ωlist), quad_func_v, xminp, xmaxp; reltol = 1e-8, abstol=0, maxevals=evals)       
    angstroms_to_nm = 1/10
    spin_dof = 2
    if method == :SHIFT || method == :INJECTION
        cnst =  spin_dof * C * ħ_ev_s * angstroms_to_nm * (1/(2*pi*p.a0))^2
        # same constant for shift and injection methods
        #UNITS: μA /(s*V^2) nm (eV s) / eV = nm μA/V^2 
    else 
        cnst = C_cd * spin_dof * angstroms_to_nm * (1/(2*pi*p.a0))^2
    end
        # the units after this are (nm μA/V s) remember that this is multiplied by  [E_j q_k]  = V/nm [h ω c] =  V/nm [2π ϵ c] = V * Jules/ s #CHECK
    return val * cnst, err # the prefactor is the jacobian     
    
end

function integral_nonlinear_substitution_slow(ωlist::Array, a, b, c, p, η, evals, method, evaluation = :slow; kws...)
    val = similar(ωlist)
    err = similar(val)
    lock = ReentrantLock()
    i = Threads.Atomic{Int}(0);
    @threads for i in 1:length(ωlist)
        # lock_threads(lock) do
            val[i], err[i] = integral_nonlinear_substitution_slow(ωlist[i], a, b, c, p, η, evals, method; kws...)
        # end
    end
    return val, err 
end

function integral_nonlinear_substitution_slow(ω::Float64, a, b, c, p, η, evals, method, evaluation = :fast; kws...)  
    check_method(method)
    integral_method = ifelse(method == :SHIFT, shift_current_TRS_ω, ifelse(method == :INJECTION, injection_current_TRS_ω, cd_ω))
    part = real
    integrand_fast(q) = part(integral_method(ω, a, b, c, q, p, η; kws...)[1])     
    _ , xmin, xmax = int_boundaries(p)
    val, err = hcubature(integrand_fast, xmin, [xmax[1]/2, xmax[2]], 
        reltol = 1e-8, abstol = 1e-6, maxevals = evals) 
    # [val]: Units: Å^3 /eV 
    angstroms_to_nm = 1/10
    spin_dof = 2
    cnst =  spin_dof * C * ħ_ev_s * angstroms_to_nm * (1/(2*pi*p.a0))^2 # same constant C for shift and injection methods
    return val * cnst, err # the prefactor is the jacobian
    # UNITS: μA /(s*V^2) nm (eV s) / eV = nm μA/V^2 
end

function check_method(method)
    try
        if method == :SHIFT
            nothing # println("Computing the Shift Current for valley ", p.ν, " ....")
        elseif method == :INJECTION
            nothing # println("Computing the Injection Current for valley ", p.ν, " ....")
        else
            throw(ArgumentError("Variable is neither :SHIFT nor :INJECTION"))
        end
    catch e
        return "Error: $(e)"
    end
end

function whichBPGE(method, r_cov, r, r_cov2, r2)
    if method == :SHIFT
        return symmetrized_imagsum(r_cov, r, r_cov2, r2)
    else
        return antisymmetrize_realsum(r_cov, r, r_cov2, r2)
    end
end

function whichBPGE(method, r_cov, r)
    if method == :SHIFT
        return symmetrized_imagsum(r_cov, r)
    else 
        return antysymmetrize_realsum(r_cov, r)
    end
end

symmetrized_imagsum(mat1, mat2)  =  imag(mat1 .* t(mat2) .+ mat2 .* t(mat1))
symmetrized_imagsum(mat1, mat2, mat3, mat4) = imag((mat1 .* t(mat2)) .+ (mat3 .* t(mat4))) 

function antisymmetrized_realsum(mat1, mat2) 
    auxmat = similar(mat1)
    dim = size(mat1, 1)
    for n in 1:dim
        for m in 1:dim
            auxmat[n,m] = real(mat1[n,m] * mat2[m,n] - mat2[n,m] * mat1[m,n])
        end
    end
    return auxmat
end

antisymmetrized_realsum(mat1, mat2, mat3, mat4) = real(mat1 .* t(mat2) .- mat3 .* t(mat4))


#_________________________________________________________________________________________
# position operator
#_________________________________________________________________________________________
r(ϵs, ψs, dh) = -1im .* vel(ψs, dh) ./ Ω(ϵs)                                    # Units: Å  
rz(ψs, dim) = ψs' * (zop(dim) .* ψs)                                            # Units: Å
function zop(dim)                                                                           
    dim2 = div(dim, 2)  
    return  3.3/2 .*  vcat(ones(dim2), -1 .* ones(dim2))       # 3.3 Å = interlayer distance                                            
end

function r_covariant(omega, ra, Δa, rb, Δb)                                    # Units: Å^2
    raa_cov = zeros(ComplexF64, size(ra,1), size(ra,1))#similar(ra)
    dim = size(omega, 1) 
    for n in 1:dim
        for m in 1:dim
            if m != n && n > m 
                aux = 0
                for p in 1:dim
                    if p != n && p != m
                        aux += 1im * (omega[n,p]*ra[n,p]*rb[p,m] - omega[p,m]*rb[n,p]*ra[p,m])
                        # Units: eV * Å * Å      #!! WARN DIVISIONS BY DIFFERENT OMEGAnm
                    else
                        aux += 0
                    end
                end
                raa_cov[n,m] = -((ra[n,m] * Δb[n,m] + rb[n,m] * Δa[n,m]) + aux)/omega[n,m]
                # Units: (Å * Å * eV) / eV
            else 
                nothing
            end
        end
    end
    return raa_cov
end

function z_covariant(zmat, rmat)
    dim = size(zmat, 1)
    mat = similar(rmat)
    @inbounds begin
        for n in 1:dim
            for m in 1:dim
                if m != n && n > m 
                    for p in 1:dim
                        if  p == n 
                            val = zmat[n,p] * rmat[p,m] 
                        elseif p == m 
                            val = - rmat[n,p] * zmat[p, m]
                        else
                            val = zmat[n,p] * rmat[p,m] - rmat[n,p] * zmat[p, m]
                        end
                        mat[n,m] += val    
                    end
                else nothing end
            end
        end
    end
    return -1im .* mat
end

#_________________________________________________________________________________________
# observables
#_________________________________________________________________________________________

lorentz(ϵs, ω, η) = 1/π .* η ./ (δ(ϵs, ω).^2 .+ η^2) 
δ(ϵs, ω) = ω .- Ω(ϵs)
Ω(ϵs) = ϵs .- t(ϵs) + 1e-7 * Diagonal(ϵs) # the last term is just to avoid infinities it'll cancel out later
vel(ψs, dh) = ψs' * dh * ψs # It really is hbar v  # Units: eV Å  

function omega_r_Δ(ϵs, ψs, dh)
    vmat = vel(ψs, dh)
    mat = similar(vmat)
    omega = Ω(ϵs)
    for i in 1:size(vmat, 1)
        for j in 1:size(vmat, 1)
            mat[i, j] = vmat[i, i] - vmat[j, j]
        end
    end
    return omega, -1im .* vmat ./ omega, mat 
end

function Δ(ψs, dh)
    vmat = vel(ψs, dh)
    mat = similar(vmat)
    for i in 1:size(vmat, 1)
        for j in 1:size(vmat, 1)
            mat[i, j] = vmat[i, i] - vmat[j, j]
        end
    end
    return mat
end

"""
fn - fm where f are the Fermi Dirac distribution in the limit of T = 0.
I set p.
"""
f(ϵs, μ, T) = [fn(ϵs[i], μ, T) - fn(ϵs[j], μ, T) for i in 1:length(ϵs), j in 1:length(ϵs)]

fn(ϵn,μ::Float64) = ifelse(ϵn < μ, 1.0, 0.0)

function fn(ϵn, μ, T)
    if T == 0
        return ifelse(ϵn < μ, 1.0, 0.0)
    else
        return 1/(exp((ϵn - μ)/(k_B * T)) + 1)
    end
end
    
#-------------------------------------------------------------------------------------------
# Operators
#-------------------------------------------------------------------------------------------

t(mat) = transpose(mat)
sum_nondiag(mat) = sum(mat) - sum(Diagonal(mat))
sum_lowerdiag(mat::Array) =  sum_lowerdiag(mat, size(mat,1))
sum_lowerdiag(mat, dim) = sum(mat .* (tril(ones(dim, dim), 0) - Diagonal(ones(dim,dim))))

function return_kwargs(kws, name::Symbol)
    kwargs = NamedTuple(kws)
    if haskey(kwargs, name)
        return kwargs[name]
    else  nothing
    end
end

"""
boundaries of the integration function 
"""
function int_boundaries(p)
    ϵ = 0*1e-4
   rcp = rec_vecs(p)
   m = rcp.M
   v_length = norm(rcp.G2 - rcp.G1)
   h_length = norm(rcp.G2 + rcp.G1)
   return m, [-h_length/2 + m[1] + ϵ, -v_length/2 - ϵ], [rcp.κ1[1]*2, v_length/2 - ϵ]#h_length/2 + m[1] + ϵ , v_length/2 - ϵ]
end

function normvolumebz(p)
    M, xmin, xmax = int_boundaries(p)
    return (xmax[1] - xmin[1]) * (xmax[2] - xmin[2])
end

function k_mesh(p, steps)
    M, xmin, xmax = int_boundaries(p)
    vec = Array{Float64, 2}(undef, steps^2 , 2)
    count = 1
    for i in collect(xmin[1]:(xmax[1]-xmin[1])/(steps-1):xmax[1])
        for j in collect(xmin[2]:(xmax[2]-xmin[2])/(steps-1):xmax[2])
            vec[count, 1] = i
            vec[count, 2] = j 
            count += 1
        end
    end
    
    return vec
end

function bzvertices(p)
    M, xmin, xmax = int_boundaries(p)
    v1 = [xmin[1], 0.]
    v2 = [xmax[1]/2, xmax[2]]
    v3 = [xmax[1], 0.]
    v4 = [xmax[1]/2, xmin[2]]
    return v1,v2,v3,v4
end

# ------------------------------------------------------------------------------------------
# JOINT DOS
# -----------------------------------------------------------------------------------------
# JOINT DOS
function jdos(p, ωmin, ωmax, step)
    ωlist = collect(ωmin:step:ωmax)
    return [integral_jdos(ω, p) for ω in ωlist]
end

function integral_jdos(ω, p)
    v_length, h_length, xmin, xmax = int_boundaries(p)    
    integrand(q) = jdos_integrand(ω, q, p, v_length, h_length)
    val, err = hcubature(integrand, xmin, xmax; reltol=1e-8, abstol=0, maxevals=200)
    return val, err
end

function jdos_integrand(ω, q, p, vl, hl)
    m = rec_vecs(p).M
    if shapeBZ(q, vl, hl, m)
        jdos_ω(ω, q, p)
    else
        0.0
    end
end

function jdos_ω(ω, q, p)
    h = bistritzer_hamiltonian(p, q)
    ϵs, _ = eigen(Matrix(h))
    Γ = 1e-4 #
    return sum(Γ ./ ((Ω(ϵs) .- ω).^2 .+ Γ.^2))
end

# ------------------------------------------------------------------------------------------
# DOS
# -----------------------------------------------------------------------------------------

function dos(p, ωmin, ωmax, step)
    ωlist = collect(ωmin:step:ωmax)
    return [integral_dos(ω, p) for ω in ωlist]
end

function integral_dos(ω, p)
    v_length, h_length, xmin, xmax = int_boundaries(p)    
    integrand(q) = dos_integrand(ω, q, p, v_length, h_length)
    val, err = hcubature(integrand, xmin, xmax; reltol=1e-8, abstol=0, maxevals=2000)
    return val, err
end

function dos_integrand(ω, q, p, vl, hl)
    m = rec_vecs(p).M
    if shapeBZ(q, vl, hl, m)
        dos_ω(ω, q, p)
    else
        0.0
    end
end

function dos_ω(ω, q, p)
    h = bistritzer_hamiltonian(p, q)
    ϵs, _ = eigen(Matrix(h))
    Γ = 1e-4 #
    return sum(Γ ./ ((ϵs .- ω).^2 .+ Γ.^2))
end


function is_inside_rhombus(q, a, b, c, d)
    # Calculate the vectors from a to q and from a to b
    vector_aq = q - a
    vector_ab = b - a

    # Calculate the vectors from a to q and from a to d
    vector_ad = d - a

    # Calculate the dot products of vector_aq with vector_ab and vector_ad
    dot_product_ab = dot(vector_aq, vector_ab)
    dot_product_ad = dot(vector_aq, vector_ad)

    # Calculate the lengths of vectors vector_ab and vector_ad
    length_ab = norm(vector_ab)
    length_ad = norm(vector_ad)

    # Calculate the maximum allowed dot product value
    max_dot_product = length_ab * length_ab

    # Check if q is outside the rhombus
    if dot_product_ab < 0.0 || dot_product_ad < 0.0 || dot_product_ab > max_dot_product || dot_product_ad > max_dot_product
        return 0
    else
        return 1
    end
end

#######
function dos(p, ωmin, ωmax, step)
    ωlist = collect(ωmin:step:ωmax)
    return [integral_dos(ω, p) for ω in ωlist]
end

function integral_dos(ω, p)
    v_length, h_length, xmin, xmax = int_boundaries(p)    
    integrand(q) = dos_integrand(ω, q, p, v_length, h_length)
    val, err = hcubature(integrand, xmin, xmax; reltol=1e-8, abstol=0, maxevals=2000)
    return val, err
end

function dos_integrand(ω, q, p, vl, hl)
    m = rec_vecs(p).M
    if shapeBZ(q, vl, hl, m)
        dos_ω(ω, q, p)
    else
        0.0
    end
end

function dos_ω(ω, q, p)
    h = bistritzer_hamiltonian(p, q)
    ϵs, _ = eigen(Matrix(h))
    Γ = 1e-4 #
    return sum(Γ ./ ((ϵs .- ω).^2 .+ Γ.^2))
end