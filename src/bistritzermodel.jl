
"""
This module contains the continuum (BM) model projected onto a basis of plane waves with layer 
and sublattice degrees of freedom.
    exported functions:
        - `bistritzer_hamiltonian(p, q; intra_only = false)`
            computes the Hamiltonian of the BM model
        -  `paramsBM(θ::Float64, nmax, ν = 1; kw...) `
            returns a struct object of type ::ParamsBM) that contains the presets of the 
            model; namely the twist angle `theta`, the cutoff in the momentum expansion
            `nmax`, and the valley index `ν = ±1`.
                see also `paramsBM(obj::ParamsBM, ν::Int)`.
        - `rec_vecs(p::ParamsBM) `
            returns a struct object of type ::Rec_vecs with reciprocal space directives
        - `Gindices(p)`
            returns the k-mesh of the plane-wave expansion up to cutoff p.nmax
For more details type ?functionname
To use the exported function type `using Main.bistritzer`
    e.g.
        > using Main.bistritzer
        > p = paramsBM(1.05, 1,  1); 
        > h = bistritzer_hamiltonian(p, [0,0]; intra_only = false);
"""
# module bistritzer
using LinearAlgebra
using SparseArrays
using Parameters
using StaticArrays

const σ0 = SA[1 0; 0 1]
const σx = SA[0 1; 1 0]
const σy = SA[0 -1im; 1im 0]
const σz = SA[1 0; 0 -1]
############################ Parameters
@with_kw struct ParamsBM # [in eV]
    a0::Float64
    ν::Int
    θ::Float64
    nmax::Int
    pointsk::Int
    mu::Float64
    t::Float64
    tAA::Float64
    tAB::Float64
    t3::Float64          # Kang Vafek mass  t1*sin(θ/2) aprox t1/100 # Phys. Rev. Lett. 125, 257602
    mass::Float64
    δ1t::Vector{Float64} # real space vector of the K points in the top layer
    δ2t::Vector{Float64}
    δ3t::Vector{Float64}
    δ1b::Vector{Float64} # real space vector of the bottom layer
    δ2b::Vector{Float64}
    δ3b::Vector{Float64}
end

struct Rec_vecs{T} 
    b1::T # unrotated primitive vector of monolayer graphene
    b2::T # unrotated primitive vector of monolayer graphene
    G1::T # bravais MBZ vector
    G2::T # bravais MBZ vector
    Kb::T # K point in the rot bottom layer
    Kt::T # K point in the rot top layer
    κ1::T # 
    κ2::T # 
    M::T # M point
end

function paramsBM!(p::ParamsBM, field_name::Symbol, new_value)
    field_names = fieldnames(typeof(p))
    if field_name in field_names
        println(field_name)
        p.field_name = new_value
    else
        error("Field '$field_name' not found in struct.")
    end
    return p
end

function paramsBM(obj::ParamsBM, ν::Int; nmax = missing)
    ParamsBM(obj.a0, ν, obj.θ, ifelse(isa(nmax, Missing), obj.nmax, nmax), obj.pointsk,
        obj.mu, obj.t, obj.tAA, obj.tAB, obj.t3, obj.mass, obj.δ1t, obj.δ2t, obj.δ3t, 
        obj.δ1b, obj.δ2b, obj.δ3b) 
end

function paramsBM(θ::Float64, nmax, ν = 1; kpoints = 30) 
    δ1 = [0, 1/√3]
    δ2 = [-1/2, -1/2/√3]
    δ3 =  [1/2, -1/2/√3]
    tAB = 0.0975
    return ParamsBM(2.46, ν, θ, nmax, kpoints, 0*0.00164448, -2.46575, 0.81 * tAB , tAB,  tAB * sin(θ*2pi/360/2), 0.017, # remove - before tAB * sin... (now is correct)
        rot_top(θ)*δ1, rot_top(θ)*δ2, rot_top(θ)*δ3, rot_bot(θ)* δ1, 
        rot_bot(θ)*δ2, rot_bot(θ)*δ3)
end

function paramsBMKAPLAN(θ::Float64, nmax, ν = 1; kpoints = 30)
    p = paramsBM(θ, nmax, ν; kpoints = kpoints)
    return ParamsBM(p, tAA = p.tAA, mass = 0.017)
end

function rec_vecs(p::ParamsBM) 
    b1 = 2π * [1, 1/√3]
    b2 = 2π * [-1, 1/√3]
    Kb = p.ν * rot_bot(p.θ) * [4π/3, 0]
    Kt = p.ν * rot_top(p.θ) * [4π/3, 0]
    G1, G2, κ1, κ2, M = _rcpvecs(p.θ, b1, b2)
    return Rec_vecs(b1, b2, G1, G2, Kb, Kt, κ1, κ2, M)
end

function _rcpvecs(θ, b1, b2) # return the bravais vectors and the momentum transfers q's
    G1 = g1(θ, b1)
    G2 = g2(θ, b2)
    return G1, G2, κ(G1, G2), κ(G2, G1), m(G1, G2)
end

g1(θ, b1) = rot_bot(θ) * b1 - rot_top(θ) * b1
g2(θ, b2) = rot_bot(θ) * b2 - rot_top(θ) * b2
κ(G1, G2) =  1/3 * (2*G1+G2) 
m(G1, G2) = 1/2 * (G1+G2)

rot_top(θ) = [cos(θ*π/360) -sin(θ*π/360); sin(θ*π/360) cos(θ*π/360)]
rot_bot(θ) = rot_top(-θ)

# _________________________________________________________________________________________

"""
        `bistritzer_hamiltonian(p, q; intra_only = false)`
computes the Hamiltonian of the BM model. `p::Params` contains the system information
set by `p = params(θ, nmax)` at momentum `q`. The `intra_only = true` optinal kwarg yields
the Hamiltonian of the two decoupled graphene layers.
    mass_term = :none
    mass_term = :kang_vafek_mass  # breaks ph but keeps the kpoints since they are protected
                                  # by C2z*T
    mass_term = :normal           # a σz term that breaks C2z*T, C2x, C2z
    mass_term = :layered          # a σz term only in the top layer that breaks C2z*T
    mass_term = :sigmaztauz       # a σzτz term that breaks C2z*T, C2y, C2z
"""

function bistritzer_hamiltonian(p, q; intra_only = false, kw...)
    rcp = rec_vecs(p) # Computes the K vectors of the MBZ
    indices = Gindices(p) # Generates a mesh of momenta up to a radius set by p.nmax
    η = ifelse(intra_only == true, 0,1)
    return intralayer_hamiltonian(p, rcp, q, indices; kw...) + interlayer_hamiltonian(p, indices; kw...)
end

""" 
    `intralayer_hamiltonian(p, rcp, q)`
computes the intralayer Hamiltonian for given presets `p::Params` and `rcp::Rec_vecs` at 
momenta `q` for a basis expanded in the MBZ reciprocal wavevectors.

    `intralayer_hamiltonian(p, rcp, q, indices)`
computes the intralayer Hamiltonian for given presets `p::Params` and `rcp::Rec_vecs` at 
momenta `q` for a basis expanded in the MBZ reciprocal wavevectors labeled by 
`indices = Gindices(p)`.
"""
intralayer_hamiltonian(p, rcp, q; kw...) = intralayer_hamiltonian(p, rcp, q, Gindices(p); kw...)

function intralayer_hamiltonian(p, rcp, q, indices; mass_term = :none, ph = true, TRS = true, TRSbreaking_mass = 0.00001, kw...)
    massmethods(mass_term)
    n_mass = ifelse(mass_term == :normal || mass_term == :szvafek, p.mass * σz, 0 * σ0)
    n_layer_mass = ifelse(mass_term == :layered, p.mass * σ0, 0 * σ0)
    n_kaplan_mass = ifelse(mass_term == :KAPLAN, p.mass * σz, 0 * σ0) #  sigmaz(1-tauz)/2  mass as in # https://journals.aps.org/prresearch/pdf/10.1103/PhysRevResearch.4.013209
    interlayer_bias_types = [:staggered, :staggered_normal_pert, :tauzvafek]
    n_interlayer_bias = (mass_term in interlayer_bias_types) ? p.mass/2 * σ0 : 0 * σ0
    n_sigmazpert = ifelse(mass_term == :staggered_normal_pert, 1e-5 * σz, 0 * σ0 )
    n_sigmaztauz_mass = ifelse(mass_term == :sigmaztauz || mass_term == :sztauzvafek, p.mass * σz, 0 * σ0)
    n_vafek_kaplan_mass = ifelse(mass_term == :vafekkaplan, p.mass * σz, 0 * σ0)
    n_noTRS = ifelse(TRS == false, σ0, 0*σ0) .* p.ν  * TRSbreaking_mass# break TRS by means of a ρz
    n_noTRS_normal = ifelse(TRS == false, p.ν * p.mass * σz, 0*σ0) # this neutralize the normal mass for one valley breaking valley symmetry

    dim_klat = 4 * Int(1 + 3p.nmax * (1 + p.nmax)) # dimension of the plane-wave basis
    mat = spzeros(ComplexF64, dim_klat, dim_klat) 
    count = 1
    halfdim = Int(dim_klat/2)
    μ = p.mu
    for i in 1:2:halfdim-1
        st = _sublattice_coupling(p, rcp, indices[count], q, layer =:top, ph)
        mat[i:i+1,i:i+1] = [-μ st; conj(st) -μ] + n_mass + n_layer_mass/2 + n_sigmaztauz_mass + n_interlayer_bias + n_sigmazpert + n_noTRS + 0n_noTRS_normal 
        sb = _sublattice_coupling(p, rcp, indices[count], q, layer = :bot, ph)
        mat[i+halfdim:i+halfdim+1,i+halfdim:i+halfdim+1] = 
            [-μ sb; conj(sb) -μ] - n_layer_mass/2+ n_mass - n_sigmaztauz_mass + n_kaplan_mass + n_vafek_kaplan_mass - n_interlayer_bias + n_sigmazpert + n_noTRS + 0n_noTRS_normal 
        count += 1
    end
    return mat
end

"intralayer coupling of graphene proyected into the layer and sublat
space:  `t * 1im * Σ_i δ^l_i'* (q + n1*G1 + n2*G2 - κ^l) exp(1im*δ_i^l*K_l)`. 
`l`` refers to the layer index and `δ^l_i` are the real space vectors connecting first 
neighbours in the rotated layer"

function _sublattice_coupling(p::ParamsBM, rcp, indices, q, ph; layer = :bot)
    if layer == :bot
        δs = [p.δ1b, p.δ2b, p.δ3b]
        qvecs = rcp.G1 * indices[1] + rcp.G2 * indices[2] + q  - p.ν * rcp.κ1
        K = rcp.Kb # The valley index has been already introduced in Kb def 
    else
        δs = [p.δ1t, p.δ2t, p.δ3t]
        qvecs = rcp.G1 * indices[1] + rcp.G2 * indices[2] + q - p.ν * rcp.κ2
        K = rcp.Kt # The valley index has been already introduced in Kb def
    end
    if ph == true
        return - p.t * 3/sqrt(3)/2 * (qvecs[1] * p.ν - 1im * qvecs[2]) # CONTINUUM MODEL - PH SYMMETRIC
    else
        return p.t * 1im * sum([δs[j]' * qvecs * exp(1im*δs[j]'*K) for j in 1:3]) #KOSHINO MODEL -> ROT PAULI MATRIX (BROKEN PH) # with a minus to reproduce kaplan
    end
    # return p.t * 1im * sum([δs[j]' * qvecs * exp(1im*δs[j]'*K) for j in 1:3])
end

# _________________________________________________________________________________________

"""
`interlayer_hamiltonian(p)`
Creates the k independent Hamiltonian associated to the hopping processes between the plane 
waves at different layers at momenta `p`.
"""
interlayer_hamiltonian(p; kw...) = interlayer_hamiltonian(p, Gindices(p); kw...)

function interlayer_hamiltonian(p, indices; mass_term = :none, kw...)
    dim_klat = Int(1 + 3p.nmax * (1 + p.nmax))
    ham_dim = 4 * dim_klat
    halfdim = Int(2*dim_klat) # dimension of the top and bottom blocks
    mat = spzeros(ComplexF64, ham_dim, ham_dim)
        # We start by coupling plane waves with the same momenta
    aux = p.tAA * σ0 + p.tAB * σx
    kang_vafek_types = [:kang_vafek_mass, :vafekkaplan, :szvafek, :tauzvafek, :sztauzvafek]
    kv_mass = (mass_term in kang_vafek_types) ?  p.t3 * p.ν * σz : 0 * σ0
    # n_kivc = (mass_term == :KIVC ? mass_term .* σx : 0 .* σ0)

    for i in 1:2:halfdim
        mat[i:i+1, halfdim+i:halfdim+i+1] = aux - 1im * kv_mass 
        mat[halfdim+i:halfdim+i+1, i:i+1] = aux + 1im * kv_mass
    end
    # Now we couple those momenta that differ by G1 or G2:
    # `(n1'-n1'') G1 + (n2'-n2'') G2 = G1 or G2`
    for i in 1:dim_klat
        for j in 1:dim_klat
            dif_ind = indices[i] .- indices[j]
            if (dif_ind[1] == 0 && dif_ind[2] == p.ν) || (dif_ind[1] == -p.ν && dif_ind[2] == 0)
                # ϕ = -p.ν * ifelse(dif_ind[1] == 0, 1, -1) * 2π/3 
                # val = p.tAA * σ0 + p.tAB * σx .* [0 exp(1im*ϕ); exp(-1im*ϕ) 0] 
                ϕ = - ifelse(dif_ind[1] == 0, 1, -1) * 2π/3 
                val = p.tAA * σ0 + p.tAB * (σx .* cos(ϕ) .- sin(ϕ*p.ν) .* σy)

                mat[2*(i-1)+1:2*(i-1)+2, 2*(j-1)+halfdim+1:2*(j-1)+halfdim+2] = val - 1im * kv_mass
                mat[2*(j-1)+halfdim+1:2*(j-1)+halfdim+2, 2*(i-1)+1:2*(i-1)+2] = val + 1im * kv_mass
            else nothing end
        end
    end
    return mat
end


# ------------------------------------------------------------------------------------------
# MOMENTA DERIVATIVES OF H
# ------------------

"""
`d_bistritzer_hamiltonian(p::ParamsBM,  dir)`
Note that the q derivative makes dH momentum independent
returns the derivative of the Bistritzer-MacDonald Hamiltonian along direction 
`dir = :x, :y, respectively, built with presets `p`.
Units = [eV * Angstrom]
Note that the interlayer_hamiltonian does not contribute to the derivative
"""
function d_bistritzer_hamiltonian(p, dir::Symbol; kw...)
    rcp = rec_vecs(p) # Computes the K vectors of the MBZ
    indices = Gindices(p)
    return d_intralayer_hamiltonian(p, rcp, indices, mapdir(dir); kw...)
end

function mapdir(dir) 
if dir == :x
    1
elseif dir == :y
    2
else throw(ArgumentError("dir must be :x or :y")) end
end

""" 
`d_intralayer_hamiltonian(p, rcp, dir)`
computes the derivative of the intralayer Hamiltonian along direction `dir`
for given presets `p::Params` and `rcp::Rec_vecs` for a basis expanded in
the MBZ reciprocal wavevectors.

`d_intralayer_hamiltonian(p, rcp, indices, dir)`
computes the derivative of the intralayer Hamiltonian along direction `dir`for given 
presets `p::Params` and `rcp::Rec_vecs`  for a basis expanded in the MBZ
reciprocal wavevectors labeled by `indices = Gindices(p)`.
"""
d_intralayer_hamiltonian(p, rcp,  dir; kw...) = 
d_intralayer_hamiltonian(p, rcp,  Gindices(p), dir; kw...)

function d_intralayer_hamiltonian(p, rcp, indices, dir; kw...)
    dim_klat = 4 * Int(1 + 3p.nmax * (1 + p.nmax)) # dimension of the plane-wave basis
    mat = spzeros(ComplexF64, dim_klat, dim_klat) 
    count = 1
    halfdim = Int(dim_klat/2)
    for i in 1:2:halfdim-1
        st = _d_sublattice_coupling(p, rcp, indices[count], dir, layer=:top; kw...)
        mat[i:i+1,i:i+1] = [0 st; conj(st) 0]    
        sb = _d_sublattice_coupling(p, rcp, indices[count], dir, layer=:bot; kw...)
        mat[i+halfdim:i+halfdim+1,i+halfdim:i+halfdim+1] = [0 sb; conj(sb) 0]    
        count += 1
    end
    return mat
end

"""
derivative of the intralayer coupling of graphene proyected into the layer and sublat
space:  `t * 1im * Σ_i δ^l_i' exp(1im*δ_i^l*K_l)`. q
`l`` refers to the layer index and `δ^l_i` are the real space vectors connecting first 
neighbours in the rotated layer. Important since the δs do not multiply the gs we must 
give them length units a0
"""
function _d_sublattice_coupling(p::ParamsBM, rcp, indices, dir; layer = :bot, ph = true,  kws...)
    if layer == :bot
        δs = [p.δ1b, p.δ2b, p.δ3b]
        K = rcp.Kb # The valley index has been already introduced in Kb def 
    else
        δs = [p.δ1t, p.δ2t, p.δ3t]
        K = rcp.Kt # The valley index has been already introduced in Kb def
    end
    if ph == true
        if dir == 1
            return - p.t * p.a0 * 3/sqrt(3)/2 * p.ν         # In Å # PH SYMMETRIC
        else
            return - p.t * p.a0 * 3/sqrt(3)/2 * - 1im       # In Å # PH SYMMETRIC
        end
    else
        return  p.t * 1im * sum([p.a0 * δs[j][dir] * exp(1im*δs[j]'*K) for j in 1:3])  # In Å # ROT PAULI MATRICES with a - to reproduce Kaplan et al
    end    
    # if ph == true
    #     if dir == 1
    #         return - p.t * p.a0 * 3/sqrt(3)/2 * p.ν         # In Å # PH SYMMETRIC
    #     else
    #         return - p.t * p.a0 * 3/sqrt(3)/2 * - 1im       # In Å # PH SYMMETRIC
    #     end
    # else
    #     return  p.t * 1im * sum([p.a0 * δs[j][dir] * exp(1im*δs[j]'*K) for j in 1:3])  # In Å # ROT PAULI MATRICES with a - to reproduce Kaplan et al
    # end    
end

##### K MESH STARS

"""
    `Gindices(p::ParamsBM)`
Creates the hexagonal lattice of momenta transfer vectors up to shell `nmax` of the form 
`Gs = n1 G1 + n2 G2` with n1, n2 = Gindices(p)[:].
See: `interlayer_hamiltonian`
"""
function Gindices(p::ParamsBM)
    Gindex = [(0, 0)]
    for i = 1:p.nmax
        n1 = i
        n2 = 0
        for k = 1:3
            for j = 1:i
                n1 -= ifelse(k > 1, 1, 0)
                n2 += 2 - k
                push!(Gindex, (n1, n2))
                j += 1
            end
            k += 1
        end
        for k = 1:3
            for j = 1:i
                n1 += ifelse(k > 1, 1, 0)
                n2 += k - 2
                push!(Gindex, (n1, n2))
                j += 1
            end
            k += 1
        end
        i += 1
    end
    return Gindex
end

Gmesh(p::ParamsBM, rcp) = Gmesh(Gindices(p), rcp)
function Gmesh(indices::Array, rcp)
    Gvecs = []
    for i in 1:length(indices)
        push!(Gvecs, rcp.G1 .* indices[i][1] + rcp.G2 .* indices[i][2])
    end
    return Gvecs
end

function massmethods(mass_term)
    masses = [:none, :kang_vafek_mass, :normal, :tauz, :layered, :sigmaztauz, :KAPLAN, 
        :szvafek, :sztauzvafek, :tauzvafek, :staggered, :staggered_normal_pert, :vafekkaplan,
        :KIVC]
    if mass_term in masses
        nothing
    else throw(ArgumentError("Provided `mass_term` is not implemented. Available masses are: ", masses))
    end
end
