using SparseArrays
include("bistritzermodel.jl")
include("observables_bistritzer.jl")

"""
creates the matrix representation of the C2 symmetry BM model for valley ν at the Γ point
this symmetry does not mix valleys it mixes sublattices and layers
"""
c2x(p) = c2x(p , rec_vecs(p))

function c2x(p, rcp) # valid for the high symmetry points My and Γ
    links = c2xlinks(p, rcp)
    dim = length(links[:,1]) # 4 sublattice and sublayer
    mat = spzeros(ComplexF64, 2dim, 2dim) # it is diagonal in layer
    count = 1
    for i in 1:2:2dim
        mat[2Int(links[count])-1:2Int(links[count]), i:i+1] = [0 1;1 0] 
        count += 1
    end
    nul = spzeros(size(mat,1), size(mat,1))
    # c2sub = c2_continuum(2dim)
    return  [nul mat; mat nul]# the C2y is σx in layer space
end

function c2xlinks(p, rcp) #M =[0.06426234562410761  , 0.0]
    gvecs = Gmesh(p , rcp)
    # find G vectors that are related by C2x
    aux = zeros(length(gvecs),1)  
    count = 1
    for i in 1:length(gvecs)
        aux[count] = findall(x -> isapprox(x, c2xop(gvecs[i])), gvecs)[1] 
        count += 1
    end
    return aux
end

c2xop(k) = [k[1], -k[2]]
_permute_flatten(v1,v2) = reduce(vcat, [[v1[i],v2[i]] for i in 1:length(v1)])
_reduce(vec) = reduce(vcat, vec)
clean(mat) = round.(sparse(mat), digits = 5)
invm(mat) = inv(Matrix(mat))

""" compute the C2x character of a selection of eigenvectors of the bistritzer model.
Since C2x does not mix valleys we can compute everything within a valley in the BZ model.
For C3 we must consider both valleys.
"""

function c2xcharacterkpath(s::System, nband; kws...)
    p = ParamsBM(paramsBM(s.θ, s.nmax), mu = s.μ, mass = s.mass_value)
    ks = kpath(p)
    aux = []
    for i in 1:length(ks)
    push!(aux, c2xcharacter(p, ks[1,:]; mass_term =  s.mass_term, ph  = s.ph, kws...)[nband])
    end
    aux
end

c2xcharacter(s::System, kpoint::Symbol = :Γ; kws...) = 
    c2xcharacter(ParamsBM(paramsBM(s.θ, s.nmax), mu = s.μ, mass = s.mass_value), kpoint; mass_term =  s.mass_term, ph  = s.ph, kws...)

c2xcharacter(p::ParamsBM, kpoint::Symbol = :Γ; kws...) = 
    c2xcharacter(p, symbol_to_momenta(kpoint); kws...)

c2xcharacter(s::System, kpoint::Array; kws...) = 
    c2xcharacter(ParamsBM(paramsBM(s.θ, s.nmax), mu = s.μ, mass = s.mass_value), kpoint::Array;  mass_term =  s.mass_term, ph  = s.ph, kws...)

function c2xcharacter(p::ParamsBM, kpoint::Array; kws...)
    nev = 6
    espk, vecpk = bistritzer_eigs(p, kpoint, :LinearAlgebra; kws...)
    indices_nev_central_eigenvalues = sortperm(abs.(espk .- p.mu))[1:nev]
    # Returns the smallest eigenvalues in absolute value around the Fermi level ordered
    sorted_indices_pk = Int.(sortperm(real.(espk[real.(sortperm(abs.(espk .- p.mu))[1:nev])])))
    indx = indices_nev_central_eigenvalues[sorted_indices_pk]
    println("Energies: ", 1000 .* espk[indx])
    #println(" ", round.(real.(espk[sorted_indices_pk]), digits = 10))
    ψ = Array{ComplexF64}(undef, size(vecpk,1), nev)
    i = 1                                      
    for ind in indx
        ψ[:, i] = vecpk[:, ind] 
        i += 1
    end
    return round.(c2xcharacter(ψ, p), digits = 4)
end


function c2xcharactertwoangles(s, θ) #angles θ = 1.08, 1.09
    mass = collect(0:0.001:0.05)
    f = Figure()
    ax = Axis(f[1,1], xlabel = "Δ (meV)", ylabel = "C2x", title =  "Angle: $(round(θ, digits = 3))")
    valsbelow = diag(real.(c2xcharacter(System(s, θ = θ,  mass_value = mass[1]))))
    scatter!(ax, mass[1], valsbelow[1], color = :green, label = "3V")
    scatter!(ax,  mass[1], valsbelow[2], color = :lightblue, label = "2V" )
    scatter!(ax,  mass[1], valsbelow[3], color = :blue,label = "1V")
    scatter!(ax,  mass[1], valsbelow[4], color = :red, label = "1C")
    scatter!(ax,  mass[1], valsbelow[5], color = :orange, label = "2C")
    scatter!(ax,  mass[1], valsbelow[6], color = :yellow, label = "3C")
    for m in mass
        valsbelow = diag(real.(c2xcharacter(System(s, θ = θ,  mass_value = m))))
        scatter!(ax, m, valsbelow[1], color = :green)
        scatter!(ax, m, valsbelow[2], color = :lightblue)
        scatter!(ax, m, valsbelow[3], color = :blue)
        scatter!(ax, m, valsbelow[4], color = :red)
        scatter!(ax, m, valsbelow[5], color = :orange)
        scatter!(ax, m, valsbelow[6], color = :yellow)

    end
    f[1, 2] = Legend(f, ax, "Transitions")
    return f
end

function c2xcharacterall(p::ParamsBM, kpoint::Array; kws...) # long way to compute it
    espk, vecpk = eigen(Matrix(bistritzer_hamiltonian(p, kpoint; kws...)))
    m = div(length(espk),2)    
    ψ = Array{ComplexF64}(undef, size(vecpk,1), 2)                                      
    ψ[:, 1] = vecpk[:, m]
    ψ[:, 2] = vecpk[:, m+1]
    #return  ψ
    return c2xcharacter(ψ, p)
end


function symbol_to_momenta(ksymb::Symbol)
    if ksymb == :Γ
        return [0.,0.]
    else throw(ArgumentError("kpoints other than the Γ point are not yet implemented")) end  
end

c2xcharacter(ψ::Array, p) = Diagonal(ψ' * c2x(p) * ψ)

