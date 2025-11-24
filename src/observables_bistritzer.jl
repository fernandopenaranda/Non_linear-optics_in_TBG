    include("bistritzermodel.jl")


    using Arpack
    using CairoMakie
    using LinearAlgebra

    function test_operators_uuu(s::System, a; valley = 1, trans = Transitions(2,1))
        p = ParamsBM(paramsBM(s.θ, s.nmax), mu = s.μ, mass = s.mass_value)
        energies(valley) = bands_bistritzer(ParamsBM(p, ν = valley, pointsk = 50), eigvecs = true, 
           ph  = s.ph, mass_term =  s.mass_term,)
        ks = kpath(p)
        
        es, phis = bistritzer_eigs(ParamsBM(p, ν = valley, pointsk = 50), [0.,0.], ph  = s.ph, mass_term =  s.mass_term,)
        v = vel(phis,  d_bistritzer_hamiltonian(p, a))
        dh = d_bistritzer_hamiltonian(p, a)
        p = ParamsBM(paramsBM(1.05, s.nmax), mu = s.μ, mass = s.mass_value)
        i, j = index_convertor(trans, es)
        omega, ra, Δa = omega_r_Δ(real.(es), phis, d_bistritzer_hamiltonian(p, a))
        rca = r_covariant(omega, ra, Δa, ra, Δa)   
        # aux = 0im
        dh = zeros(ComplexF64, size(ra,1), size(ra,1))
        for it in 1:2:size(ra, 1)
            dh[it:it+1, it:it+1] = [ 0 -im ; im  0]
        end
       
        aux = zeros(ComplexF64, size(omega, 1), size(omega, 2))
        aux01 = phis[:,j]' * dh 
        aux02 =  dh * phis[:,i]
        for p in 1:size(omega, 1)
            aux += ( phis[:, p]) * (phis[:,p]') ./ es[p]#(omega[p,i])# works
        end
        return  phis[:,j]' * dh  * aux *  dh * phis[:,i] * phis[:,i]' *  dh * phis[:,j] # works
    end 

    """given a k path it computes the bands of the `bistritzer_hamiltonian()
    along the linecut specified by `kpath` (i.e., K1-Γ-M-K2 line)"""
    function bands_bistritzer(p::ParamsBM; eigvecs = false, kws...)
        kmesh = kpath(p)
        sys_dim = 4 * Int(1 + 3p.nmax * (1 + p.nmax))
        println(p.mass)
        es = zeros(Float64, sys_dim, size(kmesh,1))
        if eigvecs == false
            for i in 1:size(kmesh,1)
                es[:,i] = real.(bistritzer_eigs(p, kmesh[i,:]; kws...)[1])
            end
            return (es, )
        elseif eigvecs == true
            eigvec_mat = zeros(ComplexF64, sys_dim, sys_dim, length(kmesh))
            for i in 1:size(kmesh,1)
                (es[:,i], eigvec_mat[:,:,i]) = bistritzer_eigs(p, kmesh[i,:]; kws...)
                return (es, eigvec_mat)
            end
        else throw(ArgumentError("`eigvecs` kwarg must be either `false` or `true`"))
        end
    end

    function bistritzersqrd_eigs(p::ParamsBM, qi, method = :LinAlg; nev = 4, kws...)
        #if method == :LinAlg
            h = bistritzer_hamiltonian(p, qi; kws...)
            l = eigen(Matrix(h^2))
            return l.values, l.vectors
        # else method == :Arpack
        #     l = eigs(bistritzer_hamiltonian(p, qi), nev = nev, sigma = 1e-6) 
        #     # shift invert method for central eigenvalues. Efficient
        #     return l[1], l[2]
        #end
    end


    function bistritzer_eigs(p::ParamsBM, qi, method = :LinAlg; nev = 4, kws...)
        #if method == :LinAlg
            l = eigen(Matrix(bistritzer_hamiltonian(p, qi; kws...)))
            return l.values, l.vectors
        # else method == :Arpack
        #     l = eigs(bistritzer_hamiltonian(p, qi), nev = nev, sigma = 1e-6) 
        #     # shift invert method for central eigenvalues. Efficient
        #     return l[1], l[2]
        #end
    end

    """given a mesh density with `pointsk` it computes the line K1, Γ, M, K2, where K1 and K2 
    are the non-equiv point of the MBZ"""
    function kpath(p::ParamsBM)
        rcp = rec_vecs(p)
        κ1, κ2, pointsk = rcp.κ1, rcp.κ2, p.pointsk
        meshdim = Int(5*pointsk/2 + 1)
        # meshdim = Int(5*pointsk + 1)

        kpoints = zeros(Float64, meshdim, 2)
    
        for q = 1:meshdim
            if q <= p.pointsk
                qx = κ1[1] - κ1[1]*(q - 1)/pointsk
                qy = κ1[2] - κ1[2]*(q - 1)/pointsk
            elseif q > pointsk && q <= 2*p.pointsk
                qx = (κ2[1] + κ1[1])/2*(q - 1 - pointsk)/pointsk
                qy = (κ2[2] + κ1[2])/2*(q - 1 - pointsk)/pointsk
            else
                qx = (κ2[1] + κ1[1])/2 + (κ2[1] - κ1[1])*(q - 1 - 2*pointsk)/pointsk
                qy = (κ2[2] + κ1[2])/2 + (κ2[2] - κ1[2])*(q - 1 - 2*pointsk)/pointsk
            end
            kpoints[q,:] = [qx, qy]
        end
        return kpoints
    end
    
    ######## Plots
    plotbands(mat; kw...) = plotbands!(Figure(resolution =(800,1000)), mat; kw...)
    function plotbands!(f::Figure, mat; dots = false, color = missing, ylimits = missing, xlimits = missing)
        ax = Axis(f[1, 1]; xlabel = "k", ylabel = "E [meV]",xlabelsize= 27, ylabelsize= 27, xticklabelsize = 27, yticklabelsize = 27)
        xarr = collect(1:size(mat,2))
        pointsk = 2/5 * (length(xarr)-1)
        if dots == false
            for i in 1:size(mat, 1)
                lines!(ax, xarr , 1e3 .* mat[i,:], color = ifelse(isa(color,Missing), :lightgray, color))
            end
        else
            for i in 1:size(mat, 1)
                scatter!(ax,  1e3 .*collect(1:size(mat,2)) , mat[i,:], markersize = 5, markeralpha = 0.8)
            end
        end   
        ax.xticks = ([1, pointsk+1, 2*pointsk+1, 5*pointsk/2 + 1], ["K1", "Γ", "M", "K2"])
        if isa(ylimits,Missing)
            nothing
        else
            ylims!(ax, ylimits[1], ylimits[2])
        end
        if isa(xlimits,Missing)
            nothing
        else
            xlims!(ax, ylimits[1], ylimits[2])
        end
        return f
    end

    function plotbands!(ax::Axis, mat; dashed = false, dots = false, color = missing, ylimits = missing, xlimits = missing)
        xarr = collect(1:size(mat,2))
        pointsk = 2/5 * (length(xarr)-1)
        if dots == false
            for i in 1:size(mat, 1)
                lines!(ax, xarr , 1e3 .* mat[i,:], color = ifelse(isa(color,Missing), :lightgray, color))#, linestyle = (dashed == true ? nothing : :dash))
            end
        else
            for i in 1:size(mat, 1)
                scatter!(ax, collect(1:size(mat,2)) ,  1e3 .*mat[i,:], markersize = 5, markeralpha = 0.8)
            end
        end   
        ax.xticks = ([1, pointsk+1, 2*pointsk+1, 5*pointsk/2 + 1], ["K1", "Γ", "M", "K2"])

        if isa(ylimits,Missing)
            nothing
        else
            ylims!(ax, ylimits[1], ylimits[2])
        end
        if isa(xlimits,Missing)
            nothing
        else
            xlims!(ax, ylimits[1], ylimits[2])
        end
        xlims!(ax, [1,  5*pointsk/2 + 1])
        return f
    end
