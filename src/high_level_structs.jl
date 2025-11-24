@with_kw struct System
    T::Float64
    μ::Float64
    θ::Float64
    ph::Bool
    mass_term::Symbol
    mass_value::Float64
    nmax::Int
    TRS::Bool
    TRSbreaking_mass::Float64
end

@with_kw struct Computation
    iterations::Int
    ωs::Vector{Float64}
    η::Float64
    valley_resolved::Bool
    a::Symbol
    b::Symbol
    c::Symbol
    jdos::Bool
end


system_constructor(; T = 0., μ = 0., θ = 1.05, ph = true,  mass_term = :none, mass_value = 0., nmax = 3, TRS = true, TRSbreaking_mass = 0.0) = 
    System(T, μ, θ, ph, mass_term, mass_value, nmax, TRS, TRSbreaking_mass)

computation_constructor(; iterations = 1, ωs = [0.], η = 5e-4, valley_resolved = false, jdos = false, a = :x, b = :x, c = :x) = 
    Computation(iterations, ωs, η, valley_resolved, a, b, c, jdos)

shift(s::System, c::Computation; kws...) = response(s, c, shift_current_TRS; kws...)
injection(s::System, c::Computation; kws...) = response(s, c, injection_current_TRS; kws...)

# shift_noTRS(s::System, c::Computation) = response(s, c, shift_current_TRS) #not yet
injection_noTRS(s::System, c::Computation; kws...) = response(s, c, injection_current_noTRS; kws...)
# HAR~TREEE
hartree_injection_noTRS(s::System, c::Computation, p_hartree::ParamsHartree; kws...) = response(s, c, injection_current_noTRS; kws...)


injection_noTRS_map(s::System, c::Computation, p_hartree::ParamsHartree, μlist; kws...) = 
    response_map(s, c, p_hartree::ParamsHartree,  μlist, injection_current_noTRS_map ; kws...)

function response_map(s::System, c::Computation, p_hartree::ParamsHartree, μlist, func; kws...) 
    println("s.μ: ", s.μ)
    start_time = time_ns()
    mat1 = zeros(Float64, length(c.ωs), length(μlist))
    mat2 = zeros(Float64, length(c.ωs), length(μlist))
    for (i,μ) in enumerate(μlist)
        print(" ", i/length(μlist))
        val = func(c.a, c.b, c.c, ParamsBM(paramsBM(s.θ, s.nmax), mass = s.mass_value, mu = μ),  p_hartree,
        c.ωs, 1, η = c.η, mass_term = s.mass_term, TRS = s.TRS, TRSbreaking_mass = s.TRSbreaking_mass, evals = c.iterations, jdos = c.jdos, 
        valley_resolved = c.valley_resolved, ph = s.ph, T = s.T; kws...)
        mat1[:,i] .= val
    end
    save_data(μlist, c.ωs, mat1)
    print(" neg valley")
    for i in 1:length(μlist)
        print(" ", i/length(μlist))
        mat2[:,i] = func(c.a, c.b, c.c, ParamsBM(paramsBM(s.θ, s.nmax), mass = s.mass_value, mu = μlist[i]),  p_hartree,
            c.ωs, -1, η = c.η, mass_term = s.mass_term, TRS = s.TRS, TRSbreaking_mass = s.TRSbreaking_mass, evals = c.iterations, jdos = c.jdos, 
            valley_resolved = c.valley_resolved, ph = s.ph, T = s.T; kws...)
    end

    end_time = time_ns()
    elapsed_time = (end_time - start_time) / 10^9 # convert to seconds
    println("Elapsed time: ", elapsed_time, " seconds")
    save_data(μlist, c.ωs, mat1, mat2)
    return μlist, c.ωs, mat1, mat2
end

# TRANSITIONS
shift(s::System, c::Computation, t::Union{Transitions, Transitions_Array}) = response(s, c, t, shift_current_TRS)
injection(s::System, c::Computation, t::Union{Transitions, Transitions_Array}) = response(s, c, t, injection_current_TRS)
injection_noTRS(s::System, c::Computation, t::Union{Transitions, Transitions_Array}) = response(s, c, t, injection_current_TRS_noTRS)

function response(s::System, c::Computation, func; kws...) 
    start_time = time_ns()
    func(c.a, c.b, c.c, ParamsBM(paramsBM(s.θ, s.nmax), mass = s.mass_value, mu = s.μ), 
        c.ωs, η = c.η, mass_term = s.mass_term, TRS = s.TRS, TRSbreaking_mass = s.TRSbreaking_mass, evals = c.iterations, jdos = c.jdos, 
        valley_resolved = c.valley_resolved, ph = s.ph, T = s.T; kws...)
    end_time = time_ns()
    elapsed_time = (end_time - start_time) / 10^9 # convert to seconds
    println("Elapsed time: ", elapsed_time, " seconds")
end

function response(s::System, c::Computation, p_hartree::ParamsHartree, func; kws...)
    start_time = time_ns()
    func(c.a, c.b, c.c, ParamsBM(paramsBM(s.θ, s.nmax), mass = s.mass_value, mu = s.μ), p_hartree, 
        c.ωs, η = c.η, mass_term = s.mass_term, TRS = s.TRS, TRSbreaking_mass = s.TRSbreaking_mass, evals = c.iterations, jdos = c.jdos, 
        valley_resolved = c.valley_resolved, ph = s.ph, T = s.T; kws...)
    end_time = time_ns()
    elapsed_time = (end_time - start_time) / 10^9 # convert to seconds
    println("Elapsed time: ", elapsed_time, " seconds")
end

function response(s::System, c::Computation, t::Union{Transitions, Transitions_Array}, func)
    start_time = time_ns()
    func(t, c.a, c.b, c.c, ParamsBM(paramsBM(s.θ, s.nmax), mass = s.mass_value, mu = s.μ), 
        c.ωs, η = c.η, mass_term = s.mass_term, TRS = s.TRS, TRSbreaking_mass = s.TRSbreaking_mass, evals = c.iterations, jdos = c.jdos, 
        valley_resolved = c.valley_resolved, ph = s.ph, T = s.T)
    end_time = time_ns()
    elapsed_time = (end_time - start_time) / 10^9 # convert to seconds
    println("Elapsed time: ", elapsed_time, " seconds")
end

# plot the integrand of the shift current formular for an ω value resolved in momentum
function k_resolved_shift(ω, steps, s::System, c::Computation, trans::Transitions; focus = :none)
    start_time = time_ns()
    p = ParamsBM(paramsBM(s.θ, s.nmax), mass = s.mass_value, mu = s.μ)    
    kres = k_resolved_shift(trans, p, c.a, c.b, c.c, steps, ω, c.η, focus = focus, mass_term = s.mass_term, 
    evals = c.iterations, jdos = c.jdos, valley_resolved = c.valley_resolved, ph = s.ph)
    end_time = time_ns()
    elapsed_time = (end_time - start_time) / 10^9 # convert to seconds
    println("Elapsed time: ", elapsed_time, " seconds")
    return kres
end

function k_resolved_injection(ω, steps, s::System, c::Computation, trans::Transitions; focus = :none)
    start_time = time_ns()
    p = ParamsBM(paramsBM(s.θ, s.nmax), mass = s.mass_value, mu = s.μ)    
    kres = k_resolved_injection(trans, p, c.a, c.b, c.c, steps, ω, c.η; focus = focus, mass_term = s.mass_term, 
    evals = c.iterations, jdos = c.jdos, valley_resolved = c.valley_resolved, ph = s.ph)
    end_time = time_ns()
    elapsed_time = (end_time - start_time) / 10^9 # convert to seconds
    println("Elapsed time: ", elapsed_time, " seconds")
    return kres
end
