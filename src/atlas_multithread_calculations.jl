include("transision_resolved_nl_photocurrents.jl")
using Base.Threads

@with_kw struct Nonlinear_params
    a::Symbol
    b::Symbol
    c::Symbol
    θ::Float64
    μ::Float64
    nmax::Int64
    η::Float64
    evals::Int64
    mass_term::Symbol
    mass::Float64
    jdos::Bool
    valley_resolved::Bool
    ph::Bool
    T::Float64 = 0
end

nonlinear_params(a, b, c, θ, μ, nmax, η, evals, mass_term, mass, jdos, valley_resolved, ph) = 
    Nonlinear_params(a, b, c, θ, μ, nmax, η, evals, mass_term, mass, evals, jdos, valley_resolved, ph)

nonlinear_params() = Nonlinear_params(:x, :x, :x, 1.05, 3, 0.0, 1e-4, 10, :none, 0.01, false, false, true)

function shift(ωlist::Array, nl::Nonlinear_params, tasks::Int, task_id::Int) 
    ωvecs = split_array(ωlist, tasks)[task_id+1]
    println("Computing the ", task_id/tasks, " partition of ωvecs ", ωvecs[1], " to ", ωvecs[end])
    shift(ωvecs, nl)
end

function split_array(arr::Array, n::Int)
    len = length(arr)
    sub_array_size = div(len, n)
    remainder = rem(len, n)
    result = []
    start_idx = 1
    for i in 1:n
        current_size = sub_array_size + (i <= remainder ? 1 : 0)
        end_idx = start_idx + current_size - 1
        push!(result, arr[start_idx:end_idx])
        start_idx = end_idx + 1
    end
    return result
end



shift() = shift(collect(0.:0.0005:0.04))
shift(ωlist::Array) = shift(ωlist, calculation_params())
shift(nl::Nonlinear_params) = shift(collect(0.:0.0005:0.04), nl::Nonlinear_params)

function shift(ωlist, nl::Nonlinear_params) 
    start_time = time_ns()
    shift_current_TRS(nl.a, nl.b, nl.c, ParamsBM(paramsBM(1., nl.nmax), mass = nl.mass,
         θ = nl.θ, mu = nl.μ), ωlist, η = nl.η, mass_term = nl.mass_term, 
        evals = nl.evals, jdos = nl.jdos, valley_resolved = nl.valley_resolved, ph = nl.ph, T = nl.T)
    end_time = time_ns()
    elapsed_time = (end_time - start_time) / 10^9 # convert to seconds
    println("Elapsed time: ", elapsed_time, " seconds")
end



function shift(transition::Union{Transitions, Transitions_Array}, ωlist, nl::Nonlinear_params) 
    start_time = time_ns()
    shift_current_TRS(transition, nl.a, nl.b, nl.c, ParamsBM(paramsBM(1., nl.nmax), mass = nl.mass,
         θ = nl.θ, mu = nl.μ), ωlist, η = nl.η, mass_term = nl.mass_term, 
        evals = nl.evals, jdos = nl.jdos, valley_resolved = nl.valley_resolved, ph = nl.ph, T = nl.T)
    end_time = time_ns()
    elapsed_time = (end_time - start_time) / 10^9 # convert to seconds
    println("Elapsed time: ", elapsed_time, " seconds")
end

injection() = injection(collect(0.:0.0005:0.04))
injection(ωlist::Array) = injection(ωlist, calculation_params())
injection(nl::Nonlinear_params) = injection(collect(0.:0.0005:0.04), nl::Nonlinear_params)
function injection(ωlist, nl::Nonlinear_params) 
    start_time = time_ns()
    injection_current_TRS(nl.a, nl.b, nl.c, ParamsBM(paramsBM(1., 3), mass = nl.mass,
        θ = nl.θ, mu = nl.μ), ωlist, η = nl.η, mass_term = nl.mass_term, 
        evals = nl.evals, jdos = nl.jdos, valley_resolved = nl.valley_resolved, ph = nl.ph, T = nl.T)
        end_time = time_ns()
        elapsed_time = (end_time - start_time) / 10^9 # convert to seconds
        println("Elapsed time: ", elapsed_time, " seconds")
end
    

c_d() = c_d(collect(0.:0.0005:0.04))
c_d(ωlist::Array) = c_d(ωlist, calculation_params())
c_d(nl::Nonlinear_params) = c_d(collect(0.:0.0005:0.04), nl::Nonlinear_params)

function c_d(ωlist, nl::Nonlinear_params) 
    start_time = time_ns()
    circular_dichroinism(nl.a, nl.b, nl.c, ParamsBM(paramsBM(1., nl.nmax), mass = nl.mass,
         θ = nl.θ, mu = nl.μ), ωlist, η = nl.η, mass_term = nl.mass_term, 
        evals = nl.evals, jdos = nl.jdos, valley_resolved = nl.valley_resolved, ph = nl.ph, T = nl.T)
    end_time = time_ns()
    elapsed_time = (end_time - start_time) / 10^9 # convert to seconds
    println("Elapsed time: ", elapsed_time, " seconds")
end

#include("observables_bistritzer.jl")
#include("irreps.jl")

function calculation_params()
    println("Enter the calculation parameters: ")
    println("a: (x, y or z)")
    a = Symbol(readline())
    println("b: (x, y or z)")
    b = Symbol(readline())
    println("c: (x, y or z)")
    c = Symbol(readline())
    println("θ: (Float64)")
    θ = parse(Float64, readline()) 
    println("μ: (Float64)")
    μ = parse(Float64, readline()) 
    println("η: (Float64) (Press Enter to use default value 1e-4)")
    η_input = parse(Float64, readline()) 
    η = !isempty(η_input) ? η_input : 1e-4
    println("mass_term: (Press Enter to use default value :none)")
    mass_term_input = Symbol(readline()) 
    mass_term = !isempty(string(mass_term_input)) ? mass_term_input : :none
    println("mass: ")
    mass = parse(Float64, readline()) 
    println("evaluations: (Press Enter to use default value 10)")
    evals_input = parse(Int64, readline())
    evals = !isempty(evals_input) ? evals_input : 10
    println("jdos: (Press Enter to use default value false)")
    jdos_input = parse(Bool, readline())
    jdos = !isempty(string(jdos_input)) ? jdos_input : false
    println("valley_resolved: (Press Enter to use default value false)")
    valley_resolved_input = parse(Bool, readline())
    valley_resolved = !isempty(string(valley_resolved_input)) ? valley_resolved_input : false
    println("ph: (Press Enter to use default value true)")
    ph_input = parse(Bool, readline())
    ph = !isempty(string(ph_input)) ? ph_input : true
    return nonlinear_params(a, b, c, θ, μ, η, evals, mass_term, mass, jdos, valley_resolved, ph)
end