include(dirname(pwd())* "ShiftCurrents_in_tBLG.jl")

#Fig 1
cfig0 = computation_constructor(
    iterations = 1000,
    ωs = collect(0.:0.00025:0.04), 
    η = 1e-3, # for pictorical purposes 
    valley_resolved = false, 
    jdos = true, 
    a = :x, 
    b = :y, 
    c = :z)

fig0a1 = system_constructor(
    T = 0., 
    μ = 0.0, 
    θ = 1.05, 
    ph = true, # no rot of pauli matrix  
    mass_term = :none, 
    mass_value = 0.0, 
    nmax = 3,
)

fig0a2 = System(fig0a1, μ = 0.01)

# Fig 2
cfig1 = computation_constructor(
    iterations = 20000,
    ωs = collect(0.:0.00025:0.04), 
    η = 5e-4, 
    valley_resolved = false, 
    jdos = false, 
    a = :x, 
    b = :y, 
    c = :z)

fig1a = system_constructor(
    T = 0., 
    μ = 0.0, 
    θ = 1.08, 
    ph = true, # no rot of pauli matrix  
    mass_term = :none, 
    mass_value = 0.0, 
    nmax = 3
)
fig1b = system_constructor(
    T = 0., 
    μ = -0.01,#0.01, with a - sign is stable
    θ = 1.08, 
    ph = true, # no rot of pauli matrix  
    mass_term = :none, 
    mass_value = 0.0, 
    nmax = 3
)

###############
# Fig 3

cfig2a = computation_constructor(
    iterations = 20000,
    ωs = collect(0.:0.00025:0.04), 
    η = 5e-4, 
    valley_resolved = false, 
    jdos = false, 
    a = :y, 
    b = :y, 
    c = :y)

fig2a1 = system_constructor(
    T = 0., 
    μ = 0.0, 
    θ = 1.08, 
    ph = true, # no rot of pauli matrix  
    mass_term = :normal, 
    mass_value = 0.01, 
    nmax = 3
)
fig2a2 = System(fig2a1, μ = 0.01)
fig2a3 = System(fig2a1, μ = -0.01)

cfig2b = computation_constructor(
    iterations = 20000,
    ωs = collect(0.:0.00025:0.04), 
    η = 5e-4, 
    valley_resolved = false, 
    jdos = false, 
    a = :x, 
    b = :x, 
    c = :x)

fig2b1 = system_constructor(
    T = 0., 
    μ = 0.0, 
    θ = 1.08, 
    ph = true, # no rot of pauli matrix  
    mass_term = :sigmaztauz, 
    mass_value = 0.01, 
    nmax = 3
)
fig2b2 = System(fig2b1, μ = 0.01)
fig2b3 = System(fig2b1, μ = -0.01)

cfig2c = computation_constructor( # uncoment
    iterations = 20000, 
    ωs = collect(0.:0.00025:0.04), 
    η = 5e-4, 
    valley_resolved = false, 
    jdos = false, 
    a = :x, 
    b = :x, 
    c = :z)

fig2c1 = system_constructor(
    T = 0., 
    μ = 0.0, 
    θ = 1.08, 
    ph = true, # no rot of pauli matrix  
    mass_term = :staggered, 
    mass_value = 0.03, 
    nmax = 3
)
fig2c2 = System(fig2c1, μ = 0.01)
fig2c3 = System(fig2c1, μ = -0.01)

# _________________________________________________________________________________________
# Figure 4

cfig3 = computation_constructor(
    iterations = 20000,
    ωs = collect(0.:0.00025:0.04), 
    η = 5e-4, 
    valley_resolved = false, 
    jdos = false, 
    a = :y, 
    b = :y, 
    c = :y)

fig3a1 = system_constructor(
    T = 0., 
    μ = 0.0, 
    θ = 1.05, 
    ph = true, # no rot of pauli matrix  
    mass_term = :normal, 
    mass_value = 0.01, 
    nmax = 3
)

fig3a2 = System(fig3a1, mass_term = :KAPLAN)
fig3a3 = System(fig3a1, mass_term = :szvafek)
fig3a4 = System(fig3a1, mass_term = :szvafek) 
# Change sign to t3 in bistritzer model by hand. kw not implemented


################################################################################################################
# Common functions

function anglesweep(s, c; method = :SHIFT, steps = 5)
    if steps % 2 != 1 # I force the calculation to pick the magic angle
        steps += 1
    else nothing end
    θM = 1.085
    δθ = 0.03
    θlist = collect(θM-δθ:2δθ/(steps-1):θM+δθ)
    if method == :SHIFT
        [shift(System(s, θ = θ), c) for  θ in θlist ]
    else 
        [injection(System(s, θ = θ), c) for  θ in θlist]
    end
end

function anglesweep(s, c, t; method = :SHIFT,  steps = 6)
    if steps % 2 != 1 
        steps += 1
    else nothing end
    θM = 1.085
    δθ = 0.03
    θlist = collect(θM-δθ:2δθ/(steps-1):θM+δθ)
    if method == :SHIFT
       [shift(System(s, θ = θ), c, t) for  θ in θlist]
    else 
        [injection(System(s, θ = θ), c, t) for  θ in θlist]
    end
end
