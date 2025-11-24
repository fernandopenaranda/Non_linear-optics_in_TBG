include(dirname(pwd())* "ShiftCurrents_in_tBLG.jl")

#=
_________________________________________________________________________________________

Supplementary 1
_________________________________________________________________________________________
=#

cfigk = computation_constructor(
    iterations = 80000,
    ωs = collect(0.:0.01:0.04),
    η = 5e-4, 
    valley_resolved = false, 
    jdos = false, 
    a = :x, 
    b = :x, 
    c = :x)

    
figka = system_constructor(
    T = 0., 
    μ = -0.00, 
    θ = 1.0, 
    ph = true, # no rot of pauli matrix
    mass_term = :normal, 
    mass_value = 0.000000,  
    TRS = false,
    TRSbreaking_mass = 0.00,
    nmax = 3
)

figkb = system_constructor(
    T = 0., 
    μ = -0.00, 
    θ = 1.0, 
    ph = true, # no rot of pauli matrix  c
    mass_term = :normal, 
    mass_value = 0.000000,  
    TRS = false,
    TRSbreaking_mass = 0.001,
    nmax = 3
)

figkc = system_constructor(
    T = 0., 
    μ = -0.005, 
    θ = 1.0, 
    ph = true, # no rot of pauli matrix  c
    mass_term = :normal, 
    mass_value = 0.000000,  
    TRS = false,
    TRSbreaking_mass = 0.00,
    nmax = 3
)

figkd = system_constructor(
    T = 0., 
    μ = -0.005, 
    θ = 1.0, 
    ph = true, # no rot of pauli matrix  c
    mass_term = :normal, 
    mass_value = 0.000000,  
    TRS = false,
    TRSbreaking_mass = 0.001,
    nmax = 3
)

figke = system_constructor(
    T = 0., 
    μ = -0.008, 
    θ = 1.0, 
    ph = true, # no rot of pauli matrix  c
    mass_term = :normal, 
    mass_value = 0.000000,  
    TRS = false,
    TRSbreaking_mass = 0.00,
    nmax = 3
)

figkf = system_constructor(
    T = 0., 
    μ = -0.008, 
    θ = 1.0, 
    ph = true, # no rot of pauli matrix  c
    mass_term = :normal, 
    mass_value = 0.000000,  
    TRS = false,
    TRSbreaking_mass = 0.001,
    nmax = 3
)

figkg = system_constructor(
    T = 0., 
    μ = -0.008, 
    θ = 1.0, 
    ph = true, # no rot of pauli matrix  c
    mass_term = :normal, 
    mass_value = 0.000000,  
    TRS = false,
    TRSbreaking_mass = 0.00,
    nmax = 3
)

figkh = system_constructor(
    T = 0., 
    μ = -0.008, 
    θ = 1.0, 
    ph = true, # no rot of pauli matrix  c
    mass_term = :normal, 
    mass_value = 0.000000,  
    TRS = false,
    TRSbreaking_mass = 0.001,
    nmax = 3
)


function computation_S1()
    injection_noTRS(figka, cfigk)
    injection_noTRS(figkb, cfigk)
    injection_noTRS(figkc, cfigk)
    injection_noTRS(figkd, cfigk)
    injection_noTRS(figke, cfigk)
    injection_noTRS(figkf, cfigk)
    injection_noTRS(figkg, cfigk)
    injection_noTRS(figkh, cfigk)
end

############

function plotbandsTRS(mass_value = 0.0, TRSbreaking_mass = 0.0, mu = 0.0)
    fig = Figure(); ax = Axis(fig[1,1]);
    esn = bands_bistritzer(ParamsBM(p, θ = 1.05,  mu = mu, ν = -p.ν, mass = mass_value), 
        eigvecs = false, mass_term = :normal, TRS = false, TRSbreaking_mass = TRSbreaking_mass );
    esp = bands_bistritzer(ParamsBM(p, θ = 1.05,  mu = mu, ν = p.ν, mass = mass_value), eigvecs = false, mass_term = :normal, 
        TRS = false, TRSbreaking_mass = TRSbreaking_mass, mass = mass_value)
    plotbands!(ax, esp[1], color = :orange)
    plotbands!(ax, esn[1], color = :gray) 
    ylims!(ax, [-0.05,0.05])
    fig
end
 #=
_________________________________________________________________________________________
                                          
                                        S4
_________________________________________________________________________________________
 =#

cfigS4 = computation_constructor(
    iterations = 20000,
    ωs = [0.0025, 0.015, 0.027],
    η = 5e-4, 
    valley_resolved = false, 
    jdos = false, 
    a = :x,
    b = :y,
    c = :z)

figS4a = system_constructor(
    T = 0., 
    μ = 0.0, 
    θ = 1.08, 
    ph = true, # no rot of pauli matrix  
    mass_term = :none, 
    mass_value = 0.0, 
    nmax = 3
)

function anglesweep(s, c; method = :SHIFT, steps = 5)
    θlist = sort(vcat(collect(1.055:0.01:1.12), [1.0845]))
    if method == :SHIFT
        [shift(System(s, θ = θ), c) for  θ in θlist ]
    else 
        [injection(System(s, θ = θ), c) for  θ in θlist]
    end
end

