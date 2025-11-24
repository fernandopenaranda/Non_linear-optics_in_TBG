using CairoMakie, Colors, ColorSchemes, LaTeXStrings, Formatting

#   BANDS
function compute_bands(s::System; ylimits = [-50,50])
    p = ParamsBM(paramsBM(s.θ, s.nmax), mu = s.μ, mass = s.mass_value)
    energies(valley) = bands_bistritzer(ParamsBM(p, ν = valley, pointsk = 50), eigvecs = false, 
       ph  = s.ph, mass_term =  s.mass_term,)[1]
    fig = plotbands(real.(1000*energies(1)), color = :gray, ylimits = ylimits)
    plotbands!(fig, real.(1000*energies(-1)), color = :orange, ylimits = ylimits)
    fig
end

function compute_bands!(fig, s::System; ylimits = [-50,50])
    p = ParamsBM(paramsBM(s.θ, s.nmax), mu = s.μ, mass = s.mass_value)
    energies(valley) = bands_bistritzer(ParamsBM(p, ν = valley, pointsk = 50), eigvecs = false, 
       ph  = s.ph, mass_term =  s.mass_term,)[1]
    plotbands!(fig, real.(1000*energies(1)), color = :gray, ylimits = ylimits)
    plotbands!(fig, real.(1000*energies(-1)), color = :orange, ylimits = ylimits)
    fig
end

plotbands(mat; kws...) = plotbands!(Figure(resolution = (300*1.5,300*1.5)), mat; kws...)
function plotbands!(f, mat; dots = false, color = missing, ylimits = missing, xlimits = missing)
    ax = Axis(f[1, 1]; xlabel = "k", ylabel = "E [meV]", yticks =  [-60,-40,-20,0,20,40,60])
    xarr = collect(1:size(mat,2))
    pointsk = 2/5 * (length(xarr)-1)
    if dots == false
        for i in 1:size(mat, 1)
            lines!(ax, xarr ,  mat[i,:], color = ifelse(isa(color,Missing), :lightgray, color))
        end
    else
        for i in 1:size(mat, 1)
            scatter!(ax, collect(1:size(mat,2)), mat[i,:], markersize = 5, markeralpha = 0.8)
        end
    end   
    ax.xticks = ([1, pointsk+1, 2*pointsk+1, 5*pointsk/2 + 1], ["K1", "Γ", "M", "K2"])
    ax.yticks = ( [-60,-40,-20,0,20,40,60])
    if isa(ylimits,Missing)
        nothing
    else
        ylims!(ax, ylimits[1], ylimits[2])
    end
    if isa(xlimits, Missing)
        nothing
    else
        xlims!(ax, ylimits[1], ylimits[2])
    end
    
    return f
end

shiftsummedvalleys(fig, a, b, c, str::String, ylimits; kws...) = 
    shiftsummedvalleys!(fig, a, b, c, Matrix(CSV.read(str, DataFrame)), ylimits; kws...)

shiftsummedvalleys(a, b, c, str::String, ylimits; kws...) = 
    shiftsummedvalleys(a, b, c, Matrix(CSV.read(str, DataFrame)), ylimits; kws...)

shiftsummedvalleys(a, b, c, mat::Array,ylimits; kws...) = 
    shiftsummedvalleys!(Figure(resolution = (700, 500)), a, b, c, mat, ylimits; kws...)

function shiftsummedvalleys!(fig, a, b, c, mat, ylimits; color = :black, style = :solid, hidex = false, hidey = false, save = true, ylabel = "", scale = 1)
    with_theme(theme_latexfonts()) do
        ax = Axis(fig[1,1], xlabel = L"$\omega$ [meV]", ylabel = ylabel, xlabelsize=27, ylabelsize=27)
            #ytickformat = v -> format.(v, commas=false, precision=1))# 
        lines!(ax, mat[:, 1]*1000, mat[:,2]./scale, color = color, linestyle = style)
        ylims!(ax, ylimits./scale)
        xlims!(ax,[0.,40])
        if hidex == true
            hidexdecorations!(ax; grid = true, minorgrid = false, minorticks = false)
        else nothing end
        if hidey == true
            hideydecorations!(ax; grid = true, minorgrid = false, minorticks = false)
        else nothing end
    return fig
    end
end

protoplot!(ax, a, b, c, str::String, ylimits; kws...) = 
    protoplot!(ax,  Matrix(CSV.read(str, DataFrame)), ylimits; kws...)

function protoplot!(ax,  mat, ylimits; color = :black, style = :solid, hidex = false, hidey = false, save = true, ylabel = "", scale = 1.0, linewidth = 1.5)
        lines!(ax, mat[:, 1]*1000, mat[:,2]./scale, color = color, linestyle = style, linewidth = linewidth )
        ylims!(ax, ylimits./scale)
    return fig
end

ylab(a::Symbol,b::Symbol,c::Symbol) = string(a)*string(b)*string(c)
    