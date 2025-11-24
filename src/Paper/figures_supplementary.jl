include(dirname(pwd())* "ShiftCurrents_in_tBLG.jl")
include(dirname(pwd())*"plot_paper_functions.jl")
# --------------------------------------------------------------------------------------------------------------
#                                          SUPPLEMENTARY FIGURES
# --------------------------------------------------------------------------------------------------------------

#=
_________________________________________________________________________________________

Supplementary 1
_________________________________________________________________________________________
=#


function lat_suppl1(str_mat)
    with_theme(merge(theme_latexfonts(), Theme(fontsize = 24))) do #
        suppl1(str_mat)     
    end
end

function suppl1(str_mat) 
    fontsizelab = 22
    fontsizetheme = 22
    labsize = 23
    scale = 1e-15
    fig = Figure(resolution = (1600,880))
    aspect = nothing
    ylab = L"$\text{Re}[\eta_{xxx}]/\tau$ [nm $\mu$A/$(V^2 fs)$]"
    
    ga = fig[1, 1] = GridLayout()
    ax1a = Axis(ga[1,1], xlabel = L"$\omega$ [meV]", ylabel = ylab, xlabelsize=labsize, ylabelsize=labsize, aspect = aspect)
    _suppl1(ax1a, str_mat, 1, scale)
    
    gb = fig[1, 2] = GridLayout()
    ax1b = Axis(gb[1,1], xlabel = L"$\omega$ [meV]", ylabel = ylab, xlabelsize=labsize, ylabelsize=labsize, aspect = aspect)
    _suppl1(ax1b, str_mat, 2, scale)
    
    gc = fig[1, 3] = GridLayout()
    ax1c = Axis(gc[1,1], xlabel = L"$\omega$ [meV]", ylabel = ylab, xlabelsize=labsize, ylabelsize=labsize, aspect = aspect)
    _suppl1(ax1c, str_mat, 3, scale)
    
    gd = fig[1, 4] = GridLayout()
    ax1d = Axis(gd[1,1], xlabel = L"$\omega$ [meV]", ylabel = ylab, xlabelsize=labsize, ylabelsize=labsize, aspect = aspect)
    _suppl1(ax1d, str_mat, 4, scale)
    
    ge = fig[2, 1] = GridLayout()
    ax1e = Axis(ge[1,1], xlabel = L"$\omega$ [meV]", ylabel = ylab, xlabelsize=labsize, ylabelsize=labsize, aspect = aspect)
    _suppl1(ax1e, str_mat, 5, scale)
    
    gf = fig[2, 2] = GridLayout()
    ax1f = Axis(gf[1,1], xlabel = L"$\omega$ [meV]", ylabel = ylab, xlabelsize=labsize, ylabelsize=labsize, aspect = aspect)
    _suppl1(ax1f, str_mat, 6, scale)
    
    gg = fig[2, 3] = GridLayout()
    ax1g = Axis(gg[1,1], xlabel = L"$\omega$ [meV]", ylabel = ylab, xlabelsize=labsize, ylabelsize=labsize, aspect = aspect)
    _suppl1(ax1g, str_mat, 7, scale)
    
    gh = fig[2, 4] = GridLayout()
    ax1h = Axis(gh[1,1], xlabel = L"$\omega$ [meV]", ylabel = ylab, xlabelsize=labsize, ylabelsize=labsize, aspect = aspect)
    _suppl1(ax1h, str_mat, 8, scale)
    
    g01 = fig[1, 0] = GridLayout()
    ax01 = Axis(g01[1,1], xlabel = "k", ylabel = "E [meV]")
    plotbandssupll1!(ax01, 0)
    
    
    g02 = fig[2, 0] = GridLayout()
    ax02 = Axis(g02[1,1], xlabel = "k", ylabel = "E [meV]")
    plotbandssupll1!(ax02, 0.001)
    
    hidexdecorations!(ax01, grid = false)
    hidexdecorations!(ax1a, grid = false)
    hidexdecorations!(ax1b, grid = false)
    hidexdecorations!(ax1c, grid = false)
    hidexdecorations!(ax1d, grid = false)
    hideydecorations!(ax1b, grid = false)
    hideydecorations!(ax1c, grid = false)
    hideydecorations!(ax1d, grid = false)
    hideydecorations!(ax1f, grid = false)
    hideydecorations!(ax1g, grid = false)
    hideydecorations!(ax1h, grid = false)

    ylims!(ax1a, [-250,250])
    ylims!(ax1b, [-250,250])
    ylims!(ax1c, [-250,250])
    ylims!(ax1d, [-250,250])
    ylims!(ax1e, [-250,250])
    ylims!(ax1f, [-250,250])
    ylims!(ax1g, [-250,250])
    ylims!(ax1h, [-250,250])
    
    cmap = ColorSchemes.viridis
    legenditems = [MarkerElement(color = :orange, marker=:rect, markersize = 15, strokecolor = :black),
    MarkerElement(color = :gray, marker=:rect, markersize = 15, strokecolor = :black)]
    terms = [L"$\nu = +", L"$\nu = -"]
    gn = fig[0,0:4] = GridLayout()
    Legend(gn[1, 1], legenditems, terms, orientation = :horizontal, tellwidth = true, tellheight =true, framewidth = .5)#, padding = (0,0,0,0))

    inset_ax1 = add_box_inset(fig; backgroundcolor=:snow2,left=510, right=700, bottom=480, top=580)
    _suppl1(inset_ax1, str_mat, 1, scale)
    xlims!(inset_ax1, [0,5.001])

    Label(ga[1, 1, Top()], L"$ \mu = 0 meV", fontsize=fontsizelab, color = :black,  padding=(0, 0, -40, 0))
    Label(gb[1, 1, Top()], L"$ \mu = -5 meV", fontsize=fontsizelab, color = cmap[0.25],  padding=(0, 0, -40, 0))
    Label(gc[1, 1, Top()], L"$ \mu = -8 meV", fontsize=fontsizelab, color = cmap[0.5],  padding=(0, 0, -40, 0))
    Label(gd[1, 1, Top()], L"$ \mu = -10 meV", fontsize=fontsizelab, color = cmap[0.75],  padding=(0, 0, -40, 0))
    Label(ge[1, 1, Top()], L"$ \mu = 0 meV", fontsize=fontsizelab, color = :black,  padding=(0, 0, -40, 0))
    Label(gf[1, 1, Top()], L"$ \mu = -5 meV", fontsize=fontsizelab, color = cmap[0.25],  padding=(0, 0, -40, 0))
    Label(gg[1, 1, Top()], L"$ \mu = -8 meV", fontsize=fontsizelab, color = cmap[0.5],  padding=(0, 0, -40, 0))
    Label(gh[1, 1, Top()], L"$ \mu = -10 meV", fontsize=fontsizelab, color = cmap[0.75],  padding=(0, 0, -40, 0))
    Label(g01[1, 1, Top()], "(a)", fontsize=fontsizelab, color = :black,  padding=(-220, -0, -40, 0))
    Label(g02[1, 1, Top()], "(f)", fontsize=fontsizelab, color = :black,  padding=(-220, -0, -40, 0))
    Label(ga[1, 1, Top()], "(b)", fontsize=fontsizelab, color = :black,  padding=(-220, -0, -40, 0))
    Label(gb[1, 1, Top()], "(c)", fontsize=fontsizelab, color = :black,  padding=(-220, -0, -40, 0))
    Label(gc[1, 1, Top()], "(d)", fontsize=fontsizelab, color = :black,  padding=(-220, -0, -40, 0))
    Label(gd[1, 1, Top()], "(e)", fontsize=fontsizelab, color = :black,  padding=(-220, -0, -40, 0))
    Label(ge[1, 1, Top()], "(g)", fontsize=fontsizelab, color = :black,  padding=(-220, -0, -40, 0))
    Label(gf[1, 1, Top()], "(h)", fontsize=fontsizelab, color = :black,  padding=(-220, -0, -40, 0))
    Label(gg[1, 1, Top()], "(i)", fontsize=fontsizelab, color = :black,  padding=(-220, -0, -40, 0))
    Label(gh[1, 1, Top()], "(j)", fontsize=fontsizelab, color = :black,  padding=(-220, -0, -40, 0))
    display(fig)
    save("suppl1.pdf", fig)
end

function _suppl1(ax, str_mat, i, scale)
    mat = Matrix(CSV.read(str_mat[i][1], DataFrame))
    lines!(ax,  1e3 .* mat[:,1], mat[:, 2]./scale, color = :gray)
    mat2 = Matrix(CSV.read(str_mat[i][2], DataFrame))
    lines!(ax,  1e3 .* mat2[:,1], mat2[:, 2]./scale, color = :orange)
    lines!(ax, 1e3 .* mat[:,1], (mat[:,2]+mat2[:,2])./scale, color = :black)
end

function plotbandssupll1!(ax, trsmass)
    p = paramsBM(1.05, 1,  1);
    p = ParamsBM(p, nmax = 3, pointsk = 100)
    esp = bands_bistritzer(ParamsBM(p, θ = 1,  mu = -0.00, ν = 1, mass =0.000), 
            eigvecs = false, mass_term = :normal, TRS = false, TRSbreaking_mass = trsmass);
    esn = bands_bistritzer(ParamsBM(p, θ = 1,  mu = -0.00, ν = -1, mass =0.000), 
            eigvecs = false, mass_term = :normal, TRS = false, TRSbreaking_mass = trsmass);
    plotbands!(ax, esp[1], color = :orange)
    plotbands!(ax, esn[1], color = :gray)
    cmap = ColorSchemes.viridis
    hlines!(ax, [0],  linestyle=:dash, color = :black, linewidth = 1)
    hlines!(ax, [-5],  linestyle=:dash, color = cmap[0.25], linewidth = 1)
    hlines!(ax, [-8],  linestyle=:dash, color = cmap[0.5], linewidth = 1)
    hlines!(ax, [-10],  linestyle=:dash, color = cmap[0.75], linewidth = 1)
    ylims!(ax, [-20, 20])
end

function add_box_inset(fig; backgroundcolor=:snow2, left=100, right=250, bottom=200, top=300)
    inset_box = Axis(fig, bbox= BBox(left, right, bottom, top), 
        xticklabelsize=18, yticklabelsize=18, backgroundcolor=:white)
    translate!(inset_box.scene, 0, 0, 10)  # bring content upfront
    return inset_box
end

#=
_________________________________________________________________________________________

Supplementary 2
_________________________________________________________________________________________
=#


function suppl2() 
    ylims = [-0.8e3,0.9e3]
    fig = Figure(resolution = (1600,450))
    aspect = nothing
    ################ PANEL A
    ylab = L"$\sigma_{xyz}$ [nm $\mu$A/$V^2$]"
    ga = fig[1, 1] = GridLayout()
    ax1a = Axis(ga[1,1], xlabel = L"$\omega$ [meV]", ylabel = ylab, xlabelsize=labsize, ylabelsize=labsize, aspect = aspect)
    scale = 1e2
    protoplot!(ax1a, :x, :y, :z, shiftxyz_mu0_btheta_ref, ylims, color = :red, hidex = true, hidey = true, ylabel = ylab, scale = scale)
    protoplot!(ax1a, :x, :y, :z, shiftxyz_mu0_btheta, ylims, color = :red, style = :dash, hidex = true, hidey = true, ylabel = ylab, scale = scale)
    protoplot!(ax1a, :x, :y, :z, shiftxyz_mu0_atheta_ref, ylims, color = :blue, hidex = true, hidey = true, ylabel = ylab, scale = scale)
    protoplot!(ax1a, :x, :y, :z, shiftxyz_mu0_atheta, ylims, color = :blue, style = :dash, hidex = false, hidey = false, ylabel = ylab, scale = scale)
    xlims!(ax1a, [0,40])
       

    gb = fig[1, 2] = GridLayout()
    ylims = [-0.29e3,0.15e3]
    ax1b = Axis(gb[1, 1], xlabel = L"$\omega$ [meV]", ylabel = ylab, xlabelsize=labsize, ylabelsize=labsize, aspect = aspect)
    scale = 1e2
    protoplot!(ax1b, :x, :y, :z, shiftxyz_mu_btheta_ref, ylims, color = :red, hidex = true, hidey = true, ylabel = ylab, scale = scale)
    protoplot!(ax1b, :x, :y, :z, shiftxyz_mu_btheta, ylims, color = :red, style = :dash, hidex = true, hidey = true, ylabel = ylab, scale = scale)
    protoplot!(ax1b, :x, :y, :z, shiftxyz_mu_atheta_ref, ylims, color = :blue, hidex = true, hidey = true, ylabel = ylab, scale = scale)
    protoplot!(ax1b, :x, :y, :z, shiftxyz_mu_atheta, ylims, color = :blue, style = :dash, hidex = false, hidey = false, ylabel = ylab, scale = scale)
    xlims!(ax1b, [0,40])
    
    ylab = L"$i\eta_{xyz}/\tau$ [nm $\mu$A/$(V^2 fs)$]"
    gc = fig[1, 3] = GridLayout()
    ylims = [-0.7e-15,0.8e-15]
    ax1c = Axis(gc[1, 1], xlabel = L"$\omega$ [meV]", ylabel = ylab, xlabelsize=labsize, ylabelsize=labsize, aspect = aspect)
    scale = 1e-16
    protoplot!(ax1c, :x, :y, :z, injectionxyz_mu0_btheta_ref, ylims, color = :red, hidex = true, hidey = true, ylabel = ylab, scale = scale)
    protoplot!(ax1c, :x, :y, :z, injectionxyz_mu0_btheta, ylims, color = :red, style = :dash, hidex = true, hidey = true, ylabel = ylab, scale = scale)
    protoplot!(ax1c, :x, :y, :z, injectionxyz_mu0_atheta_ref, ylims, color = :blue, hidex = true, hidey = true, ylabel = ylab, scale = scale)
    protoplot!(ax1c, :x, :y, :z, injectionxyz_mu0_atheta, ylims, color = :blue, style = :dash, hidex = false, hidey = false, ylabel = ylab, scale = scale)
    xlims!(ax1c, [0,40])
    

    gd = fig[1, 4] = GridLayout()
    ylims = [-1.98e-16,2.98e-16]
    ax1d = Axis(gd[1, 1], xlabel = L"$\omega$ [meV]", ylabel = ylab, xlabelsize=labsize, ylabelsize=labsize, aspect = aspect)
    scale = 1e-16
    protoplot!(ax1d, :x, :y, :z, injectionxyz_mu_btheta_ref, ylims, color = :red, hidex = true, hidey = true, ylabel = ylab, scale = scale)
    protoplot!(ax1d, :x, :y, :z, injectionxyz_mu_btheta, ylims, color = :red, style = :dash, hidex = true, hidey = true, ylabel = ylab, scale = scale)
    protoplot!(ax1d, :x, :y, :z, injectionxyz_mu_atheta_ref, ylims, color = :blue, hidex = true, hidey = true, ylabel = ylab, scale = scale)
    protoplot!(ax1d, :x, :y, :z, injectionxyz_mu_atheta, ylims, color = :blue, style = :dash, hidex = false, hidey = false, ylabel = ylab, scale = scale)
    xlims!(ax1d, [0,40])

    terms = [L"$\theta<\theta_M\ $",  L"$\theta>\theta_M\ $", L"$\omega_3 = 0\ $",  L"$\omega_3 \neq 0$"]
    ccolors = [:black, :lightblue, :red]
    elem1 = MarkerElement(color = :red, marker=:rect, markersize = 15, strokecolor = :black)
    elem2 = MarkerElement(color = :blue, marker=:rect, markersize = 15, strokecolor = :black)
    elem3 = LineElement(color = :black,  linewidth = 2)
    elem4 = LineElement(color = :black, linestyle = :dash, linewidth = 2)
    legenditems = [elem1, elem2, elem3, elem4]

    g0 = fig[2,1:4] = GridLayout()
    Legend(g0[1, 1], legenditems, terms, orientation = :horizontal, tellwidth = false, tellheight =true, framewidth = .5)#, padding = (0,0,0,0))
    rowgap!(fig.layout, 10)
    resize_to_layout!(fig)

    Label(fig[1, 1], halign = :left, L"\times 10^{2}", tellwidth=false, tellheight=false, padding = (0,0,265,0)) #
    Label(fig[1, 2], halign = :left, L"\times 10^{2}", tellwidth=false, tellheight=false, padding = (0,0,265,0)) #
    Label(fig[1, 3], halign = :left, L"\times 10^{-1}", tellwidth=false, tellheight=false, padding = (0,0,265,0)) #
    Label(fig[1, 4], halign = :left, L"\times 10^{-1}", tellwidth=false, tellheight=false, padding = (0,0,265,0)) #

    Label(ga[1, 1, TopLeft()], "(a)", fontsize=fontsizelab,  padding=(0, 30, -15, 0))
    Label(gb[1, 1, TopLeft()], "(b)", fontsize=fontsizelab,  padding=(0, 30, -15, 0))
    Label(gc[1, 1, TopLeft()], "(c)", fontsize=fontsizelab,  padding=(0, 30, -15, 0))
    Label(gd[1, 1, TopLeft()], "(d)", fontsize=fontsizelab,  padding=(0, 30, -15, 0))
    display(fig)
    save("s2.pdf", fig)
end

function lat_suppl2()
    with_theme(merge(theme_latexfonts(), Theme(fontsize = 24))) do #
        suppl2()     
    end
end

#=
_________________________________________________________________________________________

Supplementary 3
_________________________________________________________________________________________
=#


function suppl3() #
    fig = Figure(resolution = (1600,440))
    aspect = nothing
    ################ PANEL A
    ylims = [-1e5,1.7e5]
    scale = 1e5
    ylab = L"$\sigma_{yyy}$ [nm $\mu$A/$V^2$]"
    ga = fig[1, 1] = GridLayout()
    ax1a = Axis(ga[1,1], xlabel = L"$\omega$ [meV]", ylabel = ylab, xlabelsize=labsize, ylabelsize=labsize, aspect = aspect)
    protoplot!(ax1a, :x, :y, :z, shiftyyy_mu_btheta_ref, ylims, color = :red, hidex = true, hidey = true, ylabel = ylab, scale = scale)
    protoplot!(ax1a, :x, :y, :z, shiftyyy_mu_btheta, ylims, color = :red, style = :dash, hidex = true, hidey = true, ylabel = ylab, scale = scale)
    protoplot!(ax1a, :x, :y, :z, shiftyyy_mu_atheta_ref, ylims, color = :blue, hidex = true, hidey = true, ylabel = ylab, scale = scale)
    protoplot!(ax1a, :x, :y, :z, shiftyyy_mu_atheta, ylims, color = :blue, style = :dash, hidex = false, hidey = false, ylabel = ylab, scale = scale)
    xlims!(ax1a, [0,40])
    # ################ PANEL B  
    ylab = L"$\sigma_{xxx}$ [nm $\mu$A/$V^2$]"
    gb = fig[1, 2] = GridLayout()
    ylims = [-1e4,1.8e4]
    ax1b = Axis(gb[1, 1], xlabel = L"$\omega$ [meV]", ylabel = ylab, xlabelsize=labsize, ylabelsize=labsize, aspect = aspect)
    scale = 1e4
    protoplot!(ax1b, :x, :y, :z, shiftxxx_mu_btheta_ref, ylims, color = :red, hidex = true, hidey = true, ylabel = ylab, scale = scale)
    protoplot!(ax1b, :x, :y, :z, shiftxxx_mu_btheta, ylims, color = :red, style = :dash, hidex = true, hidey = true, ylabel = ylab, scale = scale)
    protoplot!(ax1b, :x, :y, :z, shiftxxx_mu_atheta_ref, ylims, color = :blue, hidex = true, hidey = true, ylabel = ylab, scale = scale)
    protoplot!(ax1b, :x, :y, :z, shiftxxx_mu_atheta, ylims, color = :blue, style = :dash, hidex = false, hidey = false, ylabel = ylab, scale = scale)
    xlims!(ax1b, [0,40])

    # ################ PANEL C
    ylab = L"$\sigma_{xxz}$ [nm $\mu$A/$V^2$]"
    gc = fig[1, 3] = GridLayout()
    ylims = [-2.8e2,0.5e2]
    ax1c = Axis(gc[1, 1], xlabel = L"$\omega$ [meV]", ylabel = ylab, xlabelsize=labsize, ylabelsize=labsize, aspect = aspect)
    scale = 1e2
    protoplot!(ax1c, :x, :y, :z, shiftxxz_mu_btheta_ref, ylims, color = :red, hidex = true, hidey = true, ylabel = ylab, scale = scale)
    protoplot!(ax1c, :x, :y, :z, shiftxxz_mu_btheta, ylims, color = :red, style = :dash, hidex = true, hidey = true, ylabel = ylab, scale = scale)
    protoplot!(ax1c, :x, :y, :z, shiftxxz_mu_atheta_ref, ylims, color = :blue, hidex = true, hidey = true, ylabel = ylab, scale = scale)
    protoplot!(ax1c, :x, :y, :z, shiftxxz_mu_atheta, ylims, color = :blue, style = :dash, hidex = false, hidey = false, ylabel = ylab, scale = scale)
    xlims!(ax1c, [0,40])

    # ################ PANEL D
    ylab = L"$i\eta_{xxz}/\tau$ [nm $\mu$A/$(V^2 fs)$]"
    gd = fig[1, 4] = GridLayout()
    ylims = [-0.5e-16,3.7e-16]
    ax1d = Axis(gd[1, 1], xlabel = L"$\omega$ [meV]", ylabel = ylab, xlabelsize=labsize, ylabelsize=labsize, aspect = aspect)
    scale = 1e-16
    protoplot!(ax1d, :x, :y, :z, injectionxxz_mu_btheta_ref, ylims, color = :red, hidex = true, hidey = true, ylabel = ylab, scale = scale)
    protoplot!(ax1d, :x, :y, :z, injectionxxz_mu_btheta, ylims, color = :red, style = :dash, hidex = true, hidey = true, ylabel = ylab, scale = scale)
    protoplot!(ax1d, :x, :y, :z, injectionxxz_mu_atheta_ref, ylims, color = :blue, hidex = true, hidey = true, ylabel = ylab, scale = scale)
    protoplot!(ax1d, :x, :y, :z, injectionxxz_mu_atheta, ylims, color = :blue, style = :dash, hidex = false, hidey = false, ylabel = ylab, scale = scale)
    xlims!(ax1d, [0,40])
    
    Label(fig[1, 1], halign = :left, L"\times 10^{5}", tellwidth=false, tellheight=false, padding = (0,0,265,0)) #
    Label(fig[1, 2], halign = :left, L"\times 10^{4}", tellwidth=false, tellheight=false, padding = (0,0,265,0)) #
    Label(fig[1, 3], halign = :left, L"\times 10^{2}", tellwidth=false, tellheight=false, padding = (0,0,265,0)) #
    Label(fig[1, 4], halign = :left, L"\times 10^{-1}", tellwidth=false, tellheight=false, padding = (0,0,265,0)) #

    terms = [L"$\theta<\theta_M\ $",  L"$\theta>\theta_M\ $", L"$\omega_3 = 0\ $",  L"$\omega_3 \neq 0$"]
    elem1 = MarkerElement(color = :red, marker=:rect, markersize = 15, strokecolor = :black)
    elem2 = MarkerElement(color = :blue, marker=:rect, markersize = 15, strokecolor = :black)
    elem3 = LineElement(color = :black,  linewidth = 2)
    elem4 = LineElement(color = :black, linestyle = :dash, linewidth = 2)
    legenditems = [elem1, elem2, elem3, elem4]

    g0 = fig[2,1:4] = GridLayout()
    Legend(g0[1, 1], legenditems, terms, orientation = :horizontal, tellwidth = false, tellheight =true, framewidth = .5)#, padding = (0,0,0,0))
    rowgap!(fig.layout, 10)
    resize_to_layout!(fig)

    Label(ga[1, 1, TopLeft()], "(a)", fontsize=fontsizelab,  padding=(0, 30, -15, 0))
    Label(gb[1, 1, TopLeft()], "(b)", fontsize=fontsizelab,  padding=(0, 30, -15, 0))
    Label(gc[1, 1, TopLeft()], "(c)", fontsize=fontsizelab,  padding=(0, 30, -15, 0))
    Label(gd[1, 1, TopLeft()], "(d)", fontsize=fontsizelab,  padding=(0, 30, -15, 0))

    display(fig)
    save("s3.pdf", fig)
end

function lat_suppl3()
    with_theme(merge(theme_latexfonts(), Theme(fontsize = 24))) do #
        suppl3()     
    end
end

#=
_________________________________________________________________________________________

Supplementary 4
_________________________________________________________________________________________
=#

##### PLOT FUNCTION


path1 = "concatenated_shift.csv"
path2 = "concatenated_injection.csv"

function suppl4(path1, path2)
    mat = Matrix(CSV.read(path1))
    mat2 = Matrix(CSV.read(path2))
    fig = Figure()
    ylab1 = L"$\sigma_{xyz}$ [nm $\mu$A/$V^2$]" 
    θlist =  sort(vcat(collect(1.05:0.005:1.12),1.0845))
    with_theme(theme_latexfonts()) do
    ax = Axis(fig[1,1], xlabel = L"$\theta (^\circ)", ylabel = ylab1,  xlabelsize=24, ylabelsize= 24 )
    lines!(ax,  θlist, mat[:, 1])
    lines!(ax,  θlist, mat[:, 2])
    lines!(ax,  θlist, mat[:, 3])
    scatter!(ax,  θlist, mat[:, 1], label = "ω = 2.5 meV")
    scatter!(ax,  θlist, mat[:, 2],label = "ω = 15 meV")
    scatter!(ax,  θlist, mat[:, 3],label = "ω = 27 meV")
    hidexdecorations!(ax, grid = true)
    scale = 1e-15
    ylab2 = L"$i\eta_{xyz}/\tau$ [nm $\mu$A/$(V^2 fs)$]" 
    ax2 = Axis(fig[2,1], xlabel = L"$\theta (^\circ)", ylabel = ylab2,  xlabelsize=24, ylabelsize=24)
    lines!(ax2,  θlist, mat2[:, 1]./scale)
    lines!(ax2,   θlist,  mat2[:, 2]./scale)
    lines!(ax2,   θlist, mat2[:, 3]./scale)
    scatter!(ax2, θlist, mat2[:, 1]./scale, label = "ω = 2.5 meV")
    scatter!(ax2,   θlist, mat2[:, 2]./scale, label = "ω = 15 meV")
    scatter!(ax2,   θlist, mat2[:, 3]./scale, label = "ω = 27 meV")
        
    vlines!(ax, [1.08475], color=:red, linewidth=0.5)
    vlines!(ax2, [1.08475], color=:red, linewidth=0.5)
    ylims!(ax, (-900,1100))
    axislegend(ax, position = (0,1))#, merge = merge, unique = unique)
    axislegend(ax2, position = (1,1))#, fontfamily="Times New Roman")
    save("referee_fig.pdf", fig)
    return fig
end

function lat_suppl4(path1, path2)
    with_theme(merge(theme_latexfonts(), Theme(fontsize = 24))) do #
        suppl4(path1, path2)     
    end
end