include(dirname(pwd())*"plot_paper_functions.jl")

fontsizelab = 22
fontsizetheme = 22
labsize = 23


# Figure 1
function dosplot(str1, str2)
    d1 = Matrix(CSV.read(str1, DataFrame))
    d2 = Matrix(CSV.read(str2, DataFrame)) 
    fig = Figure(resolution = (997, 810))
    ax = Axis(fig[1,1], xlabel = L"$\omega$ [meV]", ylabel = ylabel, xlabelsize=27, ylabelsize=27)
    lines!(ax, d1[:, 1]*1000, d2[:,2])
    lines!(ax, d1[:, 1]*1000, d1[:,2])
    ylims!(ax, [0,.7])
    xlims!(ax,[0.,40])
    display(fig)
    save("jdosfig1.pdf", fig)
end

# Figure 2

function figure2(strpath2a, strpath2b, strpath2c, strpath2d)
    outer_padding = 30
    n = 600
    color = [RGBf(v, 0, 1-v) for v in 1:-1/(length(strpath2a)-1):0]
    fig = Figure(resolution=(1.4*n, n), figure_padding = (0,0,0,0))
    # ######### PANEL 
    ylab = L"$\sigma_{xyz}$ [nm $\mu$A/$V^2$]" # question if there is TRS this should read as the absortive part there is no imag contribution to the absortive.
    ax2a = Axis(fig[1, 1], xlabel = L"$\omega$ [meV]", ylabel = ylab, xlabelsize=labsize, ylabelsize=labsize)
    ylims = [-0.99e3,0.99e3]
    f2b = protoplot!(ax2a, :x, :y, :z, strpath2a[1],  ylims, color = color[1], hidex = true, hidey = true, ylabel = ylab)
    for i in 1:length(strpath2a)
        if i>1 && i<length(strpath2a)
            f2b =  protoplot!(ax2a, :x, :y, :z, strpath2a[i], ylims, color = color[i],  hidex = true, hidey = true, ylabel = ylab)
        else 
            f2b =  protoplot!(ax2a, :x, :y, :z, strpath2a[i],  ylims, color = color[i], hidex = false, hidey = false, ylabel = ylab)
        end
    end
    Label(fig[1, 1, Top()], halign = :center,  L"$\mu = 0$", fontsize = labsize, padding = (0,0,10,0))
    hidexdecorations!(ax2a, grid = false)
    ax2a.ylabel = ylab
    xlims!(ax2a, [0.2,40])
    
    ### PANEL 
    color = [RGBf(v, 0, 1-v) for v in 1:-1/(length(strpath2b)-1):0]
    Label(fig[1, 2, Top()], halign = :center,  L"$\mu \neq 0$", fontsize = labsize, padding = (0,0,10,0))
    ax2c = Axis(fig[1, 2], xlabel = L"$\omega$ [meV]", yaxisposition = :right, xlabelsize=21, ylabelsize=21)
    ylims = [-2.9e2,1.8e2]
    ylab = L"$\sigma_{xyz}$ [nm $\mu$A/$V^2$]" # question if there is TRS this should read as the absortive part there is no imag contribution to the absortive.
    f2b = protoplot!(ax2c, :x, :y, :z, strpath2b[1],  ylims, color = color[1], hidex = true, hidey = true, ylabel = ylab)
    for i in 1:length(strpath2b)
        if i>1 && i<length(strpath2b)
            println(color[i])
            f2b =  protoplot!(ax2c, :x, :y, :z, strpath2b[i], ylims, color = color[i],  hidex = true, hidey = true, ylabel = ylab)
        else 
            f2b =  protoplot!(ax2c, :x, :y, :z, strpath2b[i],  ylims, color = color[i], hidex = false, hidey = false, ylabel = ylab)
        end
    end
    hidexdecorations!(ax2c, grid = false)
    xlims!(ax2c, [0.,40])

    # ######### PANEL 
    println("length ", length(strpath2c))
    color = [RGBf(v, 0, 1-v) for v in 1:-1/(length(strpath2c)-1):0]
    ax2b = Axis(fig[2, 1], xlabel = L"$\omega$ [meV]",  ylabel = ylab, xlabelsize=labsize, ylabelsize=labsize)
    ylims = [-0.9e-15,0.9e-15]
    scale = 1e-15
    ylab = L"$i\eta_{xyz}/\tau$ [nm $\mu$A/$(V^2 fs)$]"
    f2b = protoplot!(ax2b, :x, :y, :z, strpath2c[1],  ylims, color = color[1], hidex = true, hidey = true, ylabel = ylab, scale = scale)
    for i in 1:length(strpath2c)
        if i>1 && i<length(strpath2c)
            f2b =  protoplot!(ax2b, :x, :y, :z, strpath2c[i], ylims, color = color[i],  hidex = true, hidey = true, ylabel = ylab, scale = scale)
        else 
            f2b =  protoplot!(ax2b, :x, :y, :z, strpath2c[i],  ylims, color = color[i], hidex = false, hidey = false, ylabel = ylab, scale = scale)
        end
    end
    # Label(fig[2, 1, Top()], halign = :left, L"\times 10^{-15}")
    ax2b.ylabel = ylab
    xlims!(ax2b, [0.,40])


    ### PANEL 
    color = [RGBf(v, 0, 1-v) for v in 1:-1/(length(strpath2d)-1):0]
    ax2d = Axis(fig[2, 2], xlabel = L"$\omega$ [meV]", yaxisposition = :right, xlabelsize=labsize, ylabelsize=labsize)
    ylims = [-2e-16,3e-16]
    ylab = L"$i\eta_{xyz}/\tau$ [nm $\mu$A/$(V^2 fs)$]"
    f2b = protoplot!(ax2d, :x, :y, :z, strpath2d[1],  ylims, color = color[1], hidex = true, hidey = true, ylabel = ylab, scale = scale)
    for i in 1:length(strpath2d)
        if i>1 && i<length(strpath2d)
            f2b =  protoplot!(ax2d, :x, :y, :z, strpath2d[i], ylims, color = color[i],  hidex = true, hidey = true, ylabel = ylab, scale = scale)
        else 
            f2b =  protoplot!(ax2d, :x, :y, :z, strpath2d[i],  ylims, color = color[i], hidex = false, hidey = false, ylabel = ylab, scale = scale)
        end
    end
    xlims!(ax2d, [0.,40])
    cb = Colorbar(fig[3, 1:2], label = "", colormap = ["red", "blue"], vertical = false, flipaxis = false, 
        limits = (1.05,1.12), ticks = ([1.05, 1.0845,1.12], [L"$1.05^\circ$",L"$\theta_M",L"$1.12^\circ$"])) 
    cb.width = Relative(1/3)
    rowgap!(fig.layout, 5)
    Label(fig[1, 1, TopLeft()], "(a)", fontsize=fontsizelab, padding=(0, -30, 8, 0))
    Label(fig[1, 2, TopLeft()], "(b)", fontsize=fontsizelab, padding=(0, -20, 8, 0))
    Label(fig[2, 1, TopLeft()], "(c)", fontsize=fontsizelab, padding=(0, -30, 8, 0))
    Label(fig[2, 2, TopLeft()], "(d)", fontsize=fontsizelab, padding=(0, -20, 8, 0))
    display(fig)
    save("fig2.pdf", fig)
end

function lat_fig2(strpath2a, strpath2b, strpath2c, strpath2d)
    with_theme(merge(theme_latexfonts(), Theme(fontsize = fontsizetheme))) do #
        figure2(strpath2a, strpath2b, strpath2c, strpath2d)     
    end
end



# Figure 3


function figure3(strpath3a, strpath3b, strpath3c, strpath3d)
    outer_padding = 30
    n = 420
    fig = Figure(resolution = (3.85*n,n), figure_padding = (0,12,0,0))
    aspect = nothing
    # ######### PANEL 
    ylab = L"$\sigma_{yyy}$ [nm $\mu$A/$V^2$]" 
    ga = fig[1, 1] = GridLayout()
    ax1a = Axis(ga[1,1], xlabel = L"$\omega$ [meV]", ylabel = ylab, xlabelsize=labsize, ylabelsize=labsize, aspect = aspect)#, title = L"Point group: $C_3$, $C_{2y}")
    ylims = [-1e5,1e5]
    scale = 1e5
    f2b = protoplot!(ax1a, :y, :y, :y, strpath3a[1],  ylims, color = :red, hidex = true, hidey = true, ylabel = ylab, scale = scale)
    f2b = protoplot!(ax1a, :y, :y, :y, strpath3a[2],  ylims, color = :blue, hidex = true, hidey = true, ylabel = ylab, scale = scale, linewidth = 1)
    f2b = protoplot!(ax1a, :y, :y, :y, strpath3a[3],  ylims, color = :red, style = :dash, hidex = true, hidey = true, ylabel = ylab, scale = scale)
    f2b = protoplot!(ax1a, :y, :y, :y, strpath3a[4],  ylims, color = :blue, style = :dash, hidex = true, hidey = true, ylabel = ylab, scale = scale)
    f2b = protoplot!(ax1a, :y, :y, :y, strpath3a[5],  ylims, color = :red,  style = :dot, hidex = true, hidey = true, ylabel = ylab, scale = scale)
    f2b = protoplot!(ax1a, :y, :y, :y, strpath3a[6],  ylims, color = :blue, style = :dot, hidex = true, hidey = true, ylabel = ylab, scale = scale)
    ax1a.ylabel = ylab
    xlims!(ax1a, (0,40))
    Label(ga[1, 1, Top()], halign = :left, L"\times 10^{5}",  padding = (0, 0, 0, 0))
    Label(ga[1, 1, Top()], halign = :center, L"$\Delta_1 \sigma_z$", fontsize = 25, padding = (0, 0, 0, 0))
    
    ### PANEL 
    gb = fig[1, 2] = GridLayout()
    ylims = [-2e4,2e4]
    ylab = L"$\sigma_{xxx}$ [nm $\mu$A/$V^2$]" 
    ax1b = Axis(gb[1, 1], xlabel = L"$\omega$ [meV]", ylabel = ylab, xlabelsize=labsize, ylabelsize=labsize, aspect = aspect)#, title = L"Point group: $C_3$, $C_{2x}")#L"$\Delta_2 \sigma_z \tau_z$")
    scale = 1e4
    f2b = protoplot!(ax1b, :y, :y, :y, strpath3b[1],  ylims, color = :red, hidex = true, hidey = true, ylabel = ylab, scale = scale)
    f2b = protoplot!(ax1b, :y, :y, :y, strpath3b[2],  ylims, color = :blue, hidex = true, hidey = true, ylabel = ylab, scale = scale)
    f2b = protoplot!(ax1b, :y, :y, :y, strpath3b[3],  ylims, color = :red, style = :dash, hidex = true, hidey = true, ylabel = ylab, scale = scale)
    f2b = protoplot!(ax1b, :y, :y, :y, strpath3b[4],  ylims, color = :blue, style = :dash, hidex = true, hidey = true, ylabel = ylab, scale = scale, linewidth = 1)
    f2b = protoplot!(ax1b, :y, :y, :y, strpath3b[5],  ylims, color = :red,  style = :dot, hidex = true, hidey = true, ylabel = ylab, scale = scale, linewidth = 1)
    f2b = protoplot!(ax1b, :y, :y, :y, strpath3b[6],  ylims, color = :blue, style = :dot, hidex = true, hidey = true, ylabel = ylab, scale = scale)
    ax1b.ylabel = ylab
    xlims!(ax1b, (0,40))
    Label(gb[1, 1, Top()], halign = :left, L"\times 10^{4}",  padding = (0, 0, 0, 0))
    Label(gb[1, 1, Top()], halign = :center, L"$\Delta_2 \sigma_z \tau_z$", fontsize = 25, padding = (0, 0, 0, 0))

    # ######### PANEL 
    gc = fig[1, 3] = GridLayout()
    ylims = [-2.4e2,1.5e2]
    ylab = L"$\sigma_{xxz}$ [nm $\mu$A/$V^2$]" 
    ax1c = Axis(gc[1, 1], xlabel = L"$\omega$ [meV]", ylabel = ylab, xlabelsize=labsize, ylabelsize=labsize, aspect = aspect)#, title = L"Point group: $C_3$, $C_{2z}")#L"$\Delta_2 \sigma_z \tau_z$")
    scale = 1e2
    f2b = protoplot!(ax1c, :y, :y, :y, strpath3c[1],  ylims, color = :red, hidex = true, hidey = true, ylabel = ylab, scale = scale)
    f2b = protoplot!(ax1c, :y, :y, :y, strpath3c[2],  ylims, color = :blue,  hidex = true, hidey = true, ylabel = ylab, scale = scale)
    f2b = protoplot!(ax1c, :y, :y, :y, strpath3c[3],  ylims, color = :red, style = :dash, hidex = true, hidey = true, ylabel = ylab, scale = scale)
    f2b = protoplot!(ax1c, :y, :y, :y, strpath3c[4],  ylims, color = :blue, style = :dash, hidex = true, hidey = true, ylabel = ylab, scale = scale)
    f2b = protoplot!(ax1c, :y, :y, :y, strpath3c[5],  ylims, color = :red,  style = :dot,hidex = true, hidey = true, ylabel = ylab, scale = scale)
    f2b = protoplot!(ax1c, :y, :y, :y, strpath3c[6],  ylims, color = :blue, style = :dot, hidex = true, hidey = true, ylabel = ylab, scale = scale)
    ax1c.ylabel = ylab
    xlims!(ax1c, (0, 40))
    Label(gc[1, 1, Top()], halign = :left, L"\times 10^{2}",  padding = (0, 0, 0, 0))
    Label(gc[1, 1, Top()], halign = :center, L"$\Delta_3 \tau_z$", fontsize = 25, padding = (0, 0, 0, 0))
    gd = fig[1, 4] = GridLayout()
    ylims = [-1e-16,2.4e-16]
    ylab = L"$ i \eta_{xxz}/\tau$ [nm $\mu$A/$(V^2 fs)$]" 
    ax1d = Axis(gd[1, 1], xlabel = L"$\omega$ [meV]", ylabel = ylab, yaxisposition = :left, xlabelsize=labsize, ylabelsize=labsize, aspect = aspect)#, title = L"Point group: $C_3$, $C_{2z}")#L"$\Delta_2 \sigma_z \tau_z$")
    scale = 1e-16
    f2b = protoplot!(ax1d, :y, :y, :y, strpath3d[1],  ylims, color = :red, hidex = true, hidey = true, ylabel = ylab, scale = scale)
    f2b = protoplot!(ax1d, :y, :y, :y, strpath3d[2],  ylims, color = :blue,  hidex = true, hidey = true, ylabel = ylab, scale = scale)
    f2b = protoplot!(ax1d, :y, :y, :y, strpath3d[3],  ylims, color = :red, style = :dash,hidex = true, hidey = true, ylabel = ylab, scale = scale)
    f2b = protoplot!(ax1d, :y, :y, :y, strpath3d[4],  ylims, color = :blue, style = :dash, hidex = true, hidey = true, ylabel = ylab, scale = scale, linewidth = 1)
    f2b = protoplot!(ax1d, :y, :y, :y, strpath3d[5],  ylims, color = :red, style = :dot,hidex = true, hidey = true, ylabel = ylab, scale = scale, linewidth = 1)
    f2b = protoplot!(ax1d, :y, :y, :y, strpath3d[6],  ylims, color = :blue, style = :dot, hidex = true, hidey = true, ylabel = ylab, scale = scale)
    ax1d.ylabel = ylab
    xlims!(ax1d, (0,40))
    Label(gd[1, 1, Top()], halign = :left, L"\times 10^{-1}", tellwidth=false, tellheight=false, padding = (0, 0, 0, 0))
    Label(gd[1, 1, Top()], halign = :center, L"$\Delta_3 \tau_z$", fontsize = 25, padding = (0, 0, 0, 0))
    Label(ga[1, 1, TopLeft()], "(a)", fontsize=fontsizelab,  padding=(0, 40, 0, 0))
    Label(gb[1, 1, TopLeft()], "(b)", fontsize=fontsizelab,  padding=(0, 35, 0, 0))
    Label(gc[1, 1, TopLeft()], "(c)", fontsize=fontsizelab,  padding=(0, 35, 0, 0))
    Label(gd[1, 1, TopLeft()], "(d)", fontsize=fontsizelab,  padding=(0, 35, 0, 0))
    legenditems = [LineElement(color = :gray,  linewidth = 2),
    LineElement(color = :gray, linestyle = :dash, linewidth = 2),
    LineElement(color = :gray, linestyle = :dot, linewidth = 2),
    MarkerElement(color = :red, marker=:rect, markersize = 15, strokecolor = :black),
    MarkerElement(color = :blue, marker=:rect, markersize = 15, strokecolor = :black)]
    g0 = fig[2,1:4] = GridLayout()
    Legend(g0[1, 1], legenditems, [L"$\mu = 0$", L"$\mu > 0$", L"$\mu<0",L"$\theta<\theta_M$",L"$\theta>\theta_M"], 
    orientation = :horizontal, tellwidth = false, tellheight =true, framewidth = .5)#, padding = (0,0,0,0))
    rowgap!(fig.layout, 10)
    resize_to_layout!(fig)
    display(fig)
    save("fig3.pdf", fig)
end

function lat_fig3(strpath3a, strpath3b, strpath3c, strpath3d)
    with_theme(merge(theme_latexfonts(), Theme(fontsize = fontsizetheme))) do #
        figure3(strpath3a, strpath3b, strpath3c, strpath3d)     
    end
end

# Figure 4

function figure4(normalmass, kaplanmass, szvafekposmass, szvafeknegmass)
    ylims = [-2.51e4,2.51e4]
    ylab = L"$\sigma_{yyy}$ [nm $\mu$A/$V^2$]"
    scale = 1e4
    f4 = Figure(resolution = (700,450), figure_padding = (0,20,1,1))
    f4 = shiftsummedvalleys(f4, :y, :y, :y, normalmass,  ylims, color = :black, hidex = true, hidey = true, ylabel = ylab, scale = scale)
    f4 = shiftsummedvalleys(f4, :y, :y, :y, kaplanmass, ylims, color = :gray,  hidex = true, hidey = true, ylabel = ylab, scale = scale)
    f4 = shiftsummedvalleys(f4, :y, :y, :y, szvafekposmass, ylims, color = :red, hidex = true, hidey = true, ylabel = ylab, scale = scale)
    f4 = shiftsummedvalleys(f4, :y, :y, :y, szvafeknegmass,  ylims, color = :red, style = :dash, hidex = false, hidey = false, ylabel = ylab, scale = scale)
    terms = [L"$\Delta_1$",  L"$\Delta_1$, $\omega_3$",L"$\Delta_1$, -$\omega_3$",  L"$\Delta_1$, $\Delta_2$"]
    ccolors = [:black, :lightblue, :red] # :blue, :blue];
    elem1 = LineElement(color = :black, linewidth = 2)
    elem4 = LineElement(color = :gray, linewidth = 2)
    elem2 = LineElement(color =  :red, linewidth = 2)
    elem3 = LineElement(color = :red, linestyle = :dash, linewidth = 2)
    ax = Axis(f4[1,1],xlabelsize=34, ylabelsize=34)
    hidedecorations!(ax;)
    axislegend(ax, [elem1, elem2, elem3, elem4], terms; position= (0.35,1))
    Label(f4[1, 1], halign = :left, L"\times 10^{4}", tellwidth=false, tellheight=false, padding = (0,0,520,0))
    display(f4)
    save("fig4.pdf", f4)
end

function lat_fig4(normalmass, kaplanmass, szvafekposmass, szvafeknegmass)
    with_theme(merge(theme_latexfonts(), Theme(fontsize = 24))) do
        figure4(normalmass, kaplanmass, szvafekposmass, szvafeknegmass)
    end
end
