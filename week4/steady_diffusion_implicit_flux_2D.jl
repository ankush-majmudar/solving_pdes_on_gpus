using Plots, Plots.Measures, Printf
default(size=(1200, 800), framestyle=:box, label=false, grid=false, margin=10mm, lw=6, labelfontsize=20, tickfontsize=20, titlefontsize=24)

@views function steady_diffusion_implicit_flux_1D()
    # physics
    lx, ly     = 20.0, 20.0
    dc      = 1.0
    re      = 2π 
    ρ       = (lx / (dc * re))^2
    # numerics
    nx, ny      = 100, 101
    ϵtol    = 1e-8
    maxiter = 100nx
    ncheck  = ceil(Int, 0.25nx)
    # derived numerics
    dx, dy  = lx / nx, ly / ny
    xc, yc  = LinRange(dx / 2, lx - dx / 2, nx), LinRange(dy / 2, ly - dy / 2, ny)
    dτ      = min(dx, dy) / sqrt(1 / ρ) / sqrt(2)
    # array initialisation
    # array initialisation
    C       = @. 1.0 + exp(-(xc - lx / 4)^2 - (yc' - ly / 4)^2) - xc / lx
    qx, qy  = zeros(nx-1, ny), zeros(nx, ny-1)
    C_i     = copy(C)
    # iteration loop
    iter = 1; err = 2ϵtol; iter_evo = Float64[]; err_evo = Float64[]
    while err >= ϵtol && iter <= maxiter
        # see derivation in lecture
        qx .-= dτ ./ (ρ * dc .+ dτ) .* (qx .+ dc .* diff(C, dims=1) ./ dx)
        qy .-= dτ ./ (ρ * dc .+ dτ) .* (qy .+ dc .* diff(C, dims=2) ./ dy)
        # qx .= (qx .* ρ/dτ .- diff(C) ./ dx) ./ (ρ/dτ + 1.0/dc)
        C[2:end-1, 2:end-1] .-= dτ .* (diff(qx[:, 2:end-1], dims=1) ./ dx + diff(qy[2:end-1, :], dims=2) ./ dy)
        if iter % ncheck == 0
            err = maximum(abs.(diff(dc .* diff(C[:, 2:end-1], dims=1) ./ dx, dims=1) ./ dx .+
                               diff(dc .* diff(C[2:end-1, :], dims=2) ./ dx, dims=2) ./ dy))
            push!(iter_evo, iter / nx); push!(err_evo, err)

            p1 = heatmap(xc, yc, C'; xlims=(0, lx), ylims=(0, ly), clims=(0,1), aspect_ratio=1,
                      xlabel="lx", ylabel="Concentration", title="iter/nx=$(round(iter/nx,sigdigits=3))")
            p2 = plot(iter_evo, err_evo; xlabel="iter/nx", ylabel="err",
                      yscale=:log10, grid=true, markershape=:circle, markersize=10)
            display(plot(p1, p2; layout=(2, 1)))
        end
        iter += 1
    end
end

steady_diffusion_implicit_flux_1D()