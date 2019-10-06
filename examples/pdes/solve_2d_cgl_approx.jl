using ScatteredCollocation
using Makie
using LinearAlgebra
using Arpack
using SparseArrays, NearestNeighbors
using DifferentialEquations, OrdinaryDiffEq, DiffEqOperators, DiffEqCallbacks
using LabelledArrays
using Printf
using Sundials

function main()
    SC = ScatteredCollocation

    accuracy_order = 1
    Δt = 0.5
    Δx = 0.0075
    tspan = (0.0, 2*96.0)
    T1 = Float64
    solver = Tsit5()#CVODE_BDF(linear_solver=:GMRES)
    timekwargs  = Dict(:save_everystep=>false, :save_start=>false)
    tolkwargs   = Dict(:reltol=>1e-6, :abstol=>1e-6)

    # generate collocation points
    x, is_boundary = SC.griddedSquarePoints(Δx); x = x*100 .- Ref([25.0, 25.0])
    npts = length(x)
    tree = KDTree(reduce(hcat,x))
    idxs = collect(1:npts)
    bidxs = idxs[is_boundary]
    iidxs = setdiff(idxs, bidxs)
    nint = length(iidxs)

    # assemble differential operator matrices
    ∂²x, ∂²y = build(T1, tree, idxs, [[2,0],[0,2]], accuracy_order)
    #id,      = build(T1, tree, idxs, [[0,0]],       accuracy_order)
    id       = spzeros(T1, npts, npts) + I
    ∇² = ∂²x + ∂²y

    # define the restriction and closure operators
    # - closure means interior and boundary
    R = spzeros(T1, nint, npts)
    R[:,iidxs] += I
    B = AffineOperator(id[bidxs,:])
    D = projbc(B, bidxs)
    Q = AffineOperator(spzeros(T1, npts, nint))
    Q.A[iidxs,:] += I
    Q.A[bidxs,:] = D.A
    Q.b[bidxs] = D.b

    # construct the interior pde
    A = AffineOperator(R*∇²)
    println("Sparsity = $(@sprintf("%d",100*(1-fnz(A.A))))%")
    pde(u) = begin
        du = A(Q(u)) + u - (1+1.5im) * u .* (abs.(u).^2)
        return du
    end

    # define initial condition
    p0 = [1im*y[1]+y[2] for y in x] .* exp.(-0.03*norm.(x).^2) # pressure
    u0 = R*p0

    # initialise the plot and define the update function
    p = Node(real.(p0))
    clims = 0.3*maximum(abs.(p0))*[-1, 1]
    scene = AbstractPlotting.Scene(scale_plot=false)
    scplot!(scene, tree, p; colorrange=clims, colormap=:pu_or)
    display(scene)
    plotter(integrator) = begin
        push!(p, real.(Q(integrator.u)))
        sleep(1/60)
    end

    # define the ODE problem and solve
    cb = PeriodicCallback(plotter, Δt; save_positions=(false,false))
    ode = ODEFunction((u, p, t) -> pde(u))
    prob = ODEProblem(ode, u0, tspan)
    u = solve(prob, solver; callback=cb, timekwargs..., tolkwargs...)

    return nothing
end

main()