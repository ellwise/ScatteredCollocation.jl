using ScatteredCollocation, NearestNeighbors, SparseArrays
using DifferentialEquations, OrdinaryDiffEq, DiffEqOperators, DiffEqCallbacks
using Makie
using Sundials
using LinearAlgebra
using LabelledArrays
#using RecursiveArrayTools

function main()

    npts = 200
    xspan = (-2, 2)
    tspan = (0.0, 4.0)
    Δt = 0.05 # for plotting
    tol = 1e-6
    k = 3
    accuracy_order = 3
    T1 = Float64

    # generate collocation points
    x = collect(range(xspan[1], stop=xspan[2], length=npts))
    tree = KDTree(reduce(hcat,x))
    idxs = collect(1:npts)
    bidxs = [1, npts]
    iidxs = setdiff(idxs, bidxs)
    nint = length(iidxs)

    # assemble interior/boundary operators and define RHSs
    ∂ₓ, = build(T1, tree, idxs, [[1]], accuracy_order)
    id, = build(T1, tree, idxs, [[0]], accuracy_order)
    A = AffineOperator(∂ₓ)
    B = AffineOperator([id[bidxs[1],:]'; ∂ₓ[bidxs[2],:]'])
    D = projbc(B, bidxs)

    # define the restriction and closure operators
    # - closure means interior and boundary
    R = spzeros(T1, nint, npts)
    R[:,iidxs] += I
    Q = AffineOperator(spzeros(T1, npts, nint))
    Q.A[iidxs,:] += I
    Q.A[bidxs,:] = D.A
    Q.b[bidxs] = D.b

    # define initial condition
    p0 = exp.(-(x.+1).^2 ./ (0.2^2)) # pressure
    v0 = zero(p0)                    # velocity
    #u = @LArray [p0[iidxs]; v0] (p=1:npts-2, v=(1:npts).+(npts-2))
    #pde!(du, u) = begin
    #    du.p .= -R*A(u.v)
    #    du.v .= -A(Q(u.p))
    #    return nothing
    #end
    u = ArrayPartition(p0[iidxs], v0)
    pde!(du, u) = begin
        du.x[1] .= -R*A(u.x[2])
        du.x[2] .= -A(Q(u.x[1]))
        return nothing
    end

    # initialise the plot and define the update function
    p = Node(p0)
    v = Node(v0)
    scene = AbstractPlotting.Scene()
    scplot!(scene, tree, p; color=:blue)
    scplot!(scene, tree, v; color=:red)
    display(scene)
    function plotter(integrator)
        #pnew = Q(integrator.u.p)
        #vnew = integrator.u.v
        pnew = Q(integrator.u.x[1])
        vnew = integrator.u.x[2]
        push!(p, pnew)
        push!(v, vnew)
        sleep(1/60)
    end
    cb = PeriodicCallback(plotter, Δt; save_positions=(false,false))

    # define the ODE problem and solve
    ode! = ODEFunction((du,u,p,t)->pde!(du,u))
    prob = ODEProblem(ode!, u, tspan)
    timekwargs  = Dict(:save_everystep=>false, :save_start=>false)
    tolkwargs   = Dict(:reltol=>tol, :abstol=>tol)
    #otherkwargs = Dict(:alg_hints=>[:memorybound, :auto], :callback=>cb)
    otherkwargs = Dict(:callback=>cb)
    u = solve(prob, Tsit5(); timekwargs..., tolkwargs..., otherkwargs...)
    #pnew = Q(integrator.u.p)
    #vnew = integrator.u.v
    pnew = Q(u.x[1])
    vnew = u.x[2]
    push!(p, pnew)
    push!(v, vnew)

    return nothing
end

main()