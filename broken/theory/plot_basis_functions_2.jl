using ScatteredCollocation, Makie
SC = ScatteredCollocation

Δx = 0.05
x, bnd = SC.concentricDiscPoints(Δx)

maxorder = 3
numdims = 2
npoly = SC.countmonomials(maxorder; numdims=numdims)
m = SC.MonomialExponents{numdims,Int}(npoly)

# accorder = 1, k=3

# identity
ks = [3, 3, 3, 3, 3, 3]
islogs = [false, false, false, false, false, false]
exponents = [[0,0], [0,1], [1,0], [0,2], [1,1], [2,0]]

# first-derivative
ks = [1, 1, 3, 1, 1, 3, 1, 3, 1]
islogs = [false, false, false, false, false, false, false, false, false]
exponents = [[1,0], [1,1], [0,0], [2,0], [1,2], [0,1], [2,1], [1,0], [3,0]]

# second_derivative
ks = [1, -1, 1, -1, 1, 1, -1, 1, -1, 1, 1, -1, 3, 1, 1, -1, 1, -1, 1, 1, -1, 3, 1, 1, -1, 3, 1, 1, -1]
islogs = falses(length(ks))
exponents = [[0,0], [2,0], [0,1], [2,1], [1,0], [1,0], [3,0], [0,2], [2,2], [1,1], [1,1], [3,1], [0,0], [2,0], [2,0], [4,0],
             [0,3], [2,3], [1,2], [1,2], [3,2], [0,1], [2,1], [2,1], [4,1], [1,0], [3,0], [3,0], [5,0]]

ks = [0, 0, 0]
islogs = falses(length(ks))
exponents = [[0,0], [1,0], [0,1]]

numterms = length(ks)
diff_order = [0,0]

numrows = floor(Int, sqrt(numterms))
numcols = ceil(Int, numterms/numrows)

scenes = []
idxs = LinearIndices((numrows, numcols))
p = 0
for m = 1:numrows
    row = []
    for n = 1:numcols
        p = idxs[m, n]
        if p > numterms
            break
        else
            k = ks[p]
            islog = islogs[p]
            exponent = exponents[p]
            u = SC.phi.(x, k, Ref(exponent), islog, Ref(diff_order))
            clims = 1.0*maximum(abs.(u))*[-1, 1]
            scene = AbstractPlotting.Scene(scale_plot=false)
            scatter!(scene, map(y->y[1], x), map(y->y[2], x), markersize = 1.1*Δx, color=u,
                colorrange=clims, colormap=:pu_or, show_axis=false)
            push!(row, scene)
        end
    end
    push!(scenes, row)
    if p > numterms
        break
    end
end
tot = hbox([vbox([scene for scene in row]) for row in scenes])
display(tot)
