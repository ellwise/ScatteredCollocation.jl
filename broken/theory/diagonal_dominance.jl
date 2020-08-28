using ScatteredCollocation
SC = ScatteredCollocation
using Makie

numdims = 1
k = 3
num_poly = 8
diff_order = [2]
max_nodes = num_poly*2*2#*2#*2*2*2*2*2

C = fill(NaN, max_nodes-num_poly+1, max_nodes)
C1 = fill(NaN, max_nodes-num_poly+1, max_nodes)
C2 = fill(NaN, max_nodes-num_poly+1, max_nodes)
C3 = fill(NaN, 2*num_poly)
for j = num_poly:max_nodes
    x0 = [0.0]
    x = hcat(Float64.(collect(-1:num_poly)))
    x = Float64.(collect(-1:j-2))
    #x = vcat(x[1:num_poly], x[end-num_poly+1:end])
    x = [[y] for y in x]
    #if j==2*num_poly
    #    c = SC.expandingfd(x0, x, diff_order; k=k, npoly=num_poly)
    #    C3 .= c[:]
    #end
    c = SC.wlspoly(x0, x, diff_order; k=k, npoly=num_poly)
    c1, c2 = SC.splitfd(x0, x, diff_order; k=k, npoly=num_poly)
    C[j-num_poly+1,1:j] = c
    C1[j-num_poly+1,1:j] = c1
    C2[j-num_poly+1,1:j] = c2
    #C[j-num_poly+1,1:num_poly] = c[1:num_poly]
    #C[j-num_poly+1,j-num_poly+1:j] = c[num_poly+1:end]
    #C1[j-num_poly+1,1:num_poly] = c1[1:num_poly]
    #C1[j-num_poly+1,j-num_poly+1:j] = c1[num_poly+1:end]
    #C2[j-num_poly+1,1:num_poly] = c2[1:num_poly]
    #C2[j-num_poly+1,j-num_poly+1:j] = c2[num_poly+1:end]
end

println(C[2*num_poly,1:num_poly])
println(C[end,1:num_poly])
println(C3[1:num_poly])

xs = -1:max_nodes-2
ys = num_poly:max_nodes

clims = maximum(abs.(filter(!isnan, C)))*[-1, 1]
scene0 = AbstractPlotting.Scene(); heatmap!(scene0, C', colorrange=clims, colormap=:pu_or)
scene1 = AbstractPlotting.Scene(); heatmap!(scene1, C1', colorrange=clims, colormap=:pu_or)
scene2 = AbstractPlotting.Scene(); heatmap!(scene2, C2', colorrange=clims, colormap=:pu_or)
scene = vbox(scene0, hbox(scene1, scene2))
display(scene)

#scene = AbstractPlotting.Scene()
#plot!(scene, C[:, 2])
#display(scene)