using HDF5
using Plots
using LaTeXStrings

fid = h5open("/pscratch/sd/j/jywu/WTe2-xy-relaxed/1-epsilon/epsmat.h5")
eps_dataset = fid["/mats/matrix"]

q_idx = 12 
G_idx_range = 1 : 100 
heatmap(G_idx_range, G_idx_range, 
    abs.(eps_dataset[1, G_idx_range, G_idx_range, 1, 1, q_idx])',
    aspect_ratio = :equal, 
    xlabel = L"\mathbf{G}",
    ylabel = L"\mathbf{G}\prime",
    xlims = (minimum(G_idx_range) , maximum(G_idx_range) ), 
    ylims = (minimum(G_idx_range) , maximum(G_idx_range) ), 
)

# So indeed χ_GG′ is quite diagonal when G is large