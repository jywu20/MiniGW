using Plots
using LinearAlgebra
using Base.Iterators

q_list = map(t -> [t...], collect(product(0.0:0.05:0.5, 0.0:0.05:0.5))) 
q_list = reshape(q_list, length(q_list))

v(q, G) = 1 / (norm(q + G)^2 + 0.001^2)
weight_factor_list = (q -> v(q, [10, 4])).(q_list)

phase_factor_list = exp.(2π * im * rand(length(q_list))) 

sum(weight_factor_list .* phase_factor_list) / sum(weight_factor_list)
