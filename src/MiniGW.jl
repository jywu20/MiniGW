module MiniGW
    using HDF5
    using StaticArrays
    using LinearAlgebra
    
    export AbstractGWWaveFunction

    abstract type AbstractGWWaveFunction end
    include("bgw.jl")
    
    include("commons.jl")
end