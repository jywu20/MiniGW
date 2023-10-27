module MiniGW
    using HDF5
    using StaticArrays
    using LinearAlgebra
    using Base.Threads
    
    export AbstractGWWaveFunction

    abstract type AbstractGWWaveFunction end
    include("bgw.jl")
    
    include("commons.jl")
end