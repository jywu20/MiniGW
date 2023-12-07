module MiniGW
    using HDF5
    using StaticArrays
    using LinearAlgebra
    using Base.Threads
    using Statistics
    
    export AbstractGWWaveFunction

    abstract type AbstractGWWaveFunction end
    include("bgw.jl")
    
    include("commons.jl")
end