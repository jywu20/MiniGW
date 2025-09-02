module MiniGW
    using HDF5
    using StaticArrays
    using LinearAlgebra
    using Base.Threads
    using Statistics
    
    export AbstractGWWaveFunction,
        AbstractBSEWaveFunction

    abstract type AbstractGWWaveFunction end
    abstract type AbstractBSEWaveFunction end

    const Ry_BGW =  13.6056925

    include("bgw-gw.jl")
    include("bgw-bse.jl")
    
    include("commons.jl")
end