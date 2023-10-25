module MiniGW
    using HDF5

    export AbstractGWWaveFunction

    abstract type AbstractGWWaveFunction end
    include("bgw.jl")
    
    include("commons.jl")
end