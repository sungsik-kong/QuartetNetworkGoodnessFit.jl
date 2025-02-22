module QuartetNetworkGoodnessFit

using DataFrames
using Distributed
using NLopt
using PhyloCoalSimulations: simulatecoalescent
using PhyloNetworks
using Random: seed!
using SharedArrays
using SpecialFunctions: loggamma
using StaticArrays
using Statistics: mean, median
using StatsFuns: normccdf, chisqccdf, betacdf, betaccdf
using CSV
using Combinatorics
#using Distributions

const PN = PhyloNetworks # for easier use of internals

export
quarnetGoFtest!,
network_expectedCF,
ticr,
ticr!,
ultrametrize!,
readTopologyrand,
test,
plot_ntwk_with_Symbolic_Names,
makeEdgeLabel

include("utils.jl")
include("quarnetconcordancefactors.jl")
include("ticr.jl")
include("quarnetGoF.jl")

end # module
