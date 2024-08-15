[![doc: stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaPhylo.github.io/QuartetNetworkGoodnessFit.jl/stable)
[![doc: dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaPhylo.github.io/QuartetNetworkGoodnessFit.jl/dev)
[![Build status](https://github.com/JuliaPhylo/QuartetNetworkGoodnessFit.jl/workflows/CI/badge.svg?branch=master)](https://github.com/JuliaPhylo/QuartetNetworkGoodnessFit.jl/actions)
[![Coverage](https://codecov.io/gh/JuliaPhylo/QuartetNetworkGoodnessFit.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaPhylo/QuartetNetworkGoodnessFit.jl)
[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![PkgEval](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/Q/QuartetNetworkGoodnessFit.svg)](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/report.html)

## overview

`QuartetNetworkGoodnessFit`, aka "Quarnet GoF" or simply "QGoF",
is a Julia package for phylogenetic networks analyses using four-taxon subsets.
It includes tools to measure the
goodness of fit of a phylogenetic network to data on subsets of 4 tips.
It depends on the [PhyloNetworks](https://github.com/JuliaPhylo/PhyloNetworks.jl)
package.

## printing and writing qCFs equations
### Input
A tree or network topology with edge lengths (in coalescent unit) and inheritance probability specified written in Newick format. This topology can be read in using `readTopology` from `PhyloNetworks`. In case no parameter information is available, a newick without $\tau$ and $\gamma$ can be read using the function `readTopologyrand`. See below for example.
```
#parameters available
net=readTopology("((a:4.20458,(b:4.88952,((d:2.63556,(c:3.6052,e:3.42833):4.80495):4.4157)#H5:3.19142::0.36252):1.0612):2.64775,#H5:1.33688::0.63748);")
```
or
```
#parameters unavailable
net=readTopologyrand("((a,(b,((d,(c,e)))#H5)),#H5);")
HybridNetwork, Rooted Network
11 edges
11 nodes: 5 tips, 1 hybrid nodes, 5 internal tree nodes.
tip labels: a, b, d, c, ...
((a:3.283,(b:1.945,((d:2.583,(c:1.909,e:3.283):2.008):2.445)#H5:4.3::0.037):1.013):2.636,#H5:4.824::0.963);
```
### Output
network_expectedCF(net,savenet=false,savecsv=false,printCFs=true,symbolic=true)
Some updates were made on the function `network_expectedCF` and fice additional options were added, namely `filename`, `savenet`, `savecsv`, `printCFs`, and `symbolic`. By default the latter four options are set `false`, which results in the original `network_expectedCF`. 
```
julia> network_expectedCF(net)
(PhyloNetworks.QuartetT{StaticArraysCore.MVector{3, Float64}}[4-taxon set number 1; taxon numbers: 1,2,3,4
data: [0.994571087456838, 0.0027144562715809588, 0.0027144562715809588], 4-taxon set number 2; taxon numbers: 1,2,3,5
data: [0.9992749063953199, 0.00036254680234001707, 0.00036254680234001707], 4-taxon set number 3; taxon numbers: 1,2,4,5
data: [0.994571087456838, 0.0027144562715809588, 0.0027144562715809588], 4-taxon set number 4; taxon numbers: 1,3,4,5
data: [0.0447603628649856, 0.9104792742700288, 0.0447603628649856], 4-taxon set number 5; taxon numbers: 2,3,4,5
data: [0.0447603628649856, 0.9104792742700288, 0.0447603628649856]], ["a", "b", "c", "d", "e"])

julia> network_expectedCF(net,savenet=false,savecsv=false,printCFs=false,symbolic=false)
(PhyloNetworks.QuartetT{StaticArraysCore.MVector{3, Float64}}[4-taxon set number 1; taxon numbers: 1,2,3,4
data: [0.994571087456838, 0.0027144562715809588, 0.0027144562715809588], 4-taxon set number 2; taxon numbers: 1,2,3,5
data: [0.9992749063953199, 0.00036254680234001707, 0.00036254680234001707], 4-taxon set number 3; taxon numbers: 1,2,4,5
data: [0.994571087456838, 0.0027144562715809588, 0.0027144562715809588], 4-taxon set number 4; taxon numbers: 1,3,4,5
data: [0.0447603628649856, 0.9104792742700288, 0.0447603628649856], 4-taxon set number 5; taxon numbers: 2,3,4,5
data: [0.0447603628649856, 0.9104792742700288, 0.0447603628649856]], ["a", "b", "c", "d", "e"])
```


## citing

See [`CITATION.bib`](CITATION.bib) for the relevant reference(s).
