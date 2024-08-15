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

## printing and writing qCFs
- Input: A tree or network topology with edge lengths ($\tau$) (in coalescent unit) and inheritance probability ($\gamma$) specified written in Newick format. This topology can be read in using `readTopology` from `PhyloNetworks`. In case no parameter information is available, a newick without $\tau$ and $\gamma$ can be read using the function `readTopologyrand`. See below for example.
```
#parameters available
net=readTopology("((a:4.20458,(b:4.88952,((d:2.63556,(c:3.6052,e:3.42833):4.80495):4.4157)#H5:3.19142::0.36252):1.0612):2.64775,#H5:1.33688::0.63748);")
#parameters unavailable
net=readTopologyrand("((a,(b,((d,(c,e)))#H5)),#H5);")
```

## citing

See [`CITATION.bib`](CITATION.bib) for the relevant reference(s).
