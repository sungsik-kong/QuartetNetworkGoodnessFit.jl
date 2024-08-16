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
A tree or network topology with edge lengths (in coalescent unit) and inheritance probability specified written in Newick format. This topology can be read in using `readTopology` from `PhyloNetworks`.
```
using PhyloNetworks
net=readTopology("((a:4.201,(b:4.828,((d:2.633,(c:3.641,e:3.412):4.804):4.417)#H5:3.191::0.36):1.061:2.647,#H5:1.33::0.64));")

HybridNetwork, Rooted Network
10 edges
10 nodes: 5 tips, 1 hybrid nodes, 4 internal tree nodes.
tip labels: a, b, d, c, ...
(a:4.201,(b:4.828,#H5:3.191::0.36):1.061,((d:2.633,(c:3.641,e:3.412):4.804):4.417)#H5:1.33::0.64);
```

In case no parameter information is available, a newick without $\tau$ and $\gamma$ can be read using the function `readTopologyrand`. This function will assign random values to the parameters.
```
net=readTopologyrand("((a,(b,((d,(c,e)))#H5)),#H5);")
HybridNetwork, Rooted Network
11 edges
11 nodes: 5 tips, 1 hybrid nodes, 5 internal tree nodes.
tip labels: a, b, d, c, ...
((a:3.283,(b:1.945,((d:2.583,(c:1.909,e:3.283):2.008):2.445)#H5:4.3::0.037):1.013):2.636,#H5:4.824::0.963);
```
### Output
Some updates were made on the function `network_expectedCF` and fice additional options were added, namely `filename`, `savenet`, `savecsv`, `printCFs`, and `symbolic`. By default the latter four options are set `false` behaving identical to the original `network_expectedCF`. 
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

Set `printCFs=true` to see the equation that computes the quartet concordance factor in `network_expectedCF` as a DataFrame.
```
julia> network_expectedCF(net,savenet=false,savecsv=false,printCFs=true,symbolic=false)
Topology: ((a:3.28261,(b:1.94464,((d:2.58257,(c:1.9092,e:3.28323):2.00782):2.44487)#H5:4.29999::0.03679):1.0132):2.6364,#H5:4.82378::0.96321);
(PhyloNetworks.QuartetT{StaticArraysCore.MVector{3, Float64}}[4-taxon set number 1; taxon numbers: 1,2,3,4
data: [0.994571087456838, 0.0027144562715809588, 0.0027144562715809588], 4-taxon set number 2; taxon numbers: 1,2,3,5
data: [0.9992749063953199, 0.00036254680234001707, 0.00036254680234001707], 4-taxon set number 3; taxon numbers: 1,2,4,5
data: [0.994571087456838, 0.0027144562715809588, 0.0027144562715809588], 4-taxon set number 4; taxon numbers: 1,3,4,5
data: [0.0447603628649856, 0.9104792742700288, 0.0447603628649856], 4-taxon set number 5; taxon numbers: 2,3,4,5
data: [0.0447603628649856, 0.9104792742700288, 0.0447603628649856]], ["a", "b", "c", "d", "e"], 15×2 DataFrame
 Row │ Split   CF
     │ String  String
─────┼───────────────────────────────────────────
   1 │ ab|cd   ((1-exp(-2.44487))+(((exp(-2.444…
   2 │ ac|bd   ((((exp(-2.44487)*0.03679)*(0.03…
   3 │ ad|bc   ((((exp(-2.44487)*0.03679)*(0.03…
   4 │ ab|ce   ((1-exp(-4.45269))+(((exp(-4.452…
   5 │ ac|be   ((((exp(-4.45269)*0.03679)*(0.03…
   6 │ ae|bc   ((((exp(-4.45269)*0.03679)*(0.03…
   7 │ ab|de   ((1-exp(-2.44487))+(((exp(-2.444…
   8 │ ad|be   ((((exp(-2.44487)*0.03679)*(0.03…
   9 │ ae|bd   ((((exp(-2.44487)*0.03679)*(0.03…
  10 │ ac|de   (exp(-2.00782)/3)
  11 │ ad|ce   (1-2*exp(-2.00782)/3)
  12 │ ae|cd   (exp(-2.00782)/3)
  13 │ bc|de   (exp(-2.00782)/3)
  14 │ bd|ce   (1-2*exp(-2.00782)/3)
  15 │ be|cd   (exp(-2.00782)/3))
```
To have the equations that does not contain numerical values but parameter names (i.e., t_{#} for edge lengths, r_{#} for inheritance probabilities, and rho for the inheritance correlation, set option `symbolic=true`.
```
julia> network_expectedCF(net,savenet=false,savecsv=false,printCFs=true,symbolic=true)
Topology: ((a:t_{1},(b:t_{2},((d:t_{3},(c:t_{4},e:t_{5}):t_{6}):t_{7})#H5:t_{8}::r_{1}):t_{9}):t_{10},#H5:t_{11}::(1-r_{1}));
(PhyloNetworks.QuartetT{StaticArraysCore.MVector{3, Float64}}[4-taxon set number 1; taxon numbers: 1,2,3,4
data: [0.994571087456838, 0.0027144562715809588, 0.0027144562715809588], 4-taxon set number 2; taxon numbers: 1,2,3,5
data: [0.9992749063953199, 0.00036254680234001707, 0.00036254680234001707], 4-taxon set number 3; taxon numbers: 1,2,4,5
data: [0.994571087456838, 0.0027144562715809588, 0.0027144562715809588], 4-taxon set number 4; taxon numbers: 1,3,4,5
data: [0.0447603628649856, 0.9104792742700288, 0.0447603628649856], 4-taxon set number 5; taxon numbers: 2,3,4,5
data: [0.0447603628649856, 0.9104792742700288, 0.0447603628649856]], ["a", "b", "c", "d", "e"], 15×2 DataFrame
 Row │ Split   CF
     │ String  String
─────┼───────────────────────────────────────────
   1 │ ab|cd   ((1-exp(-t_{7}))+(((exp(-t_{7})*…
   2 │ ac|bd   ((((exp(-t_{7})*r_{1})*(r_{1}*1-…
   3 │ ad|bc   ((((exp(-t_{7})*r_{1})*(r_{1}*1-…
   4 │ ab|ce   ((1-exp(-t_{7}-t_{6}))+(((exp(-t…
   5 │ ac|be   ((((exp(-t_{7}-t_{6})*r_{1})*(r_…
   6 │ ae|bc   ((((exp(-t_{7}-t_{6})*r_{1})*(r_…
   7 │ ab|de   ((1-exp(-t_{7}))+(((exp(-t_{7})*…
   8 │ ad|be   ((((exp(-t_{7})*r_{1})*(r_{1}*1-…
   9 │ ae|bd   ((((exp(-t_{7})*r_{1})*(r_{1}*1-…
  10 │ ac|de   (exp(-t_{6})/3)
  11 │ ad|ce   (1-2*exp(-t_{6})/3)
  12 │ ae|cd   (exp(-t_{6})/3)
  13 │ bc|de   (exp(-t_{6})/3)
  14 │ bd|ce   (1-2*exp(-t_{6})/3)
  15 │ be|cd   (exp(-t_{6})/3))
```
Note the network topology with parameters are written on the top of the output to see what t_{#} or r_{#} corresponds to in the topology. The network topology can be stored at the working directory by setting `savenet=true`. Default network name is `$filename.net.txt`. $filename can be changed by setting the option `filename=$desiredfilename`. By default `$filename=result`.

To access the entire equation, set `savecsv=true` to store the DataFrame in a `csv` file. This `csv` file is stored at the working directory with the default name `$filename.csv`.

## citing

See [`CITATION.bib`](CITATION.bib) for the relevant reference(s).
