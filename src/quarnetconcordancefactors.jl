const global dpoints=5 #decimal points for all parameters when randomly generated

"""
function network_expectedCF(net0::HybridNetwork; 
                            showprogressbar=false, 
                            inheritancecorrelation=0, 
                            printCFs=false::Bool,
                            savenet=false::Bool,
                            symbolic=false::Bool,
                            savecsv=false::Bool,
                            macaulay=false::Bool,
                            matlab=false::Bool,
                            suffix=""::AbstractString)
"""
function network_expectedCF(net0::HybridNetwork; 
                            showprogressbar=false, 
                            inheritancecorrelation=0, 
                            printCFs=false::Bool,
                            savenet=false::Bool,
                            symbolic=false::Bool,
                            savecsv=false::Bool,
                            macaulay=false::Bool,
                            matlab=false::Bool,
                            suffix=""::AbstractString)
    
    #----------filename----------#
    nleaf=length(net0.leaf)
    nreticulate=length(net0.hybrid)
    filename="n$(nleaf)h$(nreticulate)"
    if !isempty(suffix) filename=filename*".$(suffix)" end

    #---------forcing things to work---------#
    if(symbolic) 
        net=deepcopy(net0)
        net=readTopologyrand(net)
    else
        net=net0 
    end
    if(matlab) macaulay=true end
    if(macaulay) symbolic=true end
    
    #-----------setting up desk-----------#
    df = DataFrame(Split=String[], CF=String[]) 
    dict=dictionary(net,inheritancecorrelation)

    #--------(almost) original stuff--------#
    net.node[net.root].leaf && error("The root can't be a leaf.")
    PN.check_nonmissing_nonnegative_edgelengths(net,
        "Edge lengths are needed in coalescent units to calcualte expected CFs.")
    all(e.gamma >= 0.0 for e in net.edge) || error("some γ's are missing for hybrid edges: can't calculate expected CFs.")
    inheritancecorrelation >= 0 || error("the inheritance correlation should be non-negative")
    inheritancecorrelation <= 1 || error("the inheritance correlation should be <= 1")
    taxa = sort!(tipLabels(net))
    taxonnumber = Dict(taxa[i] => i for i in eachindex(taxa))
    ntax = length(taxa)
    nCk = PN.nchoose1234(ntax) # matrix to rank 4-taxon sets
    qtype = MVector{3,Float64} # 3 floats: CF12_34, CF13_24, CF14_23; initialized at 0.0
    numq = nCk[ntax+1,4]
    quartet = Vector{PN.QuartetT{qtype}}(undef, numq)
    ts = [1,2,3,4]
    for qi in 1:numq
        quartet[qi] = PN.QuartetT(qi, SVector{4}(ts), MVector(0.,0.,0.))
        # next: find the 4-taxon set with the next rank,
        #       faster than using the direct mapping function
        ind = findfirst(x -> x>1, diff(ts))
        if ind === nothing ind = 4; end
        ts[ind] += 1
        for j in 1:(ind-1)
            ts[j] = j
        end
    end
    if showprogressbar
        nstars = (numq < 50 ? numq : 50)
        nquarnets_perstar = (numq/nstars)
        println("Calculation quartet CFs for $numq quartets...")
        print("0+" * "-"^nstars * "+100%\n  ")
        stars = 0
        nextstar = Integer(ceil(nquarnets_perstar))
    end
    for qi in 1:numq
        network_expectedCF!(quartet[qi], net, taxa, taxonnumber, inheritancecorrelation, df, symbolic, dict)
        if showprogressbar && qi >= nextstar
            print("*")
            stars += 1
            nextstar = Integer(ceil((stars+1) * nquarnets_perstar))
        end
    end
    showprogressbar && print("\n")


    #--------output--------#
    #extended newick string with parameters replaced with t_{1} etc. stuff or not; 
    #always print the topology at each analysis
    if(symbolic)
        eNewick=gettingSymbolicTopology(net,dict)
    else
        eNewick=PN.writeTopology(net)
    end
    #show and/or store topology as a text file if requested (printCFs=savenet=true)
    if(printCFs) display(eNewick) end
    if(savenet) open("$filename.net.txt", "w") do file write(file, eNewick) end end

    #displaying and saving csv spreadsheet if requested (savecsv=true)
    if(printCFs) display(df) end
    if(savecsv) CSV.write("$filename.csv", df, header=false) end 

    #macaulay output
    numCFs=size(df)[1]
    if(macaulay||matlab) 
        dataframe=deepcopy(df)
        params=gettingSymbolicInput(net, dataframe, inheritancecorrelation) 
    end
    if(macaulay)
        open("$filename.m2.txt", "w") do file
        str="R = QQ["
        for par in params str=str*par*"," end
        str=str*"C_1..C_$numCFs]\n"
        str=str*"I = ideal(\n"
        i=1
        while i<numCFs
            str=str*"$(dataframe[i,2])-C_$i,\n"
            i+=1
        end
        str=str*"$(dataframe[numCFs,2])-C_$numCFs);\n"
        str=str*"G = eliminate (I, {"
        numparams=length(params)   
        #for par in params str=str*par*"," end
        for par in 1:(length(params)-1) str=str*params[par]*"," end
        str=str*"$(params[numparams])})\n"
        str=str*"S = QQ[C_1..C_$numCFs]\nJ = sub(G,S)\ndim J"
        write(file, str)
        end
    end

    #matlab output
    if(matlab)
        open("$filename.matlab.txt", "w") do file 
            str="% Declare variables\n"
            str=str*"syms "
            for par in params str=str*par*" " end
            for i in 1:numCFs str=str*"C_$i " end
            str=str*"\n\n% matrix of generating polynomials\n"
            str=str*"F=["
            i=1
            while i<numCFs
                str=str*"$(dataframe[i,2])-C_$i,\n"
                i+=1
            end
            str=str*"$(dataframe[numCFs,2])-C_$numCFs];\n"
            str=str*"\n% matrix of generating polynomials\n"
            str=str*"\n% Array of all variables\n"
            str=str*"V=["
            for par in params str=str*par*" " end
            for i in 1:numCFs-1 str=str*"C_$i " end
            str=str*"C_$numCFs]\n"
            str=str*"\n% Compute dimension\nCoalDim(F,V)"
        write(file, str)
        end
    end
   
    #if !(symbolic) 
    return quartet, taxa, df #end
end


"""
    network_expectedCF!(quartet::QuartetT, net::HybridNetwork, taxa, taxonnumber,
            inheritancecorrelation)

Update `quartet.data` to contain the quartet concordance factors expected from
the multispecies coalescent along network `net` for the 4-taxon set `taxa[quartet.taxonnumber]`.
`taxa` should contain the tip labels in `net`. `quartet.taxonnumber` gives the
indices in `taxa` of the 4 taxa of interest. `taxonnumber` should be a dictionary
mapping taxon labels in to their indices in `taxa`, for easier lookup.

`net` is not modified.

For `inheritancecorrelation` see [`network_expectedCF`](@ref).
Its value should be between 0 and 1 (not checked by this internal function).
"""
function network_expectedCF!(quartet::PN.QuartetT{MVector{3,Float64}},
                                net::HybridNetwork, 
                                taxa, 
                                taxonnumber, 
                                inheritancecorrelation, 
                                df,
                                symbolic,
                                dict)
    #kong: create an array that stores the CF formulas for ab|cd, ac|bd, ad|bc
    qCFp=String["","",""] #*-_-*#

    net = deepcopy(net)
    PN.removedegree2nodes!(net)
    # delete all taxa except for the 4 in the quartet
    for taxon in taxa
        taxonnumber[taxon] in quartet.taxonnumber && continue
        deleteleaf!(net, taxon, simplify=false, unroot=false)
        # would like unroot=true but deleteleaf! throws an error when the root is connected to 2 outgoing hybrid edges
    end
    #println("[144]Input network = $(PN.writeTopologyLevel1(net))")
    
    q,qCFp=network_expectedCF_4taxa!(net, taxa[quartet.taxonnumber], inheritancecorrelation, qCFp, dict, symbolic)
    quartet.data .= q

    #kong: storing the equations to DataFrames
    for i in 1:3
        qCFp[i]=replace(qCFp[i], "&"=>"")
        qCFp[i]=filter(x -> !isspace(x), qCFp[i])
    end

    f=taxa[quartet.taxonnumber]
    push!(df, ("$(f[1])$(f[2])|$(f[3])$(f[4])", "$(qCFp[1])"))
    push!(df, ("$(f[1])$(f[3])|$(f[2])$(f[4])", "$(qCFp[2])"))         
    push!(df, ("$(f[1])$(f[4])|$(f[2])$(f[3])", "$(qCFp[3])"))
    
    return quartet
end

"""
    network_expectedCF_4taxa!(net::HybridNetwork, fourtaxa, inheritancecorrelation)

Return the quartet concordance factors expected from the multispecies coalescent
along network `net`, where the 3 quartet topologies are ordered following the
ordering of taxon names in `fourtaxa`, that is: if `fourtaxa` is a,b,c,d,
then the concordance factors are listed in this order:

    (qCF(ab|cd), qCF(ac|bd), qCF(ad,bc))

Assumptions about `net`:
- has 4 taxa, and those are the same as `fourtaxa`
- no degree-2 nodes, except perhaps for the root
- edge lengths are non-missing
- hybrid edge γ's are non-missing

The network is modified as follows: what's above the LSA is removed,
the 2 edges incident to the root are fused (if the root is of degree 2),
and external degree-2 blobs are removed. `net` is then simplified recursively
by removing hybrid edges for the recursive calculation of qCFs.

For `inheritancecorrelation` see [`network_expectedCF`](@ref).
Its value should be between 0 and 1 (not checked by this internal function).
"""
function network_expectedCF_4taxa!(net::HybridNetwork, fourtaxa, inheritancecorrelation, qCFp, dict, symbolic)
    
    #kong: begin writing qCF equations  with an opening bracket
    qCFp .*= "(" 
    
    deleteaboveLSA!(net)
    # make sure the root is of degree 3+
    if length(net.node[net.root].edge) <= 2
        PN.fuseedgesat!(net.root, net)
    end
    # find and delete degree-2 blobs along external edges
    bcc = biconnectedComponents(net, true) # true: ignore trivial blobs
    entry = PN.biconnectedcomponent_entrynodes(net, bcc)
    entryindex = indexin(entry, net.nodes_changed)
    exitnodes = PN.biconnectedcomponent_exitnodes(net, bcc, false) # don't redo the preordering
    bloborder = sortperm(entryindex) # pre-ordering for blobs in their own blob tree
    function isexternal(ib) # is bcc[ib] of degree 2 and adjacent to an external edge?
        # yes if: 1 single exit adjacent to a leaf
        length(exitnodes[ib]) != 1 && return false
        ch = PN.getChildren(exitnodes[ib][1])
        return length(ch) == 1 && ch[1].leaf
    end
    for ib in reverse(bloborder)
        isexternal(ib) || continue # keep bcc[ib] if not external of degree 2
        for he in bcc[ib]
            he.isMajor && continue
            # deletion of a hybrid can hide the deletion of another: check that he is still in net
            any(e -> e===he, net.edge) || continue
            # delete minor hybrid edge with options unroot=true: to make sure the
            # root remains of degree 3+, in case a degree-2 blob starts at the root
            # simplify=true: bc external blob
            PN.deletehybridedge!(net,he, false,true,false,true,false)
        end
    end
    ndes = 4 # number of taxa descendant from lowest hybrid node
    if net.numHybrids > 0
        preorder!(net)
        # find a lowest hybrid node and # of taxa below it
        hyb = net.nodes_changed[findlast(n -> n.hybrid, net.nodes_changed)]
        #funneledge = [e for e in hyb.edge if getparent(e) === hyb]
        funneledge = [e for e in hyb.edge if PhyloNetworks.getParent(e) === hyb]
        ispolytomy = length(funneledge) > 1
        funneldescendants = union([PN.descendants(e) for e in funneledge]...)
        ndes = length(funneldescendants)
        #n2 = (ispolytomy ? hyb : getchild(funneledge[1]))
        n2 = (ispolytomy ? hyb : PhyloNetworks.getChild(funneledge[1]))
        ndes > 2 && n2.leaf && error("2+ descendants below the lowest hybrid, yet n2 is a leaf. taxa: $(fourtaxa)")
    end
    if ndes > 2 # simple formula for qCF: find cut edge and its length
        # inheritance correlation has no impact
        # pool of cut edges below. contains NO external edge, bc n2 not leaf (if reticulation), nice tree ow
        cutpool = (net.numHybrids == 0 ? net.edge :
                    [e for e in n2.edge if PN.getParent(e) === n2])
        #filter!(e -> !getchild(e).leaf, cutpool)
        filter!(e -> !PhyloNetworks.getChild(e).leaf, cutpool)
        net.numHybrids > 0 || length(cutpool) <= 1 ||
            error("2+ cut edges, yet 4-taxon tree, degree-3 root and no degree-2 nodes. taxa: $(fourtaxa)")
        sistertofirst = 2    # arbitrarily correct if 3-way polytomy (no cut edge)
        internallength = 0.0 # correct if polytomy
        for e in cutpool
            length(cutpool) < 3 || println("more than 2 edged merged")
            internallength += e.length
            hwc = hardwiredCluster(e, fourtaxa)
            sistertofirst = findnext(x -> x == hwc[1], hwc, 2)
        end

        #internallength=round(internallength, digits = dpoints)
        minorcf = exp(-internallength)/3
        majorcf = 1.0 - 2 * minorcf
        qCF = (sistertofirst == 2 ? MVector{3,Float64}(majorcf,minorcf,minorcf) :
              (sistertofirst == 3 ? MVector{3,Float64}(minorcf,majorcf,minorcf) :
                                    MVector{3,Float64}(minorcf,minorcf,majorcf) ))
        
        #kong: writing out the equations
        #println(internallength)
        #println(dict[internallength])
        internallength=round(internallength, digits = dpoints)
        if symbolic
            minorcfp = "exp(-$(dict[internallength]))/3"
            majorcfp = "1-2*$minorcfp"
        else
            minorcfp = "exp(-$internallength)/3"
            majorcfp = "1-2*$minorcfp"
        end
        #=to be added to symbolic=#
        (sistertofirst == 2 ? (qCFp[1]*="$majorcfp",qCFp[2]*="$minorcfp",qCFp[3]*="$minorcfp") :
        (sistertofirst == 3 ? (qCFp[1]*="$minorcfp",qCFp[2]*="$majorcfp",qCFp[3]*="$minorcfp") :
                              (qCFp[1]*="$minorcfp",qCFp[2]*="$minorcfp",qCFp[3]*="$majorcfp") ))                      
        qCFp .*= ")" #kong: end qCF with an closing bracket
        
        return qCF, qCFp
    end

    ndes > 0 || error("weird: hybrid node has no descendant taxa")
    # by now, there are 1 or 2 taxa below the lowest hybrid
    qCF = MVector{3,Float64}(0.0,0.0,0.0) # mutated later
    #parenthedge = [e for e in hyb.edge if getchild(e) === hyb]
    parenthedge = [e for e in hyb.edge if PhyloNetworks.getChild(e) === hyb]
    all(h.hybrid for h in parenthedge) || error("hybrid $(hyb.number) has a parent edge that's a tree edge")
    parenthnumber = [p.number for p in parenthedge]
    nhe = length(parenthedge)
    if ndes == 1 # weighted qCFs average of the nhe (often = 2) displayed networks
        # inheritance correlation has no impact
        for i in 1:nhe # keep parenthedge[i], remove all others
            # kong: add + if there are more than a single element
            if i>1 qCFp .*= "+" end
            
            gamma = parenthedge[i].gamma
            simplernet = ( i < nhe ? deepcopy(net) : net ) # last case: to save memory allocation
            for j in 1:nhe
                j == i && continue # don't delete hybrid edge i!
                pe_index = findfirst(e -> e.number == parenthnumber[j], simplernet.edge)
                PN.deletehybridedge!(simplernet, simplernet.edge[pe_index],
                    false,true,false,false,false) # ., unroot=true, ., simplify=false,.
            end
            qCF0,qCFp = network_expectedCF_4taxa!(simplernet, fourtaxa, inheritancecorrelation,qCFp,dict,symbolic)
            qCF .+= gamma .* qCF0
            
            #kong: writing gamma into qCF equations
            if symbolic qCFp .*= "*$(dict[gamma])" 
            else qCFp .*= "*$gamma"end
        end
        
        #kong: closing qCFq with a bracket
        qCFp .*= ")"

        return qCF, qCFp
    end

    # by now: 2 descendant below the lowest hybrid node: hardest case
    # weighted qCFs average of 3 networks: 2 displayed, 1 "parental" (unless same parents)
    sameparents = (inheritancecorrelation == 1)
    #inheritancecorrelation=round(inheritancecorrelation, digits = dpoints)
    #dict[inheritancecorrelation]="&rho"
    #dict[oneminusrho] = "1-&rho"

    oneminusrho = 1 - inheritancecorrelation
    #oneminusrho=round(oneminusrho, digits = dpoints)
    if symbolic 
        oneminusrhop="1-$(dict[inheritancecorrelation])"
    else
        oneminusrhop="1-$inheritancecorrelation"
    end

    hwc = hardwiredCluster(parenthedge[1], fourtaxa)
    sistertofirst = findnext(x -> x == hwc[1], hwc, 2)
    internallength = ( ispolytomy ? 0.0 : funneledge[1].length)
    deepcoalprob = exp(-internallength)
    deepcoalprob=round(deepcoalprob, digits = dpoints)
    internallength=round(internallength, digits = dpoints)
    #dict[deepcoalprob]="e^{-$(dict[internallength])}"
    if symbolic deepcoalprobp="exp(-$internallength)" end
    # initialize qCF: when the 2 descendants coalesce before reaching the hybrid node
    qCF = (sistertofirst == 2 ? MVector{3,Float64}(1.0-deepcoalprob,0.0,0.0) :
          (sistertofirst == 3 ? MVector{3,Float64}(0.0,1.0-deepcoalprob,0.0) :
                                MVector{3,Float64}(0.0,0.0,1.0-deepcoalprob) ))

    #kong: deepcoalprobability
    if symbolic deepcoalprobp = "exp(-$(dict[internallength]))" 
    else deepcoalprobp = "exp(-$internallength)" end
    
    (sistertofirst == 2 ? (qCFp[1]*="(1-$deepcoalprobp)+",qCFp[2]*="",qCFp[3]*="") :
    (sistertofirst == 3 ? (qCFp[1]*="",qCFp[2]*="(1-$deepcoalprobp)+",qCFp[3]*="") :
                        (qCFp[1]*="",qCFp[2]*="",qCFp[3]*="(1-$deepcoalprobp)+") ))                                     
    qCFp .*= ""

    # no coalescence on cut-edge: delete it and extract parental networks
    ispolytomy || PN.shrinkedge!(net, funneledge[1])
    # shrinkedge! requires PhyloNetworks v0.15.2
    childedge = [e for e in hyb.edge if PN.getParent(e) === hyb]
    length(childedge) == 2 ||
      error("2-taxon subtree, but not 2 child edges after shrinking the cut edge.")
    all(PN.getChild(e).leaf for e in childedge) ||
      error("2-taxon subtree, cut-edge shrunk, but the 2 edges aren't both external")
    childnumber = [e.number for e in childedge]
    for i in 1:nhe
      pgam=parenthedge[i].gamma
      weighti = deepcoalprob * pgam
      #kong: round weighti to n digits
      weighti=round(weighti, digits = dpoints)
      #dict[weighti]="($(dict[deepcoalprob])*($(dict[pgam])))"
      
      for j in (sameparents ? i : 1):i # if inheritancecorrelation=1 then i!=j has probability 0
        gammaj = parenthedge[j].gamma
        #kong: round gammaj to n digits
        gammaj=round(gammaj, digits = dpoints)

        simplernet = ( i < nhe || j < nhe ? deepcopy(net) : net )
        # delete all hybedges other than i & j
        for k in 1:nhe
            (k == i || k ==j) && continue # don't delete hybrid edges i or j
            pe_index = findfirst(e -> e.number == parenthnumber[k], simplernet.edge)
            PN.deletehybridedge!(simplernet, simplernet.edge[pe_index],false,true,false,false,false) # ., unroot=true,., simplify=false,.
        end
        if i != j
            # detach childedge[2] from hyb and attach it to hyb's parent j
            pej_index = findfirst(e -> e.number == parenthnumber[j], simplernet.edge)
            pej = simplernet.edge[pej_index]
            pn = PN.getParent(pej)
            hn = PN.getChild(pej) # hyb node, but in simplernet
            ce2_index = findfirst(e -> e.number == childnumber[2], simplernet.edge)
            ce2 = simplernet.edge[ce2_index]
            PN.removeEdge!(hn,ce2)
            hn_index = findfirst(x -> x === hn, ce2.node)
            ce2.node[hn_index] = pn # ce2.isChild1 remains synchronized
            push!(pn.edge, ce2)
            # then delete hybedge j
            PN.deletehybridedge!(simplernet, pej, false,true,false,false,false) # ., unroot=true,., simplify=false,.)
            for e in simplernet.edge
                e.length=round(e.length, digits = dpoints)
            end
        end
        #kong: add "+" if there are more than a single element in the equation
        if i>1 qCFp .*= "+" end
        #initialize qCFp for qCF_subnet
        qCFps=["","",""]
        qCF_subnet, qCFps = network_expectedCF_4taxa!(simplernet, fourtaxa, inheritancecorrelation, qCFps, dict,symbolic)
        if i == j
            prob = weighti * (gammaj * oneminusrho + inheritancecorrelation)
            prob=round(prob, digits = dpoints)

            qCF .+= prob .* qCF_subnet
           
        if symbolic
            for i in 1:3
                qCFp[i] *= "((($deepcoalprobp * $(dict[pgam])) * ($(dict[gammaj]) * $(dict[oneminusrho]) + $(dict[inheritancecorrelation]))) * ($(qCFps[i])))"
            end
        else
            for i in 1:3
                qCFp[i] *= "((($deepcoalprobp * $pgam) * ($gammaj * $oneminusrhop + $inheritancecorrelation)) * ($(qCFps[i])))"
            end
        end
        else # add subnetwork with flipped assignment of the 2 taxa to parents i & j
            flipped_ij = (sistertofirst == 2 ? [1,3,2] :
                         (sistertofirst == 3 ? [3,2,1] : [2,1,3] ))
            prob = weighti * gammaj * oneminusrho
            prob=round(prob, digits = dpoints)

            qCF .+= prob .* (qCF_subnet .+ qCF_subnet[flipped_ij])
    
        if symbolic
            for i in 1:3
                qCFp[i] *= "((($deepcoalprobp * $(dict[pgam])) * $(dict[gammaj]) * $(dict[oneminusrho])) * ($(qCFps[i])+$(qCFps[flipped_ij[i]])))"
            end
        else
            for i in 1:3
                qCFp[i] *= "((($deepcoalprobp * $pgam) * $gammaj * $oneminusrhop) * ($(qCFps[i])+$(qCFps[flipped_ij[i]])))"
            end
        end
        end
      end
    end
    
    qCFp .*= ")"
    
    return qCF, qCFp
end


"""
    readTopologyrand(net; defaultval=1.1::Float64)

    Generate random parameter values for edge lengths with >`defaultval` 
    and inheritance probabilities=[0,1] and assign them to the provided network.
    Input network can be either Newick/extended Newick-formatted String or
    `HybridNetwork` object read in by function `readTopology` in `PhyloNetworks` pkg.
"""
function readTopologyrand(net; defaultval=1.1::Float64)
    #--------some warnings--------#
    defaultval>=1.1 || @warn "Setting deafultval ≥ 1.1 is recommended for various reasons."
    
    #---------read in topology if input is newick string or HybridNetwork object---------#
    if typeof(net)==HybridNetwork net=net
    else net=PN.readTopology(net) end

    #--------generaete arbitrary edge lengths--------#
    for e in net.edge
        e.length=round((rand()+defaultval), digits = dpoints)
    end

    #--------generaete arbitrary inheritance probabilities--------#
    #----preambles----#
    reticulatenodeindex=Int[]
    nreticulate=length(net.hybrid)
    gammavec=zeros(nreticulate)
    
    #getting hybrid node index numbers
    for n in net.node
        n.hybrid && push!(reticulatenodeindex,n.number)
    end
    #check the number of hybrid nodes are counted correctly
    length(reticulatenodeindex)==nreticulate || @error "Generation of gamma was incomplete. Rerun the function (error 1)."

    #generate arbitrary gamma values n=number of reticulation nodes (that will be assigned to one of the incoming edges)
    for gamma in 1:nreticulate 
        val=rand()
        while !(0<val<1) val=rand() end #check if generated gamma = [0,1]
        gammavec[gamma]=round(val, digits = dpoints)
    end
    #check if all gamma is properly generated (ie., no zero in gammavec)
    all(!=(0.0), gammavec) && all(!=(1.0), gammavec) || @error "Generation of gamma was incomplete. Rerun the function (error 2)."
    
    #assign inheritance probabilities to reticulate edges
    for i in 1:nreticulate
        nthvisit=1
        for e in net.edge
            if e.hybrid
                child=PhyloNetworks.getChild(e)
                if child.number==reticulatenodeindex[i]
                    if nthvisit==1
                        e.gamma=round(gammavec[i], digits = dpoints)
                        nthvisit+=1
                    elseif nthvisit==2 e.gamma=round((1-gammavec[i]), digits = dpoints)
                    elseif nthvisit>2 error("Three edges enters hybrid node number $(reticulatenodeindex[i]). Use binary network as input.")
                    end
                end 
            else continue
            end
        end
    end

    return net
end

function gettingSymbolicTopology(net::HybridNetwork,dict)
    eNewick=PN.writeTopology(net)
    for e in net.edge 
        eNewick=replace(eNewick,"$(e.length)"=>"t_{$(e.number)}")
        eNewick=replace(eNewick,"$(e.gamma)"=>"$(dict[e.gamma])")
    end
    return eNewick
end


function gettingSymbolicInput(net::HybridNetwork, df, inheritancecorrelation)
    edgenumber=length(net.edge)
    retnumber=length(net.hybrid)
    numCFs=size(df)[1]
    params=String[]

    for i in 1:numCFs
        df[i,2]=replace(df[i,2],"rho"=>"$inheritancecorrelation")
        string=df[i,2]
        record=false
        expressions=String[]
        
        
        for e in 1:length(net.edge)
            if(occursin("-t_{$e}","$(df[i,2])"))
                #println("-t_{$e},$(df[i,2])"#
                push!(expressions,"X$e")
            end
        end

        for e in 1:length(net.hybrid)
            if(occursin("r_{$e}","$(df[i,2])"))
                #println("r_{$e},$(df[i,2])")
                push!(expressions,"R$e")
            end
        end
        
        append!(params,expressions)

        
    end

    params=unique(params)
    ring=String[]

    for i in params
        x=split(i, "}-")
        for j in x
            push!(ring,j)
        end
    end

    params=sort!(params)
    
    for cf in 1:numCFs
        df[cf,2]=replace(df[cf,2],"r_{" => "R" )#Ropen
        df[cf,2]=replace(df[cf,2],"exp(-t_{" => "(X" )
        df[cf,2]=replace(df[cf,2],"-t_{" => "*X" )
        df[cf,2]=replace(df[cf,2],"})" => ")" )
        df[cf,2]=replace(df[cf,2],"}" => "" )#close
    end
    
    return params
end

function dictionary1(net,inheritancecorrelation)
    nmerge=2 #number of edges that are merged during CF calculation
    dict=Dict()
    edgelengths=zeros(length(net.edge))
    for i in 1:length(net.edge)
        edgelengths[i]=round(net.edge[i].length, digits=dpoints)
    end
    numedges=length(edgelengths)
    println(edgelengths)
    for eindex in 1:numedges 
        dict[edgelengths[eindex]] = "t_{$eindex}"
    end
    println(dict)


end

function dictionary(net,inheritancecorrelation; convert=false)
##########BEGINNING OF DISCTIONARY CONSTRUCTION##########E
    #kong: create a dictionary of parameter labels to values (tau and gamma)
    nmerge=2 #number of edges that are merged during CF calculation
    dict=Dict()
    edgelengths=zeros(length(net.edge))
    for i in 1:length(net.edge)
        edgelengths[i]=round(net.edge[i].length, digits=dpoints)
    end
    #println(edgelengths)

    #branch length => t_{branch id}
    for e in net.edge
        enum=e.number
        e.length=round(e.length,digits=dpoints)
        dict[e.length] = "t_{$enum}" #*-_-*#       
        if(convert) dict["t_{$enum}"] = e.length end
    end
    
    for i in 1:length(net.edge)
        for j in 2:length(net.edge)
            e1=net.edge[i]
            e2=net.edge[j]
            if e1==e2 break 
            else 
            length=round(sum([e1.length,e2.length]), digits = dpoints)
            dict[length] = "t_{$i}-t_{$j}" #*-_-*#            
            end
        end
    end
    for i in 1:length(net.edge)
        for j in 2:length(net.edge)
            for k in 3:length(net.edge)
                e1=net.edge[i]
                e1.length=round(e1.length, digits = dpoints)
                e2=net.edge[j]
                e2.length=round(e2.length, digits = dpoints)
                e3=net.edge[k]
                e3.length=round(e3.length, digits = dpoints)
                length=round(sum([e1.length,e2.length,e3.length]), digits = dpoints)
                dict[length] = "t_{$i}-t_{$j}-t_{$k}" #*-_-*#
            end
        end
    end
    for i in 1:length(net.edge)
        for j in 2:length(net.edge)
            for k in 3:length(net.edge)
                for l in 4:length(net.edge)
                    e1=net.edge[i]
                    e2=net.edge[j]
                    e3=net.edge[k]
                    e4=net.edge[l]
                    length=round(sum([e1.length,e2.length,e3.length,e4.length]), digits = dpoints)
                    dict[length] = "t_{$i}-t_{$j}-t_{$k}-t_{$l}" #*-_-*#
                end
            end
        end
    end
    for i in 1:length(net.edge)
        for j in 2:length(net.edge)
            for k in 3:length(net.edge)
                for l in 4:length(net.edge)
                    for a in 5:length(net.edge)
                        e1=net.edge[i]
                        e2=net.edge[j]
                        e3=net.edge[k]
                        e4=net.edge[l]
                        e5=net.edge[a]
                        length=round(sum([e1.length,e2.length,e3.length,e4.length,e5.length]), digits = dpoints)
                        dict[length] = "t_{$i}-t_{$j}-t_{$k}-t_{$l}-t_{$a}" #*-_-*#
                    end
                end
            end
        end
    end   
    for i in 1:length(net.edge)
        for j in 2:length(net.edge)
            for k in 3:length(net.edge)
                for l in 4:length(net.edge)
                    for a in 5:length(net.edge)
                        for b in 6:length(net.edge)
                            e1=net.edge[i]
                            e2=net.edge[j]
                            e3=net.edge[k]
                            e4=net.edge[l]
                            e5=net.edge[a]
                            e6=net.edge[b]
                            length=round(sum([e1.length,e2.length,e3.length,e4.length,e5.length,e6.length]), digits = dpoints)
                            dict[length] = "t_{$i}-t_{$j}-t_{$k}-t_{$l}-t_{$a}-t_{$b}" #*-_-*#
                        end
                    end
                end
            end
        end
    end   
    
    #dictionary for gamma
    hybnodenum=[]
    for n in net.node
        if n.hybrid push!(hybnodenum,n.number) end
    end
    lhybridnodenum=length(hybnodenum)
    
    for j in 1:lhybridnodenum
        hybnode=hybnodenum[j]
        i=1
        for e in net.edge
            child=PhyloNetworks.getChild(e)
            if child.number==hybnode
                e.gamma=round(e.gamma, digits = dpoints)
                if i==1 dict[e.gamma] = "r_{$j}"
                    if(convert) dict["r_{$j}"] = e.gamma end
                    i+=1
                elseif i==2
                    dict[e.gamma] = "(1-r_{$j})"
                else
                    println("there are tree incoming edge at hybrid node number $(child.num)")
                end
            end 
        end
    end

    inheritancecorrelation=round(inheritancecorrelation, digits = dpoints)
    oneminusrho = 1 - inheritancecorrelation
    oneminusrho=round(oneminusrho, digits = dpoints)

    dict[inheritancecorrelation]="&rho"
    dict[oneminusrho] = "1-&rho"
  
    #kong: check if there are redundant keys in dictionary. If so, regenerate parameter sets using readTopologyrand()
    x=collect(keys(dict)) 
    length(unique(x))==length(dict) || error("multiple keys present in the dictionary")

    ##########END OF DISCTIONARY CONSTRUCTION##########E

    return dict
end




##For Macaulay2
function dictionarymaca(net)
    dict=Dict()
    edgenumber=length(net.edge)
    retnumber=length(net.hybrid)
    for i in 1:edgenumber
        dict["exp(-t_{$i})"] = "X$i" #*-_-*#
    end
    for i in 1:edgenumber
        for j in 1:edgenumber
            dict["exp(-t_{$i}-t_{$j})"] = "(X$i*X$j)" #*-_-*#
        end
    end
    for i in 1:edgenumber
        for j in 1:edgenumber
            for k in 1:edgenumber
                dict["exp(-t_{$i}-t_{$j}-t_{$k})"] = "(X$i*X$j*X$k)" #*-_-*#
            end
        end
    end
    for i in 1:edgenumber
        for j in 1:edgenumber
            for k in 1:edgenumber
                for l in 1:edgenumber
                    dict["exp(-t_{$i}-t_{$j}-t_{$k}-t_{$l})"] = "(X$i*X$j*X$k*X$l)" #*-_-*#
                end
            end
        end
    end
    for i in 1:edgenumber
        for j in 1:edgenumber
            for k in 1:edgenumber
                for l in 1:edgenumber
                    for x in 1:edgenumber
                        dict["exp(-t_{$i}-t_{$j}-t_{$k}-t_{$l}-t_{$x})"] = "(X$i*X$j*X$k*X$l*X$x)" #*-_-*#
                    end
                end
            end
        end
    end
    for i in 1:edgenumber
        for j in 1:edgenumber
            for k in 1:edgenumber
                for l in 1:edgenumber
                    for x in 1:edgenumber
                        for y in 1:edgenumber
                            dict["exp(-t_{$i}-t_{$j}-t_{$k}-t_{$l}-t_{$x}-t_{$y})"] = "(X$i*X$j*X$k*X$l*X$x*X$y)" #*-_-*#
                        end
                    end
                end
            end
        end
    end
    
    for i in 1:retnumber
        dict["r_{$i}"] = "R$i" #*-_-*#
    end
    
    return dict
end




function generate_combinations(net, depth)
    # Initialize the result with single elements
    nedge=length(net.edge)
    input_array=(1:nedge)
    result = ["exp(-t{$(i)})" for i in input_array]
    
    # Generate combinations for the given depth
    function recurse(current, level)
        if level == depth
            return
        end
        
        for i in input_array
            new_combination = "" * current * "}-t_{" * string(i) * ""
            push!(result, "exp(-t{" * new_combination * "})")
            recurse(new_combination, level + 1)
        end
    end
    
    # Start recursion from each element
    for i in input_array
        recurse(string(i), 1)
    end

    #display(result)
    
    return result
end

#input_array = [1, 2, 3]
#depth = 3  # Adjust the depth as needed
#output = generate_combinations(input_array, depth)

#println(output)










function test(net;inheritancecorrelation=0,threshold=0.01::Float64,savecsv=false::Bool)
    
    df0=DataFrame(Split=String[], CF0=Float64[])
    df1=DataFrame(Split=String[], CF1=String[], CF11=Float64[])
    df2=DataFrame(Split=String[], CF2=String[], CF21=String[], CF22=Float64[])

    #1. df=split, val1, val2, val3
    print("Evaluating df0...")
    
    quartet,taxa,df=network_expectedCF(net,printCFs=false)
    numquartet=length(quartet)
    for i in 1:numquartet
        l1=taxa[quartet[i].taxonnumber[1]]
        l2=taxa[quartet[i].taxonnumber[2]]
        l3=taxa[quartet[i].taxonnumber[3]]
        l4=taxa[quartet[i].taxonnumber[4]]
        splitt1="$l1$l2|$l3$l4"
        splitt2="$l1$l3|$l2$l4"
        splitt3="$l1$l4|$l2$l3"
        push!(df0, (splitt1,quartet[i].data[1]))
        push!(df0, (splitt2,quartet[i].data[2]))
        push!(df0, (splitt3,quartet[i].data[3]))
    end

    println("Complete")
    
    
    #2. see if symbolic=false returns the same value as 1    
    print("Evaluating df1...")

    quartet,taxa,df=network_expectedCF(net,symbolic=false, printCFs=false)
    for i in 1:(size(df)[1])
        val=eval(Meta.parse(df[i,2]))
        val=round(val, digits = dpoints)
        push!(df1,(df[i,1],df[i,2],val))
    end

    println("Complete")


    #3. get symbolic one
    print("Evaluating df2, this may take a while...")

    quartet,taxa,df=network_expectedCF(net,symbolic=true, printCFs=false)
    for i in 1:(size(df)[1])
        push!(df2,(df[i,1],df[i,2],df[i,2],0.0))
    end

    dict=Dict()
    for i in net.edge
        dict["t_{$(i.number)}"]=String("$(i.length)")
    end

    dict1=dictionary(net,inheritancecorrelation, convert=true)
    for i in 1:length(net.edge)
        for j in 1:(size(df)[1])
            df2[j,3]=(replace(df2[j,3],"t_{$i}"=>dict["t_{$i}"]))
        end
    end
    
    for i in 1:length(net.hybrid)
        for j in 1:(size(df)[1])
            df2[j,3]=(replace(df2[j,3],"r_{$i}"=>dict1["r_{$i}"]))
        end
    end

    for j in 1:(size(df)[1])
        df2[j,3]=(replace(df2[j,3],"rho"=>inheritancecorrelation))
    end

    for i in 1:(size(df2)[1])
        val=eval(Meta.parse(df2[i,3]))
        val=round(val, digits = dpoints)
        df2[i,4]=val
    end

    println("Complete")


    finaldf = outerjoin(df0, df1, df2, on = :Split)

    for i in 1:(size(finaldf)[1])
        abs(finaldf[i,2]-finaldf[i,4]) < threshold || error("values do not match")
        abs(finaldf[i,2]-finaldf[i,7]) < threshold || error("values do not match")
        abs(finaldf[i,4]-finaldf[i,7]) < threshold || error("values do not match")        
    end
    
    if(savecsv) CSV.write("test.csv", finaldf, header=true) end

    return finaldf
end

#network non binary - having rho transformed for exp(-0)