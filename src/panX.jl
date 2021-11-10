module PanX

using GZip, JSON
using Rematch
using ProgressMeter

import ..Graphs

# ------------------------------------------------------------------------
# Phylogenetic tree manipulation

module Phylo

using Rematch

ENV["PYTHON"] = ""; import PyCall # use internal python so we can build a conda environment
# XXX: might need this to make it automatically work...
# import Pkg
# Pkg.build("PyCall")

# ignore syntax warnings
PyCall.pyimport("warnings").filterwarnings("ignore")
Ete3 = PyCall.pyimport_conda("ete3", "ete3", "etetoolkit")

function tree(newick)
    return Ete3.Tree(newick)
end

function binary!(tree)
    tree.resolve_polytomy(recursive=true)
    return tree
end

function root!(tree)
    midpoint = tree.get_midpoint_outgroup()
    tree.set_outgroup(midpoint)
    tree.ladderize()
    return tree
end

function rescale!(tree, by::T) where T <: Real
    for node in tree.traverse("preorder")
        node.dist *= by
    end
end

write(io::IO, tree::PyCall.PyObject) = Base.write(io, tree.write(format=1))

function dictionary(tree::PyCall.PyObject, level::Int, i::Int; mutations=false)
    node = Dict(
        "name" => if tree.name == ""
            i += 1
            suffix = string(i-1)
            suffix = @match length(suffix) begin
                1 => "0000$(suffix)"
                2 => "000$(suffix)"
                3 => "00$(suffix)"
                4 => "0$(suffix)"
                _ => suffix
            end
            "NODE_$(suffix)"
        else
            tree.name
        end,
        "branch_length" => tree.dist,
    )

    if mutations
        # TODO: actually write mutations?
        node["muts"]    = ""
        node["aa_muts"] = ""
    else
        node["clade"] = level
    end

    if tree.name != ""
        if mutations
            node["accession"]  = split(tree.name,'#')[1]
            node["annotation"] = "pan-contig"
        else
            node["attr"] = Dict(
                "host"   => tree.name,
                "strain" => tree.name,
            )
        end
        return node, i
    end

    node["children"] = Dict[]
    for child in tree.children
        c, i = dictionary(child, level+1, i; mutations=mutations)
        push!(node["children"], c)
    end

    return node, i
end

end # Phylo

function gzip(from::AbstractString, to::AbstractString; clean=false)
    isfile(from) || error("$(from) not found")

    Base.open(from, "r") do rdr; GZip.open(to, "w") do wtr
        write(wtr, read(rdr))
    end; end

    clean && rm(from)
end

gzip(path::AbstractString; clean=true) = gzip(path,"$path.gz"; clean=clean)

function polymorphisms(alignment)
    alleles = [ length(Set(row)) != 1 for row in eachrow(alignment) ]
    return alignment[alleles,:]
end

function identifiers(graph::Graphs.Graph)
    name = Dict{Graphs.Node, String}()
    for path in values(graph.sequence)
        counter = Dict{Graphs.Block, Int}()
        for node in path.node
            if node.block ∈ keys(counter)
                counter[node.block] += 1
                name[node] = "$(path.name)#$(counter[node.block])"
            else
                counter[node.block] = 1
                name[node] = "$(path.name)#1"
            end
        end
    end

    return function(node,_...)
        name[node]
    end
end

function emitblock(block::Graphs.Block, root, prefix, identifier; reduced=true)
    path = "$(root)/$(prefix)_na_aln.fa"
    open(path, "w") do io
        Graphs.marshal(io, block;
            fmt=:fasta,
            opt=(gaps=true, name=identifier)
        )
    end

    # grab reduced alignment
    alignment, scale = if reduced
        gzip(path; clean=true) # compress our previous alignment
        path = "$(root)/$(prefix)_na_aln_alleles.fa"

        align, nodes, _ = Graphs.alignment(block)
        numbp = size(align,1)

        align = polymorphisms(align)
        names = [identifier(node,i) for (i,node) in enumerate(nodes)]

        open(path, "w") do io
            for (name, seq) in zip(names, eachcol(align))
                Graphs.Utility.write_fasta(io, name, seq)
            end
        end

        path, size(align,1)/numbp
    else
        path, 1
    end

    out = IOBuffer()
    run(pipeline(
        `fasttree -nt -gtr $(alignment)`, stdout=out, stderr=devnull);
        wait=true
    )
    tree = Phylo.tree(String(take!(out)))

    Phylo.binary!(tree)
    Phylo.root!(tree)
    Phylo.rescale!(tree, scale)

    open("$root/$prefix.nwk", "w") do io
        Phylo.write(io, tree)
    end

    open("$root/$(prefix)_tree.json", "w") do io
        out, _ = Phylo.dictionary(tree,0,0;mutations=true)
        JSON.print(io, out)
    end

    gzip(path; clean=true)

    # emit a reduced alignment (panX by default looks for this)
    align, nodes, consensus = Graphs.alignment(block)
    reduce = [ length(Set(row)) == 1 for row in eachrow(align) ]
    align[reduce,:] .= '.'

    names = [identifier(node,i) for (i,node) in enumerate(nodes)]
    GZip.open("$root/$(prefix)_na_aln_reduced.fa.gz", "w") do io
        Graphs.Utility.write_fasta(io, "consensus", consensus)

        for (name, seq) in zip(names, eachcol(align))
            seq[seq.==consensus] .= '.'
            Graphs.Utility.write_fasta(io, name, seq)
        end
    end

    return tree
end

function coreblocks(G::Graphs.Graph, identifier)
    core  = Graphs.Block[]
    depth = length(G.sequence)
    for block in values(G.block)
        Graphs.depth(block) == depth || continue
        unique = reduce(&,
            map(collect(keys(block))) do node
                split(identifier(node),'#')[end] == "1"
            end
        )
        unique || continue

        push!(core, block)
    end

    return core
end

function emitcore(genes::Array{Graphs.Block}, root::String, identifier)
    length(genes) ≥ 1 || return

    isolate = (node) -> split(identifier(node),'#')[1]

    one, nodes, _ = Graphs.alignment(genes[1])
    order = Dict( isolate(node)=>i for (i,node) in enumerate(nodes) )

    aln = if length(genes) > 1
        reduce(vcat,
            map(genes[2:end]) do blk
                seq, ns, _ = Graphs.alignment(blk)
                permute = [order[isolate(n)] for n in ns]
                seq[:,permute]
            end;
        init=one)
    else
        one
    end

    numbp = size(aln,2)
    aln = polymorphisms(aln)
    scale = size(aln,2)/numbp

    names = map(isolate, nodes)
    path  = "$(root)/core_na_reduced.fa"
    open(path, "w") do io
        for (name, seq) in zip(names, eachcol(aln))
            Graphs.Utility.write_fasta(io, name, seq)
        end
    end

    out = IOBuffer()
    run(pipeline(
        `fasttree -nt -gtr $(path)`, stdout=out, stderr=devnull);
        wait=true
    )
    tree = Phylo.tree(String(take!(out)))
    gzip(path)

    Phylo.binary!(tree)
    Phylo.root!(tree)
    Phylo.rescale!(tree, scale)

    # emit as newick
    open("$root/strain_tree.nwk", "w") do io
        Phylo.write(io, tree)
    end

    # emit as bespoke json
    open("$root/coreGenomeTree.json", "w") do io
        out, _ = Phylo.dictionary(tree, 0, 0,)
        JSON.print(io, out)
    end
end

function fmt(i::Int)
    s = string(i)
    @match length(s) begin
        1 => "0000000$s"
        2 => "000000$s"
        3 => "00000$s"
        4 => "0000$s"
        5 => "000$s"
        6 => "00$s"
        7 => "0$s"
        _ => s
    end
end

function emit(G::Graphs.Graph, root::String)
    isdir(root) || mkdir(root)
    isdir("$root/geneCluster") || mkdir("$root/geneCluster")

    idents = identifiers(G)
    blocks = sort(collect(values(G.block)); by=(blk)->Graphs.depth(blk), rev=true)

    # emit all blocks
    meter = Progress(length(blocks), 1, "Aligning blocks. Building trees...")
    metadata = Array{Dict}(undef,length(blocks))
    for (i,b) in enumerate(blocks)
        prefix = "RC$(fmt(i))"
        metadata[i] = Dict(
            "geneId"     => i,
            "geneLen"    => length(b),
            "count"      => Graphs.depth(b),
            "dupli"      => "no",
            "dup_detail" => "",
            "ann"        => b.uuid,
            "msa"        => prefix,
            "divers"     => Graphs.diversity(b),
            "allAnn"     => b.uuid,
            "GName"      => "none",
            "allGName"   => "none",
            "locus"      => [ idents(node) for node in keys(b) ],
        )
        emitblock(b, "$(root)/geneCluster", "$(prefix)", idents; reduced=length(b) ≥ 5e4)
        next!(meter)
    end

    open("$root/geneCluster.json","w") do io
        JSON.print(io, metadata)
    end

    # emit core genome
    core = coreblocks(G, idents)
    emitcore(core, root, idents)

    # emit dummy strain meta info
    open("$root/metainfo.tsv", "w") do io
        write(io, join(
            ["accession",
             "strain",
             "collection_date",
             "country",
             "host",
             "organism"
            ], "\t"),
        "\n")

        for isolate in keys(G.sequence)
            write(io, join([
                isolate,
                "unknown",
                "unknown",
                "unknown",
                isolate,
                "unknown",
            ],"\t"), "\n")
        end
    end

    metainfo = [
        Dict("strain" => "unknown", "host"   => isolate, "accession" => isolate)
    for isolate  in keys(G.sequence)]

    open("$root/strainMetainfo.json","w") do io
        JSON.print(io, Dict("data"=>metainfo))
    end

    open("$root/metaConfiguration.js","w") do io
        write(io, """
        var meta_details={
            "host": $(collect(keys(G.sequence)))
        },
        meta_display={
            "meta_display_order": ["host"],
            "color_options": {"host": {"menuItem": "host", "type": "discrete"}},
        };
        var association_columns=[];
        """)
    end
end

end
