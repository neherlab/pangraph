module PanX

using GZip
using JSON
using Logging
using Rematch
using ProgressMeter
using TreeTools

import ..Graphs

# ------------------------------------------------------------------------
# Phylogenetic tree manipulation

function rescale!(tree, by::Real)
	for n in nodes(tree)
        isroot(n) && continue
		ismissing(n.tau) && error("Can't rescale tree with missing branch length.")
		n.tau *= by
	end
	return nothing
end

dictionary(tree::TreeTools.Tree; mutations=false) = dictionary(tree.root, 0, 0; mutations)
function dictionary(tt_node::TreeTools.TreeNode, level::Int, i::Int; mutations=false)

    node_lab = label(tt_node)

    node = Dict(
        "name" => if !isleaf(tt_node)
            i += 1
            suffix = lpad(i-1, 5, "0")
            "NODE_$(suffix)"
        else
            node_lab
        end,
        "branch_length" => ismissing(branch_length(tt_node)) ? 0 : branch_length(tt_node),
    )

    if mutations
        # TODO: actually write mutations?
        node["muts"]    = ""
        node["aa_muts"] = ""
    else
        node["clade"] = level
    end

    if isleaf(tt_node)
        if mutations
            node["accession"]  = split(node_lab,'#')[1]
            node["annotation"] = "pan-contig"
        else
            node["attr"] = Dict(
                "host"   => node_lab,
                "strain" => node_lab,
            )
        end
        return node, i
    end

    node["children"] = Dict[]
    for child in children(tt_node)
        c, i = dictionary(child, level+1, i; mutations=mutations)
        push!(node["children"], c)
    end

    return node, i
end

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

function produce_tree(alignment, scale)
    # run treetime and capture output
    out = IOBuffer()
    run(pipeline(
        `fasttree -nt -gtr $(alignment)`, stdout=out, stderr=devnull);
        wait=true
    )
    tree_string = strip(String(take!(out)))
    close(out)

    # build and process tree
    tree = parse_newick_string(tree_string)
    TreeTools.binarize!(tree)
    lgger = ConsoleLogger(Logging.Error) # Custom logger to filter out warnings from `root!`
    with_logger(lgger) do
    	TreeTools.root!(tree; method=:midpoint) # tree remains binary
    end
    rescale!(tree, scale)
    return tree
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

    tree = produce_tree(alignment, scale)

    write("$root/$prefix.nwk", tree, "w"; internal_labels=false)

    open("$root/$(prefix)_tree.json", "w") do io
        out, _ = dictionary(tree; mutations=true)
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

    # return tree
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


    tree = produce_tree(path, scale)
    gzip(path)

    # emit as newick
    write("$root/strain_tree.nwk", tree, "w"; internal_labels=false)

    # emit as bespoke json
    open("$root/coreGenomeTree.json", "w") do io
        out, _ = dictionary(tree)
        JSON.print(io, out)
    end
end

fmt(i::Int) = lpad(i, 8, "0")

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
        node_ids = [idents(node) for node in keys(b)]
        not_dupl = all([endswith(id, "#1") for id in node_ids])
        metadata[i] = Dict(
            "geneId"     => i,
            "geneLen"    => length(b),
            "count"      => Graphs.depth(b),
            "dupli"      => not_dupl ? "no" : "yes",
            "dup_detail" => "",
            "ann"        => b.uuid,
            "msa"        => prefix,
            "divers"     => Graphs.diversity(b),
            "allAnn"     => b.uuid,
            "GName"      => "none",
            "allGName"   => "none",
            "locus"      => node_ids,
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
