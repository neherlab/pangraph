module MMseqs

import ..PanGraph: PanContigs, Alignment
import ..PanGraph.Graphs.Utility: read_paf, write_fasta, uncigar, read_mmseqs2
import ..PanGraph.Graphs.Shell: execute

export align


"""
    align(ref::PanContigs, qry::PanContigs, klen::Int64)
Align homologous regions of `qry` and `ref` using mmseqs easy-search.
`klen` tunes the kmer length. If `klen`=0 then mmseqs default is used.
Returns the list of intervals between pancontigs.
"""
function align(ref::PanContigs, qry::PanContigs, klen::Int64)

    hits = mktempdir() do dir
        # hits = let dir = mktempdir(; cleanup = false)

        qryfa, reffa = "$dir/qry.fa", "$dir/ref.fa"

        open("$qryfa", "w") do io
            for (name, seq) in zip(qry.name, qry.sequence)
                if length(seq) ≥ 95
                    write_fasta(io, name, seq)
                end
            end
        end

        if ref != qry
            open("$reffa", "w") do io
                for (name, seq) in zip(ref.name, ref.sequence)
                    if length(seq) ≥ 95
                        write_fasta(io, name, seq)
                    end
                end
            end
        else
            reffa = qryfa
        end

        # kmer length
        klen_opt = klen == 0 ? [] : ["-k", "$klen"]
        # format of output file
        fmat_outp = "query,qlen,qstart,qend,empty,target,tlen,tstart,tend,nident,alnlen,bits,cigar,fident,raw"

        run(
            pipeline(
                `mmseqs easy-search
                $qryfa $reffa $dir/res.paf $dir/tmp
                --threads 1
                --max-seq-len 10000
                -a
                --search-type 3
                --format-output $fmat_outp
                $klen_opt`,
                stdout = devnull,
                stderr = devnull,
            ),
            wait = true,
        )

        open(read_mmseqs2, "$dir/res.paf")
    end

    for hit in hits
        # transform the cigar string in tuples of (len, char)
        hit.cigar = collect(uncigar(hit.cigar))
    end

    return hits
end

end
