module MMseqs

import ..PanGraph: PanContigs, Alignment
import ..PanGraph.Graphs.Utility: read_paf, write_fasta, uncigar, read_mmseqs2
import ..PanGraph.Graphs.Shell: execute

export align

function log(msg...)
    println(stderr, msg...)
    flush(stderr)
end

"""
    align(ref::PanContigs, qry::PanContigs)
Align homologous regions of `qry` and `ref`.
Returns the list of intervals between pancontigs.
"""
function align(ref::PanContigs, qry::PanContigs)

    hits = mktempdir() do dir
        # hits = let dir = mktempdir(; cleanup = false)
        # log("starting to aling ref and qry in $dir")

        qrydb, refdb = "$dir/qry", "$dir/ref"
        if ref == qry
            refdb = qrydb
        end

        open("$qrydb.fa", "w") do io
            for (name, seq) in zip(qry.name, qry.sequence)
                if length(seq) ≥ 95
                    write_fasta(io, name, seq)
                end
            end
        end

        # log("mmseqs createdb")
        run(
            pipeline(
                `mmseqs createdb $qrydb.fa $qrydb`,
                stdout = "$dir/out_createdb_qry.log",
                stderr = "$dir/err_createdb_qry.log",
            ),
        )

        if ref != qry

            open("$refdb.fa", "w") do io
                for (name, seq) in zip(ref.name, ref.sequence)
                    if length(seq) ≥ 95
                        write_fasta(io, name, seq)
                    end
                end
            end

            run(
                pipeline(
                    `mmseqs createdb $refdb.fa $refdb`,
                    stdout = "$dir/out_createdb_ref.log",
                    stderr = "$dir/err_createdb_ref.log",
                ),
            )

        end

        # log("mmseqs search")
        run(
            pipeline(
                `mmseqs search 
                --threads 1 
                -a 
                --max-seq-len 10000
                --search-type 3
                $qrydb $refdb $dir/res $dir/tmp`,
                stdout = "$dir/out_search.log",
                stderr = "$dir/err_search.log",
            ),
        )

        # log("mmseqs convertalis")
        run(
            pipeline(
                `mmseqs convertalis
                --threads 1 
                --search-type 3
                $qrydb $refdb $dir/res $dir/res.paf
                --format-output query,qlen,qstart,qend,empty,target,tlen,tstart,tend,nident,alnlen,bits,cigar,fident,raw`,
                stdout = "$dir/out_convert.log",
                stderr = "$dir/err_convert.log",
            ),
        )

        # log("parse paf file")
        open(read_mmseqs2, "$dir/res.paf")
    end

    for hit in hits
        # transform the cigar string in tuples of (len, char)
        hit.cigar = collect(uncigar(hit.cigar))
    end

    return hits
end

end
