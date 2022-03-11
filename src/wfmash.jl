module WFMash

import ..PanGraph: PanContigs, Alignment
import ..PanGraph.Graphs.Utility: read_paf, write_fasta
import ..PanGraph.Graphs.Shell: execute

export align

"""
    recigar!(hit::Alignment)

Transform the detailed cigar string returned from wfmash into the more conventional form returned by minimap2.
Wfmash returns detailed match/mismatch information that we do not need.
Merges them into one match category.
"""
function recigar!(hit::Alignment)
    buffer = IOBuffer()
    n = 0
    i₁, i₂ = 1, 1
    while i₁ ≤ length(hit.cigar)
        while i₂ ≤ length(hit.cigar) && isdigit(hit.cigar[i₂])
            i₂ += 1
        end

        if hit.cigar[i₂] == '=' || hit.cigar[i₂] == 'X'
            n += parse(Int, hit.cigar[i₁:i₂-1])
        else
            if n > 0
                print(buffer, "$(n)M")
                n = 0
            end
            print(buffer, hit.cigar[i₁:i₂])
        end

        i₁ = i₂ + 1
        i₂ = i₁
    end

    if n > 0
        print(buffer, "$(n)M")
        n = 0
    end

    hit.cigar = String(take!(buffer))
    return hit
end

"""
    align(ref::PanContigs, qry::PanContigs)

Align homologous regions of `qry` and `ref`.
Returns the list of intervals between pancontigs.
"""
function align(ref::PanContigs, qry::PanContigs)
    hits = mktempdir() do dir
        open("$dir/qry.fa","w") do qryio; open("$dir/ref.fa","w") do refio
            for (name, seq) in zip(qry.name, qry.sequence)
                write_fasta(qryio, name, seq)
            end

            for (name, seq) in zip(ref.name, ref.sequence)
                write_fasta(refio, name, seq)
            end

            run(`samtools faidx $dir/ref.fa`)
            result = execute(`wfmash $dir/ref.fa $dir/qry.fa`)

            result.out |> IOBuffer |> read_paf
        end; end
    end

    map(recigar!, hits)
end

end
