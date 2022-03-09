module WFMash

import ..PanGraph: PanContig
import ..PanGraph.Graph.Utility: read_paf, write_fasta
import ..PanGraph.Graph.Shell: execute

export align

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
            print(buffer, hit.cigar[i₁:i₂-1])
        end

        i₁ = i₂ + 1
        i₂ = i₁
    end

    hit.cigar = take!(buffer)
end

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
            result = execute(`wfmash $dir/ref.fa $dir/$qry.fa`)

            result.out |> read_paf |> collect
        end; end
    end

    map(recigar!, hits)
end

end
