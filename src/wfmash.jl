module WFMash

import ..PanGraph: PanContigs, Alignment
import ..PanGraph.Graphs.Utility: read_paf, write_fasta, uncigar
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
    hit.cigar = collect(uncigar(hit.cigar))
    hit = trim_flanking_mismatches!(hit)
    return hit
end

"""
    trim_flanking_mismatches!(hit::Alignment)
Remove leading and trailing deletions and insertions from the cigar,
so that the first and last entry is a match.
"""
function trim_flanking_mismatches!(hit::Alignment)

    function qry_trim!(hit, len, side)
        # side == true => trim from left side
        if side ⊻ (! hit.orientation)
            hit.qry.start += len
        else
            hit.qry.stop -= len
        end
    end

    # trim leading deletions/insertions
    trim_threshold_length = 15
    while (hit.cigar[1][2] != 'M') || (hit.cigar[1][1] <= trim_threshold_length)
        # println(stderr, "cleaning cigar start ", hit.cigar)
        len, type = popfirst!(hit.cigar)
        # println(stderr, "into", hit.cigar)
        if type == 'I'
            qry_trim!(hit, len, true)
        elseif type == 'D'
            hit.ref.start += len
        elseif type == 'M'
            qry_trim!(hit, len, true)
            hit.ref.start += len
            hit.matches -= len
        else
            error("unrecognized cigar type")
        end

        # catch the edge-case in which all the cigar string has been cleaned
        if length(hit.cigar) == 0
            hit.length = 0
            return hit
        end
    end
    
    # trim trailing deletions/insertions
    while (hit.cigar[end][2] != 'M') || (hit.cigar[end][1] <= trim_threshold_length)
        # println(stderr, "cleaning cigar end", hit.cigar)
        len, type = pop!(hit.cigar)
        # println(stderr, "len $len and type $type into -> ", hit.cigar)
        if type == 'I'
            qry_trim!(hit, len, false)
        elseif type == 'D'
            hit.ref.stop -= len
        elseif type == 'M'
            qry_trim!(hit, len, false)
            hit.ref.stop -= len
            hit.matches -= len
        else
            error("unrecognized cigar type")
        end
    end

    if (hit.qry.seq !== nothing) || (hit.ref.seq !== nothing)
        error("sequence trimming not yet implemented")
    end

    # re-evaluate length of the hit. It is necessary because wfmash does not report the standard length
    hit.length = sum([l for (l, t) in hit.cigar])
    # hit.matches = sum([l for (l, t) in hit.cigar if t == 'M'])

    return hit
end


"""
    align(ref::PanContigs, qry::PanContigs)
Align homologous regions of `qry` and `ref`.
Returns the list of intervals between pancontigs.
"""
function align(ref::PanContigs, qry::PanContigs)
    hits = nothing
    if ref != qry
        hits = mktempdir() do dir
            open("$dir/qry.fa","w") do io
                for (name, seq) in zip(qry.name, qry.sequence)
                    if length(seq) ≥ 95
                        write_fasta(io, name, seq)
                    end
                end
            end

            open("$dir/ref.fa","w") do io
                for (name, seq) in zip(ref.name, ref.sequence)
                    if length(seq) ≥ 95
                        write_fasta(io, name, seq)
                    end
                end
            end

            run(`samtools faidx $dir/ref.fa`)
            run(`samtools faidx $dir/qry.fa`)
            run(pipeline(`wfmash -g3,15,1 -p 85 -s 1k $dir/ref.fa $dir/qry.fa`,
                stdout="$dir/aln.paf",
                stderr="$dir/err.log"
               )
            )

            open(read_paf, "$dir/aln.paf")
        end
    else
        hits = mktempdir() do dir
            open("$dir/seq.fa","w") do io
                for (name, seq) in zip(qry.name, qry.sequence)
                    if length(seq) ≥ 95
                        write_fasta(io, name, seq)
                    end
                end
            end

            run(`samtools faidx $dir/seq.fa`)
            run(pipeline(`wfmash -X -g3,15,1 -p 85 -s 1k $dir/seq.fa $dir/seq.fa`,
                stdout="$dir/aln.paf",
                stderr="$dir/err.log"
               )
            )

            open(read_paf, "$dir/aln.paf")
        end
    end

    hits = map(recigar!, hits)
    # remove hits of length zero that might have been produced after trimming the cigar string
    hits = filter(h -> h.length > 0, hits)
    return hits
end

end
