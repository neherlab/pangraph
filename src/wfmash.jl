module WFMash

import ..PanGraph: PanContig
import ..PanGraph.Graph.Utility: read_paf, write_fasta

export align

function align(ref::PanContigs, qry::PanContigs)
    mktemp() do (qry_path, qry_io); mktemp() do (ref_path, ref_io)
        # qry output
        for (name, seq) in zip(qry.name, qry.sequence)
            write_fasta(qry_io, name, seq)
        end
        # ref output
        for (name, seq) in zip(ref.name, ref.sequence)
            write_fasta(ref_io, name, seq)
        end
        # samtools on ref
        run(`samtools faidx $ref_path`)

        cmd = `wfmash $ref_path $qry_path`
    end; end; end
end

end
