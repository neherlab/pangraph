module Minimap

using minimap2_jll

import Base: show

import ..PanGraph: PanContigs
import ..Utility: Alignment, Hit

export align

# ------------------------------------------------------------------------
# globals

const NO_DIAG = 0x001
const NO_DUAL = 0x002
const CIGAR   = 0x004

# ------------------------------------------------------------------------
# functions

const idx_str  = (:mm_idx_str, minimap2)
const idx_stat = (:mm_idx_stat, minimap2)
const idx_free = (:mm_idx_destroy, minimap2)

const tbuf_init = (:mm_tbuf_init, minimap2)
const tbuf_free = (:mm_tbuf_destroy, minimap2)

const set_opt    = (:mm_set_opt, minimap2)
const update_opt = (:mm_mapopt_update, minimap2)

const minimap  = (:mm_map, minimap2)
const map_frag = (:mm_map_frag, minimap2)

const divergence = (:mm_event_identity, minimap2)

# ------------------------------------------------------------------------
# types

"""
	struct Index
		handle :: Ptr{Cvoid}
	end

Store the untyped address to a minimap2 sequence index (set of minimizers).
"""
struct Index
    handle :: Ptr{Cvoid}
end

"""
	struct Buffer
		handle :: Ptr{Cvoid}
	end

Store the untyped address to a minimap2 sequence buffer (working space).
"""
struct Buffer
    handle :: Ptr{Cvoid}
end

# copied from minimap2.h
"""
	struct Extra
		capacity :: UInt32
		dp_score :: Int32
		dp_max   :: Int32
		dp_max2  :: Int32
		packed   :: UInt32 # n_ambi (30 bits) / strand (2 bits)
		n_cigar  :: UInt32
		cigar UInt32[] (variable length array)
	end

Copied from minimap2.h. See mm_extra_t.
"""
struct Extra
    capacity :: UInt32
    dp_score :: Int32
    dp_max   :: Int32
    dp_max2  :: Int32
    packed   :: UInt32 # n_ambi (30 bits) / strand (2 bits)
    n_cigar  :: UInt32
    # cigar UInt32[] (variable length array)
end

"""
	struct Record
		id :: Int32; cnt :: Int32; rid :: Int32; score :: Int32
		qs :: Int32; qe  :: Int32; rs  :: Int32; re :: Int32

		parent :: Int32; subsc :: Int32

		as     :: Int32
		mlen   :: Int32; blen :: Int32
		nsub   :: Int32;
		score0 :: Int32;

		pack :: UInt32 # mapq(8) ; split(2); rev(1); inv(1); sam_pri(1); proper_frag(1); pe_thru(1); seg_split(1); seg_id(8); split_inv(1); is_alt(1); dummy(6)
		hash :: UInt32
		div  :: Cfloat
		p    :: Ptr{Extra}
	end

Copied from minimap2.h. See mm_reg1_t.
"""
struct Record
    id :: Int32; cnt :: Int32; rid :: Int32; score :: Int32
    qs :: Int32; qe  :: Int32; rs  :: Int32; re :: Int32

    parent :: Int32; subsc :: Int32

    as     :: Int32
    mlen   :: Int32; blen :: Int32
    nsub   :: Int32;
    score0 :: Int32;

    pack :: UInt32 # mapq(8) ; split(2); rev(1); inv(1); sam_pri(1); proper_frag(1); pe_thru(1); seg_split(1); seg_id(8); split_inv(1); is_alt(1); dummy(6)
    hash :: UInt32
    div  :: Cfloat
    p    :: Ptr{Extra}
end

"""
	struct IndexOptions
		k    :: Cshort
		w    :: Cshort
		flag :: Cshort
		bucket_bits :: Cshort
		mini_batch_size :: Int64
		batch_size :: UInt64

		IndexOptions() = new()
	end

Copied from minimap2.h. See mm_idxopt_t.
"""
struct IndexOptions
    k    :: Cshort
    w    :: Cshort
    flag :: Cshort
    bucket_bits :: Cshort
    mini_batch_size :: Int64
    batch_size :: UInt64

    IndexOptions() = new()
end

"""
	mutable struct MapOptions
		flag :: Int64
		seed :: Cint
		sdust_thres :: Cint

		max_qlen :: Cint

		bw :: Cint; bw_long :: Cint
		max_gap :: Cint; max_gap_ref :: Cint 
		max_frag_len :: Cint
		max_chain_skip :: Cint; max_chain_iter :: Cint
		min_cnt :: Cint
		min_chain_score :: Cint
		chain_gap_scale :: Cfloat
		rmq_size_cap :: Cint; rmq_inner_dist :: Cint
		rmq_rescue_size :: Cint
		rmq_rescue_ratio :: Cfloat

		mask_level :: Cfloat
		mask_len :: Cint
		pri_ratio :: Cfloat
		best_n :: Cint

		alt_drop :: Cfloat

		a :: Cint; b :: Cint; q :: Cint; e :: Cint; q2 :: Cint; e2 :: Cint
		sc_ambi :: Cint
		noncan :: Cint 
		junc_bonus :: Cint
		zdrop :: Cint; zdrop_inv :: Cint
		end_bonus :: Cint
		min_dp_max :: Cint
		min_ksw_len :: Cint
		anchor_ext_len :: Cint; anchor_ext_shift :: Cint
		max_clip_ratio :: Cfloat

		rank_min_len :: Cint
		rank_frac :: Cfloat

		pe_ori :: Cint; pe_bonus :: Cint

		mid_occ_frac :: Cfloat
		min_mid_occ :: Int32; max_mid_occ :: Int32
		mid_occ :: Int32
		max_occ :: Int32; max_max_occ :: Int32; occ_dist :: Int32;
		mini_batch_size :: Int64
		max_sw_mat :: Int64
		cap_kalloc :: Int64

		split_prefix :: Cstring
	end

Copied from minimap2.h. See mm_mapopt_t.
"""
mutable struct MapOptions
	flag :: Int64
	seed :: Cint
	sdust_thres :: Cint

	max_qlen :: Cint

	bw :: Cint; bw_long :: Cint
	max_gap :: Cint; max_gap_ref :: Cint 
	max_frag_len :: Cint
	max_chain_skip :: Cint; max_chain_iter :: Cint
	min_cnt :: Cint
	min_chain_score :: Cint
	chain_gap_scale :: Cfloat
	rmq_size_cap :: Cint; rmq_inner_dist :: Cint
	rmq_rescue_size :: Cint
	rmq_rescue_ratio :: Cfloat

	mask_level :: Cfloat
	mask_len :: Cint
	pri_ratio :: Cfloat
	best_n :: Cint

	alt_drop :: Cfloat

	a :: Cint; b :: Cint; q :: Cint; e :: Cint; q2 :: Cint; e2 :: Cint
	sc_ambi :: Cint
	noncan :: Cint 
	junc_bonus :: Cint
	zdrop :: Cint; zdrop_inv :: Cint
	end_bonus :: Cint
	min_dp_max :: Cint
	min_ksw_len :: Cint
	anchor_ext_len :: Cint; anchor_ext_shift :: Cint
	max_clip_ratio :: Cfloat

	rank_min_len :: Cint
	rank_frac :: Cfloat

	pe_ori :: Cint; pe_bonus :: Cint

	mid_occ_frac :: Cfloat
	min_mid_occ :: Int32; max_mid_occ :: Int32
	mid_occ :: Int32
	max_occ :: Int32; max_max_occ :: Int32; occ_dist :: Int32;
	mini_batch_size :: Int64
	max_sw_mat :: Int64
	cap_kalloc :: Int64

	split_prefix :: Cstring

    function MapOptions(idx::Ref{IndexOptions}, minblock::Int, preset::String)
        opt = Ref{MapOptions}(new())

        # TODO: error handling
        err = ccall(set_opt, Cint, (Cstring, Ptr{Cvoid}, Ptr{Cvoid}),
              C_NULL, idx, opt
        )
        err = ccall(set_opt, Cint, (Cstring, Ptr{Cvoid}, Ptr{Cvoid}),
              preset, idx, opt
        )

        # manual overrides
        opt[].min_dp_max = minblock-10
        opt[].flag |= (NO_DIAG|CIGAR)

        return opt
    end
end

# ------------------------------------------------------------------------
# wrapper

# reference index
"""
	makeindex(w, k, names, sequence; bucketbits::Int=14)

Given a window size `w` and kmer length `k`, return a handle to a minimizer index
for sequences `sequence`.
"""
function makeindex(w, k, names, sequence; bucketbits::Int=14)
    @assert length(names) == length(sequence)

    return Index(
        ccall(idx_str, Ptr{Cvoid}, (Cint, Cint, Cint, Cint, Cint, Ptr{Cstring}, Ptr{Cstring}),
          w, k, 0, bucketbits, length(names), sequence, names
        )
    )
end

function show(ref::Index)
    ccall(idx_stat, Cvoid, (Ptr{Cvoid},), ref.handle)
end

"""
	freeindex()

Free memory associated to an opaque handle to a sequence index.
"""
function freeindex(ref::Index)
    ccall(idx_free, Cvoid, (Ptr{Cvoid},), ref.handle)
end

# thread buffer
"""
	makebuffer()

Return an opaque handle thread buffer
"""
function makebuffer()
    return Buffer(ccall(tbuf_init, Ptr{Cvoid}, ()))
end

"""
	freebuffer()

Free memory associated to an opaque handle thread buffer
"""
function freebuffer(buf::Buffer)
    return ccall(tbuf_free, Cvoid, (Ptr{Cvoid},), buf.handle)
end

# cigar string
const cigartypes = ['M','I','D','N','S','H','P','E','X']
unpackcigar(base::Ptr{UInt32}, n) = [let
    cg = unsafe_load(base,i)
    (Int(cg>>4), cigartypes[cg&0xf+1])
end for i in 1:n]

# alignment
"""
	align(ref::PanContigs, qry::PanContigs, minblock::Int, preset::String)

Call into minimap to align the set of blocks `qry` to blocks `ref`.
Preset should be a string ∈ ["asm5","asm10","asm20"]. See minimap2 manual for details.
This is probably the function you want.
If you call into the function specifically, all memory management is taken care of for you.
"""
function align(ref::PanContigs, qry::PanContigs, minblock::Int, preset::String)
    iopt = Ref{IndexOptions}(IndexOptions())
    mopt = MapOptions(iopt, minblock, preset)
    len  = Ref{Cint}(0)

    index  = makeindex(iopt[].w, iopt[].k, ref.name, ref.sequence)

    buffer = makebuffer()

    ccall(update_opt, Cvoid, (Ptr{MapOptions}, Ptr{Cvoid}), mopt, index.handle)

    homologs = Alignment[]
    for (i,seq) in enumerate(qry.sequence)
        hits = ccall(minimap, Ptr{Record},
                (Ptr{Cvoid}, Cint, Ptr{Cchar}, Ptr{Cint}, Ptr{Cvoid}, Ptr{Cvoid}, Cstring),
                index.handle, length(seq), seq, len, buffer.handle, mopt, qry.name[i]
        )

        for n in 1:len[]
            hit = unsafe_load(hits, n)
            aln = unsafe_load(hit.p, 1)

            j = hit.rid+1 # minimap2 0 indexes
            push!(homologs, Alignment(
                Hit(qry.name[i], length(seq), hit.qs+1, hit.qe, nothing),
                Hit(ref.name[j], length(ref.sequence[j]), hit.rs+1, hit.re, nothing),
                Int(hit.mlen),
                Int(hit.blen),
                Int(hit.pack & 0xff),
                ((hit.pack >> 10) & 1) == 0,
                unpackcigar(Ptr{UInt32}(hit.p+Core.sizeof(Extra)), aln.n_cigar),
                1 - ccall(divergence, Cdouble, (Ptr{Record},), Ref(hit)),
                nothing
            ))

            Libc.free(hit.p)
        end
        Libc.free(hits)
    end

    freebuffer(buffer)
    freeindex(index)

    return homologs
end

end
