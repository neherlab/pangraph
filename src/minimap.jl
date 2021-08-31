module Minimap

import Libdl
import Base: show

import ..Utility: Alignment, Hit

export align

# ------------------------------------------------------------------------
# globals

const NO_DIAG = 0x001
const NO_DUAL = 0x002
const CIGAR   = 0x004

# ------------------------------------------------------------------------
# functions

minimap2 = Libdl.dlopen("vendor/minimap2/libminimap2.so")

idx_str  = Libdl.dlsym(minimap2, :mm_idx_str)
idx_stat = Libdl.dlsym(minimap2, :mm_idx_stat)
idx_free = Libdl.dlsym(minimap2, :mm_idx_destroy)

tbuf_init = Libdl.dlsym(minimap2, :mm_tbuf_init)
tbuf_free = Libdl.dlsym(minimap2, :mm_tbuf_destroy)

set_opt    = Libdl.dlsym(minimap2, :mm_set_opt)
update_opt = Libdl.dlsym(minimap2, :mm_mapopt_update)

minimap  = Libdl.dlsym(minimap2, :mm_map)
map_frag = Libdl.dlsym(minimap2, :mm_map_frag)

divergence = Libdl.dlsym(minimap2, :mm_event_identity)

# ------------------------------------------------------------------------
# types

struct Index
    handle :: Ptr{Cvoid}
end

struct Buffer
    handle :: Ptr{Cvoid}
end

# copied from minimap2.h

struct Extra
    capacity :: UInt32
    dp_score :: Int32
    dp_max   :: Int32
    dp_max2  :: Int32
    packed   :: UInt32 # n_ambi (30 bits) / strand (2 bits)
    n_cigar  :: UInt32
    # cigar UInt32[] (variable length array)
end

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

struct IndexOptions
    k    :: Cshort
    w    :: Cshort
    flag :: Cshort
    bucket_bits :: Cshort
    mini_batch_size :: Int64
    batch_size :: UInt64

    IndexOptions() = new()
end

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
	q_occ_frac :: Cfloat
	min_mid_occ :: Int32; max_mid_occ :: Int32
	mid_occ :: Int32
	max_occ :: Int32; max_max_occ :: Int32; occ_dist :: Int32;
	mini_batch_size :: Int64
	max_sw_mat :: Int64
	cap_kalloc :: Int64

	split_prefix :: Cstring

    function MapOptions(idx::Ref{IndexOptions})
        opt = Ref{MapOptions}(new())

        # TODO: error handling
        err = ccall(set_opt, Cint, (Cstring, Ptr{Cvoid}, Ptr{Cvoid}),
              C_NULL, idx, opt
        )
        err = ccall(set_opt, Cint, (Cstring, Ptr{Cvoid}, Ptr{Cvoid}),
              "asm10", idx, opt
        )

        opt[].flag |= (NO_DIAG|CIGAR)

        return opt
    end
end

struct PanContigs{T <: AbstractArray{S} where S <: AbstractString}
    name     :: T
    sequence :: T
end

# ------------------------------------------------------------------------
# wrapper

# reference index
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

function freeindex(ref::Index)
    ccall(idx_free, Cvoid, (Ptr{Cvoid},), ref.handle)
end

# thread buffer
function makebuffer()
    return Buffer(ccall(tbuf_init, Ptr{Cvoid}, ()))
end

function freebuffer(buf::Buffer)
    return ccall(tbuf_free, Cvoid, (Ptr{Cvoid},), buf.handle)
end

# cigar string
const cigartypes = ['M','I','D','N','S','H','P','E','X']
function unpackcigar(base::Ptr{UInt32}, n)
    return [let
        cg = unsafe_load(base,i)
        Int(cg >> 4), cigartypes[cg & 0xf + 1]
    end for i in 1:n]
end

# alignment
function align(ref::PanContigs, qry::PanContigs)
    iopt = Ref{IndexOptions}(IndexOptions()) 
    mopt = MapOptions(iopt) 
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
                Hit(qry.name[i], length(seq), hit.qs, hit.qe, nothing),
                Hit(ref.name[j], length(ref.sequence[j]), hit.rs, hit.re, nothing),
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

    return homologs
end

end
