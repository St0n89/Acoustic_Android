# NetCDF3.jl — minimal NetCDF Classic (CDF-1) reader, pure Julia stdlib.
# Covers the GEBCO download-service crop files (non-record variable,
# typically "elevation" with scale_factor/add_offset/_FillValue attributes).
#
#   nc  = nc3_read("gebco.nc")              # first variable
#   nc  = nc3_read("gebco.nc", "elevation")
#   nc.elevation            -> Nc3Var (data::Array{Float64}, rows = dim 1)
#   nc3_georef(nc)          -> (lon::Vector, lat::Vector)
# Missing values become NaN; scale_factor/add_offset are applied.
#
# Format: NetCDF Classic (magic "CDF\x01"), big-endian, sections
# dim_list (tag 10) / gatt_list (tag 12) / var_list (tag 11).
# Variable record: name, type(u32), ndims(u32), dimids(u32 each),
# attr_list, vsize(u32), begin(u32). Non-record data sits at `begin`,
# length = prod(dims) * type size, big-endian.

struct Nc3Var
    name::String
    dims::Vector{String}
    attrs::Dict{String,Any}
    data::Array{Float64}   # dim order as stored; rows = first dim
end

struct Nc3File
    vars::Dict{String,Nc3Var}
    dims::Vector{Pair{String,Int}}
    gattrs::Dict{String,Any}
end
Base.getindex(f::Nc3File, k::AbstractString) = f.vars[k]

const NC3_TYPE_SIZE = Dict(0x01 => 1, 0x02 => 1, 0x03 => 2, 0x04 => 4, 0x05 => 4, 0x06 => 8)

# ---------- low-level reader over a byte buffer ----------
mutable struct Nc3Buf
    b::Vector{UInt8}
    p::Int
end

_nc3_u32(r) = begin
    v = UInt32(0)
    for q in 0:3
        v = (v << 8) | UInt32(r.b[r.p + q + 1])
    end
    r.p += 4
    return v
end
_nc3_i32(r) = reinterpret(Int32, _nc3_u32(r))
_nc3_name(r) = begin
    n = Int(_nc3_u32(r))
    s = String(r.b[r.p+1:r.p+n])
    r.p += n + ((4 - n % 4) % 4)
    return s
end
function _nc3_pad!(r, nbytes)
    r.p += (4 - nbytes % 4) % 4
end
function _nc3_values(r, typ::UInt32, n::Int)
    sz = NC3_TYPE_SIZE[typ]
    o0 = r.p
    out = Vector{Any}(undef, n)
    for k in 1:n
        o = o0 + (k - 1) * sz
        if typ == 0x01
            out[k] = Float64(reinterpret(Int8, r.b[o+1]))
        elseif typ == 0x02
            out[k] = String([r.b[o+1]])
        elseif typ == 0x03
            out[k] = Float64(reinterpret(Int16, reinterpret(UInt16, UInt16(r.b[o+1])) << 8 | UInt16(r.b[o+2])))
        elseif typ == 0x04
            v = (UInt32(r.b[o+1]) << 24) | (UInt32(r.b[o+2]) << 16) | (UInt32(r.b[o+3]) << 8) | UInt32(r.b[o+4])
            out[k] = Float64(reinterpret(Int32, v))
        elseif typ == 0x05
            bits = (UInt32(r.b[o+1]) << 24) | (UInt32(r.b[o+2]) << 16) | (UInt32(r.b[o+3]) << 8) | UInt32(r.b[o+4])
            out[k] = Float64(reinterpret(Float32, bits))
        else
            v = UInt64(0)
            for q in 1:8
                v = (v << 8) | UInt64(r.b[o+q])
            end
            out[k] = reinterpret(Float64, v)
        end
    end
    r.p = o0 + sz * n
    _nc3_pad!(r, sz * n)
    return out
end
function _nc3_attrs(r)
    d = Dict{String,Any}()
    _nc3_u32(r); cnt = Int(_nc3_u32(r))
    cnt == 0 && return d
    for _ in 1:cnt
        nm = _nc3_name(r)
        typ = UInt32(_nc3_u32(r))
        n = Int(_nc3_u32(r))
        vals = _nc3_values(r, typ, n)
        d[nm] = length(vals) == 1 ? vals[1] : vals
    end
    return d
end

function _nc3_scalar_raw(b::Vector{UInt8}, typ::UInt32, o::Int)
    if typ == 0x01
        return Float64(reinterpret(Int8, b[o+1]))
    elseif typ == 0x03
        return Float64(reinterpret(Int16, reinterpret(UInt16, UInt16(b[o+1])) << 8 | UInt16(b[o+2])))
    elseif typ == 0x04
        v = (UInt32(b[o+1]) << 24) | (UInt32(b[o+2]) << 16) | (UInt32(b[o+3]) << 8) | UInt32(b[o+4])
        return Float64(reinterpret(Int32, v))
    elseif typ == 0x05
        bits = (UInt32(b[o+1]) << 24) | (UInt32(b[o+2]) << 16) | (UInt32(b[o+3]) << 8) | UInt32(b[o+4])
        return Float64(reinterpret(Float32, bits))
    elseif typ == 0x06
        v = UInt64(0)
        for q in 1:8
            v = (v << 8) | UInt64(b[o+q])
        end
        return reinterpret(Float64, v)
    end
    error("nc3: unsupported type $typ")
end

function nc3_read(path::AbstractString, varname::AbstractString = "")
    b = read(path)
    (length(b) >= 4 && b[1:4] == UInt8[0x43, 0x44, 0x46, 0x01]) ||
        error("nc3: not a CDF-1 (classic) file: $path")
    r = Nc3Buf(b, 4)
    _nc3_u32(r)  # numrecs (unused)
    dims = Vector{Pair{String,Int}}()
    _nc3_u32(r); nd = Int(_nc3_u32(r))
    if nd != 0
        for _ in 1:nd
            nm = _nc3_name(r)
            push!(dims, nm => Int(_nc3_i32(r)))
        end
    end
    gattrs = _nc3_attrs(r)
    _nc3_u32(r); nv = Int(_nc3_u32(r))
    nv != 0 || error("nc3: file has no variables")
    decls = []  # (name, type, dimnames, dimlens, attrs, begin)
    for _ in 1:nv
        nm = _nc3_name(r)
        typ = UInt32(_nc3_u32(r))
        ndv = Int(_nc3_u32(r))
        dimidx = [Int(_nc3_i32(r)) for _ in 1:ndv]
        vattrs = _nc3_attrs(r)
        _nc3_u32(r)               # vsize (padded size)
        beg_ = Int(_nc3_u32(r))   # CDF-1: 32-bit offset
        dimnames = [dims[i+1][1] for i in dimidx]
        dimlens = [dims[i+1][2] for i in dimidx]
        push!(decls, (nm, typ, dimnames, dimlens, vattrs, beg_))
    end
    target = varname == "" ? decls[1] : begin
        hits = [d for d in decls if d[1] == varname]
        isempty(hits) && error("nc3: variable '$varname' not found")
        hits[1]
    end
    vars = Dict{String,Nc3Var}()
    for (nm, typ, dimnames, dimlens, vattrs, beg_) in decls
        sz = NC3_TYPE_SIZE[typ]
        n = prod(max.(dimlens, 1))
        dims_vl = length(dimlens) == 0 ? (1,) : Tuple(dimlens)
        data = Array{Float64}(undef, reverse(dims_vl)...)  # Julia column-major
        # read all non-record variables (GEBCO crops are small enough,
        # and the loader needs lon/lat alongside the data variable)
        raw = Float64[_nc3_scalar_raw(b, typ, beg_ + (k - 1) * sz) for k in 1:n]
        # NetCDF fastest dim is dim1 -> Julia column-major column = dim1
        data .= reshape(raw, reverse(dims_vl))
        # optional scale/offset/fill attributes
        has_sf = haskey(vattrs, "scale_factor")
        sf = has_sf ? Float64(vattrs["scale_factor"]) : 1.0
        ao = haskey(vattrs, "add_offset") ? Float64(vattrs["add_offset"]) : 0.0
        haskey(vattrs, "_FillValue") || (vattrs["_FillValue"] = nothing)
        for k in eachindex(data)
            rawv = data[k] * sf + ao
            if vattrs["_FillValue"] !== nothing && data[k] == Float64(vattrs["_FillValue"])
                rawv = NaN
            end
            data[k] = rawv
        end
        vars[nm] = Nc3Var(nm, dimnames, vattrs, data)
    end
    return Nc3File(vars, dims, gattrs)
end

function nc3_georef(nc::Nc3File; lon_name = "lon", lat_name = "lat")
    lon = get(nc, lon_name).data[:, 1]
    lat = get(nc, lat_name).data[:, 1]
    return lon, lat
end
