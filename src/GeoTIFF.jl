# GeoTIFF.jl — minimal uncompressed Float32 GeoTIFF reader, pure Julia stdlib.
# Covers what BathyApp exports and the EMODnet flavour (ModelTransformation,
# tag 34264) as well as standard tie-point + pixel-scale georeferencing.
#
#   g = gtif_read("profile.tif")
#   g.lon, g.lat (Vectors), g.z (lat x lon, Float64), g.missing
# Assumes north-up, single band, uncompressed (compression tag must be 1 =
# none; BathyApp exports uncompressed). Handles II and MM byte orders.

struct GTifGrid
    lon::Vector{Float64}
    lat::Vector{Float64}
    z::Matrix{Float64}   # [row, col]; row 1 = lat[1] (top if north-up)
    note::String
end

const _GT_SIZES = Dict(0x0001 => 1, 0x0003 => 2, 0x0004 => 4, 0x0005 => 8, 0x000c => 8, 0x0011 => 4, 0x0012 => 8)

function _gt_u16(b, o, le)
    if le
        return UInt16(b[o+1]) | (UInt16(b[o+2]) << 8)
    end
    return (UInt16(b[o+1]) << 8) | UInt16(b[o+2])
end
function _gt_u32(b, o, le)
    if le
        return UInt32(b[o+1]) | UInt32(b[o+2]) << 8 | UInt32(b[o+3]) << 16 | UInt32(b[o+4]) << 24
    end
    return UInt32(b[o+1]) << 24 | UInt32(b[o+2]) << 16 | UInt32(b[o+3]) << 8 | UInt32(b[o+4])
end
function _gt_f64(b, o, le)
    bits = UInt64(0)
    if le
        for q in 8:-1:1
            bits = (bits << 8) | UInt64(b[o+q])
        end
    else
        for q in 1:8
            bits = (bits << 8) | UInt64(b[o+q])
        end
    end
    return reinterpret(Float64, bits)
end
_gt_f32(b, o, le) = Float64(reinterpret(Float32, _gt_u32(b, o, le)))

function gtif_read(path::AbstractString)
    b = read(path)
    length(b) >= 8 || error("gtif: file too small: $path")
    le = b[1] == 0x49 && b[2] == 0x49   # "II"
    (le || (b[1] == 0x4D && b[2] == 0x4D)) || error("gtif: not a TIFF: $path")
    magic = _gt_u16(b, 2, le)
    magic == 42 || error("gtif: bad TIFF magic")
    ifd_off = _gt_u32(b, 4, le)

    tags = Dict{Int,Vector{Float64}}()
    tags_raw = Dict{Int,Int}()          # for count-dependent reading
    tags_count = Dict{Int,Int}()
    o = ifd_off
    n_entries = _gt_u16(b, o, le)
    o += 2
    local width, height, comp, spp, rps, samplefmt = 0, 0, 0, 0, 0, 0
    strip_offsets = Int[]
    strip_counts = Int[]
    local t_scale = Float64[];
    local t_tie = Float64[];
    local t_trans = Float64[]
    for _ in 1:n_entries
        tag = Int(_gt_u16(b, o, le))
        typ = Int(_gt_u16(b, o + 2, le))
        cnt = Int(_gt_u32(b, o + 4, le))
        valfield = o + 8
        sz = get(_GT_SIZES, typ, 0)
        valoff = sz * cnt <= 4 ? valfield : Int(_gt_u32(b, valfield, le))
        vals = Float64[]
        if sz > 0
            for k in 1:cnt
                q = valoff + (k - 1) * sz
                if typ == 0x0003
                    push!(vals, Float64(_gt_u16(b, q, le)))
                elseif typ == 0x0001
                    push!(vals, Float64(b[q+1]))
                elseif typ == 0x0004
                    push!(vals, Float64(_gt_u32(b, q, le)))
                elseif typ == 0x000c
                    push!(vals, _gt_f64(b, q, le))
                elseif typ == 0x0011
                    push!(vals, _gt_f32(b, q, le))
                else
                    error("gtif: unsupported tag type $typ (tag $tag)")
                end
            end
        end
        tags[tag] = vals
        tags_count[tag] = cnt
        o += 12
    end
    haskey(tags, 256) || error("gtif: no ImageWidth")
    haskey(tags, 257) || error("gtif: no ImageLength")
    width = Int(tags[256][1]); height = Int(tags[257][1])
    comp = get(tags, 259, [1.0])[1]
    comp == 1 || error("gtif: compressed TIFF not supported (tag 259 = $comp); export uncompressed")
    spp = Int(get(tags, 277, [1.0])[1])
    spp == 1 || error("gtif: SamplesPerPixel = $spp, expected 1")
    bps = get(tags, 258, [32.0])[1]
    bps == 32 || error("gtif: BitsPerSample = $bps, expected 32 (Float32)")
    rps = Int(get(tags, 278, [Float64(height)])[1])
    strip_offsets = Int.(tags[273])
    strip_counts = Int.(tags[279])
    t_scale = get(tags, 33550, Float64[])
    t_tie = get(tags, 42112, Float64[])
    t_trans = get(tags, 34264, Float64[])

    z = Matrix{Float64}(undef, height, width)
    row0 = 1
    for (so, sc) in zip(strip_offsets, strip_counts)
        nrows = min(rps, height - row0 + 1)
        ncols = Int(sc ÷ sizeof(Float32) ÷ nrows)
        ncols == width || error("gtif: strip width mismatch ($ncols != $width)")
        for rr in 1:nrows, cc in 1:width
            z[row0 + rr - 1, cc] = _gt_f32(b, so + ((rr - 1) * width + (cc - 1)) * 4, le)
        end
        row0 += nrows
    end

    # georeference
    local lon, lat, note
    if !isempty(t_trans)
        cnt = tags_count[34264]
        cnt == 16 || error("gtif: ModelTransformation count $cnt != 16")
        T = reshape(t_trans, (4, 4))   # column-major of stored row-major 4x4
        # stored row-major: row r of the file = T[r, :] here after reshape? TIFF
        # stores the matrix VALUE-wise in row-major; reshape(.,4,4) column-major
        # flips it, so use the transpose:
        T = permutedims(T)
        lon = [T[1, 4] + T[1, 1] * c + T[1, 2] * r for c in 0:(width-1)]
        lat = [T[2, 4] + T[2, 1] * c + T[2, 2] * r for r in 0:(height-1)]
        if length(lon) > 1 && lon[2] < lon[1]
            lon = reverse(lon); z = reverse(z, dims = 2)
        end
        if length(lat) > 1 && lat[2] > lat[1]
            lat = reverse(lat); z = reverse(z, dims = 1)
        end
        note = "georef: ModelTransformation (34264)"
    elseif !isempty(t_scale) && !isempty(t_tie)
        length(t_scale) >= 2 && length(t_tie) >= 6 || error("gtif: bad pixel-scale/tiepoint")
        sx, sy = t_scale[1], t_scale[2]
        tx, ty = t_tie[4], t_tie[5]
        if sy > 0
            lat = [ty - (r - t_tie[2]) * sy for r in 0:(height-1)]
        else
            sy = -sy
            lat = [ty + (r - t_tie[2]) * sy for r in 0:(height-1)]
        end
        lon = [tx + (c - t_tie[1]) * sx for c in 0:(width-1)]
        if length(lat) > 1 && lat[2] > lat[1]
            lat = reverse(lat); z = reverse(z, dims = 1)
        end
        note = "georef: tiepoint (42112) + pixel scale (33550)"
    else
        lon = Float64.(1:width)
        lat = Float64.(1:height)
        note = "no georef tags: pixel coordinates"
    end
    z[z .== -32768.0] .= NaN
    z[z .== 32767.0] .= NaN
    z[abs.(z) .> 1e6] .= NaN
    return GTifGrid(lon, lat, z, note)
end
