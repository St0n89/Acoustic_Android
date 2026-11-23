# QuickPlot.jl — Plots-free result rendering for the notebook (stdlib only).
# * tiny PNG encoder (deflate STORED blocks — no compressor needed, still a
#   fully valid PNG stream per RFC 1950/1951)
# * SVG line/profile rendering
# * reads the suite's meta.jl, traces.csv and frame_*.f32 outputs

using Printf: @sprintf
import Base64: base64encode

# ---------------- PNG (stored-block deflate) ----------------
function _qp_crc32(buf)
    # reflected CRC-32, standard table-less implementation
    crc = 0xFFFFFFFF
    for b in buf
        crc ⊻= UInt32(b)
        for _ in 1:8
            crc = (crc >>> 1) ⊻ (0xEDB88320 & (0 - (crc & 1)))
        end
    end
    return crc ⊻ 0xFFFFFFFF
end
function _qp_chunk(tag::String, payload::Vector{UInt8})
    out = UInt8[]
    len = reinterpret(UInt8, [hton(UInt32(length(payload)))])
    append!(out, len)
    body = UInt8[ codeunits(tag)..., payload... ]
    append!(out, body)
    append!(out, reinterpret(UInt8, [hton(UInt32(_qp_crc32(body)))]))
    return out
end
function qp_png(rgb::Array{UInt8,3}; flipy = false)
    h, w = size(rgb, 1), size(rgb, 2)
    @assert size(rgb, 3) == 3
    raw = UInt8[]
    for r in 1:h
        push!(raw, 0x00)  # filter: None
        row = flipy ? (h - r + 1) : r
        append!(raw, rgb[row, :, :]...)
    end
    idat = UInt8[]
    # zlib header (CMF=0x78, FLG=0x01: no compression, fastest check bits)
    push!(idat, 0x78); push!(idat, 0x01)
    o = 1
    n = length(raw)
    while true
        blk = min(65535, n - o + 1)
        final = (o + blk - 1) == n
        push!(idat, final ? 0x01 : 0x00)
        len = UInt16(blk - 1); nlen = UInt16(~len)
        push!(idat, reinterpret(UInt8, [len])...)
        push!(idat, reinterpret(UInt8, [nlen])...)
        append!(idat, raw[o:o+blk-1])
        (o + blk - 1) == n && break
        o += blk
    end
    append!(idat, reinterpret(UInt8, [hton(UInt32(_qp_adler32(raw)))]))
    # IHDR
    ihdr = UInt8[]
    append!(ihdr, reinterpret(UInt8, [hton(UInt32(w))]))
    append!(ihdr, reinterpret(UInt8, [hton(UInt32(h))]))
    append!(ihdr, [0x08, 0x02, 0x00, 0x00, 0x00])  # 8-bit, RGB
    png = UInt8[0x89, 0x50, 0x4E, 0x47, 0x0D, 0x0A, 0x1A, 0x0A]
    append!(png, _qp_chunk("IHDR", ihdr))
    append!(png, _qp_chunk("IDAT", idat))
    append!(png, _qp_chunk("IEND", UInt8[]))
    return png
end
function _qp_adler32(buf)
    a, b = UInt32(1), UInt32(0)
    for x in buf
        a = (a + UInt32(x)) % 65521
        b = (b + a) % 65521
    end
    return (b << 16) | a
end
qp_uri(bytes::Vector{UInt8}) = "data:image/png;base64," * base64encode(bytes)

# ---------------- palettes ----------------
# turbo-like 8-stop ramp (loud = red for TL fields)
const _QP_TURBO = [
    (0.19, 0.07, 0.42), (0.27, 0.35, 0.73), (0.25, 0.60, 0.72),
    (0.30, 0.82, 0.54), (0.65, 0.89, 0.30), (0.99, 0.73, 0.22),
    (0.90, 0.36, 0.14), (0.62, 0.09, 0.03)]
function qp_turbo(t)
    t = clamp(t, 0.0, 1.0)
    x = t * (length(_QP_TURBO) - 1)
    i = clamp(Int(floor(x)) + 1, 1, length(_QP_TURBO) - 1)
    f = x - (i - 1)
    a, b = _QP_TURBO[i], _QP_TURBO[i+1]
    return (a[1] * (1 - f) + b[1] * f, a[2] * (1 - f) + b[2] * f, a[3] * (1 - f) + b[3] * f)
end
const _QP_WATER = [(0.02,0.09,0.22), (0.05,0.22,0.45), (0.10,0.38,0.62)]

"""
    qp_heatmap_png(mat; W = 360, cmap = qp_turbo, invert = true,
                   lo = nothing, hi = nothing, water_mask = nothing)
Downsample mat to <= W columns and render a PNG; returns (uri, (lo, hi)).
`invert = true` maps high values to RED (TL convention loud = red).
"""
function qp_heatmap_png(mat; W = 360, cmap = qp_turbo, invert = true,
                        lo = nothing, hi = nothing, water_mask = nothing)
    H0, W0 = size(mat, 1), size(mat, 2)
    wcols = min(W, W0)
    hrows = max(Int(round(H0 * (wcols / W0))), 1)
    sx, sy = max(W0 ÷ wcols, 1), max(H0 ÷ hrows, 1)
    wcols = min(wcols, cld(W0, sx)); hrows = min(hrows, cld(H0, sy))
    sub = [begin
            vals = Float64[]
            for jj in (1 + (i - 1) * sy):min(i * sy, H0),
                       kk in (1 + (j - 1) * sx):min(j * sx, W0)
                push!(vals, mat[jj, kk])
            end
            isempty(vals) ? NaN : sum(vals) / length(vals)
        end for i in 1:hrows, j in 1:wcols]
    m = wm -> wm === nothing ? nothing :
        [begin
            vals = Float64[]
            for jj in (1 + (i - 1) * sy):min(i * sy, H0),
                       kk in (1 + (j - 1) * sx):min(j * sx, W0)
                push!(vals, wm[jj, kk])
            end
            isempty(vals) ? 0.0 : sum(vals) / length(vals)
        end for i in 1:hrows, j in 1:wcols]
    wm = m(water_mask)
    vals = [sub[i,j] for i in 1:hrows, j in 1:wcols if !isnan(sub[i,j])]
    isempty(vals) && return qp_uri(qp_png(fill(UInt8(0), hrows, wcols, 3))), (0.0, 1.0)
    if lo === nothing || hi === nothing
        q = sort(vals)
        lo = q[max(1, clamp(Int(round(0.05 * length(q))), 1, length(q)))]
        hi = q[max(1, clamp(Int(round(0.95 * length(q))), 1, length(q)))]
    end
    rgb = Array{UInt8}(undef, hrows, wcols, 3)
    for i in 1:hrows, j in 1:wcols
        v = sub[i,j]
        t = isnan(v) ? 0.0 : clamp((v - lo) / max(hi - lo, eps()), 0.0, 1.0)
        invert && (t = 1.0 - t)
        c = cmap(t)
        land = wm !== nothing && wm[i,j] > 0.5
        rgb[i,j,1] = UInt8(round(c[1] * 255))
        rgb[i,j,2] = UInt8(round(c[2] * 255))
        rgb[i,j,3] = UInt8(round(c[3] * 255))
        land && (rgb[i,j,1] = 0x36; rgb[i,j,2] = 0x64; rgb[i,j,3] = 0x4E)
    end
    return qp_uri(qp_png(rgb)), (lo, hi)
end

function qp_html_img(uri; width = "100%")
    return HTML("<img src=\"$uri\" style=\"width:$width; image-rendering:pixelated; border-radius:6px\">")
end

# ---------------- SVG ----------------
function qp_profile_svg(x, y; W = 640, H = 220, title = "", color = "#0af", land = :nothing)
    length(x) < 2 && return HTML("<em>profile too short</em>")
    xmin, xmax_ = extrema(x)
    ymax = maximum(y); ymin = min(0.0, minimum(y))
    X(v) = 30 + (v - xmin) / max(xmax_ - xmin, eps()) * (W - 45)
    Y(v) = H - 25 - (v - ymin) / max(ymax - ymin, eps()) * (H - 45)
    pts = join([@sprintf("%.1f,%.1f", X(x[k]), Y(y[k])) for k in eachindex(x)], " ")
    t = title == "" ? "" : "<text x='10' y='16' fill='#8ab' font-size='12'>$title</text>"
    return HTML("""
    <svg viewBox="0 0 $W $H" style="width:100%%; background:#0d1117; border-radius:6px">
      $t
      <line x1='30' y1='$(Y(0))' x2='$(W-15)' y2='$(Y(0))' stroke='#345' stroke-width='1'/>
      <polyline points='$pts' fill='none' stroke='$color' stroke-width='2'/>
      <text x='10' y='$(H-8)' fill='#8ab' font-size='11'>x: $(round(Int,xmin))..$(round(Int,xmax_)) m</text>
      <text x='$(W-80)' y='16' fill='#8ab' font-size='11'>max $(round(ymax,digits=1)) m</text>
    </svg>""")
end

# ---------------- suite outputs ----------------
function qp_read_meta(path)
    d = Dict{String,Any}()
    for line in eachline(path)
        m = match(r"^([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.*)$", strip(line))
        m === nothing && continue
        k, v = m.captures
        if (s = match(r"^\"(.*)\"$", v)) !== nothing
            d[k] = s.captures[1]
        elseif (n = tryparse(Float64, v)) !== nothing
            d[k] = (n === round(n) && '.' ∉ v && 'e' ∉ v) ? Int(n) : n
        elseif (i = match(r"^\[(.*)\]$", v)) !== nothing
            nums = tryparse.(Float64, split(i.captures[1], [',', ' ']))
            d[k] = [x for x in nums if x !== nothing]
        else
            d[k] = v
        end
    end
    return d
end

function qp_read_traces(path)
    lines = readlines(path)
    hdr = split(lines[1], ',')
    n = length(lines) - 1
    cols = length(hdr)
    A = Matrix{Float64}(undef, n, cols)
    for r in 1:n
        parts = split(lines[r+1], ',')
        for c in 1:cols
            A[r, c] = something(tryparse(Float64, parts[c]), 0.0)
        end
    end
    return A, hdr
end

function qp_traces_svg(path; W = 640, H = 260, dec = 1)
    A, hdr = qp_read_traces(path)
    nt, nc = size(A)
    t = A[:, 1]
    out = IOBuffer()
    print(out, "<svg viewBox='0 0 $W $H' style='width:100%%; background:#0d1117; border-radius:6px'>")
    print(out, "<text x='8' y='14' fill='#8ab' font-size='12'>$(nc-1) traces, $(Int(round(t[end]*1000))) ms</text>")
    step = max(1, nt ÷ 900)
    for c in 2:nc
        wa = maximum(abs.(A[1:step:end, c])) + eps()
        ox = 30 + (c - 2) / max(nc - 2, 1) * (W - 60)
        pts = join([@sprintf("%.1f,%.1f",
                             ox + A[r, c] / wa * (W - 60) / max(nc - 2, 1) / 2.2,
                             H - 30 - (A[r, 1] - t[1]) / max(t[end] - t[1], eps()) * (H - 50))
                    for r in 1:step:nt], " ")
        print(out, "<polyline points='$pts' fill='none' stroke='#4fd' stroke-width='0.8' opacity='0.85'/>")
        print(out, "<line x1='$(ox)' y1='20' x2='$(ox)' y2='$(H-30)' stroke='#234' stroke-width='0.5'/>")
    end
    print(out, "<line x1='30' y1='$(H-30)' x2='$(W-30)' y2='$(H-30)' stroke='#567'/>")
    print(out, "<text x='$(W-120)' y='$(H-10)' fill='#8ab' font-size='11'>t = 0 .. $(round(t[end],digits=2)) s</text>")
    print(out, "</svg>")
    return HTML(String(take!(out)))
end

function qp_read_f32(path)
    b = read(path)
    n = div(length(b), 4)
    n == 0 && return Float64[]
    return Float64.(reinterpret(Float32, b))
end
