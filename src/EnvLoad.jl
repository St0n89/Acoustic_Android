# EnvLoad.jl — BathyApp/GEBCO data -> FDTD2DElastic environment configs.
#
# Generates config files with EXACTLY the variable contract of the suite's
# config/env_pekeris.jl / config/sources_receivers.jl (2-D) and
# config/env3d.jl / config/sources_receivers3d.jl (3-D), so the unmodified
# drivers run.jl / run3d.jl pick them up:
#
#   2-D env contract:  dx, dz, xmax, zmax, water_depth, c_water, rho_water,
#       c_bottom, rho_bottom, Qp_water, Qp_bottom, pml_cells, nx, nz, x, z,
#       grid = ElasticGrid(nx, nz, dx, dz, lam, mu, rho, Qp, Qs)
#   2-D sr  contract:  f0, source_kind, src_amp, source_x, source_z,
#       check_f_max, sources, monopoles, bc, receiver_ranges_km,
#       receiver_depth, receivers, SL_REF_RANGE, SL_REF_IS_RECEIVER
#   3-D equivalents add y and use grid3d / Source3D / MonopoleSource3D /
#       BC3D / Receiver3D.
#
# Profile convention INSIDE this app: `depth(x)` positive DOWN in metres
# (bathymetry from GEBCO/BathyApp arrives positive up and is flipped here).

using Printf: @sprintf

# ---------------- CSV profiles (BathyApp exports) ----------------
"""
    env_profile_csv(path; positive_down = false) -> (x_m, depth_m, header)

Generic CSV profile reader. Accepts:
  * distance,x_m / dist_m,... + depth columns (values as-is)
  * lon,lat,... + depth/elevation columns (great-circle cumulative metres)
  * plain two-column numeric CSV (first = x [m], second = depth)
Converts positive-up elevation to positive-down depth unless positive_down.
"""
function env_profile_csv(path; positive_down = false)
    rows = Any[]
    open(path, "r") do io
        for line in eachline(io)
            s = strip(line)
            (s == "" || startswith(s, "#")) && continue
            parts = split(s, [',', ';', '\t']; keepempty = false)
            vals = Any[]
            for p in parts
                pc = replace(p, "\"" => "")
                v = tryparse(Float64, pc)
                push!(vals, v === nothing ? lowercase(pc) : v)
            end
            push!(rows, vals)
        end
    end
    isempty(rows) && error("csv: no data rows in $path")
    # is the first row a header (contains non-numeric tokens)?
    is_header = any(c -> c isa String, rows[1])
    head = is_header ? [c isa String ? c : "" for c in rows[1]] : ["" for _ in rows[1]]
    body = is_header ? rows[2:end] : rows
    body = [r for r in body if length(r) == length(head)]   # tolerate ragged lines
    isempty(body) && error("csv: only a header in $path")
    width = length(head)
    col = (needle) -> findfirst(h -> occursin(needle, h), head)
    ilon = col("lon"); ilat = col("lat")
    iz = something(col("depth"), col("elev"), col("z"))
    if ilon !== nothing && ilat !== nothing
        R = 6371000.0
        x = Float64[0.0]; last = (body[1][ilon], body[1][ilat])
        for r in body[2:end]
            push!(x, x[end] + hypot((r[ilon] - last[1]) * cosd((r[ilat] + last[2]) / 2),
                                    r[ilat] - last[2]) * pi / 180.0 * R)
            last = (r[ilon], r[ilat])
        end
        dcol = iz === nothing ? (width >= 3 ? 3 : width) : iz
        d = [r[dcol] for r in body]
        flip = !positive_down && !occursin("depth", head[dcol])
        return (x, flip ? [-v for v in d] : d, head)
    end
    idist = something(col("dist"), col("x_m"), col("range"))
    ix = idist === nothing ? 1 : idist
    idp = iz === nothing ? (width >= 2 ? 2 : width) : iz
    x = [r[ix] for r in body]
    d = [r[idp] for r in body]
    # header heuristic: a "depth" column is positive down already; an
    # "elevation" column is positive up (flip unless told otherwise)
    flip = !positive_down && !occursin("depth", head[idp])
    return (x, flip ? [-v for v in d] : d, head)
end

"`clip/smooth/interpolate profile onto the FDTD x grid; returns depth(x) positive down"
function env_sample_profile(xb, db, xq)
    length(xb) >= 2 || error("env: profile needs >= 2 points")
    # monotone increasing xb
    xs = collect(xb); ds = collect(db)
    if xs[1] > xs[end]
        xs = reverse(xs); ds = reverse(ds)
    end
    out = similar(xq, Float64)
    for (k, x) in enumerate(xq)
        i = clamp(searchsortedfirst(xs, x), 2, length(xs))
        t = (x - xs[i-1]) / (xs[i] - xs[i-1])
        t = clamp(t, 0.0, 1.0)
        out[k] = ds[i-1] * (1 - t) + ds[i] * t
    end
    return max.(out, 1.0)   # at least 1 m of water (clip land / tide)
end

function env_smooth(v::Vector{Float64}, w::Int)
    w <= 1 && return copy(v)
    n = length(v)
    out = similar(v)
    for i in 1:n
        s = 0.0; c = 0
        for j in max(1, i - w ÷ 2):min(n, i + w ÷ 2)
            s += v[j]; c += 1
        end
        out[i] = s / c
    end
    return out
end

# ---------------- BathyApp GeoTIFF -> profile ----------------
"""
    env_profile_geotiff(path; axis = :x, frac = 0.5, positive_down = false)
    -> (x_m::Vector, depth_m::Vector, note::String)

Reads a BathyApp GeoTIFF export and cuts the profile along the grid
x axis (constant lat = frac through the rows) or y axis (constant lon).
"""
function env_profile_geotiff(path; axis = :x, frac = 0.5, positive_down = false)
    g = gtif_read(path)
    if axis == :x
        r = clamp(Int(round((length(g.lat) - 1) * frac)) + 1, 1, length(g.lat))
        elev = g.z[r, :]
        lats = fill(g.lat[r], length(g.lon))
        lons = g.lon
    else
        c = clamp(Int(round((length(g.lon) - 1) * frac)) + 1, 1, length(g.lon))
        elev = g.z[:, c]
        lons = fill(g.lon[c], length(g.lat))
        lats = g.lat
    end
    xs = Float64[0.0]
    for k in 2:length(lats)
        d = hypot((lons[k] - lons[k-1]) * cosd((lats[k] + lats[k-1]) / 2),
                  lats[k] - lats[k-1]) * pi / 180.0 * 6371000.0
        push!(xs, xs[end] + d)
    end
    dvals = [isnan(v) ? 0.0 : (positive_down ? v : -v) for v in elev]
    note = @sprintf("%s, axis %s, slice %.2f (lat %.3f, lon %.3f)",
                    basename(path), axis, frac, lats[1], lons[1])
    return collect(xs), dvals, note
end

# ---------------- config generators (2-D) ----------------
"""
    env_config2d(; dx, dz, xmax, zmax, depthfn_x, pml_cells,
                  c_water, c_seabed (function of z) etc.) -> String

Generates a config/env_*.jl replacement (exact contract of
config/env_pekeris.jl). `prof_x`/`prof_d` define the bathymetric profile
(depth positive down); empty vectors give a flat seabed at `flat_depth`.
"""
function env_config2d(; dx, dz, xmax, zmax, pml_cells = 40,
                       prof_x = Float64[], prof_d = Float64[],
                       flat_depth = 100.0, smooth_cells = 0,
                       c_water = 1500.0, c_water_grad = 0.0,
                       rho_water = 1000.0, c_bottom = 1600.0,
                       rho_bottom = 1500.0, Qp_water = 1000.0, Qp_bottom = 100.0)
    buf = IOBuffer()
    p = println
    p(buf, "# ===========================================================================")
    p(buf, "# AUTO-GENERATED by Acoustic_Android (PlutoSuite notebook).")
    p(buf, "# Contract of config/env_pekeris.jl: dx, dz, xmax, zmax, water_depth,")
    p(buf, "# c_water, rho_water, c_bottom, rho_bottom, Qp_water, Qp_bottom,")
    p(buf, "# pml_cells, nx, nz, x, z, grid. Bathymetric profile embedded below.")
    p(buf, "# ===========================================================================")
    p(buf, "")
    p(buf, "const _AS_PROFILE_X = Float64[$(join(string.(prof_x), ", "))]")
    p(buf, "const _AS_PROFILE_D = Float64[$(join(string.(prof_d), ", "))]")
    p(buf, "const _AS_FLAT_DEPTH = $flat_depth")
    p(buf, "dx = $dx")
    p(buf, "dz = $dz")
    p(buf, "xmax = $xmax")
    p(buf, "zmax = $zmax")
    p(buf, "nx = Int(round(xmax / dx)) + 1")
    p(buf, "nz = Int(round(zmax / dz)) + 1")
    p(buf, "x = collect(0.0:dx:xmax)")
    p(buf, "z = collect(0.0:dz:zmax)")
    p(buf, "pml_cells = $pml_cells")
    nx = Int(round(xmax / dx)) + 1
    p(buf, "# water column sound speed (linear gradient, m/s per m) & layer params")
    p(buf, "c_water0 = $c_water")
    p(buf, "c_water_grad = $c_water_grad")
    p(buf, "const c_water = c_water0  # surface sound speed (run.jl meta + SL anchor contract)")
    p(buf, "rho_water = $rho_water")
    p(buf, "c_bottom = $c_bottom")
    p(buf, "rho_bottom = $rho_bottom")
    p(buf, "Qp_water = $Qp_water")
    p(buf, "Qp_bottom = $Qp_bottom")
    p(buf, "")
    p(buf, """_depth_at(xi) = isempty(_AS_PROFILE_D) ? _AS_FLAT_DEPTH :
    max(1.0, let
        v = 0.0
        for k in 2:length(_AS_PROFILE_X)
            if _AS_PROFILE_X[k-1] <= xi <= _AS_PROFILE_X[k]
                t = (xi - _AS_PROFILE_X[k-1]) / max(_AS_PROFILE_X[k] - _AS_PROFILE_X[k-1], eps())
                v = _AS_PROFILE_D[k-1] * (1 - t) + _AS_PROFILE_D[k] * t
                break
            end
        end
        v == 0.0 ? _AS_PROFILE_D[end] : v
    end)""")
    p(buf, "")
    p(buf, "lam = zeros(nx, nz); mu = zeros(nx, nz)")
    p(buf, "rho = zeros(nx, nz); Qp = fill(Qp_water, nx, nz); Qs = zeros(nx, nz)")
    p(buf, "lam_bottom = rho_bottom * c_bottom^2")
    p(buf, "for i in 1:nx")
    p(buf, "    d = _depth_at(x[i])  # seabed depth below surface, positive down")
    p(buf, "    for j in 1:nz")
    p(buf, "        if z[j] < d - dz")
    p(buf, "            cw = c_water0 + c_water_grad * z[j]")
    p(buf, "            rho[i,j] = rho_water;  lam[i,j] = rho_water * cw^2;  Qp[i,j] = Qp_water")
    p(buf, "        elseif z[j] > d + dz")
    p(buf, "            rho[i,j] = rho_bottom; lam[i,j] = lam_bottom;      Qp[i,j] = Qp_bottom")
    p(buf, "        else")
    p(buf, "            frac = (z[j] - (d - dz)) / (2*dz)")
    p(buf, "            frac = clamp(frac, 0.0, 1.0)")
    p(buf, "            cw = c_water0 + c_water_grad * min(z[j], d)")
    p(buf, "            rho[i,j] = rho_water*(1-frac)  + rho_bottom*frac")
    p(buf, "            lam[i,j] = rho_water*cw^2*(1-frac) + lam_bottom*frac")
    p(buf, "            Qp[i,j]  = Qp_water *(1-frac)  + Qp_bottom*frac")
    p(buf, "        end")
    p(buf, "    end")
    p(buf, "end")
    p(buf, "grid = ElasticGrid(nx, nz, dx, dz, lam, mu, rho, Qp, Qs)")
    dmean = isempty(prof_d) ? flat_depth : sum(prof_d) / length(prof_d)
    p(buf, "water_depth = $(dmean)  # mean profile depth (SL anchor is approximate on bathy runs)")
    p(buf, "println(\"Environment loaded (generated): \$(_depth_at(x[div(nx,2)])) m seabed at mid-range, \$(nx) x \$(nz) cells\")")
    return String(take!(buf))
end

"""
    env_sr2d(; f0, src_amp, source_x, source_z, receiver_ranges_km,
              receiver_depth, dx, dz, check_f_max) -> String
Generates config/sources_receivers.jl (exact contract; Hann pulse like
the suite's own Dickins/Pekeris benchmark).
"""
function env_sr2d(; f0, src_amp, source_x, source_z, receiver_ranges_km,
                  receiver_depth, dx, dz, check_f_max = nothing)
    cf = check_f_max === nothing ? 1.5 * f0 : check_f_max
    buf = IOBuffer()
    p = println
    p(buf, "# AUTO-GENERATED by Acoustic_Android. Contract of")
    p(buf, "# config/sources_receivers.jl: f0, source_kind, src_amp, source_x,")
    p(buf, "# source_z, check_f_max, sources, monopoles, bc, receiver_ranges_km,")
    p(buf, "# receiver_depth, receivers, SL_REF_RANGE, SL_REF_IS_RECEIVER.")
    p(buf, "f0          = $f0")
    p(buf, "source_kind = \"$(f0) Hz Hann pulse (Acoustic_Android)\"")
    p(buf, "src_amp     = $src_amp")
    p(buf, "source_x, source_z = $source_x, $source_z")
    p(buf, "check_f_max = $cf")
    p(buf, "ix = Int(round(source_x / dx)) + 1")
    p(buf, "iz = Int(round(source_z / dz)) + 1")
    p(buf, "T_pulse = 400.0 / f0")
    p(buf, "function pulse_pressure(t)")
    p(buf, "    (t > 0.0 && t <= T_pulse) || return 0.0")
    p(buf, "    w = 0.5 * (1.0 - cos(2.0*pi*t / T_pulse))")
    p(buf, "    return sin(2.0*pi*f0*t) * w")
    p(buf, "end")
    p(buf, "function pulse_rate(t)")
    p(buf, "    (t > 0.0 && t <= T_pulse) || return 0.0")
    p(buf, "    w  = 0.5 * (1.0 - cos(2.0*pi*t / T_pulse))")
    p(buf, "    dw = (pi / T_pulse) * sin(2.0*pi*t / T_pulse)")
    p(buf, "    return src_amp * (2.0*pi*f0*cos(2.0*pi*f0*t)*w + sin(2.0*pi*f0*t)*dw)")
    p(buf, "end")
    p(buf, "sources   = Source[]")
    p(buf, "monopoles = MonopoleSource[MonopoleSource(ix, iz, pulse_rate)]")
    p(buf, "bc = BC(:pml, :pml, :free, :pml)")
    p(buf, "receiver_ranges_km = [" * join(_numlit.(receiver_ranges_km), ", ") * "]")
    p(buf, "receiver_depth = $receiver_depth")
    p(buf, """receivers = Receiver[
    Receiver(Int(round(r*1000 / dx)) + 1, Int(round(receiver_depth / dz)) + 1)
    for r in receiver_ranges_km if Int(round(r*1000)) < xmax - pml_cells*dx]""")
    p(buf, "SL_REF_RANGE = 20.0")
    p(buf, "sl_ref_ix = Int(round((source_x + SL_REF_RANGE) / dx)) + 1")
    p(buf, "push!(receivers, Receiver(sl_ref_ix, iz))")
    p(buf, "SL_REF_IS_RECEIVER = length(receivers)")
    return String(take!(buf))
end

_numlit(v) = repr(Float64(v))

# ---------------- config generators (3-D) ----------------
function env_config3d(; dx, dy, dz, xmax, ymax, zmax, pml_cells = 20,
                       prof_x = Float64[], prof_d = Float64[],
                       flat_depth = 100.0, c_water = 1500.0, c_water_grad = 0.0,
                       rho_water = 1000.0, c_bottom = 1600.0,
                       rho_bottom = 1500.0, Qp_water = 1000.0, Qp_bottom = 100.0)
    buf = IOBuffer()
    p = println
    p(buf, "# AUTO-GENERATED by Acoustic_Android. Contract of config/env3d.jl:")
    p(buf, "# dx, dy, dz, xmax, ymax, zmax, water_depth, c_water, rho_water,")
    p(buf, "# c_bottom, rho_bottom, Qp_water, Qp_bottom, pml_cells, nx, ny, nz,")
    p(buf, "# x, y, z, grid3d. Profile extends uniformly along y.")
    p(buf, "const _AS_PROFILE_X = Float64[$(join(string.(prof_x), ", "))]")
    p(buf, "const _AS_PROFILE_D = Float64[$(join(string.(prof_d), ", "))]")
    p(buf, "const _AS_FLAT_DEPTH = $flat_depth")
    p(buf, "dx = $dx; dy = $dy; dz = $dz")
    p(buf, "xmax = $xmax; ymax = $ymax; zmax = $zmax")
    p(buf, "nx = Int(round(xmax / dx)) + 1")
    p(buf, "ny = Int(round(ymax / dy)) + 1")
    p(buf, "nz = Int(round(zmax / dz)) + 1")
    p(buf, "x = collect(0.0:dx:xmax); y = collect(0.0:dy:ymax); z = collect(0.0:dz:zmax)")
    p(buf, "pml_cells = $pml_cells")
    p(buf, "c_water0 = $c_water; c_water_grad = $c_water_grad")
    p(buf, "const c_water = c_water0  # surface sound speed (run3d.jl meta contract)")
    p(buf, "rho_water = $rho_water; c_bottom = $c_bottom; rho_bottom = $rho_bottom")
    p(buf, "Qp_water = $Qp_water; Qp_bottom = $Qp_bottom")
    p(buf, """_depth_at(xi) = isempty(_AS_PROFILE_D) ? _AS_FLAT_DEPTH :
    max(1.0, let
        v = 0.0
        for k in 2:length(_AS_PROFILE_X)
            if _AS_PROFILE_X[k-1] <= xi <= _AS_PROFILE_X[k]
                t = (xi - _AS_PROFILE_X[k-1]) / max(_AS_PROFILE_X[k] - _AS_PROFILE_X[k-1], eps())
                v = _AS_PROFILE_D[k-1] * (1 - t) + _AS_PROFILE_D[k] * t
                break
            end
        end
        v == 0.0 ? _AS_PROFILE_D[end] : v
    end)""")
    p(buf, "lam = zeros(nx, ny, nz); mu = zeros(nx, ny, nz)")
    p(buf, "rho = zeros(nx, ny, nz); Qp = fill(Qp_water, nx, ny, nz); Qs = zeros(nx, ny, nz)")
    p(buf, "lam_bottom = rho_bottom * c_bottom^2")
    p(buf, "for i in 1:nx")
    p(buf, "    d = _depth_at(x[i])")
    p(buf, "    for j in 1:ny, k in 1:nz")
    p(buf, "        if z[k] < d - dz")
    p(buf, "            cw = c_water0 + c_water_grad * z[k]")
    p(buf, "            rho[i,j,k] = rho_water;  lam[i,j,k] = rho_water * cw^2;  Qp[i,j,k] = Qp_water")
    p(buf, "        elseif z[k] > d + dz")
    p(buf, "            rho[i,j,k] = rho_bottom; lam[i,j,k] = lam_bottom;      Qp[i,j,k] = Qp_bottom")
    p(buf, "        else")
    p(buf, "            frac = clamp((z[k] - (d - dz)) / (2*dz), 0.0, 1.0)")
    p(buf, "            cw = c_water0 + c_water_grad * min(z[k], d)")
    p(buf, "            rho[i,j,k] = rho_water*(1-frac)  + rho_bottom*frac")
    p(buf, "            lam[i,j,k] = rho_water*cw^2*(1-frac) + lam_bottom*frac")
    p(buf, "            Qp[i,j,k]  = Qp_water *(1-frac)  + Qp_bottom*frac")
    p(buf, "        end")
    p(buf, "    end")
    p(buf, "end")
    p(buf, "grid3d = ElasticGrid3D(nx, ny, nz, dx, dy, dz, lam, mu, rho, Qp, Qs)")
    dmean = isempty(prof_d) ? flat_depth : sum(prof_d) / length(prof_d)
    p(buf, "water_depth = $(dmean)")
    return String(take!(buf))
end

function env_sr3d(; f0, src_amp, source_x, source_y, source_z,
                  receiver_ranges_km, receiver_depth, dx, dy, dz, ymax,
                  check_f_max = nothing)
    cf = check_f_max === nothing ? 1.5 * f0 : check_f_max
    buf = IOBuffer()
    p = println
    p(buf, "# AUTO-GENERATED by Acoustic_Android. Contract of")
    p(buf, "# config/sources_receivers3d.jl")
    p(buf, "f0          = $f0")
    p(buf, "source_kind = \"$(f0) Hz Hann pulse (Acoustic_Android 3D)\"")
    p(buf, "src_amp     = $src_amp")
    p(buf, "source_x, source_y, source_z = $source_x, $source_y, $source_z")
    p(buf, "check_f_max = $cf")
    p(buf, "ix = Int(round(source_x / dx)) + 1")
    p(buf, "iy = Int(round(source_y / dy)) + 1")
    p(buf, "iz = Int(round(source_z / dz)) + 1")
    p(buf, "T_pulse = 400.0 / f0")
    p(buf, "function pulse_rate3d(t)")
    p(buf, "    (t > 0.0 && t <= T_pulse) || return 0.0")
    p(buf, "    w  = 0.5 * (1.0 - cos(2.0*pi*t / T_pulse))")
    p(buf, "    dw = (pi / T_pulse) * sin(2.0*pi*t / T_pulse)")
    p(buf, "    return src_amp * (2.0*pi*f0*cos(2.0*pi*f0*t)*w + sin(2.0*pi*f0*t)*dw)")
    p(buf, "end")
    p(buf, "sources3d   = Source3D[]")
    p(buf, "monopoles3d = MonopoleSource3D[MonopoleSource3D(ix, iy, iz, pulse_rate3d)]")
    p(buf, "bc3d = BC3D(:pml, :pml, :pml, :pml, :free, :pml)")
    p(buf, "receiver_ranges_km = [" * join(_numlit.(receiver_ranges_km), ", ") * "]")
    p(buf, "receiver_depth = $receiver_depth")
    p(buf, """receivers3d = Receiver3D[
    Receiver3D(Int(round(r*1000 / dx)) + 1, iy, Int(round(receiver_depth / dz)) + 1)
    for r in receiver_ranges_km if Int(round(r*1000)) < xmax - pml_cells*dx]""")
    p(buf, "SL_REF_RANGE3D = 20.0")
    p(buf, "sl_ref_ix = Int(round((source_x + SL_REF_RANGE3D) / dx)) + 1")
    p(buf, "push!(receivers3d, Receiver3D(sl_ref_ix, iy, iz))")
    p(buf, "SL_REF_IS_RECEIVER3D = length(receivers3d)")
    return String(take!(buf))
end
