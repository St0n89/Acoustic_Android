# test/selftest.jl — offline self-check: synthetic NetCDF/TIFF round-trips,
# config-generator contract checks, PNG encoder sanity.
#   julia test/selftest.jl
# Exit nonzero on failure. Uses only stdlib (the readers get fully exercised
# without any suite clone or network access).

const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "src", "AcousticCore.jl"))

passed, failed = 0, 0
check(name, cond) = begin
    if cond
        global passed += 1
        println("  PASS  $name")
    else
        global failed += 1
        println("  FAIL  $name")
    end
end

tmpdir = mktempdir()
println("scratch: $tmpdir")

# ---------- synthetic NetCDF CDF-1 (GEBCO-style) ----------
function write_be(x::UInt32, io)
    for q in (24, 16, 8, 0)
        write(io, UInt8((x >> q) & 0xFF))
    end
end
function write_be16(x::UInt16, io)
    write(io, UInt8((x >> 8) & 0xFF)); write(io, UInt8(x & 0xFF))
end
function nc_name(s::String, io)
    b = codeunits(s)
    write_be(UInt32(length(b)), io)
    write(io, b)
    write(io, fill(UInt8(0), (4 - length(b) % 4) % 4))
end
# Simplify: use a practical builder that records offsets in a second pass.
function build_cdf2(path, lats, lons, store::Matrix{Int16})
    io = open(path, "w")
    println(io, "")  # placeholder, replaced below
    close(io)
    io = open(path, "w")
    write(io, UInt8[0x43, 0x44, 0x46, 0x01])
    write_be(UInt32(0), io)
    write_be(UInt32(10), io); write_be(UInt32(2), io)
    nc_name("lat", io); write_be(UInt32(length(lats)), io)
    nc_name("lon", io); write_be(UInt32(length(lons)), io)
    write_be(UInt32(12), io); write_be(UInt32(0), io)
    write_be(UInt32(11), io); write_be(UInt32(3), io)   # var_list: 3 vars
    var_slots = Dict{String,Int}()
    function varhdr(name, typ, dimids)
        nc_name(name, io); write_be(UInt32(typ), io)
        write_be(UInt32(length(dimids)), io)
        for d in dimids
            write_be(UInt32(d), io)
        end
    end
    # lon
    varhdr("lon", UInt32(6), [1])
    write_be(UInt32(12), io); write_be(UInt32(0), io)
    write_be(UInt32(length(lons) * 8), io)
    var_slots["lon"] = position(io); write_be(UInt32(0), io)
    # lat
    varhdr("lat", UInt32(6), [0])
    write_be(UInt32(12), io); write_be(UInt32(0), io)
    write_be(UInt32(length(lats) * 8), io)
    var_slots["lat"] = position(io); write_be(UInt32(0), io)
    # elevation
    varhdr("elevation", UInt32(3), [0, 1])
    # attrs
    write_be(UInt32(12), io); write_be(UInt32(2), io)
    nc_name("scale_factor", io)
    write_be(UInt32(5), io); write_be(UInt32(1), io)
    write(io, reinterpret(UInt8, [hton(Float32(0.01))]))
    nc_name("_FillValue", io)
    write_be(UInt32(3), io); write_be(UInt32(1), io)
    write_be16(UInt16(0x8000), io); write(io, [UInt8(0), UInt8(0)])  # -32768 + pad
    write_be(UInt32(4096), io)                # vsize (padded)
    var_slots["elevation"] = position(io); write_be(UInt32(0), io)
    # data
    flush(io)
    data_start = position(io)
    lon_off = data_start
    for v in lons
        write(io, hton(Float64(v)))
    end
    lat_off = position(io)
    for v in lats
        write(io, hton(Float64(v)))
    end
    el_off = position(io)
    h, w = size(store)   # h = lat, w = lon
    for i in 1:h, j in 1:w
        write_be16(reinterpret(UInt16, store[i, j]), io)
    end
    pad = (4 - (h * w * 2) % 4) % 4
    write(io, fill(UInt8(0), pad))
    close(io)
    # patch offsets
    b = read(path)
    patch(o, v) = (b[o+1] = UInt8((v >> 24) & 0xFF); b[o+2] = UInt8((v >> 16) & 0xFF);
                   b[o+3] = UInt8((v >> 8) & 0xFF); b[o+4] = UInt8(v & 0xFF))
    patch(var_slots["lon"], lon_off)
    patch(var_slots["lat"], lat_off)
    patch(var_slots["elevation"], el_off)
    open(path, "w") do io
        write(io, b)
    end
    return path
end

println("\n[1] NetCDF CDF-1 reader")
lats = [54.0, 53.5, 53.0]
lons = [7.0, 7.5, 8.0, 8.5]
store = Int16[(10 * i + 100 * j) for i in 1:3, j in 1:4]  # raw counts (0.01 m each)
ncp = joinpath(tmpdir, "toy.nc")
build_cdf2(ncp, lats, lons, store)
nc3p_dbg = read(ncp)
nc = nc3_read(ncp, "elevation")
elev = nc["elevation"].data
check("dims vars present", haskey(nc.vars, "lon") && haskey(nc.vars, "lat"))
check("lon roundtrip", nc["lon"].data[:, 1] ≈ lons)
check("lat roundtrip", nc["lat"].data[:, 1] ≈ lats)
check("elevation shape (lon,lat)", size(elev) == (length(lons), length(lats)))
println("  elev[2,3]=", elev[2,3], " expected=", store[3,2]*0.01)
check("elevation scaled roundtrip", isapprox(elev[2, 3], store[3, 2] * 0.01; rtol = 1e-6))  # lon2, lat3 = store[3,2]

# ---------- synthetic GeoTIFF ----------
println("\n[2] GeoTIFF reader")
function build_tiff(path, h, w, lat0, lon0, dy, dxp)
    data = [Float32(-100.0 - 10.0 * ((i - 1) * w + (j - 1))) for i in 1:h, j in 1:w]
    nb = data .* 4
    io = open(path, "w")
    write(io, UInt8[0x49, 0x49, 0x2A, 0x00])   # II, 42
    write(io, reinterpret(UInt8, [UInt32(8)]))  # IFD at 8
    ifd_pos = 8
    ifd_size = 2 + 10 * 12 + 4
    data_off = ifd_pos + ifd_size
    extra_off = data_off + h * w * 4
    pxs_off = extra_off
    tie_off = pxs_off + 3 * 8
    local function entry(t, typ, cnt, val_or_off)
        write(io, reinterpret(UInt8, [UInt16(t), UInt16(typ)]))
        write(io, reinterpret(UInt8, [UInt32(cnt)]))
        if typ == 3 && cnt == 1
            write(io, reinterpret(UInt8, [UInt16(val_or_off), UInt16(0)]))
        else
            write(io, reinterpret(UInt8, [UInt32(val_or_off)]))
        end
    end
    write(io, reinterpret(UInt8, [UInt16(10)]))  # 10 entries
    entry(256, 4, 1, w)      # width
    entry(257, 4, 1, h)      # height
    entry(258, 3, 1, 32)     # bits per sample
    entry(259, 3, 1, 1)      # compression: none
    entry(273, 4, 1, data_off)   # strip offsets
    entry(277, 3, 1, 1)      # spp
    entry(278, 4, 1, h)      # rows per strip
    entry(279, 4, 1, h * w * 4)# strip byte counts
    entry(33550, 12, 3, pxs_off)
    entry(42112, 12, 6, tie_off)
    write(io, reinterpret(UInt8, [UInt32(0)]))   # next IFD
    @assert position(io) == data_off
    write(io, reinterpret(UInt8, vec(data')))
    @assert position(io) == pxs_off
    for q in (dxp, dy, 0.0)
        write(io, reinterpret(UInt8, [Float64(q)]))
    end
    for q in (0.0, 0.0, 0.0, lon0, lat0, 0.0)
        write(io, reinterpret(UInt8, [Float64(q)]))
    end
    close(io)
    return path
end
tifp = joinpath(tmpdir, "toy.tif")
build_tiff(tifp, 3, 4, 54.5, 7.0, 0.5, 0.25)
g = gtif_read(tifp)
check("tiff dims", size(g.z) == (3, 4))
check("tiff lon", g.lon ≈ [7.0, 7.25, 7.5, 7.75])
check("tiff lat (north-up)", g.lat ≈ [54.5, 54.0, 53.5])
check("tiff z[1,1]", g.z[1, 1] == -100.0)

# ---------- profile loading ----------
println("\n[3] Profiles + config generation")
csvp = joinpath(tmpdir, "prof.csv")
open(csvp, "w") do io
    println(io, "x_m,depth_m")
    println(io, "0,120.0")
    println(io, "1000,110.0")
    println(io, "2000,140.0")
end
x, d, h = env_profile_csv(csvp)
check("csv x", x == [0.0, 1000.0, 2000.0])
check("csv depth (no sign flip on depth cols)", d == [120.0, 110.0, 140.0])
xq = collect(0.0:100.0:2000.0)
dq = env_sample_profile(x, d, xq)
check("profile interp", dq[21] == 140.0 && dq[6] == 115.0)

cfg2d = env_config2d(; dx = 3.0, dz = 2.0, xmax = 1000.0, zmax = 200.0,
    pml_cells = 40, prof_x = x, prof_d = env_smooth(d, 2), flat_depth = 120.0)
for needle in ["grid = ElasticGrid(nx, nz, dx, dz, lam, mu, rho, Qp, Qs)",
               "water_depth =", "pml_cells = 40", "nx = Int(round(xmax / dx)) + 1"]
    check("config2d contains: $(needle[1:min(end,40)])", occursin(needle, cfg2d))
end
sr2d = env_sr2d(; f0 = 50.0, src_amp = 1.0e3, source_x = 150.0, source_z = 25.0,
    receiver_ranges_km = [0.5, 0.9], receiver_depth = 50.0, dx = 3.0, dz = 2.0)
for needle in ["monopoles = MonopoleSource[MonopoleSource(ix, iz, pulse_rate)]",
               "bc = BC(:pml, :pml, :free, :pml)", "SL_REF_IS_RECEIVER"]
    check("sr2d contains: $(needle[1:min(end,40)])", occursin(needle, sr2d))
end
cfg3d = env_config3d(; dx = 3.0, dy = 3.0, dz = 2.0, xmax = 600.0, ymax = 300.0,
    zmax = 200.0, pml_cells = 20, flat_depth = 100.0)
check("config3d grid3d", occursin("grid3d = ElasticGrid3D(nx, ny, nz, dx, dy, dz, lam, mu, rho, Qp, Qs)", cfg3d))

# ---------- PNG encoder ----------
println("\n[4] QuickPlot PNG")
img = fill(UInt8(200), 20, 30, 3)
for i in 1:20, j in 1:30
    img[i, j, 1] = UInt8(i * 10); img[i, j, 2] = UInt8(j * 8); img[i, j, 3] = UInt8(7)
end
png = qp_png(img)
check("png signature", png[1:8] == UInt8[0x89, 0x50, 0x4E, 0x47, 0x0D, 0x0A, 0x1A, 0x0A])
check("png IHDR", png[13:16] == UInt8[0x49, 0x48, 0x44, 0x52])
w  = reinterpret(UInt32, png[17:20])[1] |> ntoh
h2 = reinterpret(UInt32, png[21:24])[1] |> ntoh
check("png dims field", w == 30 && h2 == 20)
uri = qp_uri(png)
check("png data uri", startswith(uri, "data:image/png;base64,"))
mat = [Float64(sin(i / 10) * cos(j / 20)) for i in 1:60, j in 1:80]
hmuri, lims = qp_heatmap_png(mat)
check("heatmap uri", startswith(hmuri, "data:image/png;base64,"))

# ---------- meta parse ----------
println("\n[5] meta.jl reader")
metap = joinpath(tmpdir, "meta.jl")
open(metap, "w") do io
    println(io, "# AUTO")
    println(io, "out_dir = \"pekeris_run\"")
    println(io, "dx = 3.0")
    println(io, "nx = 334")
    println(io, "nframes = 51")
    println(io, "receiver_ranges_km = [1.0, 2.0]")
end
m = qp_read_meta(metap)
check("meta out_dir", m["out_dir"] == "pekeris_run")
check("meta nx int", m["nx"] == 334)
check("meta dx float", m["dx"] == 3.0)
check("meta recs", m["receiver_ranges_km"] == [1.0, 2.0])

println("\n==== $passed PASS, $failed FAIL ====")
failed == 0 || exit(1)
