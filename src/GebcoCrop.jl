# GebcoCrop.jl — GEBCO download-service crop API client (Julia, stdlib only).
# Mirrors BathyApp's Kotlin client (BathyApp :core gebco/GebcoCrop.kt):
#   POST /api/queue   -> basketId   (JSON body exactly as the web client)
#   GET  /api/queue/status/{id}    -> "finished" poll
#   GET  /api/queue/download/{id}  -> ZIP containing the NetCDF (CDF-1,
#                                     metres, POSITIVE UP)
# HTTP via the system `curl`; ZIP extraction via system `unzip`.
#
#   grid = gebco_crop(3.0, 53.0, 6.0, 55.0; cachedir = pwd())
#   grid.lon, grid.lat, grid.z      (z positive up, metres; land = positive)

using Printf: @sprintf

const GEBCO_BASE = "https://download.gebco.net/api"
const GEBCO_MAX_AREA_DEG2 = 14400.0

struct GebcoGrid
    lon::Vector{Float64}
    lat::Vector{Float64}
    z::Matrix{Float64}   # [lat, lon], positive up
    src::String
end

function _curl_str(cmd::Cmd)
    out = Pipe()
    proc = run(pipeline(cmd; stdout = out, stderr = out), wait = false)
    close(out.in)
    s = read(out, String)
    close(out)
    wait(proc)
    proc.exitcode == 0 || error("curl failed ($(proc.exitcode)): $(first(s, 300))")
    return s
end

function gebco_body(west, south, east, north)
    (east > west && north > south) || error("gebco: invalid bounding box")
    area = (east - west) * (north - south)
    area <= GEBCO_MAX_AREA_DEG2 ||
        error(@sprintf("gebco: bbox %.1f deg^2 exceeds the %.0f deg^2 server limit", area, GEBCO_MAX_AREA_DEG2))
    return @sprintf("{\"items\":[{\"id\":0,\"grid_id\":%d,\"data_source_ids\":[%d],\"formats\":[%d],\"left\":%.6f,\"right\":%.6f,\"top\":%.6f,\"bottom\":%.6f}]}",
                    1, 1, 1, west, east, north, south)
end

"""
    gebco_crop(west, south, east, north; cachedir, poll_s, timeout_s, onstatus)
    -> GebcoGrid

Downloads a GEBCO 2026 bathymetry crop as NetCDF (cache hit if the file
already exists in `cachedir`). Returns a grid with z positive up.
"""
function gebco_crop(west, south, east, north;
                    cachedir = pwd(), poll_s = 3.0, timeout_s = 300.0,
                    onstatus = println)
    mkpath(cachedir)
    key = @sprintf("gebco_%.4f_%.4f_%.4f_%.4f.nc", west, south, east, north)
    ncp = joinpath(cachedir, key)
    if isfile(ncp)
        onstatus("GEBCO: cache hit $(basename(ncp))")
        return gebco_load_nc(ncp)
    end
    body = gebco_body(west, south, east, north)
    onstatus("GEBCO: submitting crop request...")
    resp = _curl_str(`curl -sS --max-time 60 -X POST -H "Content-Type: application/json" --data $body $GEBCO_BASE/queue`)
    m = match(r"\"basketId\"\s*:\s*\"([^\"]+)\"", resp)
    m === nothing && error("gebco: no basketId in response: $(first(resp, 300))")
    basket = m.captures[1]
    onstatus("GEBCO: basket $basket submitted")
    t0 = time()
    while true
        st = _curl_str(`curl -sS --max-time 60 $GEBCO_BASE/queue/status/$basket`)
        ms = match(r"\"status\"\s*:\s*\"([^\"]*)\"", st)
        status = ms === nothing ? "?" : ms.captures[1]
        onstatus("GEBCO: server status: $status")
        lowercase(status) == "finished" && break
        lowercase(status) in ("error", "failed") && error("gebco: basket $basket failed ($status)")
        time() - t0 > timeout_s && error("gebco: timeout waiting for basket $basket")
        sleep(poll_s)
    end
    zipp = ncp * ".zip"
    run(`curl -sS --fail --max-time 600 -o $zipp $GEBCO_BASE/queue/download/$basket`)
    isfile(zipp) || error("gebco: download produced no file")
    # extract the .nc entry (system unzip)
    io = Pipe()
    proc = run(pipeline(`unzip -o -j $zipp "*.nc" -d $cachedir`; stdout = io, stderr = devnull), wait = false)
    close(io.in)
    _out = read(io, String); close(io); wait(proc)
    proc.exitcode == 0 || error("gebco: unzip failed — is `unzip` installed? (Debian: apt-get install unzip)")
    found = filter(readdir(cachedir)) do f
        endswith(lowercase(f), ".nc") && f != key && mtime(joinpath(cachedir, f)) > time() - 120
    end
    isempty(found) && error("gebco: no NetCDF extracted from the basket ZIP")
    mv(joinpath(cachedir, found[1]), ncp; force = true)
    rm(zipp; force = true)
    onstatus("GEBCO: NetCDF ready -> $(basename(ncp))")
    return gebco_load_nc(ncp)
end

function gebco_load_nc(path)
    nc = nc3_read(path, "elevation")
    lon = collect(nc["lon"].data[:, 1])
    lat = collect(nc["lat"].data[:, 1])
    z = nc["elevation"].data
    if size(z, 1) == length(lat) && size(z, 2) == length(lon)
        return GebcoGrid(lon, lat, z, basename(path))
    elseif size(z, 1) == length(lon) && size(z, 2) == length(lat)
        return GebcoGrid(lon, lat, permutedims(z), basename(path))
    end
    error("gebco: dims ($(size(z))) do not match lon/lat ($(length(lon)) x $(length(lat)))")
end

"""
    gebco_profile(g, lat0, lon0, lat1, lon1; n = 512)
    -> (x_m, elev_up_m, track)  with track = [(lat, lon), ...]

Great-circle profile through the crop; x in metres from (lat0, lon0).
"""
function gebco_profile(g::GebcoGrid, lat0, lon0, lat1, lon1; n = 512)
    R = 6371000.0
    d = hypot((lon1 - lon0) * cosd((lat0 + lat1) / 2), lat1 - lat0) * pi / 180.0 * R
    ts = range(0.0, 1.0; length = n)
    xs = collect(range(0.0, d; length = n))
    lats = [lat0 + (lat1 - lat0) * t for t in ts]
    lons = [lon0 + (lon1 - lon0) * t for t in ts]
    elev = [env_bilinear(g.lat, g.lon, g.z, la, lo) for (la, lo) in zip(lats, lons)]
    return xs, elev, collect(zip(lats, lons))
end

function env_bilinear(ys, xs, A, y, x)
    x1 = clamp(searchsortedfirst(xs, x), 2, length(xs))
    y1 = clamp(searchsortedfirst(ys, y), 2, length(ys))
    x0, y0 = x1 - 1, y1 - 1
    fx = (x - xs[x0]) / (xs[x1] - xs[x0])
    fy = (y - ys[y0]) / (ys[y1] - ys[y0])
    v = A[y0,x0]*(1-fx)*(1-fy) + A[y0,x1]*fx*(1-fy) + A[y1,x0]*(1-fx)*fy + A[y1,x1]*fx*fy
    return isnan(v) ? 0.0 : v
end
