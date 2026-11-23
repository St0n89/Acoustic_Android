### A Pluto.jl notebook ###
# v0.20.10
# editable in a browser at http://127.0.0.1:8890 — see fdtdapp.jl

using Markdown
using InteractiveUtils

# ╔═╡ a0000001-0000-4000-8000-000000000001
md"""
# 🌊 Acoustic_Android — Underwater Acoustics Suite on your phone

State-of-the-art trigger app for the
[St0n89/FDTD2DElastic](https://github.com/St0n89/FDTD2DElastic) suite:
2-D / 3-D viscoelastic FDTD, coupled modes (Kraken), PE, Bellhop, FWI.

**Environment** comes from [BathyApp](https://github.com/St0n89/BathyApp):
* a BathyApp GeoTIFF or CSV export from your /workdir, or
* a fresh GEBCO crop straight from download.gebco.net, or
* a flat textbook waveguide.

**Runs are triggered as subprocesses** (`run.jl`, `run3d.jl`, `fwi.jl`,
`toy/xvalidate.jl`) — the suite itself stays unmodified; only its
`config/*.jl` files are written (originals backed up once per clone).

Sections:
**1 · Environment** → **2 · Suite config** → **3 · Trigger** → **4 · Results** → **5 · More solvers**.
"""

# ╔═╡ a0000002-0000-4000-8000-000000000002
begin
	APP = get(ENV, "ACOUSTIC_ANDROID_HOME", "") != "" ?
		  ENV["ACOUSTIC_ANDROID_HOME"] : normpath(joinpath(@__DIR__, ".."))
	if isfile(joinpath(APP, "src", "AcousticCore.jl"))
		Base.include(Main, joinpath(APP, "src", "AcousticCore.jl"))
	else
		@warn "src/ not found under $APP — start via fdtdapp.jl" APP
	end
end

# ╔═╡ a0000003-0000-4000-8000-000000000003
md"## 1 · Environment source"

# ╔═╡ a0000004-0000-4000-8000-000000000004
# PlutoUI is installed by fdtdapp.jl
import PlutoUI

# ╔═╡ a0000005-0000-4000-8000-000000000005
@bind env_source Select([
		"flat" => "Flat waveguide (textbook)",
		"geotiff" => "BathyApp GeoTIFF export",
		"csv" => "BathyApp CSV profile",
		"gebco" => "GEBCO crop (download here)",
	]; default = "flat")

# ╔═╡ a0000006-0000-4000-8000-000000000006
begin
	md"""**Flat waveguide:** water depth and seabed
		
	$(@bind flat_depth NumberField(10.0:10.0:2000.0; default = 100.0)) m depth
	"""
end

# ╔═╡ a0000007-0000-4000-8000-000000000007
begin
	if env_source == "geotiff"
		md"""**BathyApp GeoTIFF export** (path inside your workdir or absolute)
		
		$(@bind tif_path TextField((60, 1); placeholder = "gebco_export.tif / profile_case.tif", default = ""))
		
		Cut the FDTD profile along the image
		$(@bind tif_axis Select([:x => "x axis (along lon)", :y => "y axis (along lat)"]; default = :x))
		 through slice fraction
		$(@bind tif_frac Slider(0.0:0.05:1.0; default = 0.5, show_value = true))
		"""
	else
		md"_Select GeoTIFF to see options._"
	end
end

# ╔═╡ a0000008-0000-4000-8000-000000000008
begin
	if env_source == "csv"
		md"""**BathyApp CSV profile** (columns: `x_m|dist,...,depth|elev,...` or `lon,lat,...,depth,...`)
		
		$(@bind csv_path TextField((60, 1); placeholder = "my_profile.csv", default = ""))
		
		$(@bind csv_posdown CheckBox(default = false)) values are positive **down** (untick = elevation, positive up)
		"""
	else
		md"_Select CSV to see options._"
	end
end

# ╔═╡ a0000009-0000-4000-8000-000000000009
begin
	if env_source == "gebco"
		md"""**GEBCO crop** — bounding box (decimal degrees, °N/°E positive)
		
		west $(@bind g_west NumberField(-180.0:0.1:180.0; default = 7.0))
		east $(@bind g_east NumberField(-180.0:0.1:180.0; default = 8.0))
		south $(@bind g_south NumberField(-90.0:0.1:90.0; default = 53.5))
		north $(@bind g_north NumberField(-90.0:0.1:90.0; default = 54.5))
		
		Profile through the crop: from (lat, lon) $(@bind p_lat0 NumberField(-90.0:0.1:90.0; default = 54.0)) , $(@bind p_lon0 NumberField(-180.0:0.1:180.0; default = 7.3))
		to (lat, lon) $(@bind p_lat1 NumberField(-90.0:0.1:90.0; default = 54.0)) , $(@bind p_lon1 NumberField(-180.0:0.1:180.0; default = 7.9))
		"""
	else
		md"_Select GEBCO crop to see options._"
	end
end

# ╔═╡ a000000a-0000-4000-8000-00000000000a
begin
	md"""$(@bind load_now Button("Load / refresh profile"; primary = true))"""
end

# ╔═╡ a000000b-0000-4000-8000-00000000000b
begin
	load_now
	PROFILE = try
		if env_source == "flat"
			(Float64[0.0, 1.0], Float64[flat_depth, flat_depth], "flat $(flat_depth) m")
		elseif env_source == "geotiff"
			@assert strip(tif_path) != "" "enter a GeoTIFF path"
			p = isfile(tif_path) ? tif_path : joinpath(get(ENV, "ACOUSTIC_RUNS", pwd()), tif_path)
			env_profile_geotiff(p; axis = tif_axis, frac = tif_frac)
		elseif env_source == "csv"
			@assert strip(csv_path) != "" "enter a CSV path"
			p = isfile(csv_path) ? csv_path : joinpath(get(ENV, "ACOUSTIC_RUNS", pwd()), csv_path)
			x, d, h = env_profile_csv(p; positive_down = csv_posdown)
			(x, d, "CSV: $(join(filter(!isempty, h), ", "))")
		else
			cachedir = joinpath(get(ENV, "ACOUSTIC_RUNS", pwd()), "data")
			grid = gebco_crop(g_west, g_south, g_east, g_north; cachedir = cachedir)
			x, e, track = gebco_profile(grid, p_lat0, p_lon0, p_lat1, p_lon1; n = 512)
			(x, -e, "GEBCO $(basename(grid.src)), $(round(Int, x[end]/1000)) km track")
		end
	catch err
		md"⚠️ **$(sprint(showerror, err))**"
	end
end

# ╔═╡ a000000c-0000-4000-8000-00000000000c
if PROFILE isa Tuple
	md"""**Profile loaded** — $(PROFILE[3])
	
	$(qp_profile_svg(PROFILE[1], PROFILE[2]; title = "seabed depth (positive down)"))
	"""
end

# ╔═╡ a000000d-0000-4000-8000-00000000000d
md"""## 2 · Suite config (2-D FDTD)

Grid + physics + acquisition. The generated files follow the exact
`config/env_pekeris.jl` / `config/sources_receivers.jl` contract.
"""

# ╔═╡ a000000e-0000-4000-8000-00000000000e
begin
	md"""**Grid:** dx $(@bind dx_nb NumberField(0.5:0.5:100.0; default = 3.0)) m
	| xmax $(@bind xmax_nb NumberField(500.0:500.0:60000.0; default = 5000.0)) m
	| zmax water+bottom $(@bind zmax_nb NumberField(50.0:50.0:2000.0; default = 250.0)) m
	| smoothing $(@bind smooth_nb Slider(0:2:20; default = 4, show_value = true)) cells
	"""
end

# ╔═╡ a000000f-0000-4000-8000-00000000000f
begin
	md"""**Water / seabed:** c at surface $(@bind c_w NumberField(1400.0:1.0:1600.0; default = 1500.0)) m/s
	| gradient $(@bind c_g NumberField(-0.05:0.01:0.05; default = 0.0)) 1/s
	| c seabed $(@bind c_b NumberField(1450.0:10.0:4000.0; default = 1600.0)) m/s
	"""
end

# ╔═╡ a00000010-0000-4000-8000-0000000010
begin
	md"""**Acquisition (2-D):** f0 $(@bind f0_nb NumberField(5.0:5.0:1000.0; default = 50.0)) Hz
	| source (x,z) $(@bind sx_nb NumberField(50.0:50.0:1000.0; default = 150.0)) , $(@bind sz_nb NumberField(5.0:5.0:500.0; default = 25.0)) m
	| receivers [km, comma] $(@bind recs_tf TextField((40, 1); default = "1.0, 2.0, 3.0"))
	| receiver depth $(@bind rdep_nb NumberField(5.0:5.0:500.0; default = 50.0)) m
	| tmax $(@bind tmax_nb NumberField(0.5:0.5:60.0; default = 5.0)) s
	"""
end

# ╔═╡ a00000011-0000-4000-8000-0000000011
begin
	md"""**Numerics:** spatial order O$(@bind order_sel Select(["4", "6", "8", "10", "12"]; default = "4"))
	| PML cells $(@bind pml_nb Slider(10:5:60; default = 40, show_value = true))
	"""
end

# ╔═╡ a00000012-0000-4000-8000-0000000012
begin
	ENV2D_TEXT = if PROFILE isa Tuple
		xs, ds, _ = PROFILE
		dsm = smooth_nb > 1 ? env_smooth(ds, smooth_nb) : ds
		env_config2d(; dx = dx_nb, dz = 2.0, xmax = xmax_nb, zmax = zmax_nb,
			pml_cells = pml_nb, prof_x = xs, prof_d = dsm,
			flat_depth = flat_depth,
			c_water = c_w, c_water_grad = c_g, c_bottom = c_b)
	end
	recs = [parse(Float64, s) for s in split(recs_tf, ","; keepempty = false)]
	SR2D_TEXT = if PROFILE isa Tuple
		env_sr2d(; f0 = f0_nb, src_amp = 1.0e3, source_x = sx_nb, source_z = sz_nb,
			receiver_ranges_km = recs, receiver_depth = rdep_nb, dx = dx_nb, dz = 2.0)
	end
	nx = Int(round(xmax_nb / dx_nb)) + 1
	nz = Int(round(zmax_nb / 2.0)) + 1
	md"""**Config preview** — grid `nx × nz = $(nx) × $(nz)` cells = $(round(nx*nz/1e6, digits=2)) M cells (state + CPML ≈ $(round(nx*nz*8*6/1024^2, digits=0)) MiB). dt comes from the solver's CFL (order-dependent); the **pre-flight gate** prints the dispersion verdict before the run allocates anything. Receivers kept inside the grid: $(sum(r -> r*1000 < xmax_nb - pml_nb*dx_nb, recs)) of $(length(recs))"""
end

# ╔═╡ a00000013-0000-4000-8000-0000000013
begin
	md"""Install the generated config into the suite clone:
	
	$(@bind install_btn Button("Write config files"; primary = true))"""
end

# ╔═╡ a00000014-0000-4000-8000-0000000014
begin
	install_btn
	SUITE = sc_suite_path()
	INSTALL_MSG = if SUITE == ""
		md"⚠️ **Suite not found — run fdtdapp.jl once, or set ACOUSTIC_SUITE.**"
	elseif !(ENV2D_TEXT isa String) || !(SR2D_TEXT isa String)
		md"⚠️ **Load a profile first (section 1).**"
	else
		f = sc_install_env2d(SUITE, ENV2D_TEXT, SR2D_TEXT)
		# fresh run: stale checkpoints from other configs would corrupt the run
		rm(joinpath(SUITE, "pekeris_ckpt"); force = true, recursive = true)
		rm(joinpath(SUITE, "pekeris_frames"); force = true, recursive = true)
		rm(joinpath(SUITE, "pekeris_run"); force = true, recursive = true)
		# driver derived from the suite’s run.jl, but with our tmax/frame settings
		drv = sc_make_driver(SUITE, tmax_nb; frame_every = 100)
		md"✅ config written: `$(relpath(f.env, SUITE))` + `$(relpath(f.sr, SUITE))`, driver `$(basename(drv))` (tmax = $(tmax_nb) s); ckpt/frames cleared"
	end
end

# ╔═╡ a00000015-0000-4000-8000-0000000015
md"""## 3 · Trigger the calculation"""

# ╔═╡ a00000016-0000-4000-8000-0000000016
begin
	md"""Threads (phone: 4 is a good start) $(@bind nthreads_nb Slider(1:1:16; default = 4, show_value = true))
	
	$(@bind run2d_btn Button("🚀 Run 2-D FDTD"; primary = true))
	"""
end

# ╔═╡ a00000017-0000-4000-8000-0000000017
begin
	run2d_btn
	SUITE = sc_suite_path()
	RUN2D_MSG = if SUITE == ""
		md"⚠️ no suite"
	else
		sc_start_job(SUITE, "run2d", `run_app.jl`; threads = nthreads_nb, order = order_sel)
		md"▶️ 2-D FDTD started — refresh below."
	end
end

# ╔═╡ a00000018-0000-4000-8000-0000000018
begin
	md"""$(@bind refresh2d Button("↻ Refresh status"))"""
end

# ╔═╡ a00000019-0000-4000-8000-0000000019
begin
	refresh2d
	st = sc_status("run2d")
	md"""**run2d** — $(st.running ? "🟢 running $(st.secs)s" : "⚪ finished/stopped")
	
	```
$(st.tail)
	```"""
end

# ╔═╡ a0000001a-0000-4000-8000-000000001a
md"""## 4 · Results"""

# ╔═╡ a0000001b-0000-4000-8000-000000001b
begin
	md"""$(@bind show_results Button("📸 Show traces + wavefield"))"""
end

# ╔═╡ a0000001c-0000-4000-8000-000000001c
begin
	show_results
	SUITE2 = sc_suite_path()
	RESULTS = try
		SUITE2 == "" && error("no suite")
		meta_p = joinpath(SUITE2, "meta.jl")
		isfile(meta_p) || error("no meta.jl yet — run the 2-D job first")
		meta = qp_read_meta(meta_p)
		traces_p = joinpath(SUITE2, get(meta, "out_dir", "pekeris_run"), "traces.csv")
		isfile(traces_p) || error("no traces.csv")
		format_p = get(meta, "frame_dir", "pekeris_frames")
		frames = [f for f in readdir(joinpath(SUITE2, format_p)) if occursin(r"^frame_\d+\.f32$", f)]
		(meta, traces_p, joinpath(SUITE2, format_p), sort(frames))
	catch err
		md"⚠️ **$(sprint(showerror, err))**"
	end
end

# ╔═╡ a0000001d-0000-4000-8000-000000001d
md"Receiver traces (normalised per trace):"

# ╔═╡ a0000001e-0000-4000-8000-000000001e
begin
	if RESULTS isa Tuple
		qp_traces_svg(RESULTS[2])
	else
		RESULTS
	end
end

# ╔═╡ a0000001f-0000-4000-8000-000000001f
begin
	if RESULTS isa Tuple
		md"""**Pressure wavefield** — frame $(@bind frame_idx_nb NumberField(1:1:2000; default = 1))"""
	end
end

# ╔═╡ a00000020-0000-4000-8000-0000000020
begin
	if RESULTS isa Tuple
		meta, _, fdir, frames = RESULTS
		k = clamp(frame_idx_nb, 1, length(frames))
		buf = qp_read_f32(joinpath(fdir, frames[k]))
		nx = get(meta, "nx", 0)::Int
		nz = get(meta, "nz", 0)::Int
		mat = reshape(buf, nx, nz)          # column-major p, i fastest
		mat = permutedims(mat)              # rows = z
		uri, (lo, hi) = qp_heatmap_png(mat; W = 360)
		md"""frame $(frames[k]) mapped from $(round(lo, digits=2)) .. $(round(hi, digits=2)) Pa
		$(qp_html_img(uri))
		"""
	else
		md""
	end
end

# ╔═╡ a00000021-0000-4000-8000-0000000021
md"""## 5 · More solvers
3-D FDTD uses the same profile (extruded in y). FWI / cross-solver triggers
run the suite CLIs directly.
"""

# ╔═╡ a00000022-0000-4000-8000-0000000022
begin
	md"""**3-D FDTD slab:** dy $(@bind dy_nb NumberField(1.0:1.0:50.0; default = 3.0)) m |
	ymax $(@bind ymax_nb NumberField(200.0:100.0:3000.0; default = 600.0)) m |
	xmax3d $(@bind xmax3d_nb NumberField(500.0:500.0:10000.0; default = 2000.0)) m | tmax $(@bind tmax3d_nb NumberField(0.2:0.2:10.0; default = 1.0)) s
	
	$(@bind install3d_btn Button("Write 3-D config + run (run3d_app.jl)"; primary = true))
	"""
end

# ╔═╡ a00000023-0000-4000-8000-0000000023
begin
	install3d_btn
	SUITE = sc_suite_path()
	RUN3D_MSG = if SUITE == "" || !(PROFILE isa Tuple)
		md"⚠️ need a suite + loaded profile"
	else
		xs, ds, _ = PROFILE
		dsm = smooth_nb > 1 ? env_smooth(ds, smooth_nb) : ds
		e3 = env_config3d(; dx = dx_nb, dy = dy_nb, dz = 2.0, xmax = xmax3d_nb,
			ymax = ymax_nb, zmax = zmax_nb, pml_cells = 20, prof_x = xs, prof_d = dsm,
			flat_depth = flat_depth, c_water = c_w, c_water_grad = c_g, c_bottom = c_b)
		s3 = env_sr3d(; f0 = f0_nb, src_amp = 1.0e3, source_x = sx_nb,
			source_y = ymax_nb / 2, source_z = sz_nb,
			receiver_ranges_km = [r for r in recs if r * 1000 < xmax3d_nb - 20 * dx_nb],
			receiver_depth = sz_nb, dx = dx_nb, dy = dy_nb, dz = 2.0, ymax = ymax_nb)
		sc_install_env3d(SUITE, e3, s3)
		rm(joinpath(SUITE, "pekeris3d_ckpt"); force = true, recursive = true)
		rm(joinpath(SUITE, "pekeris3d_frames"); force = true, recursive = true)
		rm(joinpath(SUITE, "pekeris3d_run"); force = true, recursive = true)
		sc_make_driver(SUITE, tmax3d_nb; three_d = true)
		sc_start_job(SUITE, "run3d", `run3d_app.jl`; threads = nthreads_nb, order = order_sel)
		md"🚀 3-D FDTD started — 3-D grids are EXPENSIVE on a phone: keep nx·ny small."
	end
end

# ╔═╡ a00000024-0000-4000-8000-0000000024
begin
	md"""$(@bind refresh3d Button("↻ 3-D status"))"""
end

# ╔═╡ a00000025-0000-4000-8000-0000000025
begin
	refresh3d
	st = sc_status("run3d")
	md"""**run3d** — $(st.running ? "🟢 running $(st.secs)s" : "⚪ finished/stopped")
	
	```
	$(st.tail)
	```"""
end

# ╔═╡ a00000026-0000-4000-8000-0000000026
begin
	md"""**Full-waveform inversion** — trigger the suite CLI:
	`julia fwi.jl <env> <sr> <obs.jls> <outdir> [opts]`
	
	env.jl $(@bind fwi_env TextField((50, 1); default = "toy/out/ms_env.jl"))
	sr.jl $(@bind fwi_sr TextField((50, 1); default = "toy/out/ms_sr.jl"))
	obs.jls $(@bind fwi_obs TextField((50, 1); default = "toy/out/obs.jls"))
	options $(@bind fwi_opts TextField((50, 1); default = "--iters 10 --opt lbfgs-pre"))
	
	$(@bind fwi_btn Button("🚀 Run fwi.jl"; primary = true))
	"""
end

# ╔═╡ a00000027-0000-4000-8000-0000000027
begin
	fwi_btn
	SUITE = sc_suite_path()
	ok, a, b, c = true, fwi_env, fwi_sr, fwi_obs
	FWI_MSG = if SUITE == ""
		md"⚠️ no suite"
	else
		outd = joinpath(get(ENV, "ACOUSTIC_RUNS", pwd()), "fwi_out")
		mkpath(outd)
		extra = [t for t in split(fwi_opts, " "; keepempty = false)]
		cmd = `fwi.jl $(abspath(SUITE, a)) $(abspath(SUITE, b)) $(abspath(SUITE, c)) $outd $extra`
		sc_start_job(SUITE, "fwi", cmd; threads = nthreads_nb)
		md"🚀 FWI started (run.jl-style resume applies per fwi.jl)."
	end
end

# ╔═╡ a00000028-0000-4000-8000-0000000028
begin
	md"""$(@bind refreshfwi Button("↻ FWI status"))"""
end

# ╔═╡ a00000029-0000-4000-8000-0000000029
begin
	refreshfwi
	st = sc_status("fwi")
	md"""**fwi** — $(st.running ? "🟢 running $(st.secs)s" : "⚪ finished/stopped")
	
	```
	$(st.tail)
	```"""
end

# ╔═╡ a0000002a-0000-4000-8000-000000002a
begin
	md"""**Cross-solver validation** — Kraken × PE × Bellhop at one frequency:
	f0 $(@bind xv_f NumberField(5.0:5.0:500.0; default = 50.0)) Hz
	
	$(@bind xv_btn Button("Run toy/xvalidate.jl"))
	"""
end

# ╔═╡ a0000002b-0000-4000-8000-000000002b
begin
	xv_btn
	SUITE = sc_suite_path()
	XV_MSG = SUITE == "" ? md"⚠️ no suite" :
		(sc_start_job(SUITE, "xvalidate", `toy/xvalidate.jl $(xv_f)`; threads = min(nthreads_nb, 2));
		 md"🚀 xvalidate started")
end

# ╔═╡ a0000002c-0000-4000-8000-000000002c
begin
	md"""$(@bind refreshxv Button("↻ xvalidate status"))"""
end

# ╔═╡ a0000002d-0000-4000-8000-000000002d
begin
	refreshxv
	st = sc_status("xvalidate")
	md"""**xvalidate** — $(st.running ? "🟢 running $(st.secs)s" : "⚪ finished/stopped")
	
	```
	$(st.tail)
	```"""
end

# ╔═╡ a0000002e-0000-4000-8000-000000002e
md"""
## Notes & gotchas
* **RAM on the phone**: a 1 M-cell grid ≈ 50 MiB of state; OOC frames are
  enabled by `run.jl` itself (frames + checkpoints on disk). Keep
  `nx·nz` under ~5 M cells on an 8 GB phone.
* **Resume**: `run.jl` resumes from `pekeris_ckpt`. The install buttons
  clear ckpt/frames so a config change can never resume stale state.
* **SL anchor**: with real bathymetry the meta `water_depth` is the mean
  profile depth, so absolute TL levels carry that approximation.
* **Frames**: `frame_%06d.f32` are column-major `p = -(txx+tzz)/2`;
  frame cadence is set by run.jl (`frame_every = 100`).
* The suite stays untouched: config files are backed up once
  (`.plutosuite_orig`).
"""

# ╔═╡ Cell order:
# ╟─ a0000001-0000-4000-8000-000000000001
# ╟─ a0000002-0000-4000-8000-000000000002
# ╟─ a0000003-0000-4000-8000-000000000003
# ╟─ a0000004-0000-4000-8000-000000000004
# ╟─ a0000005-0000-4000-8000-000000000005
# ╟─ a0000006-0000-4000-8000-000000000006
# ╟─ a0000007-0000-4000-8000-000000000007
# ╟─ a0000008-0000-4000-8000-000000000008
# ╟─ a0000009-0000-4000-8000-000000000009
# ╟─ a000000a-0000-4000-8000-00000000000a
# ╠═ a000000b-0000-4000-8000-00000000000b
# ╟─ a000000c-0000-4000-8000-00000000000c
# ╟─ a0000003-0000-4000-8000-00000000000d
# ╟─ a000000e-0000-4000-8000-00000000000e
# ╟─ a000000f-0000-4000-8000-00000000000f
# ╟─ a00000010-0000-4000-8000-0000000010
# ╟─ a00000011-0000-4000-8000-0000000011
# ╠═ a00000012-0000-4000-8000-0000000012
# ╟─ a00000013-0000-4000-8000-0000000013
# ╠═ a00000014-0000-4000-8000-0000000014
# ╠═ a00000016-0000-4000-8000-0000000016
# ╟─ a00000018-0000-4000-8000-0000000018
# ╠═ a00000019-0000-4000-8000-0000000019
# ╟─ a0000001a-0000-4000-8000-000000001a
# ╟─ a0000001b-0000-4000-8000-000000001b
# ╠═ a0000001c-0000-4000-8000-000000001c
# ╟─ a0000001d-0000-4000-8000-000000001d
# ╠═ a0000001e-0000-4000-8000-000000001e
# ╟─ a0000001f-0000-4000-8000-000000001f
# ╠═ a00000020-0000-4000-8000-0000000020
# ╟─ a00000021-0000-4000-8000-0000000021
# ╠═ a00000022-0000-4000-8000-0000000022
# ╠═ a00000023-0000-4000-8000-0000000023
# ╟─ a00000024-0000-4000-8000-0000000024
# ╠═ a00000025-0000-4000-8000-0000000025
# ╠═ a00000026-0000-4000-8000-0000000026
# ╠═ a00000027-0000-4000-8000-0000000027
# ╟─ a00000028-0000-4000-8000-0000000028
# ╠═ a00000029-0000-4000-8000-0000000029
# ╠═ a0000002a-0000-4000-8000-000000002a
# ╠═ a0000002b-0000-4000-8000-000000002b
# ╟─ a0000002c-0000-4000-8000-000000002c
# ╠═ a0000002d-0000-4000-8000-000000002d
# ╟─ a0000002e-0000-4000-8000-000000002e
