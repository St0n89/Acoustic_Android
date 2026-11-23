# Acoustic_Android

A Pluto.jl app that turns your phone (or laptop) into an underwater-acoustics
workstation: trigger 2-D / 3-D FDTD, FWI and cross-validation runs of the
[FDTD2DElastic](https://github.com/St0n89/FDTD2DElastic) suite from a browser
notebook, with the environment set up from real bathymetry (GEBCO 2024 crop
API, BathyApp GeoTIFF / CSV exports, or a flat Pekeris waveguide).

Built to run inside [JuliaTerm](https://github.com/St0n89/JuliaTerm) /
[JuliaTermux](https://github.com/St0n89/JuliaTermux) on Android; also works on
any Linux box with Julia >= 1.9.

## Quick start (any machine with Julia)

```bash
julia fdtdapp.jl
# then open http://127.0.0.1:8890 in a (phone) browser
```

`fdtdapp.jl` checks for the guest tools it needs (`curl`, `unzip`, `tar`, CA
certs), tries to bootstrap the suite (tarball of `St0n89/FDTD2DElastic` into
`~/.acoustic_android/FDTD2DElastic`), installs Pluto + PlutoUI and opens
`notebooks/Suite.jl` bound to 127.0.0.1:8890.

**The suite repo is private**: if the anonymous tarball fails, clone it once
yourself (or point `ACOUSTIC_SUITE=/path/to/FDTD2DElastic` at an existing
clone); the launcher prints the exact expected paths. Pluto starts either
way, and the notebook tells you where to put the suite.

On the phone: push this repo into the JuliaTerm(Tux) home, run the launcher
with the juliaup julia (`~/.juliaup/bin/julia fdtdapp.jl`), and browse on the
same device. On a desktop you can also point `ACOUSTIC_SUITE=/path/to/FDTD2DElastic`
at an existing clone and the app uses that, untouched-but-for-config.

## What the notebook does

1. **Environment section** — choose the seabed:
   * *flat* Pekeris waveguide (depth + sound-speed gradient),
   * *GeoTIFF* exported from BathyApp (uncompressed Float32, geotags or
     ModelTransformation),
   * *CSV* profile exported from BathyApp (`depth` positive-down; `elev` is
     flipped, great-circle distances if lon/lat columns are present),
   * *GEBCO live crop*: enter a bounding box, the app posts to
     `https://download.gebco.net/api/queue` (grid_id=1, data_source_ids=[1],
     formats=[1]), polls, downloads and unzips the ZIP, and reads the NetCDF
     CDF-1 `.nc` itself (scale_factor / add_offset / _FillValue handled).
   The loaded profile is previewed as SVG; a smoothing window can be applied.
2. **Configure** — grid (dx, xmax, zmax, PML), water column (c0, gradient,
   bottom c / rho), acquisition (f0, source position, receiver ranges / depth),
   numerics (spatial order O4..O12 via FDTD_ORDER, tmax, threads).
3. **Install** — writes `config/env_pekeris.jl` and
   `config/sources_receivers.jl` (2-D) or `config/env3d.jl` /
   `config/sources_receivers3d.jl` (3-D) into the suite, with the exact
   variable contract the suite drivers expect, backs up the originals once
   (`.plutosuite_orig`), clears stale ckpt/frames, and derives
   `run_app.jl` / `run3d_app.jl` from the suite's own run scripts with only
   the RUN SETTINGS block (tmax, frame_every, ckpt_every) replaced.
4. **Trigger** — starts the run as a background job (`FDTD_ORDER`, threads),
   logs to `~/.acoustic_runs/logs/`, shows live status, lets you kill / tail.
   FWI and the Kraken x PE x Bellhop cross-validation are exposed as further
   one-click jobs.
5. **Results** — receiver traces and wavefield frames are rendered inline
   (pure-Julia SVG / hand-written PNG encoder; no Plots dependency in the
   notebook host either).

The suite itself stays unmodified: only the generated `config/` files and the
derived `run_app.jl` / `run3d_app.jl` are added inside the suite folder.

## Layout

```
fdtdapp.jl          launcher (guest-tool check, suite bootstrap, Pluto server)
notebooks/Suite.jl  the app (Pluto notebook)
src/NetCDF3.jl      CDF-1 (classic NetCDF) reader, GEBCO tiles
src/GeoTIFF.jl      uncompressed Float32 GeoTIFF reader (BathyApp exports)
src/GebcoCrop.jl    GEBCO 2024 crop queue client (curl + unzip subprocesses)
src/EnvLoad.jl      profile loaders + env/source-receiver config generators
src/QuickPlot.jl    SVG traces/profiles, PNG heatmap (hand-written DEFLATE)
src/SuiteCtl.jl     suite discovery, config install, driver derivation, jobs
src/AcousticCore.jl include aggregator (everything lands in Main, Pluto-style)
test/selftest.jl    29-check selftest, no network needed
```

Pure Julia stdlib everywhere: no packages beyond `Pluto` / `PlutoUI` for the
notebook server itself.

## Tests / CI

```bash
julia test/selftest.jl
```

CI (Ubuntu, Julia stable) runs the selftest and a parse gate over every `.jl`
file including the notebook.

## Limitations, honestly

* **SL anchor on bathy runs is approximate**: the config sets
  `water_depth` to the mean profile depth, while `sl_reference` treats it as a
  constant water column; keep that in mind when reading absolute TL numbers
  from rough bathymetry.
* **3-D is expensive** on a phone; the notebook exposes a thin slab
  (y-stretched `dy`) rather than a cubic grid. If `ny` is smaller than the
  PML bands the suite's SL shell sampler finds an empty shell and stamps
  `p_ref = 0` (a warning is printed; use a larger `ny`).
* GEBCO crop needs network and the guest's `curl` / `unzip`; results are
  cached under `~/.acoustic_android/data/`.
* The app derives its drivers from the suite's `run.jl` / `run3d.jl` at
  install time and refuses (loudly) if that file drifts from the expected
  RUN SETTINGS anchors.

## License

MIT, (c) 2026 Dennis Kunde.
