# SuiteCtl.jl — suite discovery, config installation and job triggering.
#
# Acoustic_Android never modifies the suite's code: it only writes the two
# config files run.jl / run3d.jl include (config/env_pekeris.jl,
# config/sources_receivers.jl, config/env3d.jl, config/sources_receivers3d.jl),
# after backing the originals up once per clone (suffix .plutosuite_orig).
# A config fingerprint (.plutosuite_hash) guards checkpoint-resume across
# config changes (stale checkpoints would corrupt a new run).

using Printf: @sprintf

const SUITE_DEFAULT_DIR = joinpath(homedir(), ".acoustic_android", "FDTD2DElastic")

function sc_suite_path(; quiet = false)
    p = get(ENV, "ACOUSTIC_SUITE", "")
    (p != "" && isfile(joinpath(p, "run.jl"))) && return p
    isfile(joinpath(SUITE_DEFAULT_DIR, "run.jl")) && return SUITE_DEFAULT_DIR
    # maybe the notebook lives inside a suite checkout
    local root = get(ENV, "ACOUSTIC_ANDROID_HOME", "")
    for cand in (root, homedir(), pwd())
        cand != "" || continue
        q = joinpath(cand, "FDTD2DElastic")
        isfile(joinpath(q, "run.jl")) && return q
        isfile(joinpath(cand, "run.jl")) && return cand
    end
    return ""
end

function sc_ensure_orig(suite, relpath)
    full = joinpath(suite, relpath)
    orig = full * ".plutosuite_orig"
    if isfile(full) && !isfile(orig)
        cp(full, orig)
    end
    return full
end

function sc_write_config(suite, relpath, text)
    isdir(joinpath(suite, "config")) || error("suite config/ not found under $suite")
    full = sc_ensure_orig(suite, relpath)
    open(full, "w") do io
        print(io, text)
    end
    return full
end

function sc_install_env2d(suite, envtext, srtext)
    f1 = sc_write_config(suite, "config/env_pekeris.jl", envtext)
    f2 = sc_write_config(suite, "config/sources_receivers.jl", srtext)
    sc_flag_config(suite, "run2d", envtext * srtext)
    return (env = f1, sr = f2)
end

function sc_install_env3d(suite, envtext, srtext)
    f1 = sc_write_config(suite, "config/env3d.jl", envtext)
    f2 = sc_write_config(suite, "config/sources_receivers3d.jl", srtext)
    sc_flag_config(suite, "run3d", envtext * srtext)
    return (env = f1, sr = f2)
end

function sc_flag_config(suite, tag, text)
    h = hash(text)
    open(joinpath(suite, ".plutosuite_hash"), "w") do io
        println(io, "$tag $(UInt(h))")
    end
end
function sc_config_changed(suite, tag, text)
    p = joinpath(suite, ".plutosuite_hash")
    isfile(p) || return true
    lines = readlines(p)
    any(l -> startswith(l, "$tag ") && l != "$tag $(UInt(hash(text)))", lines)
end

# ---------------- job drivers ----------------
"""
    sc_make_driver(suite, tmax; frame_every = 100, ckpt_every = 3000,
                   three_d = false) -> path

Derives `run_app.jl` (or `run3d_app.jl`) from the suite's own run.jl /
run3d.jl at install time: identical logic, only the RUN SETTINGS block
(proj/tmax/frame_every/ckpt_every) carries the notebook's values. The
suite's own files stay untouched.
"""
function sc_make_driver(suite, tmax; frame_every = 100, ckpt_every = 3000,
                        three_d = false)
    srcf = joinpath(suite, three_d ? "run3d.jl" : "run.jl")
    isfile(srcf) || error("driver source not found: $srcf")
    txt = read(srcf, String)
    n1 = count(r"(?m)^tmax\s*="m, txt)
    n1 == 1 || error("run driver drifted: tmax anchor x$n1")
    txt = replace(txt, r"(?m)^tmax\s*=.*$" => "tmax        = $tmax")
    txt = replace(txt, r"(?m)^frame_every\s*=.*$" => "frame_every = $(frame_every)")
    txt = replace(txt, r"(?m)^ckpt_every\s*=.*$" => "ckpt_every  = $(ckpt_every)")
    dst = joinpath(suite, three_d ? "run3d_app.jl" : "run_app.jl")
    open(dst, "w") do io
        print(io, txt)
    end
    return dst
end

# ---------------- jobs ----------------
mutable struct ScJob
    name::String
    proc::Base.Process
    log::String
    started::Float64
end

const SC_JOBS = Dict{String,ScJob}()

"""
    sc_start_job(suite, name, cmd::Cmd; threads = 2, order = "4", extra_env = Dict())
Runs `julia -t N <cmd>` inside the suite dir with FDTD_ORDER set.
Jobs with the same name are killed first. Log goes to
<workdir>/logs/<name>.log (workdir: ENV["ACOUSTIC_RUNS"] || pwd()).
"""
function sc_start_job(suite, name, cmd::Cmd; threads = 2, order = "4", extra_env = Dict{String,String}())
    sc_kill_job(name)
    workdir = get(ENV, "ACOUSTIC_RUNS", pwd())
    logdir = joinpath(workdir, "logs")
    mkpath(logdir)
    logp = joinpath(logdir, name * ".log")
    jbin = sc_julia_bin()
    full = addenv(`$jbin -t $(threads) $cmd`, merge(Dict("FDTD_ORDER" => order), extra_env))
    full = setcwd(full, suite)
    io = open(logp, "w")
    proc = run(pipeline(full; stdout = io, stderr = io), wait = false)
    SC_JOBS[name] = ScJob(name, proc, logp, time())
    return SC_JOBS[name]
end

function sc_julia_bin()
    # JuliaTerm guest: the host prefix is on PATH there — NEVER rely on PATH
    # lookups; use the juliaup absolute path (the jt wrapper guarantees it).
    p = get(ENV, "ACOUSTIC_JULIA", "")
    p != "" && isfile(p) && return p
    cand = joinpath(homedir(), ".juliaup", "bin", "julia")
    isfile(cand) && return cand
    pc = something(Sys.which("julia"), "")
    pc != "" && return pc
    error("no julia binary found (set ACOUSTIC_JULIA)")
end

function sc_running(j::ScJob)
    return process_running(j.proc)
end

function sc_kill_job(name)
    j = get(SC_JOBS, name, nothing)
    j === nothing && return false
    if sc_running(j)
        try
            kill(j.proc)
            sleep(0.5)
        catch
        end
    end
    delete!(SC_JOBS, name)
    return true
end

function sc_status(name)
    j = get(SC_JOBS, name, nothing)
    j === nothing && return (running = false, secs = 0.0, tail = "(no job '$name' yet)")
    running = process_running(j.proc)
    return (running = running, secs = round(time() - j.started; digits = 1),
            tail = sc_tail(j.log))
end

function sc_tail(path; n = 30)
    isfile(path) || return "(no log yet)"
    lines = readlines(path)
    join(last(lines, n), "\n")
end

# ---------------- result discovery ----------------
sc_newest(pathpat, dir) = begin
    cand = filter(f -> occursin(pathpat, f), readdir(dir))
    isempty(cand) ? "" : joinpath(dir, last(sort!(cand)))
end
