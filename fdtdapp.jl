#!/usr/bin/env julia
# ============================================================================
# fdtdapp.jl — the Acoustic_Android launcher.
#
# Boots everything needed to run the notebook app on a phone (JuliaTerm,
# Fairphone-class device) or a desktop:
#   1. ensures guest tools (curl, ca-certificates, unzip, tar) when running
#      inside the JuliaTerm proot Debian guest (absolute-path checks only —
#      in the guest the HOST Termux prefix sits on PATH, so `command -v`
#      would find host binaries and lie);
#   2. ensures a clone of the FDTD2DElastic suite under
#      ~/.acoustic_android/FDTD2DElastic (tarball via curl, no git needed);
#   3. ensures Pluto + PlutoUI in the default Julia environment
#      (installed first, then imported via eval — imports of to-be-installed
#      packages cannot happen at top level);
#   4. starts Pluto on http://127.0.0.1:8890 and tries to open the notebook.
#
# Phone (JuliaTerm), with a workdir selected in the drawer:
#     cd /workdir && julia /path/to/fdtdapp.jl
# then open http://127.0.0.1:8890 in the phone browser.
# Desktop:
#     cd anywhere/ && julia fdtdapp.jl
# ============================================================================

import Pkg
import Printf: @printf

const APP_HOME = normpath(dirname(abspath(PROGRAM_FILE)))
const SUITE_DIR = joinpath(homedir(), ".acoustic_android", "FDTD2DElastic")
const SUITE_TARBALL = "https://codeload.github.com/St0n89/FDTD2DElastic/tar.gz/refs/heads/main"

say(_) = println("[fdtdapp] ", _)

function apt_install(name)
    say("apt-get install -y $name")
    ok = success(`apt-get update -qq`)
    ok || say("apt-get update failed (offline?) — continuing anyway")
    success(`apt-get install -y $name`) || error("could not install $name; install it manually and retry")
end

function ensure_guest_tools()
    isfile("/usr/bin/apt-get") || return   # desktop / not a Debian guest
    say("proot Debian guest detected — checking tools (absolute paths only)")
    isfile("/usr/bin/curl") || apt_install("curl")
    isfile("/etc/ssl/certs/ca-certificates.crt") || apt_install("ca-certificates")
    isfile("/usr/bin/unzip") || apt_install("unzip")
    (isfile("/bin/tar") || isfile("/usr/bin/tar")) || apt_install("tar")
    if isfile("/etc/ssl/certs/ca-certificates.crt")
        ENV["SSL_CERT_FILE"] = "/etc/ssl/certs/ca-certificates.crt"
    end
end

function ensure_suite()
    isfile(joinpath(SUITE_DIR, "run.jl")) && return SUITE_DIR
    say("suite not found — trying the public tarball")
    mkpath(dirname(SUITE_DIR))
    tmp = joinpath(dirname(SUITE_DIR), "suite.tgz")
    ok = try
        run(`curl -sS --fail -L -o $tmp $SUITE_TARBALL`)
        run(pipeline(`tar -xzf $tmp -C $(dirname(SUITE_DIR))`; stdout = devnull))
        inner = joinpath(dirname(SUITE_DIR), "FDTD2DElastic-main")
        if isdir(inner)
            isdir(SUITE_DIR) && mv(SUITE_DIR, SUITE_DIR * ".backup_$(time_ns())")
            mv(inner, SUITE_DIR)
        end
        true
    catch err
        say("tarball fetch failed ($(sprint(showerror, err)))")
        false
    end
    rm(tmp; force = true)
    if !ok || !isfile(joinpath(SUITE_DIR, "run.jl"))
        say("=================================================================")
        say("FDTD2DElastic is a PRIVATE repo — I cannot fetch it anonymously.")
        say("Provide the suite once, then rerun. Any of these works:")
        say("  * export ACOUSTIC_SUITE=/path/to/your/FDTD2DElastic clone")
        say("  * git clone it (with your credentials) to:")
        say("        $SUITE_DIR")
        say("  * unzip a suite zip there so that $SUITE_DIR/run.jl exists")
        say("Starting Pluto anyway — the notebook shows a hint until then.")
        return ""
    end
    say("suite ready: $SUITE_DIR")
    return SUITE_DIR
end

function ensure_pluto()
    say("ensuring Pluto + PlutoUI (default environment; idempotent)")
    Pkg.add(["Pluto", "PlutoUI"])
    Core.eval(Main, :(import Pluto))
    return Main.Pluto
end

function main()
    say("Acoustic_Android launcher")
    ensure_guest_tools()
    ENV["ACOUSTIC_ANDROID_HOME"] = APP_HOME
    ENV["ACOUSTIC_RUNS"] = pwd()
    haskey(ENV, "ACOUSTIC_SUITE") || (ENV["ACOUSTIC_SUITE"] = SUITE_DIR)
    suite = ensure_suite()
    isempty(suite) || (ENV["ACOUSTIC_SUITE"] = suite)
    jup = joinpath(homedir(), ".juliaup", "bin", "julia")
    isfile(jup) && (ENV["ACOUSTIC_JULIA"] = jup)
    Pluto = ensure_pluto()
    nb = joinpath(APP_HOME, "notebooks", "Suite.jl")
    say("starting Pluto on http://127.0.0.1:8890 ...")
    try
        session = Pluto.ServerSession()
        Pluto.run!(session, 8890, "127.0.0.1")
        say("browser: http://127.0.0.1:8890")
        Pluto.SessionActions.open!(session, nb)
        say("notebook open: $nb")
        while true
            sleep(10)
        end
    catch err
        @info "auto-open failed; falling back to plain Pluto.run" error = sprint(showerror, err)
        Pluto.run(8890; launch_browser = false, dismiss_update_notification = true)
    end
end

main()
