# AcousticCore.jl — include helper. Acoustic_Android's Julia-side helpers
# (plain top-level functions, no modules, so a plain Base.include into Main
# from a notebook cell or a script works everywhere).
#
#   Base.include(Main, "/path/to/src/AcousticCore.jl")

_acoustic_inc(f) = include(joinpath(@__DIR__, f))
_acoustic_inc("NetCDF3.jl")
_acoustic_inc("GeoTIFF.jl")
_acoustic_inc("GebcoCrop.jl")
_acoustic_inc("EnvLoad.jl")
_acoustic_inc("QuickPlot.jl")
_acoustic_inc("SuiteCtl.jl")
