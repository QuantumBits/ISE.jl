module ISE

using Random, LinearAlgebra

using Unitful, DynamicalSystems
using SatelliteToolbox

using BenchmarkTools, Plots, GLMakie

include("Types.jl")
using .Types

#! include("GNC.jl")
#! 
#! struct GNC_MODEL
#!     c::GNC.commands
#!     p::GNC.parameters
#!     o::GNC.output
#! end

include("LibGNC.jl")

struct GNC_MODEL
    c::LibGNC.commands
    p::LibGNC.parameters
    o::LibGNC.output
end


include("THR.jl")

struct THR_MODEL
    c::THR.commands
    p::THR.parameters
    o::THR.output
end

include("ORB.jl")

struct ORB_MODEL
    c::ORB.commands
    p::ORB.parameters
    o::ORB.output
end

include("GPS.jl")

struct GPS_MODEL
    c::GPS.commands
    p::GPS.parameters
    o::GPS.output
end


mutable struct ISE_MODEL
    thr::THR_MODEL
    orb::ORB_MODEL
    gps::GPS_MODEL
    gnc::GNC_MODEL
    Δt::Time
end

function extract_state(m::ISE_MODEL)::SVector
    return ustrip.(SVector(
        m.orb.o.t,
        m.orb.o.r...,
        m.orb.o.v...,
        m.thr.o.m_prop,
        m.gnc.o.specific_energy
    ))
end

function extract_parameters(m::ISE_MODEL)::SVector
    return ustrip.(SVector(m.thr.p.m_dry, m.Δt))
end

function reset_model!(m::ISE_MODEL, u::SVector, p::SVector)
end

function sim!(m::ISE_MODEL)

    # Actuators

    thr_i = THR.input(m.thr.o.m_prop, m.orb.o.t, m.thr.o.t)
    thr_c = THR.commands(m.gnc.o.thr_on)
    thr_o, thr_c = THR.step(thr_i, thr_c, m.thr.p)

    # Physics

    thr_a = thr_o.T .* m.gnc.o.thr_dir
    t = m.orb.o.t + m.Δt

    orb_i = ORB.input(m.orb.o.r, m.orb.o.v, thr_a, t, m.orb.o.t)
    orb_c = ORB.commands()
    orb_o, orb_c = ORB.step(orb_i, orb_c, m.orb.p)

    # Sensors

    gps_i = GPS.input(m.orb.o.r, m.orb.o.v, m.orb.o.t)
    gps_c = GPS.commands()
    gps_o, gps_c = GPS.step(gps_i, gps_c, m.gps.p)

    # GNC

    #! gnc_i = GNC.input(gps_o.r_est, gps_o.v_est, thr_o.m_prop)
    #! gnc_c = GNC.commands()
    #! gnc_o, gnc_c = GNC.step(gnc_i, gnc_c, m.gnc.p)

    gnc_i = LibGNC.input(Tuple(ustrip.(upreferred.(gps_o.r_est))), Tuple(ustrip.(upreferred.(gps_o.v_est))), ustrip(upreferred(thr_o.m_prop)))
    gnc_c = LibGNC.commands(false)
    gnc_o = LibGNC.step(Ref(gnc_i), Ref(gnc_c), Ref(m.gnc.p))

    m.thr = THR_MODEL(thr_c, m.thr.p, thr_o)
    m.orb = ORB_MODEL(orb_c, m.orb.p, orb_o)
    m.gps = GPS_MODEL(gps_c, m.gps.p, gps_o)
    m.gnc = GNC_MODEL(gnc_c, m.gnc.p, gnc_o)

end

function run_stuff(duration::Time=1.0u"hr", Δt::Time=0.1u"s")

    # Initial Conditions

    RE = 6371.009u"km"
    μE = 398600.4418u"km^3/s^2"
    R0 = RE + 400u"km"
    V0 = sqrt(μE / R0)

    r0 = Position(R0, 0.0u"m", 0.0u"m")
    v0 = Velocity(0.0u"m/s", V0, 0.0u"m/s")
    a0 = Acceleration(0.0u"m/s^2", 0.0u"m/s^2", 0.0u"m/s^2")
    t0 = 0.0u"s"
    ϵ0 = norm(v0)^2 / 2 - μE / norm(r0)

    # VASIMR
    Isp = 4900.0u"s" * 0.72 # (72% efficiency)
    g0 = 9.80665u"m/s^2"
    dmdt = 107.0u"mg/s"
    m_dry = 3500.0u"kg"
    m_prop = 310.0u"kg"

    # Bus Initialization

    thr_c = THR.commands(true)
    thr_p = THR.parameters(m_dry, Isp, g0, dmdt)
    thr_o = THR.output(0u"m/s^2", m_prop, t0)

    orb_c = ORB.commands()
    orb_p = ORB.parameters(μE) # , load_gravity_model(JGM3()))
    orb_o = ORB.output(r0, v0, a0, t0)

    gps_c = GPS.commands()
    gps_p = GPS.parameters(
        MersenneTwister(1), 10.0u"m",
        MersenneTwister(2), 1.0u"mm/s",
        MersenneTwister(3), 1.0u"μs")
    gps_o = GPS.output(r0, v0, t0)

    #! gnc_c = GNC.commands()
    #! gnc_p = GNC.parameters(μE)
    #! gnc_o = GNC.output(false, SVector(0.0, 0.0, 0.0), ϵ0)

    gnc_c = LibGNC.commands(false)
    gnc_p = LibGNC.parameters(ustrip(upreferred(μE)))
    gnc_o = LibGNC.output(false, (0.0, 0.0, 0.0), ustrip(upreferred(ϵ0)))

    m = ISE_MODEL(
        THR_MODEL(thr_c, thr_p, thr_o),
        ORB_MODEL(orb_c, orb_p, orb_o),
        GPS_MODEL(gps_c, gps_p, gps_o),
        GNC_MODEL(gnc_c, gnc_p, gnc_o),
        Δt)

    ise = DynamicalSystems.ArbitrarySteppable(m, ISE.sim!, ISE.extract_state, ISE.extract_parameters, ISE.reset_model!)

    N = upreferred(duration ./ m.Δt)

    DynamicalSystems.step!(ise, N)

    return DynamicalSystems.trajectory(ise, N)

end

function plot_stuff(trajectory)

    X, t = trajectory

    t = uconvert.(u"d", X[:, 1]u"s")
    rX = X[:, 2]u"m"
    rY = X[:, 3]u"m"
    rZ = X[:, 4]u"m"
    vX = X[:, 5]u"m/s"
    vY = X[:, 6]u"m/s"
    vZ = X[:, 7]u"m/s"
    m_prop = X[:, 8]u"kg"
    ϵ = X[:, 9]u"m^2/s^2"

    RE = 6371009.0u"m"

    plotly()

    p_ϵ = Plots.plot(t, ϵ, color=:red, label="ϵ")
    p_m_prop = Plots.plot(t, m_prop, color=:blue, label="m_prop",
        ylims=(0.0u"kg", maximum(m_prop)))
    Plots.plot(p_ϵ, p_m_prop)
end

"""
Make a 3D video of trajectory!
"""
function viz_stuff(trajectory)

    X, t = trajectory

    t = X[:, 1] # "s"
    rX = X[:, 2] # "m"
    rY = X[:, 3] # "m"
    rZ = X[:, 4] # "m"
    vX = X[:, 5] # "m/s"
    vY = X[:, 6] # "m/s"
    vZ = X[:, 7] # "m/s"
    m_prop = X[:, 8] # "kg"

    RE = 6371009.0 # u"m"


    points = Observable(Point3f[])
    colors = Observable(Int[])

    RMAX = 1.1 * maximum(sqrt.(rX .^ 2 + rY .^ 2 + rZ .^ 2)) / RE

    set_theme!(theme_black())

    fig, ax, l = lines(points, color=colors,
        colormap=:inferno, transparency=true,
        axis=(; type=Axis3, protrusions=(0, 0, 0, 0),
            viewmode=:fit, limits=RMAX .* (-1, 1, -1, 1, -1, 1)))

    # earth_img = GLMakie.load(download("https://svs.gsfc.nasa.gov/vis/a000000/a002900/a002915/bluemarble-2048.png"))
    # mesh!(fig, Sphere(Point3f(0), 1.0), color=earth_img, shading = true)

    FR = 3000 # floor(Int, min(length(t), 120))
    TI = floor(Int, length(t) / FR)

    record(fig, "orbit.mp4", 1:TI:length(t); framerate=FR) do f
        push!(points[], Point3f(rX[f] / RE, rY[f] / RE, rZ[f] / RE))
        push!(colors[], f)
        notify.((points, colors))
        l.colorrange = (0, f)
    end

end

end