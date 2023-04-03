module Types

using Unitful, StaticArrays

export Length, Speed, Time, Mass
export AccelerationMagnitude, GravitationalParameter, SpecificEnergy

export Position, Velocity, Acceleration

const Length = typeof(1.0u"m")
const Speed = typeof(1.0u"m/s")
const Time = typeof(1.0u"s")
const Mass = typeof(1.0u"kg")

const AccelerationMagnitude = typeof(1.0u"m/s^2")
const GravitationalParameter = typeof(1.0u"m^3/s^2")
const SpecificEnergy = typeof(1.0u"m^2/s^2")

const Position = SVector{3,Length}
const Velocity = SVector{3,Speed}
const Acceleration = SVector{3,typeof(1.0u"m/s^2")}

end