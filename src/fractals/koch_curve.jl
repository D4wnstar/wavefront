using Meshes
using Rotations: Angle2d
import WGLMakie

A = Point(0, 0)
B = Point(1, 0)

vec = (B - A) |> Rotate(Angle2d(π / 6))
viz(vec)