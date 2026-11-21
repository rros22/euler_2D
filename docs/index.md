# euler_2D documentation

euler_2D solves the 2D compressible Euler equations for an ideal gas on a domain made of
rectangular boxes. These pages explain how to run it, what every option does, what each built-in
case contains, how to write your own case, and how the code is organized.

| Page | Read it when you want to |
|---|---|
| [Getting started](getting-started.md) | build the code, run a first simulation and make a video |
| [Options](options.md) | know exactly what each flag of `euler_2D` and `render` does |
| [Cases](cases.md) | see the geometry, boundary conditions and defaults of each built-in case |
| [Writing a new case](new-case.md) | set up your own problem from scratch |
| [Output and rendering](output-and-rendering.md) | understand the frame files, the renderer and the Python reader |
| [Architecture](architecture.md) | understand what each part of the code does and where to change things |
| [Testing](testing.md) | know what the tests check and what to do when one fails |
| [Maintaining](maintaining.md) | follow the conventions and keep these docs up to date |

## Units

All quantities in the code are **dimensionless**. A case chooses three reference scales, and every
number in it is measured in those units:

| Quantity | Unit |
|---|---|
| length | reference length \(L\) (e.g. the domain width, or the nozzle width) |
| density | reference density \(\rho_\text{ref}\) |
| pressure | reference pressure \(p_\text{ref}\) |
| velocity | \(\sqrt{p_\text{ref} / \rho_\text{ref}}\) |
| time | \(L / \sqrt{p_\text{ref} / \rho_\text{ref}}\) |

The ratio of specific heats is \(\gamma = 1.4\) (air), set by `GAMMA` in `src/main.c`. With these
units the sound speed of the reference gas is \(\sqrt{\gamma} \approx 1.18\).

**Converting to seconds.** Pick real values for the three scales. For example, if the reference gas
is air at room conditions (\(p = 10^5\) Pa, \(\rho = 1.2\) kg/m³), the velocity unit is
\(\sqrt{10^5 / 1.2} \approx 289\) m/s. With \(L = 1\) cm, one time unit is
\(0.01 / 289 \approx 35\) µs, so `t_end = 25` is about 0.87 ms. The equations have no viscosity, so
the flow looks the same at any size; only the time scale changes.

Each case's reference scales are listed in [Cases](cases.md).
