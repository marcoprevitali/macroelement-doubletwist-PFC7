The .dat files in this folders are examples for creating a single contact interaction and doing basic tests (e.g. axial pull, simple bending etc..) to verify that you are getting the desired response with the input parameters.

All scripts use the load_properties.dat file to assign contact parameters. In general, the particles positions and rotations are fixed, displacements are prescribed (e.g. axial, rotational) and the resulting contact force is measured.

2ball_driver: axial pull (constant velocity applied to the second particle)
2ball_driver_compression: axial compression (constant negative velocity applied to the second particle)
2ball_driver_rotation: simple bending through constant spin
2ball_driver_hybrid: simple bending + axial loading applied together
2ball_non_mono: set an initial velocity, observe oscillatory behaviour (use to debug sign inversion issues).
