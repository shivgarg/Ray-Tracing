Ray Tracer
==========

### Description
Recursive ray tracer has been implemented . Simple scenes which contain plane and sphere can be rendered. The following things are supported in rendering:-
* Objects affine properties.
* Material properties of objects ( reflection coefficients etc).
* Local illumination model (diffuse, spectacular and ambient components).
* Multiple light Sources.
* Global illumination with reflection, refraction and shadows.
* Anti-aliasing using supersampling

### Usage
1. **make affine** makes the binary for ray tracing when affine transformation has been applied to object.
2. **make main** makes the binary which incorporates rest of all the features. 
3. ``**./<binary name> < <sample input file name> ``
Sample input files have been provided in ``sample_inputs`` directory.
The format of the input files is given in format file.
