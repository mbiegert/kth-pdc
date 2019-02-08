# Mandelbrot

This program utilizes the MPI platform to compute a mandelbrot set and render an image of it.

## Building

To build the program you need the mpi platform setup and use the mpi compiler.
On Archlinux simply install the package `openmpi`.

If you want to build png support you also need a version of `libpng` installed.

Don't forget to adjust the parameters in the beginning of `main()` before compiling.
You can set the resolution and viewing area of the generated set, as well as the bound.

```
mpicc mandelbrot.c -DENABLE_PNG -o mandelbrot -lm -lpng
```
or without png support
```
mpicc mandelbrot.c -o mandelbrot -lm
```

## Running

You can either run the program directly, or (preferably) utilising the mpi platform via
`mpirun`.

The program takes up to parameters.
```
--output-name,-o <filename>         Write output to <filename>.png|.txt. Overwrites existing files!
--output-mode,-o ascii|png          Output mode. Must be build with -DENABLE_PNG for png option.
--palette,-p interpolated|so|parula Only applies in png mode. Specify the color palette to use. Note, that parula is not included in the repository, since it might be copyrighted by Mathworks Inc. You can use the parula palette from a working matlab installation. You might find `convert.py` come in handy.
```
