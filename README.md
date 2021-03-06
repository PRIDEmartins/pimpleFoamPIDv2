# pimpleFoamPIDv2

```pimpleFoamPIDv2``` is an alternative version of the state observer-based data-assimilation technique ```XXXX```. This custom OF solver assimilates pressure data through a state observer
that is constructed based on a PID control law and that acts on
the pressure equation. Reference data can be experimentally-acquied pressure readings, or results from a more  computationally-intensive simulation run in a much denser grid.

The approach allows for more accurate RANS simulations in inexpensive simulation setups. 

## Installation

Simply download the content of this project into your machine and run the following lines in your OpenFOAM installation:

```
wclean
wmake
```

If the compiler does not present errors, ```pimpleFoamPIDv2``` will be ready to use.

## Usage

The folder _./SampleCase_ contains one example of application of ```pimpleFoamPIDv2```

## Known issues

- [] The code does not work in parellal processing (see [#1](https://github.com/PRIDEmartins/pimpleFoamPIDv2/issues/2#issue-962988203))
- []  The code does not work with moving mesh cases
- []  The project has only been tested with OF 8.0


## References

Neeteson, Nathan J., and David E. Rival. "State observer-based data assimilation: a PID control-inspired observer in the pressure equation." Measurement Science and Technology 31.1 (2019): 014003. ([online](https://iopscience.iop.org/article/10.1088/1361-6501/ab40d4)).