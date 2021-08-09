# pimpleFoamPIDv2

pimpleFoamPIDv2 is an alternative version of the state observer-based data-assimilation technique XXXX. This novel approach assimilates pressure data through a state observer
that is constructed based on a PID control law and that acts on
the pressure equation. Reference data can be experimentally-acquied pressure readins, or results from an expensive computational simulation run in a much denser grid.

The approach allows for more accurate RANS simulations in inexpensive setups. 

## Installation

Simply download the content of this project into your machine and run the following lines in your OpenFOAM installation:

```
wclean
wmake
```

If the compiler does not present errors, pimpleFoamPIDv2 will be ready to use.

Included in this 

## Usage

