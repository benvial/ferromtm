ferromtm
==============================

.. inclusion-marker-do-not-remove

### Coupled model and homogenization of ferroelectric-dielectric metamaterials.

This repository provides the codes to run and postprocess the data for the
results obtained in this research project.

####Requirements

- python 3
- [gmsh]
- [getdp]
- make


####Installation


First clone this repository:

```
git clone https://github.com/benvial/ferromtm.git
```

Then create, activate the environment and test it:


```
cd ferromtm
make env
source activate ferromtm
make testenv
```


Finally install the required packages:
```
make req
```


[getdp]: http://getdp.info/
[gmsh]: http://gmsh.info/
