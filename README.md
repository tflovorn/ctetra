An implementation of the tetrahedron method for Brillouin zone summation, as
described in: [Blöchl, Jepsen, and Andersen, PRB 49, 16223 (1994)](http://journals.aps.org/prb/abstract/10.1103/PhysRevB.49.16223).
This paper is referred to in the documentation of this package as 'BJA94'.

# Dependencies

Requires the GNU Scientific Library. To obtain on Debian-based distribtions:

    sudo apt-get install libgsl0ldbl libgsl0-dev

# Acknowledgement

The implementation of the tetrahedron method in [Quantum ESPRESSO](http://www.quantum-espresso.org/)
was consulted during the development of this implementation in order
to clarify understanding of the method and to help verify correctness.
