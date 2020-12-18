FC = gfortran

all:
	${FC} mo_rte_kind.F90  mo_rte_spherical_correction.F90 test_spherical_correction.f90 -o test

clean:
	rm -f *.mod* *.o	./test