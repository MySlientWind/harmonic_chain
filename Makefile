build_all:
	gfortran -fopenmp ./main.f90 -o ./bin/main_p.exe
	gfortran ./main.f90 -o ./bin/main_np.exe