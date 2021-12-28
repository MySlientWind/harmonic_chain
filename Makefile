build_all:
	gfortran -fopenmp ./main.f90 -o ./bin/main_p.exe
	gfortran ./main.f90 -o ./bin/main_np.exe
push-to-cheng:
	rsync -av --exclude-from=.gitignore ./ Cheng:/home/jingtao/work/fortran_test/harmonic_chain