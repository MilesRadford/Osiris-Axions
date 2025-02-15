write(output_unit,'(A)') 'build_info = {'
write(output_unit,'(A)') '  "version": "",'
write(output_unit,'(A)') '  "branch": "",'
write(output_unit,'(A)') '  "dimensions": "2",'
write(output_unit,'(A)') '  "precision": "DOUBLE",'
write(output_unit,'(A)') '  "build_type": "production",'
write(output_unit,'(A)') '  "system": "scarf",'
write(output_unit,'(A)') '  "build_tool": "make",'
write(output_unit,'(A)') '  "build_flags": {'
write(output_unit,'(A)') '    "FPP"  :"mpicc -C -nostdinc -E -x assembler-w&
   &ith-cpp -D__HAS_MPI_IN_PLACE__ -D_OPENMP -DHDF5",'
write(output_unit,'(A)') '    "FPPF" :"-DP_X_DIM=2 -DOS_REV=\"\" -DFORTRANS&
   &INGLEUNDERSCORE -DPRECISION_DOUBLE -DENABLE_QED -DENABLE_NEUTRAL_SPIN -&
   &DENABLE_XXFEL -DENABLE_GR",'
write(output_unit,'(A)') '    "F90"  :"mpif90 -fallow-argument-mismatch -fa&
   &llow-invalid-boz",'
write(output_unit,'(A)') '    "F90F" :"-pipe -ffree-line-length-none -fno-r&
   &ange-check -fopenmp -O3 -march=native   -I /apps20/sw/easybuilt/rocky/9&
   &/amd/zen2/software/HDF5/1.14.0-gompi-2023a/include",'
write(output_unit,'(A)') '    "CF"   :"-O3 -march=native -std=c99  -DFORTRA&
   &NSINGLEUNDERSCORE -DPRECISION_DOUBLE -D__POSIX_TIMER__",'
write(output_unit,'(A)') '    "cc"   :"mpicc -Wa,-q",'
write(output_unit,'(A)') '    "LDF"  :" -L /apps20/sw/easybuilt/rocky/9/amd&
   &/zen2/software/HDF5/1.14.0-gompi-2023a/lib -lhdf5_fortran -lhdf5 -lz -l&
   &dl -lm -Wl,-rpath -Wl,/apps20/sw/easybuilt/rocky/9/amd/zen2/software/HD&
   &F5/1.14.0-gompi-2023a/lib   ",'
write(output_unit,'(A)') '   }'
write(output_unit,'(A)') '}'

