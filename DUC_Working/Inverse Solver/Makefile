ccfile=fem.v30
exefile=fem.v30

cbexe=/home/shreyas/Desktop/Files/Cubit_LINUX64.13.0/cubit -nographics -batch -noecho -nobanner -nojournal
gf=/usr/bin/gfortran
seldon=/home/shreyas/Desktop/Files/seldon-5.1.2
minc=/home/shreyas/Desktop/Files/MUMPS_4.10.0/include
mseq=/home/shreyas/Desktop/Files/MUMPS_4.10.0/libseq
mlib=/home/shreyas/Desktop/Files/MUMPS_4.10.0/lib

fem: ${ccfile}.cc
	g++ -o ${exefile} -O3 -ffast-math ${ccfile}.cc -I${seldon} -I${minc} -I${mseq} -L${mlib} -ldmumps -lzmumps -lmumps_common -lpord -L${mseq} -lmpiseq -L${gf} -lgfortran -lfftw3 -lm -lpthread -llapack -lblas

mesh: ${cbsfile}.jou
	${cbexe} ${cbsfile}.jou

all: fem mesh

clean:
	rm -f ${exefile}
