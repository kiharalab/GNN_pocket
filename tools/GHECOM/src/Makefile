CC	      = gcc

CFLAGS	      = -O2 -g -Wall

DEST	      = .

EXTHDRS	      =

HDRS	      = \

INSTALL	      = /etc/install

LD	      = gcc

LDFLAGS	      = -O2 -g -lm -Wall

LIBS	      = 

MAKEFILE      = Makefile

OBJS	      = ghecom.o\
	        Radius.o\
	        Amino.o\
	        Grid3D.o\
	        VecFunc.o\
	        OpeMMorph.o\
	        OpeStd3DPr.o\
	        OpeClus3DPr.o\
	        MscProbe.o\
	        MscShPrb.o\
	        MscLabelAtm.o\
	        MscLabelRes.o\
	        MscPockClus.o\
	        MscSkelClus.o\
	        HeapSort.o\
	        GclsMergeSort.o\
	        GrPtMergeSort.o\
	        C3mapMergeSort.o\
	        OpnRotPrb.o\
	        TauOpnRotPrb.o\
	        MapCCP4_CHAR.o\
	        MapCCP4.o\
	        LigsitePSP.o\
	        GenSphePrb.o\
	        ClusSphePrb.o\
	        TreeSphePrb.o\
	        ReadOpt.o\
	        FldFillQue.o\
	        ClusBinPock.o\
	        Jacobi3.o\
	        ConSpheres.o\
	        io_mmCIF.o\
	        PDB_from_mmCIF.o\
	        ATOMs_ghecom_from_PDBmmCIF.o\
    	        PdbIO.o 

PRINT	      = pr

PROGRAM       = ../ghecom

SHELL	      = /bin/sh

SRCS	      = ghecom.c \
	        Radius.c \
	        Amino.c\
	        Grid3D.c \
	        VecFunc.c \
	        OpeMMorph.c \
	        OpeStd3DPr.c \
	        OpeClus3DPr.c \
	        MscProbe.c \
	        MscShPrb.c \
	        MscLabelAtm.c \
	        MscLabelRes.c \
	        MscPockClus.c \
	        MscSkelClus.c \
	        HeapSort.c \
	        GclsMergeSort.c \
	        GrPtMergeSort.c \
	        C3mapMergeSort.c \
	        OpnRotPrb.c \
	        TauOpnRotPrb.c \
	        MapCCP4_CHAR.c \
	        MapCCP4.c \
	        LigsitePSP.c \
	        GenSphePrb.c \
	        ClusSphePrb.c \
	        TreeSphePrb.c \
	        ReadOpt.c \
	        FldFillQue.c \
	        ClusBinPock.c \
	        Jacobi3.c \
	        ConSpheres.c \
	        io_mmCIF.c\
	        PDB_from_mmCIF.c\
	        ATOMs_ghecom_from_PDBmmCIF.c\
    	        PdbIO.c

all:		$(PROGRAM)

$(PROGRAM):     $(OBJS) $(LIBS)
		@echo "Linking $(PROGRAM) ..."
		@$(LD) $(OBJS) $(LDFLAGS) $(LIBS)  -o $(PROGRAM)
		@echo "done"
