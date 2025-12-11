bin_PROGRAMS = mfld
mfld_SOURCES = \
	mfld.cpp \
	mfldFRDM.cpp \
	mflddataread.cpp \
	mfldphconfig.cpp \
	mfldprint.cpp \
	mfldsd.cpp \
	mfld.h \
	$(COHSRC)

COHSRC = \
	FRDMpairing.cpp \
	FRDMpotential.cpp \
	FRDMprint.cpp \
	FRDMreadshape.cpp \
	FRDMshape.cpp \
	FRDMshell.cpp \
	FRDMsolvebcs.cpp \
	HFdensity.cpp \
	HObasis.cpp \
	HOhamiltonian.cpp \
	HOpolynomial.cpp \
	HOsorteigenvalue.cpp \
	coupling.cpp \
	eispack.cpp \
	etc.cpp \
	mt19937.cpp \
	polycalc.cpp

AM_CPPFLAGS = \
	-I../source 

AM_CXXFLAGS = \
	-Wall -Wextra -DDATADIR=\"${pkgdatadir}\"
