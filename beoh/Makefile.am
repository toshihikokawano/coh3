bin_PROGRAMS = beoh
beoh_LDADD = ../omplib/libomp.a
beoh_SOURCES = \
	beoh.cpp \
	beoh.h \
	beohbeta.h\
	beohstat.h \
	beohdataread.cpp \
	beohdummy.cpp \
	beohensdf.cpp \
	beohffrag.cpp \
	beohffrag.h \
	beohffragcalc.cpp \
	beohffragbeta.cpp \
	beohffragdecay.cpp \
	beohffragenergy.cpp \
	beohlabspectrum.cpp \
	beohmain.cpp \
	beohneutrino.cpp \
	beohomcalc.cpp \
	beohoutput.h \
	beohoutput.cpp \
	beohpopinit.cpp \
	beohprofile.cpp \
	beohsetup.cpp \
	beohspectra.cpp \
	beohspecmc.cpp \
	ffpmultichance.cpp \
	ffpparameter.cpp \
	ffpoutput.cpp \
	ffpyield.cpp \
	ffpyieldread.cpp \
	ffp.h \
	wahl.cpp \
	wahl.h \
	betachain.h \
	fobs.h \
	mchf.h \
	$(COHSRC) $(OMPLIB)

COHSRC = \
	asympt.cpp \
	beta2.cpp \
	builtin.cpp \
	coulomb.cpp \
	coupling.cpp \
	datatable.cpp \
	etc.cpp \
	extwave.cpp \
	fisdens.cpp \
	fispop.cpp \
	fpotential.cpp \
	gampop.cpp \
	gamcas.cpp \
	gdr.cpp \
	gtrans.cpp \
	gtransread.cpp \
	intwave.cpp \
	kcksyst.cpp \
	levden.cpp \
	masstable.cpp \
	mt19937.cpp \
	omsetform.cpp \
	omsetparm.cpp \
	outprep.cpp \
	parameter.cpp \
	parpop.cpp \
	setupinit.cpp \
	spline.cpp \
	smatrix.cpp \
	statcompounddecay.cpp \
	statlevdecay.cpp \
	statmaniplevel.cpp \
	statmodeltools.cpp \
	stattrans.cpp \
	statparm.cpp \
	ripl2levels.cpp

AM_CPPFLAGS = \
	-I../source -I../omplib

AM_CXXFLAGS = \
	-Wall -Wextra -DDATADIR=\"${pkgdatadir}\" -DBeoH

