#
# $Id: Makefile.am,v 1.2 2005/10/20 11:35:25 verdoolaege Exp $
#

OBJ_DIR = $(top_builddir)

ZPOLYTESTS =  \
	LatDiff1.in \
	LatDiff2.in \
	LatDiff3.in \
	LatDiff4.in \
	LatDiff5.in \
	LatHNF1.in \
	LatHNF2.in \
	LatHNF3.in \
	LatInter1.in \
	LatInter2.in \
	LatInter3.in \
	LatInter4.in \
	LatInter5.in \
	LatInc1.in \
	LatInc2.in \
	LeftHermite1.in \
	LeftHermite2.in \
	ZEmpty1.in \
	ZEmpty2.in \
	ZInter1.in \
	ZInter2.in  \
	ZInter3.in  \
	ZUnion1.in \
	ZUnion2.in \
	ZDiff1.in \
	ZDiff2.in \
	ZDiff3.in \
	ZDiff4.in \
	ZDiff5.in \
	ZImage1.in \
	ZPre1.in \
	ZPre2.in  \
	ZPre3.in \
	ZImPre1.in \
	SplitZpol1.in \
	SplitZpol2.in \
	SplitZpol3.in

EXTRA_DIST = $(ZPOLYTESTS) \
	$(ZPOLYTESTS:%.in=%.out) 

TEST=$(ZPOLYTESTS:%.in=%.tested) 

tests:
	@failedtest=0; \
	failedtestnames=""; \
	for x in $(ZPOLYTESTS) ; do \
		file=`basename $$x .in`; \
		/bin/echo -n "Verify file $$file... " ; \
		$(OBJ_DIR)/Zpolytest$(TEST_BITS) < $(srcdir)/$$x > xyz; \
		diff -w xyz $(srcdir)/$$file.out ; \
		result=$$?; \
		if [ "$$result" -eq "0" ]; then \
			/bin/echo "$$out passed"; \
		else \
			/bin/echo "Error: $$file.out is not the same"; \
			failedtest=`expr $$failedtest + 1`; \
			failedtestnames="$$failedtestnames:$$file"; \
	    fi; \
	done ; \
	if [ $$failedtest != 0 ]; then \
		/bin/echo "$$failedtest tests failed"; \
		/bin/echo "These tests failed: $$failedtestnames:"; \
	else \
		/bin/echo "Test Zpolytest works correctly"; \
	fi

CLEANFILES=xyz
