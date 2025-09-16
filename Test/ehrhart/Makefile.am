#
# $Id: Makefile.am,v 1.8 2005/10/20 11:35:25 verdoolaege Exp $
# Modified: Vincent Loechner, 2025
#

EHRHART_32 =  \
	0D.in \
	0D2.in \
	8hedron.in \
	alex11.in \
	alex1.in \
	bart6.in \
	c3.in\
	c4.in\
	e0.in\
	e10.in\
	e11.in\
	e12.in\
	e13.in\
	e14.in\
	e15.in\
	e1.in\
	e2b.in\
	e2.in\
	e3.in\
	e4.in\
	e5.in\
	e6.in\
	e7.in\
	e8.in\
	e9.in\
	g12.in\
	g13.in\
	g14.in\
	g1.in\
	g2.in\
	g4.in\
	g5.in\
	t1.in\
	t2.in\
	t31.in\
	t3.in\
	t4.in\
	tpar.in

# those extra tests need 64 bits or GMP precision
EHRHART_64 = \
	bart8.in \
	isnm2.in \
	e16.in \
	esum.in \
	t12.in

TEST_EXE = $(top_builddir)/testehrhart

# 32 bits specific:
TEST_FILES_32 = $(EHRHART_32)

# default:
TEST_FILES = $(EHRHART_32) $(EHRHART_64)

include ../Tests.mk


# Extra in Ehrhart: check the correctness of some initial polynomials
longtest:
	@if [ "$(TEST_BITS)" = "32" -a -n "$(TEST_FILES_32)" ]; then \
		testfiles="$(TEST_FILES_32)"; \
	else \
		testfiles="$(TEST_FILES)"; \
	fi ; \
	for x in $$testfiles ; do \
		echo "Verify Ehrhart for $$x"; \
		$(top_builddir)/verif_ehrhart$(TEST_BITS) < $(srcdir)/$$x ;\
	done ;
