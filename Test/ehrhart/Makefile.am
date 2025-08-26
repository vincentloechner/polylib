#
# $Id: Makefile.am,v 1.8 2005/10/20 11:35:25 verdoolaege Exp $
#

OBJ_DIR = $(top_builddir)


EHRHART_BASE =  \
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

EHRHART_EXTRA=bart8.in isnm2.in e16.in esum.in t12.in

EHRHART=$(EHRHART_BASE) $(EHRHART_EXTRA)

EXTRA_DIST = $(EHRHART)

# not using %.in GNU extension
dist-hook:
	for pl_file in $(EHRHART) ; do test -f $(distdir)/`basename $$pl_file .in`.out || cp -p $(srcdir)/`basename $$pl_file .in`.out $(distdir) ; done

# not using filter-out gnu extension
EHRHART_TESTS = `if test $(TEST_BITS) = 32 ; then echo $(EHRHART_BASE) ; else echo $(EHRHART) ; fi `

tests:
	@failedtest=0; \
	failedtestnames=""; \
	for x in $(EHRHART_TESTS) ; do \
		file=`basename $$x .in`; \
		printf "Verify file $$file... " ; \
		$(OBJ_DIR)/testehrhart$(TEST_BITS) < $(srcdir)/$$x > xyz; \
		diff -w xyz $(srcdir)/$$file.out ; \
		result=$$?; \
		if [ "$$result" -eq "0" ]; then \
			printf "passed\n"; \
		else \
			printf "\033[31mError: $$file.out is not the same\033[0m\n"; \
			failedtest=`expr $$failedtest + 1`; \
			failedtestnames="$$failedtestnames:$$file"; \
	    fi; \
	done ; \
	if [ $$failedtest != 0 ]; then \
		echo "\033[31m$$failedtest tests failed\033[0m"; \
		echo "These tests failed: $$failedtestnames:"; \
		exit 1; \
	else \
		echo "Test Ehrhart works correctly"; \
	fi


longtest:
	@for x in $(EHRHART_TESTS) ; do \
		echo "Verify Ehrhart for $$x"; \
		$(OBJ_DIR)/verif_ehrhart$(TEST_BITS) < $(srcdir)/$$x ;\
	done ;
CLEANFILES=xyz
