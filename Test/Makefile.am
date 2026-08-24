EXTRA_DIST = cases old others
SUBDIRS = Zpolytest general pp ehrhart ehrhart_union ranking
# TODO: ehrhart_union and ranking leak memory. All others pass with the -fsanitize=address option.
# SUBDIRS = Zpolytest general pp ehrhart

check:
	@for BITS in $(ALL_BITS); do \
	    printf "\n------------------ Performing tests (bits=$$BITS) ------------------\n"; \
	    for d in $(SUBDIRS); do \
				echo "-- Entering $$d --"; \
				(cd $$d; $(MAKE) TEST_BITS=$$BITS tests) || exit 1; \
	    done \
	done

longtest:
	@for BITS in $(ALL_BITS); do \
	    printf "\n----------- Performing Ehrhart longtests (bits=$$BITS) -----------\n"; \
	    (cd ehrhart; $(MAKE) TEST_BITS=$$BITS longtest) || exit 1; \
	done

dist-hook:
	rm -rf `find $(distdir) -name CVS`
