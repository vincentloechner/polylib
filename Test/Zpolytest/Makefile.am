#
# $Id: Makefile.am,v 1.2 2005/10/20 11:35:25 verdoolaege Exp $
# Modified: Vincent Loechner, 2025
#

TEST_EXE = $(top_builddir)/Zpolytest
TEST_FILES_32 =  \
	LatDiff1.in \
	LatDiff2.in \
	LatDiff3.in \
	LatDiff3b.in \
	LatDiff5.in \
	LatDiff7.in \
	LatDiff7b.in \
	LatDiff8.in \
	LatDiff8b.in \
	LatHNF1.in \
	LatHNF2.in \
	LatHNF3.in \
	LatInc1.in \
	LatInc2.in \
	LatInc3.in \
	LatInter1.in \
	LatInter1b.in \
	LatInter2.in \
	LatInter4.in \
	LatInter5.in \
	LatInter6.in \
	LatInter7.in \
	phideo2.in \
	ZAlloc1.in \
	ZAlloc1b.in \
	ZAlloc1c.in \
	ZAlloc2.in \
	ZAlloc3.in \
	ZAlloc3b.in \
	ZAlloc3c.in \
	ZComp1.in \
	ZComp3.in \
	ZComp3b.in \
	ZCompComp1.in \
	ZCompComp3.in \
	ZCompComp3b.in \
	ZDiff1.in \
	ZDiff2.in \
	ZDiff3.in \
	ZDiff4.in \
	ZDiff5.in \
	ZDiff6.in \
	ZDiff7.in \
	ZDiff8.in \
	ZDiff9.in \
	ZDiff10.in \
	ZDiff12.in \
	ZDiff13.in \
	ZDiffInter6a.in \
	ZDomain2b.in \
	ZDomain3.in \
	ZDomain4.in \
	ZDomain6.in \
	ZDomain10.in \
	ZDomain12.in \
	ZDomain14.in \
	ZEmpty2.in \
	ZImage1.in \
	ZImPre1.in \
	ZImPre2.in \
	ZImPre3.in \
	ZInter1.in \
	ZInter2.in  \
	ZInter3.in  \
	ZInter8.in  \
	ZInter11.in  \
	ZInter12.in  \
	ZPre1.in \
	ZPre2.in  \
	ZPre3.in \
	ZUnion1.in \
	ZUnion2.in \
	ZUnion3.in


# 32 bits overflow in those:
TEST_FILES = $(TEST_FILES_32) \
	ZComp11.in \
	ZComp11b.in \
	ZCompComp11.in \
	ZCompComp11b.in \
	ZDiff11.in \
	ZDiff11b.in \
	ZDomain11.in

# not checked:
# - long execution time:
# phideo1.in
# phideo3.in
# ZDisj11b.in
# ZDomain13.in  # need gmp (overflow even in 64 bits)

# - not yet included:
# ZDisj8.in

# - old (to be updated):
# SplitZpol1.in
# SplitZpol2.in
# SplitZpol3.in
# ZAlloc[4-7].in
# ZCheckEmpty1.in

# - temporary/no longer used:
# LeftHermite*.in
# ZEmpty1.in
# ZImage[2-4]*.in


include ../Tests.mk
