#ifndef _EXT_EHRHART_H_
#define _EXT_EHRHART_H_

extern Enumeration *Domain_Enumerate(Polyhedron *D, Polyhedron *C,
				     unsigned MAXRAYS, const char **pn);

extern void new_eadd (evalue *e1,evalue *res);

extern void Scalar_product(Value *p1,Value *p2,unsigned length, Value *r);

#endif
