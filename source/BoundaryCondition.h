/*
 * BoundaryCondition.h
 *
 *  Created on: Apr 15, 2014
 *      Author: david
 */

#ifndef BOUNDARYCONDITION_H_
#define BOUNDARYCONDITION_H_

#include "Basic.h"

class CBC_c;
class  BC_c;

//----------------------------------------------------------------------------------------------------
//
//                      About Partial Differential Equations' Boundary Conditions
//
//----------------------------------------------------------------------------------------------------
class CBC_c {//Cell Boundary Conditions expressed in coefficient and RHS value
public:
	BCType_e T;
	numeric_t aP; //Coefficient for itself
	numeric_t b; //Right Side Value
	std::vector<   index_t> I;
	std::vector< numeric_t> a;

	CBC_c():T(NoBC), aP(0), b(0), I(0), a(0) {}

	CBC_c(const BCType_e & t, const numeric_t & self, const numeric_t & rhs)
	: T(t), aP(self), b(rhs)
	, I(0), a(0) {}

	CBC_c(const BCType_e & t, const numeric_t & self, const numeric_t & rhs
		, const index_t & i1, const numeric_t & a1)
	: T(t), aP(self), b(rhs)
	, I(1, i1), a(1, a1) {}

	CBC_c(const BCType_e & t, const numeric_t & self, const numeric_t & rhs
		, const index_t & i1, const numeric_t & a1
		, const index_t & i2, const numeric_t & a2)
	: T(t), aP(self), b(rhs)
	, I(2, i1), a(2, a1) {
		I[1]=i2; a[1]=a2;
	}

	~CBC_c()=default;

	friend std::ostream & operator << (std::ostream & os, const CBC_c & b);
};
//----------------------------------------------------------------------------------------------------
//
//                      Boundary Conditions
//
//----------------------------------------------------------------------------------------------------
class BC_c {//Boundary Conditions
public:
	BCType_e T;
	numeric_t u, v, p, c;
	index_t UpstreamLI, DownstreamLI, PI;
	numeric_t TL;

	BC_c()
	: T(NoBC), u(0), v(0), p(0), c(0), UpstreamLI(voidIndex), DownstreamLI(voidIndex), PI(voidIndex)
	, TL(0) {}

	BC_c(const BCType_e & t, const numeric_t & uv, const numeric_t & vv, const numeric_t & pv, const numeric_t & cv) {
		T=t;
		if(T==Inlet) {
			u=uv;
			v=vv;
			p=pv;
			c=cv;
			UpstreamLI=DownstreamLI=PI=voidIndex;
			TL=0;
		}
	}

	BC_c(const BCType_e & t) {
		T=t;
		if(T==Outlet) {
			u=v=p=c=0;
			UpstreamLI=DownstreamLI=PI=voidIndex;
			TL=0;
		}
	}

	BC_c(const BCType_e & t, const numeric_t & uv, const numeric_t & vv) {
		T=t;
		if(T==Wall) {
			u=uv;
			v=vv;
			p=c=0;
			UpstreamLI=DownstreamLI=PI=voidIndex;
			TL=0;
		}
	}

	BC_c(const BCType_e & t, const numeric_t & pv, const index_t & i) {
		T=t;
		if(T==PeriodicUpstream) {
			u=v=c=0;
			p=pv;
			UpstreamLI=PI=voidIndex;
			DownstreamLI=i;
			TL=0;
		}
	}

	BC_c(const BCType_e & t, const index_t & i, const numeric_t & pv, const numeric_t & cv) {
		T=t;
		if(T==PeriodicDownstream) {
			u=v=0;
			p=pv;
			c=cv;
			UpstreamLI=i;
			DownstreamLI=PI=voidIndex;
			TL=0;
		}
	}

	BC_c(const BCType_e & t, const index_t & i, const numeric_t & l) {
		T=t;
		if(T==PorousOpen) {
			u=v=p=c=0;
			UpstreamLI=DownstreamLI=voidIndex;
			PI=i;
			TL=l;
		}
	}

	BC_c(const BCType_e & t, const numeric_t & uv, const numeric_t & vv, const numeric_t & cv) {
		T=t;
		if(T==VelocityInlet) {
			u=uv;
			v=vv;
			p=0;
			c=cv;
			UpstreamLI=DownstreamLI=PI=voidIndex;
			TL=0;
		}
	}

	BC_c(const BCType_e & t, const numeric_t & pcv) {
		T=t;
		if(T==PressureOutlet) {
			u=v=c=0;
			p=pcv;
		} else if(T==ConstConcentration) {
			u=v=p=0;
			c=pcv;
		} else {}
		UpstreamLI=DownstreamLI=PI=voidIndex;
		TL=0;
	}

	inline bool operator == (const BC_c & b) const {return (T==b.T) && (u==b.u) && (v==b.v) && (p==b.p) && (c==b.c) &&
			                                     (UpstreamLI==b.UpstreamLI) && (DownstreamLI==b.DownstreamLI) && (PI==b.PI);}
	inline bool operator != (const BC_c & b) const {return !((T==b.T) && (u==b.u) && (v==b.v) && (p==b.p) && (c==b.c) &&
				                                     (UpstreamLI==b.UpstreamLI) && (DownstreamLI==b.DownstreamLI) && (PI==b.PI));}

	friend std::istream & operator >> (std::istream & is, BC_c & b);
	friend std::ostream & operator << (std::ostream & os, const BC_c & b);
};

#endif /* BOUNDARYCONDITION_H_ */
