/*
 * BoundaryCondition.cpp
 *
 *  Created on: Apr 15, 2014
 *      Author: david
 */

#include "BoundaryCondition.h"

using namespace std;

//----------------------------------------------------------------------------------------------------
//
//                      Boundary Conditions Member Functions
//
//----------------------------------------------------------------------------------------------------
istream & operator >> (istream & is, BC_c & b) {
	string ts;
	is>>ts;

	if(ts=="NoBC") {
		b.T=NoBC;
		b.u=b.v=b.c=0; b.TL=0;
		b.UpstreamLI=b.DownstreamLI=b.PI=voidIndex;
	} else if(ts=="Inlet") {
		b.T=Inlet;
		is>>b.u>>b.v>>b.p>>b.c;
		b.TL=0;
		b.UpstreamLI=b.DownstreamLI=b.PI=voidIndex;
	} else if(ts=="Outlet") {
		b.T=Outlet;
		b.u=b.v=b.p=b.c=0; b.TL=0;
		b.UpstreamLI=b.DownstreamLI=b.PI=voidIndex;
	} else if(ts=="Wall") {
		b.T=Wall;
		is>>b.u>>b.v;
		b.p=b.c=0; b.TL=0;
		b.UpstreamLI=b.DownstreamLI=b.PI=voidIndex;
	} else if(ts=="PeriodicUpstream") {
		b.T=PeriodicUpstream;
		is>>b.p>>b.c;
		b.u=b.v=0; b.TL=0;
		b.UpstreamLI=b.DownstreamLI=b.PI=voidIndex;
	} else if(ts=="PeriodicDownstream") {
		b.T=PeriodicDownstream;
		is>>b.p;
		b.u=b.v=b.c=0; b.TL=0;
		b.UpstreamLI=b.DownstreamLI=b.PI=voidIndex;
	} else if(ts=="PorousOpen") {
		b.T=PorousOpen;
		b.u=b.v=b.p=b.c=0; b.TL=0;
		b.UpstreamLI=b.DownstreamLI=b.PI=voidIndex;
	} else if(ts=="VelocityInlet") {
		b.T=VelocityInlet;
		is>>b.u>>b.v>>b.c;
		b.TL=0;
		b.UpstreamLI=b.DownstreamLI=b.PI=voidIndex;
	} else if(ts=="PressureOutlet") {
		b.T=PressureOutlet;
		is>>b.p;
		b.u=b.v=b.c=0; b.TL=0;
		b.UpstreamLI=b.DownstreamLI=b.PI=voidIndex;
	} else if(ts=="ConstConcentration") {
		b.T=ConstConcentration;
		is>>b.c;
		b.u=b.v=b.p=0; b.TL=0;
		b.UpstreamLI=b.DownstreamLI=b.PI=voidIndex;
	}

	return is;
}

ostream & operator << (ostream & os, const BC_c & b) {
	switch (b.T) {
	case NoBC:
		os<<"NoBC";
		break;
	case Inlet:
		os<<"Inlet with U="<<b.u<<"[m/s] & V="<<b.v<<"[m/s] & P="<<b.p<<"[Pa] & C="<<b.c<<"[mol/m^3]";
		break;
	case Outlet:
		os<<"Outlet";
		break;
	case Wall:
		if(b.u==0 && b.v==0) {
			os<<"Static Wall";
		} else if(b.u==0 && b.v!=0) {
			os<<"Sliding Wall with Speed V="<<b.v<<"[m/s]";
		} else if(b.u!=0 && b.v==0) {
			os<<"Sliding Wall with Speed U="<<b.u<<"[m/s]";
		} else {
			os<<"Wall can only slide!";
		}
		break;
	case PeriodicUpstream:
		os<<"Periodic To   "<<b.DownstreamLI<<" with Pressure="<<b.p<<"[Pa] & Concentration="<<b.c<<"[mol/m^3]";
		break;
	case PeriodicDownstream:
		os<<"Periodic from "<<  b.UpstreamLI<<" with Pressure="<<b.p<<"[Pa]";
		break;
	case PorousOpen:
		os<<"PorousOpen"<<" with PI="<<b.PI;
		break;
	case VelocityInlet:
		os<<"Inlet with U="<<b.u<<"[m/s] & V="<<b.v<<"[m/s] & C="<<b.c<<"[mol/m^3]";
		break;
	case PressureOutlet:
		os<<"Outlet with P="<<b.p<<"[Pa]";
		break;
	case ConstConcentration:
		os<<"Concentration: C="<<b.c<<"[mol/m^3]";
		break;
	case FixedValueBC:
		os<<"What? This is wrong!";
		break;
	case Relation0BC:
		os<<"What? This is wrong!";
		break;
	case Relation1BC:
		os<<"What? This is wrong!";
		break;
	case Relation2BC:
		os<<"What? This is wrong!";
		break;
	}

	return os;
}
