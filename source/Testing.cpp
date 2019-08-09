/*
 * Testing.cpp
 *
 *  Created on: May 13, 2015
 *      Author: david
 */

#include "GlobalDeclaration.h"

using namespace std;
using namespace Eigen;

void testDivision();

/*----------------------------------------------------------------------------------------------------

                                        function Testing

----------------------------------------------------------------------------------------------------*/
void Testing() {
	cout<<"Please choose which part of the code to test:"<<endl;
	size_t wp(0);
	cin>>wp;

	switch (wp) {
	case  1:
		testDivision();
		break;
	case  2:
		break;
	case  3:
		break;
	case  4:
		break;
	}
}
/*----------------------------------------------------------------------------------------------------

                                        function testSegment

----------------------------------------------------------------------------------------------------*/
void testDivision() {
	segment_c s(0, 1.1, 0, 1, 0);
	vector<numeric_t> vs=s.fillBetween(0, 10);
	for(size_t i=0; i<vs.size(); ++i) {
		cout<<vs[i]<<endl;
	}

	vector<point_c> vp=s.fillBetween(point_c(0,0), point_c(0,10));
	for(size_t i=0; i<vp.size(); ++i) {
		cout<<vp[i]<<endl;
	}

	numberAxis_c na(0, 10);
	na.insertNumber(5, s, s, BC_c(), BC_c());
	na.fill();
	for(size_t i=0; i<na.num.size(); ++i) {
		cout<<na.num[i]<<endl;
	}

	numberAxis_c na1(0, 15);
	na1.insertSegment(5, 10, s, s, s, BC_c(), BC_c(), BC_c());
	na1.fill();
	for(size_t i=0; i<na1.num.size(); ++i) {
		cout<<na1.num[i]<<endl;
	}

	line_c l1(point_c(0,0), point_c(5,10));
	l1.insertPoint(point_c(3, 6), s, s, BC_c(), BC_c());
	l1.fill();
	for(size_t i=0; i<l1.Pt.size(); ++i) {
		cout<<l1.Pt[i]<<endl;
	}

	cout<<l1.Pt.size()<<endl;
	cout<<l1.PtI.size()<<endl;
	cout<<l1.S.size()<<endl;
	cout<<l1.B.size()<<endl;
	cout<<l1.PELI.size()<<endl;
	cout<<l1.SCI.size()<<endl;
	cout<<l1.XorY<<endl;

	line_c l2(point_c(0,0), point_c(0,10));
	l2.insertSegment(point_c(0, 3), point_c(0, 7), s, s, s, BC_c(), BC_c(Outlet), BC_c());
	l2.fill();
	for(size_t i=0; i<l2.Pt.size(); ++i) {
		cout<<l2.Pt[i]<<endl;
	}
	for(size_t i=0; i<l2.B.size(); ++i) {
		cout<<l2.B[i]<<endl;
	}
	for(size_t i=0; i<l2.S.size(); ++i) {
		cout<<l2.S[i].SN<<endl;
	}

	cout<<l2.Pt.size()<<endl;
	cout<<l2.PtI.size()<<endl;
	cout<<l2.S.size()<<endl;
	cout<<l2.B.size()<<endl;
	cout<<l2.PELI.size()<<endl;
	cout<<l2.SCI.size()<<endl;
	cout<<l2.XorY<<endl;
}
/*----------------------------------------------------------------------------------------------------

                                        function testNumberAxis

----------------------------------------------------------------------------------------------------*/
void testNumberAxis() {

}
//		for(size_t i=0; i<U.size(); ++i) {
//			for(size_t j=0; j<nodeNumber; ++j) {
//				if(U[i].NI[j]!=voidIndex) {
//					if(U[i].NP[j]!=&U[U[i].NI[j]]) cerr<<"U["<<i<<"].NP["<<j<<"] Changed!"<<endl;
//				}
//				if(U[i].FI[j]!=voidIndex) {
//					if((j==e || j==w) & U[i].FP[j]!=&C[U[i].FI[j]]) cerr<<"U["<<i<<"].FP["<<j<<"] Changed!"<<endl;
//					if((j==ne || j==nw || j==sw || j==se) && U[i].FP[j]!=&V[U[i].FI[j]]) cerr<<"U["<<i<<"].FP["<<j<<"] Changed!"<<endl;
//				}
//			}
//		}
//		for(size_t i=0; i<V.size(); ++i) {
//			for(size_t j=0; j<nodeNumber; ++j) {
//				if(V[i].NI[j]!=voidIndex) {
//					if(V[i].NP[j]!=&V[V[i].NI[j]]) cerr<<"V["<<i<<"].NP["<<j<<"] Changed!"<<endl;
//				}
//				if(V[i].FI[j]!=voidIndex) {
//					if((j==n || j==s) && V[i].FP[j]!=&C[V[i].FI[j]]) cerr<<"V["<<i<<"].FP["<<j<<"] Changed!"<<endl;
//					if((j==ne || j==nw || j==sw || j==se) && V[i].FP[j]!=&U[V[i].FI[j]]) cerr<<"V["<<i<<"].FP["<<j<<"] Changed!"<<endl;
//				}
//			}
//		}
//		for(size_t i=0; i<C.size(); ++i) {
//			for(size_t j=0; j<nodeNumber; ++j) {
//				if(C[i].NI[j]!=voidIndex) {
//					if(C[i].NP[j]!=&C[C[i].NI[j]]) cerr<<"C["<<i<<"].NP["<<j<<"] Changed!"<<endl;
//				}
//				if(C[i].FI[j]!=voidIndex) {
//					if((j==n || j==s) && C[i].FP[j]!=&V[C[i].FI[j]]) cerr<<"C["<<i<<"].FP["<<j<<"] Changed!"<<endl;
//					if((j==e || j==w) && C[i].FP[j]!=&U[C[i].FI[j]]) cerr<<"C["<<i<<"].FP["<<j<<"] Changed!"<<endl;
//				}
//			}
//		}
