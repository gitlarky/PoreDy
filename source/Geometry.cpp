/*
 * Geometry.cpp
 *
 *  Created on: Apr 15, 2014
 *      Author: david
 */

#include "Geometry.h"
#include "CalculationDomain.h"
//#include "GlobalDeclaration.h"

using namespace std;

//----------------------------------------------------------------------------------------------------
//
//                                        point_c member functions
//
//----------------------------------------------------------------------------------------------------
istream & operator >> (istream & is, point_c & p) {
	is>>p.x>>p.y;

	return is;
}

ostream & operator << (ostream & os, const point_c & p) {
	os<<"("<<p.x<<", "<<p.y<<")";

	return os;
}

//----------------------------------------------------------------------------------------------------
//
//                                        segment_c member functions
//
//----------------------------------------------------------------------------------------------------
vector<numeric_t> fillNumbers(const numeric_t & beg, const numeric_t & end, const size_t & sn, const numeric_t & g0, const numeric_t & g1) {
	vector<numeric_t> pos;
	numeric_t L(1), g(0), seg(0), segi(0);
	numeric_t length(end-beg);

	seg=L;
	pos.push_back(0);
	pos.push_back(seg);
	for(size_t i=2; i<=sn; ++i) {
		if(g1==0 && g0==0) g=defaultGrowth;
		else if(g0!=0 && g1==0) g=g0;
		else if(g0==0 && g1!=0) g=1/g1;
		else {
			if(sn!=2) g=g0+(i-2)*(1/g1-g0)/(sn-2);
			else g=0.5*(1/g1+g0);
		}

		segi=seg*g;
		L+=segi;
		pos.push_back(L);
		seg=segi;
	}
	for(size_t i=0; i<=sn; ++i) {
		if(i==0) pos[i]=beg;
		else if(i==sn) pos[i]=end;
		else pos[i]=beg+length*pos[i]/L;
	}

	return pos;
}

vector<numeric_t> fillNumbers(const numeric_t & length, const numeric_t & iniseg, const numeric_t & growth) {
	vector<numeric_t> pos;
	pos.push_back(0);
	numeric_t seg(iniseg);
	numeric_t newpos(seg);
	numeric_t therest(length-newpos);
	while(therest>0) {
		pos.push_back(newpos);

		seg*=growth;
		newpos+=seg;
		therest=length-newpos;
	}
	if(fabs(therest)<seg*(1-1/growth || pos.size()==1)) {
		pos.push_back(length);
	} else {
		vector<numeric_t> newpos=fillNumbers(0, length, pos.size()-1, growth, 0);
		for(size_t i=0; i<pos.size(); ++i) {
			pos[i]=newpos[i];
		}
	}

	return pos;
}

vector<numeric_t> segment_c::fillBetween(const numeric_t & beg, const numeric_t & end) {
	vector<numeric_t> pos;
	if(SN!=0) {
		pos=fillNumbers(beg, end, SN, GS, GL);
	} else if(SN==0 && GL==0 && LL==0) {
		pos=fillNumbers(end-beg, LS, GS);
		SN=pos.size()-1;
		pos[0]=beg; pos[SN]=end;
		for(size_t i=1; i!=SN; ++i) {
			pos[i]+=beg;
		}
	} else if(SN==0 && GS==0 && LS==0) {
		vector<numeric_t> rp=fillNumbers(end-beg, LL, GL);
		SN=rp.size()-1;
		pos.push_back(beg);
		for(size_t i=1; i!=SN; ++i) {
			pos.push_back(end-rp[SN-i]);
		}
		pos.push_back(end);
	}
	return pos;
}

vector<  point_c> segment_c::fillBetween(const point_c & p0, const point_c & p1) {
	vector<point_c> temp(0);
	if(p0.x==p1.x && p0.y!=p1.y) {
		vector<numeric_t> Y=fillBetween(p0.y, p1.y);

		for(size_t i=0; i<=SN; ++i) {
			temp.push_back(point_c(0.5*(p0.x+p1.x),Y[i]));
		}
	} else if(p0.x!=p1.x && p0.y==p1.y) {
		vector<numeric_t> X=fillBetween(p0.x, p1.x);

		for(size_t i=0; i<=SN; ++i) {
			temp.push_back(point_c(X[i],0.5*(p0.y+p1.y)));
		}
	} else {
		vector<numeric_t> X=fillBetween(p0.x, p1.x);
		vector<numeric_t> Y=fillBetween(p0.y, p1.y);

		for(size_t i=0; i<=SN; ++i) {
			temp.push_back(point_c(X[i],Y[i]));
		}
	}

	return temp;
}

//----------------------------------------------------------------------------------------------------
//
//                                        numberAxis_c member functions
//
//----------------------------------------------------------------------------------------------------
bool numberAxis_c::clear() {
	num.clear();
	S.clear();
	B.clear();

	return true;
}

bool numberAxis_c::clearMiddle() {
	numeric_t begn(*num.begin()), endn(*(num.end()-1));
	BC_c      b   (*B  .begin());

	num.clear(); num.push_back(begn); num.push_back(endn);
	S  .clear(); S  .push_back(segment_c());
	B  .clear(); B  .push_back(b);

	return true;
}

bool numberAxis_c::initialize() {
//	clear();
//
//	num.push_back(0);
//	num.push_back(0);
//	S.push_back(segment_c());
//	B.push_back(BC_c());

	return true;
}

bool numberAxis_c::insertNumber(const numeric_t & nb, const segment_c & s1, const segment_c & s2,
		            const BC_c & b1, const BC_c & b2) {
	if(nb>*num.begin() && nb<*(num.end()-1)) {
		bool RangeNotFound(true);
		for(size_t i=1; RangeNotFound && i!=num.size(); ++i) {
			if(nb<num[i] && nb!=num[i-1]) {
				vector< numeric_t>::iterator it=num.begin()+i;
				num.insert(it, nb);

				vector< segment_c>::iterator st=S.begin()+i;
				*(st-1)=s1;
				S.insert(st, s2);

				vector<BC_c>::iterator bt=B.begin()+i;
				*(bt-1)=b1;
				B.insert(bt, b2);

				RangeNotFound=false;
			}
		}
		if(RangeNotFound) {
			cerr<<"Number already existed!"<<endl;
			return false;
		} else {
			return true;
		}
	} else {
		cerr<<"Number is not in Range of this Number Axis!"<<endl;
		return false;
	}
}

bool numberAxis_c::insertSegment(const numeric_t & nb0, const numeric_t & nb1, //user is responsible for making nb0<=nb1
		             const segment_c & s1, const segment_c & s2, const segment_c & s3,
		             const BC_c & b1, const BC_c & b2, const BC_c & b3) {
	if(nb0>=nb1) {
		cerr<<"First number should be the smaller one, and second should be the larger one!"<<endl;
		return false;
	} else {

		if(nb0>*num.begin() && nb1<*(num.end()-1)) {

			bool RangeNotFound(true);
			for(size_t i=1; RangeNotFound && i!=num.size(); ++i) {

				if(nb1<num[i] && nb0>num[i-1]) {
					vector< numeric_t>::iterator it=num.begin()+i;
					num.insert(it, nb1);
					it=num.begin()+i;
					num.insert(it, nb0);
					vector< segment_c>::iterator st=S.begin()+i;
					*(st-1)=s1;
					S.insert(st, s3);
					st=S.begin()+i;
					S.insert(st, s2);
					vector<BC_c>::iterator bt=B.begin()+i;
					*(bt-1)=b1;
					B.insert(bt, b3);
					bt=B.begin()+i;
					B.insert(bt, b2);
					RangeNotFound=false;
				}
			}
			if(RangeNotFound) {
				cerr<<"Number already existed! or Segment overlap with existed one!"<<endl;
				return false;
			} else {
				return true;
			}
		} else {
			cerr<<"Number is not in Range of this Number Axis!"<<endl;
			return false;
		}
	}
}

bool numberAxis_c::fill() {
	for(size_t i=S.size(); i>=1; --i) {
		vector<numeric_t> subNumber=S[i-1].fillBetween(num[i-1], num[i]);

		vector<numeric_t>::iterator it=num.begin()+i;
		vector<segment_c>::iterator st=S.begin()+i;
		vector<     BC_c>::iterator bt=B.begin()+i;

		for(size_t j=subNumber.size()-2; j>=1; --j) {
			num.insert(it, subNumber[j]);
			it=num.begin()+i;
			S.insert(st, segment_c());
			st=S.begin()+i;
			B.insert(bt, *(bt-1));
			bt=B.begin()+i;
		}
		*(st-1)=segment_c();
	}
	return true;
}

bool numberAxis_c::filled() {
	for(size_t i=0; i!=S.size(); ++i) {
		if(S[i].SN!=1) return false;
	}
	return true;
}

/*----------------------------------------------------------------------------------------------------

                                       line_c member functions

----------------------------------------------------------------------------------------------------*/
bool line_c::clear() {
	point_c beg (*Pt .begin()), end (*(Pt .end()-1));
	index_t begi(*PtI.begin()), endi(*(PtI.end()-1));
	BC_c    b   (*B  .begin());

	Pt .clear(); Pt .push_back(beg ); Pt .push_back(end );
	PtI.clear(); PtI.push_back(begi); PtI.push_back(endi);

	S.clear(); S.push_back(segment_c());

	B.clear(); B.push_back(b);

	return true;
}

bool line_c::initialize() {
//	clear();
//
//	Pt.push_back(point_c());
//	Pt.push_back(point_c());
//
//	S.push_back(segment_c());
//
//	B.push_back(BC_c());

	return true;
}

numberAxis_c line_c::extractNumberAxis(const coordinate_e & ce) {
	numberAxis_c theAxis;
	theAxis.clear();
	for(size_t i=0; i<S.size(); ++i) {
		theAxis.S.push_back(S[i]);
		theAxis.B.push_back(B[i]);
	}
	for(size_t i=0; i<Pt.size(); ++i) {
		switch(ce) {
		case origin:
			break;
		case xAxis:
			theAxis.num.push_back(Pt[i].x);
			break;
		case yAxis:
			theAxis.num.push_back(Pt[i].y);
			break;
		}
	}

	return theAxis;
}

bool line_c::plus2NumberAxis(const numberAxis_c & xa, const numberAxis_c & ya) {
	Pt.clear(); S.clear(); B.clear();
	if(xa.S.size()==ya.S.size()) {
		for(size_t i=0; i!=xa.S.size(); ++i) {
			 Pt.push_back(point_c(xa.num[i],ya.num[i]));

			 if( xa.S[i]==ya.S[i] && xa.B[i]==ya.B[i] ) {
			     S.push_back(xa.S[i]);
				 B.push_back(xa.B[i]);
			 } else {
			     cerr<<"The setting of these two axis don't match!"<<endl;
				 return false;
			 }
		}
		Pt.push_back(point_c(xa.num[S.size()], ya.num[S.size()]));
		return true;
	} else {
		cerr<<"The size of these two axis don't match!"<<endl;
		return false;
	}
}

bool line_c::NumberAxisPlusNumber(const numberAxis_c & axis, const numeric_t & number) {
	Pt.clear(); S.clear(); B.clear();
	if(XorY==xAxis) {
		for(size_t i=0; i<axis.S.size(); ++i) {
			 Pt.push_back(point_c(axis.num[i], number));
			 S.push_back(axis.S[i]);
			 B.push_back(axis.B[i]);
		}
		Pt.push_back(point_c(axis.num[S.size()], number));
	} else if(XorY==yAxis) {
		for(size_t i=0; i<axis.S.size(); ++i) {
			 Pt.push_back(point_c(number, axis.num[i]));
			 S.push_back(axis.S[i]);
			 B.push_back(axis.B[i]);
		}
		Pt.push_back(point_c(number, axis.num[S.size()]));
	} else {}

	return true;
}

numeric_t line_c::distance(const line_c & l) {
	if(XorY==xAxis && l.XorY==xAxis) {
		return fabs(Pt[0].y-l.Pt[0].y);
	} else if(XorY==yAxis && l.XorY==yAxis) {
		return fabs(Pt[0].x-l.Pt[0].x);
	} else {
		return 0;
	}
}


bool line_c::insertPoint(const point_c & bp=point_c()
		               , const segment_c & s1=segment_c(), const segment_c & s2=segment_c()
					   , const BC_c & b1=BC_c(), const BC_c & b2=BC_c()) {
	index_t begi(*PtI.begin()), endi(*(PtI.end()-1));

	if(XorY==xAxis) {
		numberAxis_c XNA(extractNumberAxis(xAxis));
		XNA.insertNumber(bp.x, s1, s2, b1, b2);
		NumberAxisPlusNumber(XNA, 0.5*(Pt.begin()->y+(Pt.end()-1)->y));
	} else if(XorY==yAxis) {
		numberAxis_c YNA(extractNumberAxis(yAxis));
		YNA.insertNumber(bp.y, s1, s2, b1, b2);
		NumberAxisPlusNumber(YNA, 0.5*(Pt.begin()->x+(Pt.end()-1)->x));
	} else {
		numberAxis_c XNA(extractNumberAxis(xAxis)), YNA(extractNumberAxis(yAxis));
		XNA.insertNumber(bp.x, s1, s2, b1, b2);
		YNA.insertNumber(bp.y, s1, s2, b1, b2);
		plus2NumberAxis(XNA, YNA);
	}

	PtI.clear();
	PtI.resize(Pt.size(), voidIndex);
	*PtI.begin()=begi; *(PtI.end()-1)=endi;

	return true;
}

bool line_c::insertSegment(const point_c & bp0=point_c(), const point_c & bp1=point_c()
		                 , const segment_c & s1=segment_c(), const segment_c & s2=segment_c(), const segment_c & s3=segment_c()
						 , const BC_c & b1=BC_c(), const BC_c & b2=BC_c(), const BC_c & b3=BC_c()) {
	index_t begi(*PtI.begin()), endi(*(PtI.end()-1));

	if(XorY==xAxis) {
		numberAxis_c XNA(extractNumberAxis(xAxis));
		XNA.insertSegment(bp0.x, bp1.x, s1, s2, s3, b1, b2, b3);
		NumberAxisPlusNumber(XNA, 0.5*(Pt[0].y+Pt[Pt.size()-1].y));
	} else if(XorY==yAxis) {
		numberAxis_c YNA(extractNumberAxis(yAxis));
		YNA.insertSegment(bp0.y, bp1.y, s1, s2, s3, b1, b2, b3);
		NumberAxisPlusNumber(YNA, 0.5*(Pt[0].x+Pt[Pt.size()-1].x));
	} else {
		numberAxis_c XNA(extractNumberAxis(xAxis)), YNA(extractNumberAxis(yAxis));
		XNA.insertSegment(bp0.x, bp1.x, s1, s2, s3, b1, b2, b3);
		YNA.insertSegment(bp0.y, bp1.y, s1, s2, s3, b1, b2, b3);
		plus2NumberAxis(XNA, YNA);
	}

	PtI.clear();
	PtI.resize(Pt.size(), voidIndex);
	*PtI.begin()=begi; *(PtI.end()-1)=endi;

	return true;
}

bool line_c::fill() {
	index_t begi(*PtI.begin()), endi(*(PtI.end()-1));

	for(size_t i=S.size(); i>=1; --i) {
		vector<  point_c> subSeg=S[i-1].fillBetween(Pt[i-1], Pt[i]);

		vector<  point_c>::iterator  it=Pt .begin()+i;

		vector<segment_c>::iterator  st=S  .begin()+i;
		vector<     BC_c>::iterator  bt=B  .begin()+i;

		for(size_t j=subSeg.size()-2; j>=1; --j) {
			Pt.insert(it, subSeg[j]);
			it=Pt.begin()+i;
			S.insert(st, segment_c());
			st=S.begin()+i;
			B.insert(bt, *(bt-1));
			bt=B.begin()+i;
		}
		*(st-1)=segment_c();
	}

	PtI.clear();
	PtI.resize(Pt.size(), voidIndex);
	*PtI.begin()=begi; *(PtI.end()-1)=endi;

	SCI.clear();
	SCI.resize(B.size(), voidIndex);

	return true;
}

bool line_c::filled() {
	for(size_t i=0; i!=S.size(); ++i) {
		if(S[i].SN!=1) return false;
	}
	return true;
}

coordinate_e line_c::ParaEqTo(const line_c & L) {// parallel & equal length, this function is only suitable for vertical or horizontal lines
	if       (XorY==xAxis && L.XorY==xAxis
	       && Pt.begin()->x==L.Pt.begin()->x && (Pt.end()-1)->x==(L.Pt.end()-1)->x
	       && Pt.begin()->y!=L.Pt.begin()->y && (Pt.end()-1)->y!=(L.Pt.end()-1)->y) {
		return xAxis;
	} else if(XorY==yAxis && L.XorY==yAxis
		   && Pt.begin()->x!=L.Pt.begin()->x && (Pt.end()-1)->x!=(L.Pt.end()-1)->x
		   && Pt.begin()->y==L.Pt.begin()->y && (Pt.end()-1)->y==(L.Pt.end()-1)->y) {
		return yAxis;
	} else {
		return origin;
	}
}

bool line_c::copyParaEq(const line_c & L) {
	coordinate_e ce(ParaEqTo(L));
	if(ce) {
		S=L.S;

		BC_c b=B[0];//keep Boundary unchanged and distribute it to all segment
		B.clear();
		for(size_t i=0; i<L.B.size(); ++i) {
			B.push_back(b);
		}

		point_c beg (*Pt .begin()), end (*(Pt .end()-1));
		index_t begi(*PtI.begin()), endi(*(PtI.end()-1));
		Pt.clear(); PtI.clear();
		for(size_t i=0; i<L.Pt.size(); ++i) {
			if(i==0) {
				Pt .push_back(beg );
				PtI.push_back(begi);
			} else if(i==L.Pt.size()-1) {
				Pt .push_back(end );
				PtI.push_back(endi);
			} else {
				if(ce==xAxis) {
					 Pt.push_back(point_c(L.Pt[i].x, 0.5*(beg.y+end.y)));
				} else if(ce==yAxis) {
					 Pt.push_back(point_c(0.5*(beg.x+end.x), L.Pt[i].y));
				} else {}//Do nothing!
				PtI.push_back(voidIndex);
			}
		}

		S.clear();
		for(size_t i=0; i<L.S.size(); ++i) {
			S.push_back(L.S[i]);
		}
		return true;
	} else {
		return false;
	}
}
ostream & operator << (std::ostream & os, const line_c & eg) {
	for(size_t i=0; i!=eg.PtI.size(); ++i) {
		os<<eg.PtI[i]<<"\t";
	}
	os<<endl;
	for(size_t i=0; i!=eg.Pt.size(); ++i) {
		os<<eg.Pt[i]<<"\t";
	}
	os<<endl;
	for(size_t i=0; i!=eg.S.size(); ++i) {
		os<<eg.S[i].SN<<"\t";
	}
	os<<endl;
	for(size_t i=0; i!=eg.S.size(); ++i) {
		os<<eg.S[i].GS<<"\t";
	}
	os<<endl;
	for(size_t i=0; i!=eg.S.size(); ++i) {
		os<<"\t"<<eg.S[i].GL;
	}
	os<<endl;
	for(size_t i=0; i!=eg.S.size(); ++i) {
		os<<eg.S[i].LS<<"\t";
	}
	os<<endl;
	for(size_t i=0; i!=eg.S.size(); ++i) {
		os<<"\t"<<eg.S[i].LL;
	}
	os<<endl;
	for(size_t i=0; i!=eg.B.size(); ++i) {
		os<<eg.B[i]<<"\t";
	}
	os<<endl;
	return os;
}

//----------------------------------------------------------------------------------------------------
//
//                                        block_c member functions
//
//----------------------------------------------------------------------------------------------------

istream & operator >> (istream & is, block_c & b) {
	is>>b.CnI[0]>>b.CnI[1]>>b.CnI[2]>>b.CnI[3];

	string t;
	is>>t;
	if(t=="ghost") b.T=ghost;
	else if(t=="flowfield") b.T=flowfield;
	else if(t=="porousmedia") b.T=porousmedia;

	is>>b.B[0]>>b.B[1]>>b.B[2]>>b.B[3];

	if(b.T==ghost) {
		//doing nothing
	} else if (b.T==flowfield) {
		is>>b.Nx>>b.Ny
		  >>b.G[0]>>b.G[1]>>b.G[2]>>b.G[3];
	} else if (b.T==porousmedia) {
		string pt;
		is>>pt;
		if(pt=="singlePorosity") b.porosityType=singlePorosity;
		else if(pt=="downupPorosity") b.porosityType=downupPorosity;
		else if(pt=="leftrightPorosity") b.porosityType=leftrightPorosity;
		else if(pt=="assignPorosity") b.porosityType=assignPorosity;
		is>>b.Nx>>b.Ny;
		is>>b.PPConnectionRate>>b.PPCrossRate;

		is>>b.avgTLength>>b.TLVarianceRate;
		size_t theSize(b.porosityType==singlePorosity?1:2);
		b.PolyN.resize(theSize, 0); b.avgTDiameter.resize(theSize, 0); b.TDVarianceRate.resize(theSize, 0);
		b.TDisplayScale.resize(theSize, 1);
		for(size_t i=0; i<theSize; ++i) {
			is>>b.PolyN[i]>>b.avgTDiameter[i]>>b.TDVarianceRate[i];
		}
		is>>b.RandomRealization;

		is>>b.NinDiameter>>b.NinLength;
	}

	return is;
}

ostream & operator << (ostream & os, const block_c & b) {
	os<<"Corner Points: "<<b.CnI[0]<<"\t"<<b.CnI[1]<<"\t"<<b.CnI[2]<<"\t"<<b.CnI[3]<<endl;

	os<<"Block Type: "<<b.T<<endl;

	os<<"Boundary Conditions{East Side: "<<b.B[0]
	                   <<", North Side: "<<b.B[1]
	                    <<", West Side: "<<b.B[2]
	                   <<", South Side: "<<b.B[3]<<"}"<<endl;
	if(b.T==ghost) {
		//doing nothing
	} else if (b.T==flowfield) {
		os<<"Nx="<<b.Nx<<"\tNy="<<b.Ny<<endl;
		//os<<"\tLx="<<b.Lx<<"\tLy="<<b.Ly<<endl;
		os<<"Growth Rate {E:"<<b.G[0]<<"\tN:"<<b.G[1]<<"\tW:"<<b.G[2]<<"\tS:"<<b.G[3]<<"}"<<endl;
	} else if (b.T==porousmedia) {
		if(b.porosityType==singlePorosity) {
			os<<"singlePorosity"<<endl;
		} else if(b.porosityType==leftrightPorosity) {
			os<<"leftrightPorosity"<<endl;
		} else if(b.porosityType==downupPorosity) {
			os<<"updownPorosity"<<endl;
		} else if(b.porosityType==assignPorosity) {
			os<<"assignPorosity"<<endl;
		}
		os<<";\tTPPConnectionRate: "<<b.PPConnectionRate
		  <<";\tTPPCrossRate: "     <<b.PPCrossRate<<endl;

		os<<";\tavgTLength: "       <<b.avgTLength
		  <<";\tTLVarianceRate: "   <<b.TLVarianceRate;
		for(size_t i=0; i<b.PolyN.size(); ++i) {
			os<<"Porosity #"<<i;
			os<<": Polygon Type: "  <<b.PolyN[i];
			os<<"\tavgTDiameter: "  <<b.avgTDiameter[i]
			  <<"\tTDVarianceRate: "<<b.TDVarianceRate[i];
			os<<"\tTDisplayScale: " <<b.TDisplayScale[i];
		}
		os<<";\tRandomRealization: "<<b.RandomRealization;
		cout<<endl;

		if(b.NinDiameter!=0 && b.NinLength==0) {
			os<<"Fine Mesh: "<<b.NinDiameter<<"\tCells in each throat Diameter."<<endl;
		} else if(b.NinDiameter==0 && b.NinLength!=0) {
			os<<"coarse Mesh: "<<b.NinLength<<"\tCells in each throat Length."<<endl;
		} else {
			os<<"Something is wrong with the PorousOpen side setting."<<endl;
		}
	}

	return os;
}

void block_c::nondimensionalize(const physics_c & phy) {
	for(size_t i=0; i<B.size(); ++i) {
		B[i].u/=cd_c::phy.RefVelocity;
		B[i].v/=cd_c::phy.RefVelocity;
		B[i].p/=cd_c::phy.RefPressure;
		B[i].c/=cd_c::phy.RefConcentration;
	}

	avgTLength/=cd_c::phy.RefLength;
	for(size_t i=0; i<avgTDiameter.size(); ++i) {
		avgTDiameter[i]/=cd_c::phy.RefLength;
	}
}

