/*
 * Geometry.h
 *
 *  Created on: Apr 15, 2014
 *      Author: david
 */

#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#include "Basic.h"
#include "Physics.h"
#include "BoundaryCondition.h"

class      point_c;
class    segment_c;
class numberAxis_c;
class       line_c;
class      block_c;

///*----------------------------------------------------------------------------------------------------
//
//                                        shape_c class
//
//----------------------------------------------------------------------------------------------------*/
//class shape_c {
//public:
//	numeric_t A, R; // Area Total, Area Bulk, Radius of Corner
//	shape_c(const numeric_t & r=0): R(r) {A=0;}
//
//	virtual ~shape_c()=default;
//
//};
//
//class circle_c:shape_c {
//public:
//	circle_c(const numeric_t & r=0): R(r) {A=pi*R*R;}
//	inline void calArea() {A=pi*R*R;}
//};
//
//class polygon_c:shape_c {
//public:
//	unsigned int N;
//	numeric_t RC, AB;
//
//};

/*----------------------------------------------------------------------------------------------------

                                        point_c class

----------------------------------------------------------------------------------------------------*/
class point_c {
public:
	numeric_t x, y;

	point_c(const numeric_t & xVal=0, const numeric_t & yVal=0):x(xVal), y(yVal) {}
	virtual ~point_c()=default;

	inline void offset(const numeric_t & dx, const numeric_t & dy) {x+=dx; y+=dy;}

	inline point_c offsetCopy(const numeric_t & dx, const numeric_t & dy) {return point_c(x+dx, y+dy);}

	inline bool operator == (const point_c & otherP) {return (x==otherP.x) && (y==otherP.y);}
	inline bool operator != (const point_c & otherP) {return (x!=otherP.x) || (y!=otherP.y);}

	inline numeric_t distance(const point_c & p) {return sqrt((x-p.x)*(x-p.x)+(y-p.y)*(y-p.y));}
	inline numeric_t distanceX(const point_c & p) {return fabs(x-p.x);}
	inline numeric_t distanceY(const point_c & p) {return fabs(y-p.y);}

	inline void nondimensionalize(const numeric_t & ref) {x/=ref; y/=ref;}

	friend std::istream & operator >> (std::istream & is, point_c & p);
	friend std::ostream & operator << (std::ostream & os, const point_c & p);
};

/*----------------------------------------------------------------------------------------------------

                                        segment_c class
                                                    -the preparation for numberAxis_c & line_c
----------------------------------------------------------------------------------------------------*/
class segment_c {
public:
	size_t SN;
	numeric_t GS, GL; //Growth Small/Large Side
	numeric_t LS, LL; //Length Small/Large Side

	segment_c(const size_t & n=defaultSegment
			, const numeric_t & gss=defaultGrowth, const numeric_t & gls=defaultGrowth
			, const numeric_t & lss=0            , const numeric_t & lls=0            )
	: SN(n), GS(gss), GL(gls), LS(lss), LL(lls) {}

	std::vector<numeric_t> fillBetween(const numeric_t & beg, const numeric_t & end);
	std::vector<  point_c> fillBetween(const point_c & p0, const point_c & p1);

	inline bool operator == (const segment_c & s) const {return (SN==s.SN) && (GS==s.GS) && (GL==s.GL) && (LS==s.LS) && (LL==s.LL);}
	inline bool operator != (const segment_c & s) const {return (SN!=s.SN) || (GS!=s.GS) || (GL!=s.GL) || (LS!=s.LS) || (LL!=s.LL);}
};

/*----------------------------------------------------------------------------------------------------

                                        numberAxis_c class
                                                    -the preparation for line_c
----------------------------------------------------------------------------------------------------*/
class numberAxis_c { //a string of numbers from small to large
public:
	std::vector<numeric_t> num;

	std::vector<segment_c> S;

	std::vector<     BC_c> B; // boundary condition of the segment

	numberAxis_c(const numeric_t & beg=0, const numeric_t & end=0
			   , const segment_c & s=segment_c(), const BC_c & b=BC_c())
	: num(2, beg), S(1, s), B(1, b) {*(num.end()-1)=end;}

	bool clear(); //Clear all vectors inside object
	bool clearMiddle();

	bool initialize(); //Initialize it to the status created by default constructor

	bool insertNumber(const numeric_t & nb, const segment_c & s1, const segment_c & s2
			        , const BC_c & b1, const BC_c & b2);

	bool insertSegment(const numeric_t & nb0, const numeric_t & nb1 //user is responsible for making nb0<=nb1
			         , const segment_c & s1, const segment_c & s2, const segment_c & s3
				     , const BC_c & b1, const BC_c & b2, const BC_c & b3);

	bool fill();
	bool filled(); //Check if the number axis is filled or not
};

/*----------------------------------------------------------------------------------------------------

                                        line_c class

----------------------------------------------------------------------------------------------------*/
class line_c {
public:
	coordinate_e XorY;

	std::vector<  point_c> Pt;

	std::vector<segment_c> S;

	std::vector<     BC_c> B; // boundary condition of the segment

	std::vector<  index_t> PtI; // the corresponding point's position in global vector<point_c>

	std::vector<  index_t> PELI; //parallel equal line index

	std::vector<  index_t> SCI; //Segment Cell Index in global cells vector

	line_c(const point_c & beg=point_c(), const point_c & end=point_c()
		 , const size_t & begi=voidIndex, const size_t & endi=voidIndex
		 , const coordinate_e & xy=origin
		 , const segment_c & s=segment_c()
		 , const BC_c & b=BC_c())
	: XorY(xy), Pt(2, beg), S(1, s), B(1, b), PtI(2, begi), PELI(0), SCI(0) {
		*(Pt.end()-1)=end; *(PtI.end()-1)=endi;
		if     (Pt[0].x==Pt[Pt.size()-1].x && Pt[0].y!=Pt[Pt.size()-1].y) XorY=yAxis;
		else if(Pt[0].x!=Pt[Pt.size()-1].x && Pt[0].y==Pt[Pt.size()-1].y) XorY=xAxis;
		else;
	}

    bool clear();
    bool initialize();

    numberAxis_c extractNumberAxis(const coordinate_e & ce);
    bool plus2NumberAxis(const numberAxis_c & xa, const numberAxis_c & ya);
    bool NumberAxisPlusNumber(const numberAxis_c & axis, const numeric_t & num);

    inline numeric_t length() {return Pt[0].distance(Pt[Pt.size()-1]);}
    inline numeric_t slope() {return (Pt[Pt.size()-1].y-Pt[0].y)/(Pt[Pt.size()-1].x-Pt[0].x);}
    numeric_t distance(const line_c & l);

	bool insertPoint(const point_c & bp
                   , const segment_c & s1, const segment_c & s2
			       , const BC_c & b1, const BC_c & b2);

	bool insertSegment(const point_c & bp0, const point_c & bp1
                     , const segment_c & s1, const segment_c & s2, const segment_c & s3
			         , const BC_c & b1, const BC_c & b2, const BC_c & b3);

	bool fill();
	bool filled();

	bool operator == (const line_c & L) {return (Pt[0]==L.Pt[0] && Pt[Pt.size()-1]==L.Pt[L.Pt.size()-1]) ||
			                                       (PtI[0]==L.PtI[0] && PtI[PtI.size()-1]==L.PtI[L.PtI.size()-1]);}
	bool operator != (const line_c & L) {return !((Pt[0]==L.Pt[0] && Pt[Pt.size()-1]==L.Pt[L.Pt.size()-1]) ||
			                                       (PtI[0]==L.PtI[0] && PtI[PtI.size()-1]==L.PtI[L.PtI.size()-1]));}

	coordinate_e ParaEqTo(const line_c & L);
	bool copyParaEq(const line_c & L);

	friend std::ostream & operator << (std::ostream & os, const line_c & eg);

};

/*----------------------------------------------------------------------------------------------------

                                        block_c class

----------------------------------------------------------------------------------------------------*/
class block_c {
public:
	std::vector<index_t> CnI; //Corner Points Index
	blockType_e            T; //Block Type
	std::vector<BC_c>      B; //Boundary Conditions
	size_t            Nx, Ny; //For FlowField NxNy are about Points & For PorousMedia NxNy are about Pores
	std::vector<numeric_t> G; //Growth Rate

	porosityType_e porosityType;
	numeric_t PPConnectionRate, PPCrossRate;
    numeric_t                avgTLength, TLVarianceRate;//ThroatLengthVarianceRate, PoretoPoreConnectionRate
	std::vector<size_t>    PolyN;// 0 circle, 3 triangle, 4 square
	std::vector<numeric_t> avgTDiameter, TDVarianceRate;//ThroatDiameterVarianceRate
	std::vector<numeric_t> TDisplayScale;
	int RandomRealization;//This number can only be positive integer

	size_t NinDiameter, NinLength; // segments in diameter/length, 0 means no divide

	std::vector<index_t>  NI; //neighbor Blocks Index

	std::vector<index_t> LnI; //Line Index

	numeric_t Lx, Ly;
	numeric_t Ox, Oy; //Coordinates of block center or SouthWest Corner(For PorousMedia only)

	std::vector<std::vector<index_t> > PtI, CI, UI, VI;

	std::vector<std::vector<index_t> > PI;
	std::vector<index_t> TI;//Pore Index, Throat Index

	block_c()
	: CnI(sideNumber, voidIndex), T(ghost), B(sideNumber, BC_c()), Nx(1), Ny(1), G(sideNumber, defaultGrowth)
	, porosityType(singlePorosity)
	, PPConnectionRate(1), PPCrossRate(1)
	, avgTLength(ThroatAverageLength), TLVarianceRate(0.1)
	, PolyN(1, 0), avgTDiameter(1, ThroatAverageDiameter), TDVarianceRate(1, 0.1), TDisplayScale(1, 1)
	, RandomRealization(1)
	, NinDiameter(1), NinLength(0)
	, NI(nodeNumber, voidIndex), LnI(sideNumber, voidIndex)
	, Lx(0), Ly(0), Ox(0), Oy(0)
	, PtI(0), CI(0), UI(0), VI(0), PI(0), TI(0) {}

	friend std::istream & operator >> (std::istream & is, block_c & b);
    friend std::ostream & operator << (std::ostream & os, const block_c & b);

    void nondimensionalize(const physics_c & phy);
};



#endif /* GEOMETRY_H_ */
