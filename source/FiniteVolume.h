/*
 * FiniteVolume.h
 *
 *  Created on: Apr 15, 2014
 *      Author: david
 */

#ifndef FINITEVOLUME_H_
#define FINITEVOLUME_H_

#include "Basic.h"
#include "Physics.h"
#include "BoundaryCondition.h"
#include "Geometry.h"
#include "CalculationDomain.h"

/*----------------------------------------------------------------------------------------------------
 *
 *                                        cell_c class definition
 *
  ----------------------------------------------------------------------------------------------------*/
class cell_c:public cPoint_c {
public:
	std::vector<numeric_t> D;//distance from cell center to neighbor cell center: dE, dN, dW, dS
	std::vector<numeric_t> d;//distance from cell center to the side wall: de, dn, dw, ds

	numeric_t Ax, Ay, V;//Area & Volume
	std::vector<numeric_t> p;

	std::vector<index_t  > NI;
	std::vector<index_t  > FI;

	std::vector< cell_c *> NP;
	std::vector< cell_c *> FP;

	CBC_c FB;//FB: Flow field Boundary; EB: Evaporation Boundary
	std::bitset<sideNumber> isWall; //Neighbor Available
	bool nearBC;//close to Boundary Condition, which means have only one level neighbors

	std::vector<numeric_t> FA;//FAe, FAn, FAw, FAs
	std::vector<numeric_t> DA;//DAe, DAn, DAw, DAs

	numeric_t aP; //Coefficient for itself, left side value
	std::vector<index_t  > CI;//Index       of related
	std::vector<numeric_t  > a ;//Coefficient of related, left side value
	numeric_t b; //Right Side Value
	numeric_t r;//r=A/aP

	cell_c(const numeric_t & xV=0, const numeric_t & yV=0
		 , const cellType_e & tv=fCell
		 , const numeric_t & DEV=0, const numeric_t & DNV=0, const numeric_t & DWV=0, const numeric_t & DSV=0
         , const numeric_t & deV=0, const numeric_t & dnV=0, const numeric_t & dwV=0, const numeric_t & dsV=0
		 , const numeric_t & pV=0 , const CBC_c & fbV=CBC_c(), const std::bitset<sideNumber> isWallV=0, const bool & nearbc=false)
	: cPoint_c(xV, yV, tv)
	, D(sideNumber, 0), d(sideNumber, 0), Ax((dnV+dsV)*cd_c::dz), Ay((deV+dwV)*cd_c::dz), V((dnV+dsV)*(deV+dwV)*cd_c::dz), p(HowManyCopy, pV)
	, NI(nodeNumber, voidIndex), FI(nodeNumber, voidIndex)
	, NP(nodeNumber, nullptr), FP(nodeNumber, nullptr)
	, FB(fbV), isWall(isWallV), nearBC(nearbc)
	, FA(sideNumber, 0), DA(sideNumber, 0)
    , aP(0), CI(0), a(0), b(0), r(0) {
		d[e]=deV; d[n]=dnV; d[w]=dwV; d[s]=dsV;
		D[E]=DEV; D[N]=DNV; D[W]=DWV; D[S]=DSV;
	}

	virtual ~cell_c()=default;

	inline numeric_t dx() {return  d[e]+d[w];}
	inline numeric_t dy() {return  d[n]+d[s];}
	inline numeric_t Ae() {return  Ax;}
	inline numeric_t Aw() {return -Ax;}
	inline numeric_t An() {return  Ay;}
	inline numeric_t As() {return -Ay;}

    virtual void calculateFD(const size_t & uvli)=0;

	        void ConvectionItemUseHybridScheme (const size_t & uvli, const size_t & pli );//calculate aP, a, b, d with Hybrid Scheme
	        void ConvectionItemUseQUICKScheme  (const size_t & uvli, const size_t & pli );//calculate aP, a, b, d with QUICK Scheme
	virtual void DiffussionItemUseCentralScheme(const size_t & uvli, const size_t & pcli)=0;//calculate aP, a, b, d with Central Scheme

    virtual void resetCoefficient()=0;
    virtual void calculateCoefficient(const size_t & uvli, const size_t & pli, const std::bitset<EqComponent> & ec)=0;//From Left to Right: Convection, Diffusion, DeltaPA/DeltaVT
};

/*----------------------------------------------------------------------------------------------------
 *
 *                                        uvCell_c class definition
 *
  ----------------------------------------------------------------------------------------------------*/
class uvCell_c:public cell_c {
public:
	static numeric_t MassIn, MassOut, MassRatio;//Mass Balance Ratio=MassIn/MassOut
	static std::vector< index_t  > InletUI, OutletUI, InletVI, OutletVI;
//	static std::vector<uvCell_c *> InletUVP, OutletUVP;

	numeric_t DeltaPA;
	numeric_t NormalA;//Inlet or Outlet Area; For Outlet: + for e n, - for w s; For Inlet: + for w s, - for e n

	uvCell_c(const numeric_t & xV=0, const numeric_t & yV=0
		   , const cellType_e & tv=uvCell
		   , const numeric_t & DEV=0, const numeric_t & DNV=0, const numeric_t & DWV=0, const numeric_t & DSV=0
	       , const numeric_t & deV=0, const numeric_t & dnV=0, const numeric_t & dwV=0, const numeric_t & dsV=0
		   , const numeric_t & pV=0 , const CBC_c & fbV=CBC_c(), const std::bitset<sideNumber> isWallV=0, const bool & nearbc=false
		   , const int & normaldirection=0)
	: cell_c(xV, yV, tv, deV, dnV, dwV, dsV, DEV, DNV, DWV, DSV, pV, fbV, isWallV, nearbc)
	, DeltaPA(0), NormalA(normaldirection) {}

	virtual ~uvCell_c()=default;

    virtual void calculateFD(const size_t & uvli)=0;//calculate FA and DA

	        void DiffussionItemUseCentralScheme(const size_t & uvli, const size_t & pli);
	virtual void calculateDeltaPA(const size_t & pli)=0;
	virtual void calculateRorD()=0;

	        void resetCoefficient();
	        void calculateCoefficient(const size_t & uvli, const size_t & pli, const std::bitset<EqComponent> & ec);

};

//----------------------------------------------------------------------------------------------------
class uCell_c:public uvCell_c {
public:
	uCell_c(const numeric_t & xV=0, const numeric_t & yV=0
		  , const cellType_e & tv=uCell
		  , const numeric_t & DEV=0, const numeric_t & DNV=0, const numeric_t & DWV=0, const numeric_t & DSV=0
		  , const numeric_t & deV=0, const numeric_t & dnV=0, const numeric_t & dwV=0, const numeric_t & dsV=0
		  , const numeric_t & pV=0 , const CBC_c & fbV=CBC_c(), const std::bitset<sideNumber> isWallV=0, const bool & nearbc=false
		  , const int & normaldirection=0)
	: uvCell_c(xV, yV, tv, deV, dnV, dwV, dsV, DEV, DNV, DWV, DSV, pV, fbV, isWallV, nearbc, normaldirection) {
		NormalA*=Ax;
	}

	void calculateFD(const size_t & uvli);

	void calculateDeltaPA(const size_t & pli);
	void calculateRorD();
};

//----------------------------------------------------------------------------------------------------
class vCell_c:public uvCell_c {
public:
	vCell_c(const numeric_t & xV=0, const numeric_t & yV=0
		  , const cellType_e & tv=vCell
		  , const numeric_t & DEV=0, const numeric_t & DNV=0, const numeric_t & DWV=0, const numeric_t & DSV=0
		  , const numeric_t & deV=0, const numeric_t & dnV=0, const numeric_t & dwV=0, const numeric_t & dsV=0
		  , const numeric_t & pV=0 , const CBC_c & fbV=CBC_c(), const std::bitset<sideNumber> isWallV=0, const bool & nearbc=false
		  , const int & normaldirection=0)
	: uvCell_c(xV, yV, tv, deV, dnV, dwV, dsV, DEV, DNV, DWV, DSV, pV, fbV, isWallV, nearbc, normaldirection) {
		NormalA*=Ay;
	}

	void calculateFD(const size_t & uvli);

	void calculateDeltaPA(const size_t & pli);
	void calculateRorD();
};

/*----------------------------------------------------------------------------------------------------
 *
 *                                        cCell_c class definition
 *
  ----------------------------------------------------------------------------------------------------*/
class cCell_c:public cell_c {
public:
	static numeric_t dpdd;

	CBC_c EB;//FB: Flowfield Boundary; EB: Evaporation Boundary

//	std::vector<numeric_t  > L;

	std::vector<  index_t  > PI;
	std::vector<   pore_c *> PP;
	std::vector<  index_t  > TI;
	std::vector< throat_c *> TP;

	std::vector<numeric_t  > ap;//ap: apore: coefficient for pore
	numeric_t DeltaVT;


	cCell_c(const numeric_t & xV=0, const numeric_t & yV=0
		  , const cellType_e & tv=cCell
		  , const numeric_t & DEV=0, const numeric_t & DNV=0, const numeric_t & DWV=0, const numeric_t & DSV=0
	      , const numeric_t & deV=0, const numeric_t & dnV=0, const numeric_t & dwV=0, const numeric_t & dsV=0
		  , const numeric_t & pV=0 , const CBC_c & fbV=CBC_c(), const std::bitset<sideNumber> isWallV=0, const bool & nearbc=false
		  , const numeric_t & cV=0,  const CBC_c & ebV=CBC_c())
	: cell_c(xV, yV, tv, deV, dnV, dwV, dsV, DEV, DNV, DWV, DSV, pV, fbV, isWallV, nearbc)
	, EB(ebV), PI(0), PP(0), TI(0), TP(0), ap(0), DeltaVT(0) {
		p.resize(2*HowManyCopy, 0);
		for(size_t i=0; i<HowManyCopy; ++i) {
			p[i]=cV;
		}
		for(size_t i=HowManyCopy; i<2*HowManyCopy; ++i) {
			p[i]=pV;
		}
	}
	~cCell_c()=default;

    void resetPCoefficient();
    void calculatePCoefficient(const size_t & uvli, const size_t & pli, const std::bitset<VariableCount> & vc);

	void calculateFD(const size_t & uvli);//calculate FA and DA items in Evaporation Equation

	void DiffussionItemUseCentralScheme(const size_t & uvli, const size_t & cli);
	void calculateDeltaVT(const size_t & cli);

	void resetCoefficient();
    void calculateCoefficient(const size_t & uvli, const size_t & pli, const std::bitset<EqComponent> & ec);

};

/*----------------------------------------------------------------------------------------------------
 *
 *                                        fCell_c class definition
 *
  ----------------------------------------------------------------------------------------------------*/
class fCell_c:public cell_c {
public:
	fCell_c(const numeric_t & xV=0, const numeric_t & yV=0
		  , const cellType_e & tv=fCell
		  , const numeric_t & DEV=0, const numeric_t & DNV=0, const numeric_t & DWV=0, const numeric_t & DSV=0
	      , const numeric_t & deV=0, const numeric_t & dnV=0, const numeric_t & dwV=0, const numeric_t & dsV=0
		  , const numeric_t & pV=0 , const CBC_c & fbV=CBC_c(), const std::bitset<sideNumber> isWallV=0, const bool & nearbc=false
		  , const int & normaldirection=0)
	: cell_c(xV, yV, tv, deV, dnV, dwV, dsV, DEV, DNV, DWV, DSV, pV, fbV, isWallV, nearbc) {}

	fCell_c(const point_c & pt=point_c())
	: cell_c(pt.x, pt.y, fCell, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, CBC_c(), 0) {}

	virtual ~fCell_c()=default;

	numeric_t     U() {return p[0];}
	numeric_t     V() {return p[1];}
	numeric_t     P() {return p[2];}
	numeric_t     C() {return p[3];}

	numeric_t RealX() {return x   *cd_c::phy.RefLength;}
	numeric_t RealY() {return y   *cd_c::phy.RefLength;}
	numeric_t RealU() {return p[0]*cd_c::phy.RefVelocity;}
	numeric_t RealV() {return p[1]*cd_c::phy.RefVelocity;}
	numeric_t RealP() {return p[2]*cd_c::phy.RefPressure;}
	numeric_t RealC() {return p[3]*cd_c::phy.RefConcentration;}

    void calculateFD(const size_t & uvli) {}

	void DiffussionItemUseCentralScheme(const size_t & uvli, const size_t & cli) {}
    void resetCoefficient() {}
    void calculateCoefficient(const size_t & uvli, const size_t & pli, const std::bitset<EqComponent> & ec) {}
};


#endif /* FINITEVOLUME_H_ */
