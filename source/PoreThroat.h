/*
 * PoreThroat.h
 *
 *  Created on: Apr 15, 2014
 *      Author: david
 */

#ifndef PORETHROAT_H_
#define PORETHROAT_H_

#include "Basic.h"
#include "Physics.h"
#include "BoundaryCondition.h"
#include "Geometry.h"
#include "CalculationDomain.h"

/*----------------------------------------------------------------------------------------------------

                                        pore_t class definition

----------------------------------------------------------------------------------------------------*/
class pore_c:public cPoint_c {
public:
	static std::vector<index_t> openPI;

	numeric_t RM, RI, A;
	numeric_t PhiF, PhiG, PhiLFI, PhiFGI;

	std::vector< numeric_t > p, f, r;

	std::vector<  index_t  > CI;
	std::vector<  index_t  > PI;
	std::vector<  index_t  > TI;
	std::vector<   cell_c *> NP;
	std::vector<   pore_c *> PP;
	std::vector< throat_c *> TP;

	std::vector<numeric_t  > CA;//Cell Crosssection Area, corresponding to NP
	std::vector<numeric_t  > CD;//Cell Distance, corresponding to NP

	numeric_t aP, b;
	std::vector<numeric_t> a;
	std::vector<numeric_t> ap;



	bool labelled;
	index_t clusterLabel;

	position_e position;
	 sStatus_e sStatus;
	vpStatus_e vpStatus;

	numeric_t vaporPressure;

	std::vector<index_t  > PtI;
	std::vector<point_c *> PtP;

    pore_c(const numeric_t & xv=0, const numeric_t & yv=0
    	 , const cellType_e & tv=pore
		 , const numeric_t & cv=0, const numeric_t & fv=0
		 , const bool & lbd=false, const index_t & clv=voidIndex
		 , const position_e & posv=liquidSide, const vpStatus_e & vpsv=none, const sStatus_e & ssv=liquid, const numeric_t & vpv=Psat)
	: cPoint_c(xv, yv, tv)
	, RM(0), RI(0), A(0)
	, PhiF(0), PhiG(0), PhiLFI(0), PhiFGI(0)
	, p(HowManyCopy, cv), f(HowManyCopy, fv), r(HowManyCopy, 0)
	, CI(0), PI(0), TI(0)
	, NP(0), PP(0), TP(0)
	, aP(0), b(0), a(0), ap(0)
    , labelled(lbd), clusterLabel(clv), position(posv), sStatus(ssv), vpStatus(vpsv), vaporPressure(vpv)
	, PtI(0), PtP(0) {}

	inline bool operator == (const pore_c & otherP) {return (x==otherP.x) && (y==otherP.y);}
	inline bool operator != (const pore_c & otherP) {return (x!=otherP.x) || (y!=otherP.y);}

    bool setPosition(const position_e & pos);
    bool setLabel(const size_t & CL);
    bool removeLabel();
    bool isEmpty();
//    inline numeric_t lnp() {return log(Pa-vaporPressure);}
//    bool lnp2vp(numeric_t v);
    bool updateStatus();
    bool isMeniscus();
    numeric_t MassOutFlux(const size_t & cli);
    numeric_t MassToCell(const size_t & cli);

    numeric_t C1Phi(const size_t & cli);//Create a function: C=C1Phi*f[cli]+C2Phi, in which C1Phi=
    numeric_t C2Phi(const size_t & cli);//Create a function: C=C1Phi*f[cli]+C2Phi, in which C2Phi=

    bool updatef2pr(const std::bitset<HowManyCopy> & hmc);
    bool updatepr2f(const std::bitset<HowManyCopy> & hmc);

    void resetCoefficient();
    void calculateCoefficient(const size_t & cli, const std::bitset<EqComponent> & ec);
    void calculatePhiCoefficient(const size_t & cli, const std::bitset<EqComponent> & ec);
};

/*----------------------------------------------------------------------------------------------------

                                        throat_t class definition

----------------------------------------------------------------------------------------------------*/
class throat_c:public cd_c {
public:
	std::vector<index_t  > PI;
	std::vector< pore_c *> PP;

	int N; //Polygon, N=3, 4, 5, 6...; N=0 for circle
	numeric_t Alpha, ThetaC, Kai;
	numeric_t RI, RM, r;
	numeric_t Psi1, Psi2, Psi3, Psi4, Psi5, Beta, Kappa;
	numeric_t L, A, V, M, MFull;
	numeric_t PhiF, PhiG;

	 sStatus_e sStatus;
	position_e position;

	index_t clusterLabel;
	bool labelled;

	std::vector<index_t  > PtI;
	std::vector<point_c *> PtP;

	throat_c(const index_t & pi0=voidIndex, const index_t & pi1=voidIndex
		   , const int & nv=0
		   , const numeric_t & riv=0
		   , const numeric_t & lv=0
		   , const sStatus_e & ssv=liquid, const position_e & posv=liquidSide, const index_t & clv=voidIndex, const bool & lbdv=false)
	: PI(2, voidIndex), PP(2, nullptr)
	, N(nv)
	, Alpha(0), ThetaC(0), Kai(0)
	, RI(riv), RM(0), r(0)
	, Psi1(0), Psi2(0), Psi3(0), Psi4(0), Psi5(0), Beta(0), Kappa(0)
	, L(lv), A(0), V(0), M(0), MFull(0)
	, PhiF(0), PhiG(0)
	, sStatus(ssv), position(posv)
	, clusterLabel(clv), labelled(lbdv)
	, PtI(0), PtP(0) {
		if(N<3) {
			A=pi*RI*RI;
		} else {
			Alpha=pi/2-pi/N;
			ThetaC=pi/N;
			Psi1=cos(Alpha+cd_c::phy.Wet.contactAngle)
				*(cos(Alpha+cd_c::phy.Wet.contactAngle)+sin(Alpha+cd_c::phy.Wet.contactAngle)*tan(Alpha));
			Psi2=1-cd_c::phy.Wet.contactAngle/(pi/2-Alpha);
			Psi3=cos(Alpha+cd_c::phy.Wet.contactAngle)/cos(Alpha);
			Psi4=cos(cd_c::phy.Wet.contactAngle)*cos(Alpha+cd_c::phy.Wet.contactAngle)/sin(Alpha)-pi/2+Alpha+cd_c::phy.Wet.contactAngle;
			Psi5=(pi/2-Alpha)*tan(Alpha);
			Kai=cos(cd_c::phy.Wet.contactAngle)
			   +sqrt((ThetaC-cd_c::phy.Wet.contactAngle+sin(cd_c::phy.Wet.contactAngle)*cos(cd_c::phy.Wet.contactAngle))/tan(ThetaC));
			RM=RI/Kai;
			r=RM;
			Beta=12*(Psi1-Psi5*Psi2)*pow(sin(Alpha)*(1-Psi5)*(Psi3-(1-Psi5)*cd_c::phy.RC/r)/(1-sin(Alpha))/Psi5, 2)
				/pow(Psi1-Psi5*Psi2-(1-Psi5)*pow(cd_c::phy.RC/r, 2), 3);
			Kappa=Psi4/Beta;
			A=pi*RI*RI+N*(RI*RI-cd_c::phy.RC*cd_c::phy.RC)*(1/tan(Alpha)-pi/2+Alpha);

			PhiF=N*Kappa*cd_c::phy.LiquidConcentration/3/cd_c::phy.Ca;
			PhiG=A;
		}

		V=A*L;
		M=V*cd_c::phy.LiquidConcentration;
		MFull=M;
		PI[0]=pi0; PI[1]=pi1;
	}

	inline bool operator == (const throat_c & otherT) const {
		return ( ((( *PI.begin()==*otherT.PI.begin() ) && (*(PI.end()-1)==*(otherT.PI.end()-1))) || ((*PI.begin()==*(otherT.PI.end()-1)) && (*(PI.end()-1)==*otherT.PI.begin())))
				&& *PI.begin()!=voidIndex && *(PI.end()-1)!=voidIndex && *otherT.PI.begin()!=voidIndex && *(otherT.PI.end()-1)!=voidIndex );
	}
	inline bool operator != (const throat_c & otherT) const {
		return !(((PI[0]==otherT.PI[0]) && (PI[1]==otherT.PI[1])) ||
				 ((PI[0]==otherT.PI[1]) && (PI[1]==otherT.PI[0])));
	}

	bool setPosition(const position_e & pos);
    inline index_t otherP(const index_t & oneP) {return oneP==PI[0]?PI[1]:PI[0];}
    bool setLabel(const index_t & CN);
    bool removeLabel();

    inline numeric_t Saturation() {return M/MFull;}
    inline numeric_t Concentration() {return M/V;}

    inline bool isDry()   {return Concentration()<=(1+cd_c::phy.RelativeConvergeCriteria)*cd_c::phy.SaturatedConcentration  ?true:false;}
    inline bool isEmpty() {return Concentration()<=(1+cd_c::phy.RelativeConvergeCriteria)*cd_c::phy.EnvironmentConcentration?true:false;}

    bool isMeniscus();
    bool updateStatus();
    numeric_t MassOutFlux(const size_t & cli);

    numeric_t S(const index_t & p1, const index_t & p2);
};

/*----------------------------------------------------------------------------------------------------

                                        cluster_t class definition

----------------------------------------------------------------------------------------------------*/
class cluster_c:public cd_c {
public:
	std::vector<index_t> PI;
	std::vector<index_t> TI;
	std::vector<index_t> maxThroat;

	numeric_t MassOutFlux;
	numeric_t timeToEmptyMax;
	numeric_t MassOut;

	cluster_c()
	: PI(0), TI(0)
	, maxThroat(0), MassOutFlux(0), timeToEmptyMax(0), MassOut(0) {}

	numeric_t mass();

	bool clear();
	bool calcMassOutFlux(const size_t & cli);
	bool findMaxThroat();
	numeric_t maxMass();
	inline numeric_t DividedMassOutFlux();
	bool calctimeToEmptyMax();
	bool drying(const numeric_t & time, const bool GotEmptyItNow);
};


//bool createThroatWithPores(const pore_c &, const pore_c &, const numeric_t &);


#endif /* PORETHROAT_H_ */
