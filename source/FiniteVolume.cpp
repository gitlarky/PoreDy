/*
 * FiniteVolume.cpp
 *
 *  Created on: Apr 15, 2014
 *      Author: david
 */

#include "FiniteVolume.h"
#include "PoreThroat.h"

using namespace std;

/*----------------------------------------------------------------------------------------------------
 *
 *                                        Initialize Static Members
 *
  ----------------------------------------------------------------------------------------------------*/
numeric_t uvCell_c::MassIn   =0;
numeric_t uvCell_c::MassOut  =0;
numeric_t uvCell_c::MassRatio=1;

vector<index_t> uvCell_c:: InletUI(0);
vector<index_t> uvCell_c:: InletVI(0);
vector<index_t> uvCell_c::OutletUI(0);
vector<index_t> uvCell_c::OutletVI(0);

numeric_t cCell_c::dpdd=0;

/*====================================================================================================
 *
 *             Flowfield related cells member functions
 *
  ====================================================================================================*/
/*----------------------------------------------------------------------------------------------------
 *
 *             For calculate continuity equation, to calculate pressure correction
 *
  ----------------------------------------------------------------------------------------------------*/
void cCell_c::resetPCoefficient() {
	aP=0; b=0;
	a.clear(); CI.clear();
}

//----------------------------------------------------------------------------------------------------
void cCell_c::calculatePCoefficient(const size_t & uvli, const size_t & pli, const std::bitset<VariableCount> & vc) {
	resetPCoefficient();
	a.resize(sideNumber, 0); CI.resize(sideNumber, voidIndex);
	for(size_t i=0; i<CI.size(); ++i) {
		CI[i]=NI[i];
	}

	if(FB.T==NoBC) {
		if(FP[e]->FB.T==NoBC) {
			b   -= FP[e]->p[uvli]*Ax;
			a[E] =-FP[e]->r*Ax;
			aP  -= a[E];
		} else {
			b   -= FP[e]->p[uvli]*Ax;
			a[E] = 0;
			aP  -= a[E];
		}

		if(FP[w]->FB.T==NoBC) {
			b   += FP[w]->p[uvli]*Ax;
			a[W] =-FP[w]->r*Ax;
			aP  -= a[W];
		} else {
			b   += FP[w]->p[uvli]*Ax;
			a[W] = 0;
			aP  -= a[W];
		}

		if(FP[n]->FB.T==NoBC) {
			b   -= FP[n]->p[uvli]*Ay;
			a[N] =-FP[n]->r*Ay;
			aP  -= a[N];
		} else {
			b   -= FP[n]->p[uvli]*Ay;
			a[N] = 0;
			aP  -= a[N];
		}

		if(FP[s]->FB.T==NoBC) {
			b   += FP[s]->p[uvli]*Ay;
			a[S] =-FP[s]->r*Ay;
			aP  -= a[S];
		} else {
			b   += FP[s]->p[uvli]*Ay;
			a[S] = 0;
			aP  -= a[S];
		}
	} else {
		if(vc==4) {
			aP=FB.aP;
			b =FB.b ;
		} else {
			aP=1;
			b =0;
		}
	}

	for(size_t i=0; i<CI.size(); ++i) {
		if(CI[i]==voidIndex) {
			a.erase(a.begin()+i); CI.erase(CI.begin()+i);
			--i;
		}
	}
}

/*----------------------------------------------------------------------------------------------------
 *
 *             For calculate momentum equation
 *
  ----------------------------------------------------------------------------------------------------*/
void uvCell_c::resetCoefficient() {
	aP=0; b=0; r=0; DeltaPA=0;
	a.clear(); CI.clear();
}

//----------------------------------------------------------------------------------------------------
void uvCell_c::calculateCoefficient(const size_t & uvli, const size_t & pli, const std::bitset<EqComponent> & ecbs) {
	if(ecbs[convectionItem] && ecbs[diffussionItem] && ecbs[RorD]) {
		resetCoefficient();

		if(FB.T==NoBC) {
			calculateFD(uvli);

			switch (cd_c::phy.SchemeFlowfield) {
			case Hybrid :
				a.resize(sideNumber, 0); CI.resize(sideNumber, voidIndex);
				ConvectionItemUseHybridScheme(uvli, pli);
				break;
			case HayaseQUICK :
				a.resize(NI.size(), 0); CI.resize(NI.size(), voidIndex);
				ConvectionItemUseQUICKScheme(uvli, pli);
//				if(nearBC) 	{
//					a.resize(sideNumber, 0); CI.resize(sideNumber, voidIndex);
//					ConvectionItemUseHybridScheme(uvli, pli);
//				} else {
//					a.resize(NI.size(), 0); CI.resize(NI.size(), voidIndex);
//					ConvectionItemUseQUICKScheme(uvli, pli);
//				}
				break;
			case Upwind: case PowerLaw: case TVD: case QUICK:
				break;
			}

			DiffussionItemUseCentralScheme(uvli, pli);

			calculateRorD();

			for(size_t i=0; i<CI.size(); ++i) {
				if(CI[i]==voidIndex) {
					a.erase(a.begin()+i); CI.erase(CI.begin()+i);
					--i;
				}
			}

			b +=aP*(1-cd_c::phy.UVRelaxFactor)/cd_c::phy.UVRelaxFactor*p[uvli];
			aP/=cd_c::phy.UVRelaxFactor;

		} else {
			aP     =FB.aP;
			b      =FB.b ;
			CI     =FB.I ;
			a      =FB.a ;
			r      =0    ;
		}
	}

	if(ecbs[deltaPAVT]) {
		if(FB.T==NoBC) {
			calculateDeltaPA(pli);
		} else {
			DeltaPA=0    ;
		}
	}

//	cout<<"aP="<<aP<<"\tb="<<b<<"\tDeltaPA="<<DeltaPA<<"\tr="<<r;
//	for(size_t i=0; i<CI.size(); ++i) {
//		cout<<"\tCI="<<CI[i]<<"\ta="<<a[i];
//	}
//	cout<<endl;
}

//----------------------------------------------------------------------------------------------------
void uvCell_c::DiffussionItemUseCentralScheme(const size_t & uvli, const size_t & pli) {
	a[E]-=isWall[e]?0:DA[e];
	a[W]-=isWall[w]?0:DA[w];
	a[N]-=isWall[n]?0:DA[n];
	a[S]-=isWall[s]?0:DA[s];

	aP  +=DA[e]+DA[w]+DA[n]+DA[s];


//	if(cd_c::phy.SchemeFlowfield==Hybrid) {// For Triditional Hybrid Method, the DA sometimes can be dropped
//		if(FA[e]>=2*DA[e] || FA[e]<=-2*DA[e]) {
//			a[E]-=0;
//		} else {
//			a[E]-=isWall[e]?0:DA[e]; aP+=DA[e];
//		}
//		if(FA[w]>=2*DA[w] || FA[w]<=-2*DA[w]) {
//			a[W]-=0;
//		} else {
//			a[W]-=isWall[w]?0:DA[w]; aP+=DA[w];
//		}
//		if(FA[n]>=2*DA[n] || FA[n]<=-2*DA[n]) {
//			a[N]-=0;
//		} else {
//			a[N]-=isWall[n]?0:DA[n]; aP+=DA[n];
//		}
//		if(FA[s]>=2*DA[s] || FA[s]<=-2*DA[s]) {
//			a[S]-=0;
//		} else {
//			a[S]-=isWall[s]?0:DA[s]; aP+=DA[s];
//		}
//	}

}

/*----------------------------------------------------------------------------------------------------
 *
 *             calculateD, which in this implementation is calculate r
 *
  ----------------------------------------------------------------------------------------------------*/
void uCell_c::calculateRorD() {
	r=Ax/aP;
}

void vCell_c::calculateRorD() {
	r=Ay/aP;
}

//----------------------------------------------------------------------------------------------------
void uCell_c::calculateDeltaPA(const size_t & pli) {
	DeltaPA=(FP[w]->p[pli]-FP[e]->p[pli])*Ax;
}

void vCell_c::calculateDeltaPA(const size_t & pli) {
	DeltaPA=(FP[s]->p[pli]-FP[n]->p[pli])*Ay;
}

/*----------------------------------------------------------------------------------------------------
 *
 *             calculate FA & DA for uCell_c & vCell_c
 *
  ----------------------------------------------------------------------------------------------------*/
void uCell_c::calculateFD(const size_t & uvli) {
	FA[e]=0.5*(p[uvli]+NP[E]->p[uvli])*Ax;
	FA[w]=0.5*(p[uvli]+NP[W]->p[uvli])*Ax;

	DA[e]=Ax/(phy.Re*D[E]);
	DA[w]=Ax/(phy.Re*D[W]);

	if(isWall[n]) {
		FA[n]=0;

		DA[n]=Ay/(phy.Re*d[n]);
	} else {
		FA[n]=0.5*(FP[ne]->p[uvli]+FP[nw]->p[uvli])*Ay;

		DA[n]=Ay/(phy.Re*D[N]);
	}
	if(isWall[s]) {
		FA[s]=0;

		DA[s]=Ay/(phy.Re*d[s]);
	} else {
		FA[s]=0.5*(FP[sw]->p[uvli]+FP[se]->p[uvli])*Ay;

		DA[s]=Ay/(phy.Re*D[S]);
	}

	numeric_t rou(cd_c::phy.Dry.density/cd_c::phy.RefDensity), miu(cd_c::phy.Dry.dynamicViscosity/cd_c::phy.RefViscosity);
	for(size_t d=e; d<=s; ++d) {
		FA[d]*=rou;
		DA[d]*=miu;
	}
}

//----------------------------------------------------------------------------------------------------
void vCell_c::calculateFD(const size_t & uvli) {
	FA[n]=0.5*(p[uvli]+NP[N]->p[uvli])*Ay;
	FA[s]=0.5*(p[uvli]+NP[S]->p[uvli])*Ay;

	DA[n]=Ay/(phy.Re*D[N]);
	DA[s]=Ay/(phy.Re*D[S]);

	if(isWall[e]) {
		FA[e]=0;

		DA[e]=Ax/(phy.Re*d[e]);
	} else {
		FA[e]=0.5*(FP[ne]->p[uvli]+FP[se]->p[uvli])*Ax;

		DA[e]=Ax/(phy.Re*D[E]);
	}
	if(isWall[w]) {
		FA[w]=0;

		DA[w]=Ax/(phy.Re*d[w]);
	} else {
		FA[w]=0.5*(FP[nw]->p[uvli]+FP[sw]->p[uvli])*Ax;

		DA[w]=Ax/(phy.Re*D[W]);
	}

	numeric_t rou(cd_c::phy.Dry.density/cd_c::phy.RefDensity), miu(cd_c::phy.Dry.dynamicViscosity/cd_c::phy.RefViscosity);
	for(size_t d=e; d<=s; ++d) {
		FA[d]*=rou;
		DA[d]*=miu;
	}
}

/*====================================================================================================
 *
 *    Hybrid & QUICK & Central Scheme Definitions, Common for flow field and evaporation calculation
 *
  ====================================================================================================*/
void cell_c::ConvectionItemUseHybridScheme (const size_t & uvli, const size_t & pli) {
	for(size_t i=E; i<=S; ++i) {
		CI[i]=NI[i];
	}

	if(FA[e]>=2*DA[e]) {
		a[E]+=0        ; aP+=FA[e]    ;
	} else if(FA[e]<=-2*DA[e]) {
		a[E]+=FA[e]    ; aP+=0        ;
	} else {
		a[E]+=0.5*FA[e]; aP+=0.5*FA[e];
	}

	if(FA[n]>=2*DA[n]) {
		a[N]+=0        ; aP+=FA[n]    ;
	} else if(FA[n]<=-2*DA[n]) {
		a[N]+=FA[n]    ; aP+=0        ;
	} else {
		a[N]+=0.5*FA[n]; aP+=0.5*FA[n];
	}

	if(FA[w]>=2*DA[w]) {
		a[W]-=FA[w]    ; aP+=0        ;
	} else if(FA[w]<=-2*DA[w]) {
		a[W]+=0        ; aP-=FA[w]    ;
	} else {
		a[W]-=0.5*FA[w]; aP-=0.5*FA[w];
	}

	if(FA[s]>=2*DA[s]) {
		a[S]-=FA[s]    ; aP+=0        ;
	} else if(FA[s]<=-2*DA[s]) {
		a[S]+=0        ; aP-=FA[s]    ;
	} else {
		a[S]-=0.5*FA[s]; aP-=0.5*FA[s];
	}
}

void cell_c::ConvectionItemUseQUICKScheme  (const size_t & uvli, const size_t & pli) {
	for(size_t i=E; i<=SS; ++i) {
		CI[i]=NI[i];
	}

	if(CI[EE]==voidIndex || CI[W]==voidIndex) {
		if(FA[e]>=2*DA[e]) {
			a[E]+=0        ; aP+=FA[e]    ;
		} else if(FA[e]<=-2*DA[e]) {
			a[E]+=FA[e]    ; aP+=0        ;
		} else {
			a[E]+=0.5*FA[e]; aP+=0.5*FA[e];
		}
	} else {
		if(FA[e]>=0) {
			a[W ]-=0.125*FA[e]; aP   +=0.75*FA[e]; a[E ]+=0.375*FA[e];
		} else {
			aP   +=0.375*FA[e]; a[E ]+=0.75*FA[e]; a[EE]-=0.125*FA[e];
		}
	}

	if(CI[NN]==voidIndex || CI[S]==voidIndex) {
		if(FA[n]>=2*DA[n]) {
			a[N]+=0        ; aP+=FA[n]    ;
		} else if(FA[n]<=-2*DA[n]) {
			a[N]+=FA[n]    ; aP+=0        ;
		} else {
			a[N]+=0.5*FA[n]; aP+=0.5*FA[n];
		}
	} else {
		if(FA[n]>=0) {
			a[S ]-=0.125*FA[n]; aP   +=0.75*FA[n]; a[N ]+=0.375*FA[n];
		} else {
			aP   +=0.375*FA[n]; a[N ]+=0.75*FA[n]; a[NN]-=0.125*FA[n];
		}
	}

	if(CI[WW]==voidIndex || CI[E]==voidIndex) {
		if(FA[w]>=2*DA[w]) {
			a[W]-=FA[w]    ; aP+=0        ;
		} else if(FA[w]<=-2*DA[w]) {
			a[W]+=0        ; aP-=FA[w]    ;
		} else {
			a[W]-=0.5*FA[w]; aP-=0.5*FA[w];
		}
	} else {
		if(FA[w]>=0) {
			a[WW]+=0.125*FA[w]; a[W ]-=0.75*FA[w]; aP   -=0.375*FA[w];
		} else {
			a[W ]-=0.375*FA[w]; aP   -=0.75*FA[w]; a[E ]+=0.125*FA[w];
		}
	}

	if(CI[SS]==voidIndex || CI[N]==voidIndex) {
		if(FA[s]>=2*DA[s]) {
			a[S]-=FA[s]    ; aP+=0        ;
		} else if(FA[s]<=-2*DA[s]) {
			a[S]+=0        ; aP-=FA[s]    ;
		} else {
			a[S]-=0.5*FA[s]; aP-=0.5*FA[s];
		}
	} else {
		if(FA[s]>=0) {
			a[SS]+=0.125*FA[s]; a[S ]-=0.75*FA[s]; aP   -=0.375*FA[s];
		} else {
			a[S ]-=0.375*FA[s]; aP   -=0.75*FA[s]; a[N ]+=0.125*FA[s];
		}
	}

//
//	if(CI[EE]==voidIndex) {//almost generate same results as above
//		if(FA[e]>=2*DA[e]) {
//			a[E]+=0        ; aP+=FA[e]    ;
//		} else if(FA[e]<=-2*DA[e]) {
//			a[E]+=FA[e]    ; aP+=0        ;
//		} else {
//			a[E]+=0.5*FA[e]; aP+=0.5*FA[e];
//		}
//	} else {
//		if(FA[e]>=0) {
//			a[W ]-=0.125*FA[e]; aP   +=0.75*FA[e]; a[E ]+=0.375*FA[e];
//		} else {
//			aP   +=0.375*FA[e]; a[E ]+=0.75*FA[e]; a[EE]-=0.125*FA[e];
//		}
//	}
//
//	if(CI[NN]==voidIndex) {
//		if(FA[n]>=2*DA[n]) {
//			a[N]+=0        ; aP+=FA[n]    ;
//		} else if(FA[n]<=-2*DA[n]) {
//			a[N]+=FA[n]    ; aP+=0        ;
//		} else {
//			a[N]+=0.5*FA[n]; aP+=0.5*FA[n];
//		}
//	} else {
//		if(FA[n]>=0) {
//			a[S ]-=0.125*FA[n]; aP   +=0.75*FA[n]; a[N ]+=0.375*FA[n];
//		} else {
//			aP   +=0.375*FA[n]; a[N ]+=0.75*FA[n]; a[NN]-=0.125*FA[n];
//		}
//	}
//
//	if(CI[WW]==voidIndex) {
//		if(FA[w]>=2*DA[w]) {
//			a[W]-=FA[w]    ; aP+=0        ;
//		} else if(FA[w]<=-2*DA[w]) {
//			a[W]+=0        ; aP-=FA[w]    ;
//		} else {
//			a[W]-=0.5*FA[w]; aP-=0.5*FA[w];
//		}
//	} else {
//		if(FA[w]>=0) {
//			a[WW]+=0.125*FA[w]; a[W ]-=0.75*FA[w]; aP   -=0.375*FA[w];
//		} else {
//			a[W ]-=0.375*FA[w]; aP   -=0.75*FA[w]; a[E ]+=0.125*FA[w];
//		}
//	}
//
//	if(CI[SS]==voidIndex) {
//		if(FA[s]>=2*DA[s]) {
//			a[S]-=FA[s]    ; aP+=0        ;
//		} else if(FA[s]<=-2*DA[s]) {
//			a[S]+=0        ; aP-=FA[s]    ;
//		} else {
//			a[S]-=0.5*FA[s]; aP-=0.5*FA[s];
//		}
//	} else {
//		if(FA[s]>=0) {
//			a[SS]+=0.125*FA[s]; a[S ]-=0.75*FA[s]; aP   -=0.375*FA[s];
//		} else {
//			a[S ]-=0.375*FA[s]; aP   -=0.75*FA[s]; a[N ]+=0.125*FA[s];
//		}
//	}


//	if(FA[e]>=0) {
//		a[W ]-=0.125*FA[e]; aP   +=0.75*FA[e]; a[E ]+=0.375*FA[e];
//	} else {
//		aP   +=0.375*FA[e]; a[E ]+=0.75*FA[e]; a[EE]-=0.125*FA[e];
//	}
//
//	if(FA[w]>=0) {
//		a[WW]+=0.125*FA[w]; a[W ]-=0.75*FA[w]; aP   -=0.375*FA[w];
//	} else {
//		a[W ]-=0.375*FA[w]; aP   -=0.75*FA[w]; a[E ]+=0.125*FA[w];
//	}
//
//	if(FA[n]>=0) {
//		a[S ]-=0.125*FA[n]; aP   +=0.75*FA[n]; a[N ]+=0.375*FA[n];
//	} else {
//		aP   +=0.375*FA[n]; a[N ]+=0.75*FA[n]; a[NN]-=0.125*FA[n];
//	}
//
//	if(FA[s]>=0) {
//		a[SS]+=0.125*FA[s]; a[S ]-=0.75*FA[s]; aP   -=0.375*FA[s];
//	} else {
//		a[S ]-=0.375*FA[s]; aP   -=0.75*FA[s]; a[N ]+=0.125*FA[s];
//	}
}

/*====================================================================================================
 *
 *             Evaporation related cCells member functions
 *
  ====================================================================================================*/
void cCell_c::resetCoefficient() {
	aP=0; b=0; r=0; DeltaVT=0;
	a.clear(); CI.clear(); ap.resize(PI.size(), 0);
	for(size_t i=0; i<PI.size(); ++i) {
		ap[i]=0;
	}
}

//----------------------------------------------------------------------------------------------------
void cCell_c::calculateCoefficient(const size_t & uvli, const size_t & cli, const std::bitset<EqComponent> & ecbs) {
	resetCoefficient();

	if(EB.T==NoBC) {
//		if(cd_c::EDT!=0) {
		if(true) {
			calculateFD(uvli);

			if(ecbs[convectionItem]) {
				switch (cd_c::phy.SchemeEvaporation) {
				case Hybrid :
					a.resize(sideNumber, 0); CI.resize(sideNumber, voidIndex);
					ConvectionItemUseHybridScheme(uvli, cli);
					break;
				case HayaseQUICK :
					a.resize(NI.size(), 0); CI.resize(NI.size(), voidIndex);
					ConvectionItemUseQUICKScheme(uvli, cli);
//					if(nearBC) 	{//this won't work
//						a.resize(sideNumber, 0); CI.resize(sideNumber, voidIndex);
//						ConvectionItemUseHybridScheme(uvli, cli);
//					} else {
//						a.resize(NI.size(), 0); CI.resize(NI.size(), voidIndex);
//						ConvectionItemUseQUICKScheme(uvli, cli);
//					}
					break;
				case Upwind: case PowerLaw: case TVD: case QUICK:
					break;
				}
			}

			if(ecbs[diffussionItem]) {
				if(!ecbs[convectionItem]) {
					a.resize(sideNumber, 0); CI.resize(sideNumber, voidIndex);
					for(size_t i=E; i<=S; ++i) {
						CI[i]=NI[i];
					}
				}
				DiffussionItemUseCentralScheme(uvli, cli);
			}

			for(size_t i=0; i<CI.size(); ++i) {
				if(CI[i]==voidIndex) {
					a.erase(a.begin()+i); CI.erase(CI.begin()+i);
					--i;
				}
			}
		}

		if(ecbs[deltaPAVT]) calculateDeltaVT(cli);
	} else if(EB.T==ConstConcentration) {
		aP=EB.aP;
		b =EB.b ;
		ap.clear();
		PI.clear();
	} else {
		aP=EB.aP;
		b =EB.b ;
		CI=EB.I ;
		a =EB.a ;
	}
}

//----------------------------------------------------------------------------------------------------
void cCell_c::DiffussionItemUseCentralScheme(const size_t & uvli, const size_t & cli) {
	a[E]-=isWall[e]?0:DA[e];
	a[W]-=isWall[w]?0:DA[w];
	a[N]-=isWall[n]?0:DA[n];
	a[S]-=isWall[s]?0:DA[s];

	aP  +=DA[e]+DA[w]+DA[n]+DA[s];



//
//	if(cd_c::phy.SchemeFlowfield==Hybrid) {// For Triditional Hybrid Method, the DA sometimes can be dropped
//		if(FA[e]>=2*DA[e] || FA[e]<=-2*DA[e]) {
//			a[E]-=0;
//		} else {
//			a[E]-=isWall[e]?0:DA[e]; aP+=DA[e];
//		}
//		if(FA[w]>=2*DA[w] || FA[w]<=-2*DA[w]) {
//			a[W]-=0;
//		} else {
//			a[W]-=isWall[w]?0:DA[w]; aP+=DA[w];
//		}
//		if(FA[n]>=2*DA[n] || FA[n]<=-2*DA[n]) {
//			a[N]-=0;
//		} else {
//			a[N]-=isWall[n]?0:DA[n]; aP+=DA[n];
//		}
//		if(FA[s]>=2*DA[s] || FA[s]<=-2*DA[s]) {
//			a[S]-=0;
//		} else {
//			a[S]-=isWall[s]?0:DA[s]; aP+=DA[s];
//		}
//	}



	for(size_t i=0; i<PP.size(); ++i) {
		ap[i]-=PP[i]->CA[0]/(cd_c::phy.Pe*PP[i]->CD[0]);
		aP   -=ap[i];

		if(cd_c::phy.filmEffect) {//*******************Film Effect*****
			b    -=PP[i]->C2Phi(cli)*ap[i];
			ap[i]*=PP[i]->C1Phi(cli);

		}
	}
}

//----------------------------------------------------------------------------------------------------
void cCell_c::calculateDeltaVT(const size_t & cli) {
	DeltaVT=cd_c::EDT!=0?V/cd_c::EDT:V/cd_c::phy.RefTime;

	b +=DeltaVT*p[cli];
	aP+=DeltaVT;
}

//----------------------------------------------------------------------------------------------------
void cCell_c::calculateFD(const size_t & uvli) {
	FA[e]=Ax*FP[e]->p[uvli];
	FA[w]=Ax*FP[w]->p[uvli];
	FA[n]=Ay*FP[n]->p[uvli];
	FA[s]=Ay*FP[s]->p[uvli];

	if(isWall[e] || NI[E]==voidIndex) {
//		DA[e]=Ax/(phy.Pe*d[e]);
		DA[e]=0;
	} else {
		DA[e]=Ax/(phy.Pe*D[E]);
	}

	if(isWall[w] || NI[W]==voidIndex) {
//		DA[w]=Ax/(phy.Pe*d[w]);
		DA[w]=0;
	} else {
		DA[w]=Ax/(phy.Pe*D[W]);
	}

	if(isWall[n] || NI[N]==voidIndex) {
//		DA[n]=Ay/(phy.Pe*d[n]);
		DA[n]=0;
	} else {
		DA[n]=Ay/(phy.Pe*D[N]);
	}

	if(isWall[s] || NI[S]==voidIndex) {
//		DA[s]=Ay/(phy.Pe*d[s]);
		DA[s]=0;
	} else {
		DA[s]=Ay/(phy.Pe*D[S]);
	}
}
