/*
 * PoreThroat.cpp
 *
 *  Created on: Apr 15, 2014
 *      Author: david
 */

#include "PoreThroat.h"
#include "FiniteVolume.h"

using namespace std;

/*----------------------------------------------------------------------------------------------------
 *
 *                                        Initialize Static Members
 *
  ----------------------------------------------------------------------------------------------------*/
std::vector<index_t> pore_c::openPI(0);

/*----------------------------------------------------------------------------------------------------
 *
 *                                        class member functions definition
 *
  ----------------------------------------------------------------------------------------------------*/
bool pore_c::setPosition(const position_e & pos) {
	if(pos==open) {
		position=open;
		sStatus=empty;
		vpStatus=unknown;
		for(size_t i=0; i<HowManyCopy; ++i) {
			p[i]=cd_c::phy.EnvironmentConcentration;
		}
		aP=1;
		b =cd_c::phy.EnvironmentConcentration;
	} else if(pos==meniscus) {
		position=meniscus;
		sStatus=empty;
		vpStatus=saturated;
		for(size_t i=0; i<HowManyCopy; ++i) {
			p[i]=cd_c::phy.SaturatedConcentration;
		}
		aP=1;
		b =cd_c::phy.SaturatedConcentration;
	} else if(pos==inside) {
		position=inside;
		sStatus=empty;
		vpStatus=unknown;
	} else if(pos==liquidSide) {
		position=liquidSide;
		sStatus=liquid;
		vpStatus=none;
		for(size_t i=0; i<HowManyCopy; ++i) {
			p[i]=cd_c::phy.SaturatedConcentration;
		}
		aP=1;
		b =cd_c::phy.SaturatedConcentration;
	}
	return true;
}

bool pore_c::setLabel(const size_t & CL) {
	clusterLabel=CL;
	labelled=true;

	return true;
}

bool pore_c::removeLabel() {
	clusterLabel=voidIndex;
	labelled=false;

	return true;
}

bool pore_c::isEmpty() {
	if(CI.size()!=0) return true;
	for(size_t i=0; i<TP.size(); i++) {
		if(TP[i]->isEmpty()) return true;
	}
	return false;
}

bool pore_c::isMeniscus() {
	if(isEmpty()) {
		for(size_t j=0; j<TP.size(); j++) {
			if(!TP[j]->isEmpty()) return true;
		}
	}
	return false;
}

bool pore_c::updateStatus() {
	if(isEmpty()) {
		if(position==open) setPosition(open);
		else if(isMeniscus()) setPosition(meniscus);
		else setPosition(inside);
	} else {
		setPosition(liquidSide);
	}

	return true;
}

numeric_t pore_c::MassOutFlux(const size_t & cli) {
	numeric_t MOF(0);
	if(cd_c::phy.filmEffect) {
		for(size_t i=0; i<PP.size(); ++i) {
			MOF+=(f[cli]-PP[i]->f[cli])/TP[i]->L;
		}
	} else {
		for(size_t i=0; i<PP.size(); ++i) {
			MOF+=(TP[i]->A/TP[i]->L)*(p[cli]-PP[i]->p[cli]);
		}
	}

	MOF+=MassToCell(cli);

	return MOF;
}

numeric_t pore_c::MassToCell(const size_t & cli) {
	numeric_t MTC(0);
	if(cd_c::phy.filmEffect) {
		for(size_t i=0; i<NP.size(); ++i) {
			MTC+=(CA[i]   /CD[i]   )*((C1Phi(cli)*f[cli]+C2Phi(cli))-NP[i]->p[cli]);
		}
	} else {
		for(size_t i=0; i<NP.size(); ++i) {
			MTC+=(CA[i]   /CD[i]   )*(p[cli]-NP[i]->p[cli]);
		}
	}

	return MTC;
}

void pore_c::resetCoefficient() {
	aP=0; b=0; a.resize(CI.size(),0); ap.resize(PI.size(), 0);
	for(size_t i=0; i<CI.size(); ++i) {
		a [i]=0;
	}
	for(size_t i=0; i<PI.size(); ++i) {
		ap[i]=0;
	}
}

void pore_c::calculateCoefficient(const size_t & cli, const std::bitset<EqComponent> & ec) {
	resetCoefficient();

	if(position==inside) {
		if(ec==6 || ec==7) {
			for(size_t i=0; i<CI.size(); ++i) {
				a [i]-=CA[i]   /CD[i]   ;
				aP   -=a[i];
			}

			for(size_t i=0; i<TP.size(); ++i) {
				ap[i]-=TP[i]->A/TP[i]->L;
				aP   -=ap[i];
			}
		} else if(ec==5) {
			aP=1;
			b =p[cli];
		} else {
			//no such case
		}
	} else if(position==open) {
		aP=1;
		b =cd_c::phy.EnvironmentConcentration;
	} else if(position==meniscus) {
		aP=1;
		b =cd_c::phy.SaturatedConcentration;
	} else if(position==liquidSide) {
		aP=1;
		b =cd_c::phy.SaturatedConcentration;
	} else {
		//no such case
	}

}

void pore_c::calculatePhiCoefficient(const size_t & cli, const std::bitset<EqComponent> & ec) {
	resetCoefficient();

	if(position==inside) {
		if(ec==6 || ec==7) {
			for(size_t i=0; i<CI.size(); ++i) {
				a [i]-=CA[i]   /CD[i]   ;
				aP   +=CA[i]   /CD[i]   *C1Phi(cli);
				b    -=CA[i]   /CD[i]   *C2Phi(cli);
			}

			for(size_t i=0; i<TP.size(); ++i) {
				ap[i]-=1/TP[i]->L;
				aP   -=ap[i];
			}
		} else if(ec==5) {
			aP=1;
			b =f[cli];
		} else {
			//no such case
		}
	} else if(position==open) {
		aP=1;
		b =PhiG*cd_c::phy.EnvironmentConcentration;
	} else if(position==meniscus) {
		aP=1;
		b =PhiLFI;
	} else if(position==liquidSide) {
		aP=1;
		b =PhiLFI;
	} else {
		//no such case
	}
}

bool pore_c::updatef2pr(const std::bitset<HowManyCopy> & hmc) {
	for(size_t i=0; i<HowManyCopy; ++i) {
		if(hmc[i]) {
			switch (cd_c::phy.FilmApproximationType) {
			case 1:
				p[i]=C1Phi(i)*f[i]+C2Phi(i);
				break;
			}
		}
	}

	return true;
}

numeric_t pore_c::C1Phi(const size_t & cli) {
	numeric_t C1Phi(2/f[cli]/(1+exp(-cd_c::phy.FilmApproximationCoefficient*f[cli]/PhiG)));

	return C1Phi;
}
numeric_t pore_c::C2Phi(const size_t & cli) {
	numeric_t C2Phi(-1);

	return C2Phi;
}
/*----------------------------------------------------------------------------------------------------
 *
 *                                        class member functions definition
 *
  ----------------------------------------------------------------------------------------------------*/
bool throat_c::setPosition(const position_e & pos) {
	if(pos==open) {
		position=open;
		sStatus=empty;
	} else if(pos==meniscus) {
		position=meniscus;
		sStatus=liquid;
	} else if(pos==inside) {
		position=inside;
		sStatus=empty;
	} else if(pos==liquidSide) {
		position=liquidSide;
		sStatus=liquid;
	}

	return true;
}

bool throat_c::setLabel(const index_t & CN) {
	clusterLabel=CN;
	labelled=true;

	return true;
}

bool throat_c::removeLabel() {
	clusterLabel=voidIndex;
	labelled=false;

	return true;
}

bool throat_c::isMeniscus() {
	if(!isEmpty()) {
		for(size_t j=0; j<PP.size(); ++j) {
			if(PP[j]->isEmpty()) return true;
		}
	}

	return false;
}

bool throat_c::updateStatus() {
	if(isEmpty()) {
		setPosition(inside);
	} else {
		if(isMeniscus()) {
			setPosition(meniscus);
		} else {
			setPosition(liquidSide);
		}
	}

	return true;
}

numeric_t throat_c::MassOutFlux(const size_t & cli) {
	numeric_t MOF(0);

	for(size_t i=0; i<PI.size(); ++i) {
		if(PP[i]->isMeniscus()) {
			numeric_t totA(0);
			for(size_t j=0; j<PP[i]->TI.size(); ++j) {
				if(!((PP[i]->TP[j])->isEmpty())) totA+=(PP[i]->TP[j])->A;
			}
			MOF+=PP[i]->MassOutFlux(cli)*A/totA;
		}
	}

	return MOF;
}

