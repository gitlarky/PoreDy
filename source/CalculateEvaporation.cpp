/*
 * CalculateEvaporation.cpp
 *
 *  Created on: Apr 15, 2014
 *      Author: david
 */

#include "GlobalDeclaration.h"

using namespace std;
using namespace Eigen;

vector<size_t> unlabelledT;
vector<size_t> unlabelledP;
vector<size_t>  TtoExam;
vector<size_t>  PtoExam;
vector<size_t>  examinedT;
vector<size_t>  examinedP;
vector<size_t> unknownP;

bool initializeNetwork();
bool initializeExternalField();
bool checkNetwork();
bool clusterLabeling();
bool calculateConcentration(const bool & cal, const bool & update);
bool calcVaporPressure();
bool emptyingOneThroat();
bool frequentOutput(const size_t & sc);
bool calculateFluxNTime(const size_t & cli);

bool inClusterToEmpty(const size_t & index);//a function to determine if the index of clusters in cd_c::ClusterToEmpty or not

bool plotSaturation(ofstream & ofs);
bool plotAllField(ofstream & ofs);
bool plotExtField(ofstream & ofs);

bool solveExternalConcentrationEq(const size_t & cli, const std::bitset<EqComponent> & ec, const size_t & clo, const size_t & rp);
bool solveConcentrationEq(const size_t & cli, const std::bitset<EqComponent> & ec, const size_t & clo);
bool solvePhiEq(const size_t & cli, const std::bitset<EqComponent> & ec, const size_t & clo);

bool inExaminedT(const size_t & index);//a function to determine if the index of throats ever searched
bool inExaminedP(const size_t & index);//a function to determine if the index of pores ever searched
bool inTtoExam(const size_t & index);//a function to determine if the index of throats ever searched
bool inPtoExam(const size_t & index);//a function to determine if the index of pores ever searched
bool inUnlabelledT(const size_t & index);//a function to determine if the index of throats in unlabelledT or not
bool inUnlabelledP(const size_t & index);//a function to determine if the index of pores in unlabelledP or not
vector<size_t>::iterator ItInUnlabelledT(const size_t & index);//a function to return the position of a index of throats if in, if not return NULL
vector<size_t>::iterator ItInUnlabelledP(const size_t & index);//a function to return the position of a index of pores if in, if not return NULL

size_t unknownPi(size_t & index);//find out the gas pore with unknown vapor pressure

bool liquidExistYet();
bool wetYet();
inline numeric_t TotSaturation() {return cd_c::TotLiquid/cd_c::TotLiquid0;}

/*====================================================================================================

                                        calculateEvaporation

====================================================================================================*/
bool calculateEvaporation() {
	size_t stepCount(0);
	numeric_t EvaporationTime(0), StepTime(0);
	numeric_t TotMassFluxOut(0), TotMassOut(0), TotLiquidToGasFlux(0), TotLiquidToGas(0);

	string   seh(caseName+".eh");
	ofstream ofeh(seh.c_str());

	initializeNetwork();
	initializeExternalField();

    do {
    	clusterLabeling();
//    	checkNetwork();

    	if(stepCount%cd_c::phy.DryOutputFrequency==0) frequentOutput(stepCount);

    	cd_c::ESS=cd_c::phy.EvaporationSubStep;
    	StepTime=0; TotLiquidToGas=0; TotMassOut=0;

    	calculateFluxNTime(0); //cli=0

    	do {
			if(cd_c::ESS>1) {
				cd_c::EDT=cd_c::TimeToEmptyMax/cd_c::ESS;

				calculateConcentration(true, true); //calculate and update, copy 2->0
				calculateFluxNTime(0); // cli=0

				cd_c::TimeToEmptyMax-=cd_c::EDT;
			} else if(cd_c::ESS==1) {
    			size_t SubStep(0);
	    		do {
		    		cd_c::EDT=cd_c::TimeToEmptyMax;

	    			if(SubStep==1) cout<<"\tSubStep # "<<"0"<<";\tTimeToEmptyMax="<<cd_c::TimeToEmptyMax<<endl;

					calculateConcentration(true, false); //calculate but not update
					calculateFluxNTime(2);

	    			if(SubStep>=1) cout<<"\tSubStep # "<<SubStep<<";\tTimeToEmptyMax="<<cd_c::TimeToEmptyMax<<endl;
		    		SubStep++;
	    		} while(cd_c::EDT>cd_c::TimeToEmptyMax || (cd_c::TimeToEmptyMax-cd_c::EDT)/cd_c::TimeToEmptyMax>cd_c::phy.RelativeConvergeCriteria);//
//    		    } while(abs(cd_c::TimeToEmptyMax-cd_c::EDT)/cd_c::TimeToEmptyMax>cd_c::phy.RelativeConvergeCriteria);//
	    		calculateConcentration(false, true); //update only, just copy 2->0
			} else {
				//do nothing
			}

        	bool GotEmptyItNow=cd_c::ESS==1?true:false;//If it is the last step, then got empty the maxThroat and adjust the cd_c::EDT
        	for(size_t i=1; i<Ct.size(); ++i) {
        		if(inClusterToEmpty(i)) {
            		Ct[i].drying(cd_c::EDT, GotEmptyItNow);
        		} else {
            		Ct[i].drying(cd_c::EDT, false);
        		}

        		TotLiquidToGas+=Ct[i].MassOut;
        	}
			for(size_t i=0; i<pore_c::openPI.size(); ++i) {
				TotMassOut+=P[pore_c::openPI[i]].MassToCell(0)*cd_c::EDT;
			}

        	StepTime       +=cd_c::EDT;
        	EvaporationTime+=cd_c::EDT;
        	cd_c::ESS--;
    	}while(cd_c::ESS>0);

    	stepCount++;
    	TotMassFluxOut    =TotMassOut    /StepTime;
    	TotLiquidToGasFlux=TotLiquidToGas/StepTime;
    	cd_c::TotLiquid  -=TotMassOut;

//    	oflg<<"During Step # "<<stepCount<<" : StepTime="<<StepTime<<"; "<<"Total Evaporation Time Accumulated To: "<<EvaporationTime<<endl;
//    	oflg<<"Emptied";
//    	oflg<<"TotLiquidToGasFlux="<<TotLiquidToGasFlux<<"\tTotLiquidToGas="<<TotLiquidToGas<<endl;
//    	oflg<<"TotMassFluxOut    ="<<TotMassFluxOut    <<"\tTotMassOut    ="<<TotMassOut    <<endl;

    	if((stepCount-1)%10==0) {
    		cout<<"Step #"<<"\tStepTime"<<"\tEvaporationTime"<<"\tNetworkSaturation"<<"\tMassOutFlux"        <<"\tMassOut"
    				                                                                <<"\tLiquidToGasMassFlux"<<"\tLiquidToGasMass"<<endl;
    	}
    	cout<<stepCount<<"\t"<<StepTime<<"\t"<<EvaporationTime<<"\t"<<TotSaturation()<<"\t"<<TotMassFluxOut    <<"\t"<<TotMassOut
    			                                                                     <<"\t"<<TotLiquidToGasFlux<<"\t"<<TotLiquidToGas<<endl;
    	if((stepCount-1)   ==0) {
    		ofeh<<"Step #"<<"\tStepTime [hr]"<<"\tEvaporationTime [hr]"<<"\tNetworkSaturation"
    				      <<"\tMassOutFlux [kg/hr]"        <<"\tMassOut [kg]"
    				      <<"\tLiquidToGasMassFlux [kg/hr]"<<"\tLiquidToGasMass [kg]"<<endl;
    	}
    	ofeh<<stepCount<<"\t"<<StepTime*cd_c::phy.RefTime/3600<<"\t"<<EvaporationTime*cd_c::phy.RefTime/3600<<"\t"<<TotSaturation()
    			       <<"\t"<<3600*cd_c::phy.RefMassFlux*TotMassFluxOut    <<"\t"<<TotMassOut*cd_c::phy.RefMass
    			       <<"\t"<<3600*cd_c::phy.RefMassFlux*TotLiquidToGasFlux<<"\t"<<TotLiquidToGas*cd_c::phy.RefMass<<endl;

    	for(size_t i=0; i<P.size(); ++i) {
    		P[i].removeLabel();
    		P[i].updateStatus();
    	}
    	for(size_t i=0; i<T.size(); ++i) {
    		T[i].removeLabel();
    		T[i].updateStatus();
    	}
    } while(liquidExistYet() && stepCount<cd_c::phy.EvaporationMaxStep);
    frequentOutput(stepCount);

    ofeh.close();
    cout<<endl<<"Drying Process Finished!"<<endl;

	return true;
}

//----------------------------------------------------------------------------------------------------
bool inClusterToEmpty(const size_t & index) {//a function to determine if the index of clusters in cd_c::ClusterToEmpty or not
	for(size_t i=0; i<cd_c::ClusterToEmpty.size(); ++i) {
		if(index==cd_c::ClusterToEmpty[i]) return true;
	}

	return false;
}

/*----------------------------------------------------------------------------------------------------

                                        Function initializeNetwork

----------------------------------------------------------------------------------------------------*/
bool initializeNetwork() {
	cout<<"Initialize Pore-Throat Network......"<<endl;
	cd_c::TotLiquid=0;
	for(size_t i=0; i<T.size(); ++i) {
		T[i].setPosition(liquidSide);
		cd_c::TotLiquid+=T[i].M;
	}
	cd_c::TotLiquid0=cd_c::TotLiquid;

	for(size_t i=0; i<P.size(); ++i) {
		P[i].a .resize(P[i].CI.size(), 0);
		P[i].ap.resize(P[i].PI.size(), 0);
		if(P[i].CI.size()!=0) {
			P[i].setPosition(meniscus);
			pore_c::openPI.push_back(i);
			T[P[i].TI[0]].setPosition(meniscus);
		} else {
			P[i].setPosition(liquidSide);
		}
	}

	for(size_t i=0; i<C.size(); ++i) {
		for(size_t cli=0; cli<HowManyCopy; ++cli) {
			if(C[i].EB.T==ConstConcentration) {
				C[i].p[cli]=C[i].EB.b/C[i].EB.aP;
			} else {
				C[i].p[cli]=cd_c::phy.EnvironmentConcentration;
//				C[i].p[cli]=cd_c::phy.SaturatedConcentration;
//				if(C[i].PI.size()==0) {
//
//				} else {
//					C[i].p[cli]=cd_c::phy.SaturatedConcentration;
//				}
			}
		}
	}

	if(cd_c::phy.filmEffect) {
		for(size_t k=0; k<porousBlocks; ++k) {
			numeric_t avgRI(0), avgRM(0), avgA(0);
			numeric_t avgPhiF(0), avgPhiG(0), avgPhiLFI(0), avgPhiFGI(0);
			for(size_t l=0; l<Bk[k].avgTDiameter.size(); ++l) {
				numeric_t RI=Bk[k].avgTDiameter[l]/2;

				if(Bk[k].PolyN[l]>=3) {
					numeric_t Alpha=pi/2-pi/Bk[k].PolyN[l];
					numeric_t ThetaC=pi/Bk[k].PolyN[l];
					numeric_t Kai=cos(cd_c::phy.Wet.contactAngle)
					             +sqrt((ThetaC-cd_c::phy.Wet.contactAngle+sin(cd_c::phy.Wet.contactAngle)*cos(cd_c::phy.Wet.contactAngle))/tan(ThetaC));
					numeric_t RM=RI/Kai;
					numeric_t Psi1=cos(Alpha+cd_c::phy.Wet.contactAngle)
						          *(cos(Alpha+cd_c::phy.Wet.contactAngle)+sin(Alpha+cd_c::phy.Wet.contactAngle)*tan(Alpha));
					numeric_t Psi2=1-cd_c::phy.Wet.contactAngle/(pi/2-Alpha);
					numeric_t Psi3=cos(Alpha+cd_c::phy.Wet.contactAngle)/cos(Alpha);
					numeric_t Psi4=cos(cd_c::phy.Wet.contactAngle)*cos(Alpha+cd_c::phy.Wet.contactAngle)/sin(Alpha)-pi/2+Alpha+cd_c::phy.Wet.contactAngle;
					numeric_t Psi5=(pi/2-Alpha)*tan(Alpha);
					numeric_t Beta=12*(Psi1-Psi5*Psi2)*pow(sin(Alpha)*(1-Psi5)*(Psi3)/(1-sin(Alpha))/Psi5, 2)/pow(Psi1-Psi5*Psi2, 3);
					numeric_t Kappa=Psi4/Beta;
					numeric_t A=pow(RI, 2)*(pi+N*(1/tan(Alpha)-pi/2+Alpha));
					cout<<"Psi4="<<Psi4<<"\tBeta="<<Beta<<"\tKappa="<<Kappa<<endl;
					numeric_t PhiF=Bk[k].PolyN[l]*Kappa*cd_c::phy.LiquidConcentration/3/cd_c::phy.Ca;
					numeric_t PhiG=A;

					avgRI  +=RI;
					avgRM  +=RM;
		            avgA   +=A;
					avgPhiF+=PhiF;
					avgPhiG+=PhiG;
				} else {
					numeric_t A=pi*RI*RI;

					avgRI  +=RI;
		            avgA   +=A;
				}
			}
			avgRI  /=Bk[k].avgTDiameter.size();
			avgRM  /=Bk[k].avgTDiameter.size();
			avgA   /=Bk[k].avgTDiameter.size();
			avgPhiF/=Bk[k].avgTDiameter.size();
			avgPhiG/=Bk[k].avgTDiameter.size();

			avgPhiLFI=avgPhiF*pow(avgRM, 3)+avgPhiG*cd_c::phy.SaturatedConcentration;
			avgPhiFGI=                      avgPhiG*cd_c::phy.SaturatedConcentration;

			cout<<"Block#"<<k<<":\tRM="<<avgRM<<"\tA="<<avgA<<"\tCa="<<cd_c::phy.Ca
					          <<"\tPhiF="<<avgPhiF<<"\tPhiG="<<avgPhiG<<"\tPhiLFI="<<avgPhiLFI<<"\tPhiFGI="<<avgPhiFGI<<endl;

			for(size_t i=0; i<=Bk[k].Nx; ++i) {
				for(size_t j=0; j<=Bk[k].Ny; ++j) {
					if(Bk[k].PI[i][j]==voidIndex) {
						//do nothing
					} else {
						P[Bk[k].PI[i][j]].RI    =avgRI    ;
						P[Bk[k].PI[i][j]].RM    =avgRM    ;
						P[Bk[k].PI[i][j]].A     =avgA     ;
						P[Bk[k].PI[i][j]].PhiF  =avgPhiF  ;
						P[Bk[k].PI[i][j]].PhiG  =avgPhiG  ;

						P[Bk[k].PI[i][j]].PhiLFI=avgPhiLFI;
						P[Bk[k].PI[i][j]].PhiFGI=avgPhiFGI;

						for(size_t cli=0; cli<HowManyCopy; ++cli) {
							P[Bk[k].PI[i][j]].f[cli]=P[Bk[k].PI[i][j]].PhiLFI;
						}
						P[Bk[k].PI[i][j]].updatef2pr(7);//all copy updated
					}//jump out of if void
				}
			}//jump out of loop
		}
	}


	cout<<"Pore Throats Network Initialized!"<<endl;

	return true;
}
/*----------------------------------------------------------------------------------------------------

                                        Function initializeExternalField

----------------------------------------------------------------------------------------------------*/
bool initializeExternalField() {
	cout<<"Start initialize the external field..."<<endl;

	size_t iterCount(0);
	numeric_t avgCError(0), avgCErr0r(0);
	size_t cLevelIn(0), cLevelOut(1);//C values, which to use and where to put the results
	bool FilmEffectTurnedOff(false);
	size_t RelatedPores(0);

	cout<<"Prepare boundary conditions......"<<endl;
	if(cd_c::phy.filmEffect) {
		cd_c::phy.filmEffect=false;
		FilmEffectTurnedOff=true;
	}

	for(size_t i=0; i<P.size(); ++i) {
		if(P[i].CI.size()) {
			RelatedPores++;
		}
	}
	cout<<"RelatedPores"<<RelatedPores<<endl;

	cout<<"Start iterating......"<<endl;
	do {
		solveExternalConcentrationEq(cLevelIn, 3,cLevelOut, RelatedPores);//3=1+2: Calculate Convection+Diffusion

		avgCError=0;
		for(size_t i=0; i<C.size(); ++i) {
			avgCError+=abs(C[i].p[cLevelOut]-C[i].p[cLevelIn]);
		}
		avgCError/=C.size();

		if(iterCount==0) {
			avgCErr0r=avgCError;
			cout<<"Initial Average Concentration Error="<<avgCErr0r<<endl;
		}
		avgCError/=avgCErr0r;

		if(iterCount%10==0) cout<<"Iteration #"<<"\tavgCError"<<endl;
		iterCount++;
		cout<<iterCount<<"\t"<<avgCError<<endl;

		for(size_t i=0; i<C.size(); ++i) {
			C[i].p[cLevelIn]=C[i].p[cLevelOut];
		}

	} while(avgCError>cd_c::phy.FlowCriteria && iterCount<cd_c::phy.FlowfieldMaxStep);
	if(avgCError<=cd_c::phy.FlowCriteria) {
		cout<<"External Concentration Field Calculation Converged!"<<endl;
	} else {
		cout<<"External Concentration Calculation Reached Max Steps!"<<endl;
	}

	cout<<"Restore boundary conditions......"<<endl;
	if(FilmEffectTurnedOff) cd_c::phy.filmEffect=true;

	return true;
}

//---------------------------------------------------------------------------------------------------------------
bool solveExternalConcentrationEq(const size_t & cli, const std::bitset<EqComponent> & ec, const size_t & clo, const size_t & rp) {
	if(cell_c::phy.ImplicitEvaporation) {
		size_t newcp(C.size()+rp);
		VectorXd X(newcp), B(newcp);

		int EachSize(0);
		if(cd_c::phy.SchemeEvaporation==Hybrid     ) EachSize=5;
		if(cd_c::phy.SchemeEvaporation==HayaseQUICK) EachSize=9;
		typedef Eigen::Triplet<double> Tri;
		vector<Tri> TripletList;
		TripletList.reserve(newcp*EachSize);

		size_t rpcount(0);

		for(size_t i=0; i<C.size(); ++i) {
			C[i].calculateCoefficient(0, cli, ec);//uvli=0
			B(i)=C[i].b;
			TripletList.push_back(Tri(i, i, C[i].aP));
			for(size_t j=0; j<C[i].CI.size(); ++j) {
				TripletList.push_back(Tri(i, C[i].CI[j], C[i].a[j]));
			}
			for(size_t j=0; j<C[i].PI.size(); ++j) {
				TripletList.push_back(Tri(i, rpcount   +C.size(), C[i].ap[j]));
				rpcount++;
			}
		}

		for(size_t i=0, I=C.size(); i<rp      ; ++i, ++I) {

			B(I)=cd_c::phy.SaturatedConcentration;
			TripletList.push_back(Tri(I, I, 1.));
		}
		SparseMatrix<double, Eigen::ColMajor> A(newcp, newcp);
		A.setFromTriplets(TripletList.begin(), TripletList.end());
		SparseLU<SparseMatrix<double>> solver;
		solver.analyzePattern(A);
		solver.factorize(A);
		X= solver.solve(B);

		for(size_t i=0; i<C.size(); ++i) {
			C[i].p[clo]=X[i];
		}
	} else {

	}

	return true;
}
/*----------------------------------------------------------------------------------------------------

                                        Function checkNetwork

----------------------------------------------------------------------------------------------------*/
bool checkNetwork() {
	oflg<<"Pores information:"<<endl;
	for(size_t i=0; i<P.size(); ++i) {
		oflg<<"P#"<<i<<":"<<"\tCt#"<<P[i].clusterLabel<<":"<<"\tsStatus="<<P[i].sStatus<<"\tposition="<<P[i].position;
		if(P[i].isMeniscus()) oflg<<"\tMassOutFlux="<<P[i].MassOutFlux(0);
		for(size_t j=0; j<P[i].PI.size(); ++j) {
			oflg<<"\tP"<<P[i].PI[j];
		}
		for(size_t j=0; j<P[i].CI.size(); ++j) {
			oflg<<"\tC"<<P[i].CI[j];
		}
		oflg<<"\t|";
		for(size_t j=0; j<P[i].TI.size(); ++j) {
			oflg<<"\tT"<<P[i].TI[j];
		}
		oflg<<endl;
	}

	oflg<<"Throats information:"<<endl;
	for(size_t i=0; i<T.size(); ++i) {
		oflg<<"T#"<<i<<": [P"<<T[i].PI[0]<<"\tP"<<T[i].PI[1]<<"]\t"<<"\tCt#"<<T[i].clusterLabel<<":"
			<<"\tRI="<<T[i].RI<<"\tL="<<T[i].L<<"\tM="<<T[i].M
			<<"\tsStatus="<<T[i].sStatus<<"\tposition="<<T[i].position;
		if(T[i].isMeniscus()) oflg<<"\tMassOutFlux="<<T[i].MassOutFlux(0);
		oflg<<endl;
	}

	oflg<<"Clusters information:"<<endl;
	for(size_t i=0; i<Ct.size(); ++i) {
		oflg<<"Ct#"<<i<<":"<<"\tMaxThroat="<<Ct[i].maxThroat.size()
			<<"\tMassOutFlux="<<Ct[i].MassOutFlux<<"\tTimeToEmptyMax="<<Ct[i].timeToEmptyMax<<"\tMassOut="<<Ct[i].MassOut;

		for(size_t j=0; j<Ct[i].PI.size(); ++j) {
			oflg<<"\tP"<<Ct[i].PI[j];
		}

		oflg<<"\t|";
		for(size_t j=0; j<Ct[i].TI.size(); ++j) {
			oflg<<"\tT"<<Ct[i].TI[j];
		}

		oflg<<endl;
	}

	oflg<<"Network information: "<<"Total Liquid Mass="<<cd_c::TotLiquid<<"; Total Network Saturation="<<TotSaturation()<<endl;





//	cout<<"\t"<<P[0].PhiF<<"\t"<<P[0].PhiG<<"\t"<<P[0].PhiLFI<<"\t"<<P[0].PhiFGI<<"\t"<<P[0].RM<<endl;
//	cout<<"\t"<<T[0].PhiF<<"\t"<<T[0].PhiG<<"\t"<<T[0].Kappa<<"\t"<<T[0].Kai<<"\t"<<T[0].A<<endl;

	return true;
}

/*----------------------------------------------------------------------------------------------------

                                        Function clusterLabelling

----------------------------------------------------------------------------------------------------*/
bool clusterLabeling() {
	Ct.clear();
	unlabelledT.clear(); unlabelledP.clear();
	TtoExam.clear(); PtoExam.clear();
	examinedT.clear(); examinedP.clear();

	index_t currentCluster(Ct.size());
//	oflg<<"Create a new empty cluster: Ct# "<<currentCluster<<" to hold all empty pores and throats:"<<endl;
	Ct.push_back(cluster_c());

//	oflg<<"Add";
	for(size_t i=0; i<T.size(); ++i) {
		if(T[i].sStatus==empty) {
			Ct[currentCluster].TI.push_back(i);
			T[i].setLabel(currentCluster);
//			oflg<<" T"<<i;
		} else {
			unlabelledT.push_back(i);
		}
	}
//	oflg<<" &";
	for(size_t i=0; i<P.size(); ++i) {
		if(P[i].sStatus==empty) {
			Ct[currentCluster].PI.push_back(i);
			P[i].setLabel(currentCluster);
//			oflg<<" P"<<i;
		} else {
			unlabelledP.push_back(i);
		}
	}
//	oflg<<" into Ct"<<currentCluster<<endl;

	while(unlabelledT.size()>0 || unlabelledP.size()>0) {
		for(size_t i=0; i<Ct.size(); ++i) {
//			oflg<<"In Ct#"<<i<<":\t";
			for(size_t j=0; j<Ct[i].TI.size(); ++j) {
//				oflg<<"\tT"<<Ct[i].TI[j];
			}
//			oflg<<"\t|";
			for(size_t j=0; j<Ct[i].TI.size(); ++j) {
//				oflg<<"\tT"<<Ct[i].TI[j];
			}
		}
//		oflg<<endl<<"Unlabelled T:";
		for(size_t i=0; i<unlabelledT.size(); ++i) {
//			oflg<<"\t#"<<unlabelledT[i];
		}
//		oflg<<endl<<"Unlabelled P:";
		for(size_t i=0; i<unlabelledP.size(); ++i) {
//			oflg<<"\t#"<<unlabelledP[i];
		}
//		oflg<<endl;

		currentCluster=Ct.size();
//		oflg<<"Create a new empty cluster: Ct# "<<currentCluster<<" :"<<endl;
		Ct.push_back(cluster_c());

		TtoExam.push_back(unlabelledT[0]); //set a starting Throat and put it into TtoExam
		while(TtoExam.size()!=0 || PtoExam.size()!=0) {
//			oflg<<"Add";
			for(size_t i=0; i<TtoExam.size(); ++i) {
				Ct[currentCluster].TI.push_back(TtoExam[i]);
//				oflg<<" T"<<TtoExam[i];
				T[TtoExam[i]].setLabel(currentCluster);
				unlabelledT.erase(ItInUnlabelledT(TtoExam[i]));

				for(size_t j=0; j<T[TtoExam[i]].PI.size(); ++j) {
					size_t indexP(T[TtoExam[i]].PI[j]);
					if((P[indexP].sStatus!=empty) && (inUnlabelledP(indexP)) && (!inExaminedP(indexP)) && (!inPtoExam(indexP))) {
						PtoExam.push_back(indexP);
					}
				}
			}
//			oflg<<" into Ct"<<currentCluster<<endl;
			examinedT.clear();
			examinedT=TtoExam;
			TtoExam.clear();

//			oflg<<"Add";
			for(size_t i=0; i!=PtoExam.size(); ++i) {
				Ct[currentCluster].PI.push_back(PtoExam[i]);
//				oflg<<" P"<<PtoExam[i];
				P[PtoExam[i]].setLabel(currentCluster);
				unlabelledP.erase(ItInUnlabelledP(PtoExam[i]));

				for(vector<index_t>::size_type j=0; j<P[PtoExam[i]].TI.size(); j++) {
					index_t indexT(P[PtoExam[i]].TI[j]);
					if((T[indexT].sStatus!=empty) && (inUnlabelledT(indexT)) && (!inExaminedT(indexT)) && (!inTtoExam(indexT))) {
						TtoExam.push_back(indexT);
					}
				}
			}
//			oflg<<" into Ct"<<currentCluster<<endl;
			examinedP.clear();
			examinedP=PtoExam;
			PtoExam.clear();
		}
	}



    return true;
}

//----------------------------------------------------------------------------------------------------
bool inExaminedT(const size_t & index) {
	for(size_t i=0; i!=examinedT.size(); ++i) {
		if(index==examinedT[i]) return true;
	}
	return false;
}

bool inExaminedP(const size_t & index) {
	for(size_t i=0; i!=examinedP.size(); ++i) {
		if(index==examinedP[i]) return true;
	}
	return false;
}

bool inTtoExam(const size_t & index) {
	for(size_t i=0; i!=TtoExam.size(); ++i) {
		if(index==TtoExam[i]) return true;
	}
	return false;
}

bool inPtoExam(const size_t & index) {
	for(size_t i=0; i!=PtoExam.size(); ++i) {
		if(index==PtoExam[i]) return true;
	}
	return false;
}

bool inUnlabelledT(const size_t & index) {
	for(size_t i=0; i!=unlabelledT.size(); ++i) {
		if(index==unlabelledT[i]) return true;
	}
	return false;
}

bool inUnlabelledP(const size_t & index) {
	for(size_t i=0; i!=unlabelledP.size(); ++i) {
		if(index==unlabelledP[i]) return true;
	}
	return false;
}

vector<size_t>::iterator ItInUnlabelledT(const size_t & index) {
	for(vector<size_t>::iterator it=unlabelledT.begin(); it!=unlabelledT.end(); ++it) {
		if(index==*it) return it;
	}
	return unlabelledT.end();
}

vector<size_t>::iterator ItInUnlabelledP(const size_t & index) {
	for(vector<size_t>::iterator it=unlabelledP.begin(); it!=unlabelledP.end(); ++it) {
		if(index==*it) return it;
	}
	return unlabelledP.end();
}

/*----------------------------------------------------------------------------------------------------
 *
 *                                     Function calulateTimeToEmptyMax
 *
 ----------------------------------------------------------------------------------------------------*/
bool calculateFluxNTime(const size_t & cli) {
	cd_c::ClusterToEmpty.clear();
	for(size_t i=1; i<Ct.size(); i++) {
//		oflg<<"Cluster # "<<i<<":"<<endl;

		Ct[i].calcMassOutFlux(cli);
//		oflg<<"MassOutFlux="<<Ct[i].MassOutFlux<<"\t";

		Ct[i].findMaxThroat();
//		oflg<<"Throat Empty First: "<<Ct[i].maxThroat[0]<<"\t";

		Ct[i].calctimeToEmptyMax();
//		oflg<<"TimeToEmptyMax="<<Ct[i].timeToEmptyMax<<endl;

		if(i==1) {
			cd_c::TimeToEmptyMax=Ct[i].timeToEmptyMax;
			cd_c::ClusterToEmpty.push_back(i);
		} else if(Ct[i].timeToEmptyMax< cd_c::TimeToEmptyMax && Ct[i].timeToEmptyMax>0) {
			cd_c::TimeToEmptyMax=Ct[i].timeToEmptyMax;
			cd_c::ClusterToEmpty.clear();
			cd_c::ClusterToEmpty.push_back(i);
		} else if(Ct[i].timeToEmptyMax==cd_c::TimeToEmptyMax && Ct[i].timeToEmptyMax>0) {
			cd_c::ClusterToEmpty.push_back(i);
		} else {
			//do nothing
		}
	}

	return true;
}

/*----------------------------------------------------------------------------------------------------
 *
 *                                     Function calulateConcentration
 *
 ----------------------------------------------------------------------------------------------------*/

//----------------------------------------------------------------------------------------------------
bool calculateConcentration(const bool & cal, const bool & update) {
	if(cal) {
		if(cd_c::phy.filmEffect) {
			if(cd_c::phy.AlgorithmEvaporation==ConvectionDiffusion) {
				solvePhiEq(0, 5, 1);//cLevelIn=0; 5=1+4: Convection+DeltaVT; cLevelOut=1
				solvePhiEq(1, 6, 2);//cLevelIn=1; 6=2+4: Diffussion+DeltaVT; cLevelOut=2
			} else if(cd_c::phy.AlgorithmEvaporation==DiffusionConvection) {
				solvePhiEq(0, 6, 1);//cLevelIn=1; 6=2+4: Diffussion+DeltaVT; cLevelOut=1
				solvePhiEq(1, 5, 2);//cLevelIn=0; 5=1+4: Convection+DeltaVT; cLevelOut=2
			} else if(cd_c::phy.AlgorithmEvaporation==NonOperatorSplitting) {
				solvePhiEq(0, 7, 2);//cLevelIn=0; 7=1+2+4: Convection+Diffussion+DeltaVT; cLevelOut=2
			}

			for(size_t i=0; i<P.size(); ++i) {
				P[i].updatef2pr(6); //update the new copies
			}
		} else {
			if(cd_c::phy.AlgorithmEvaporation==ConvectionDiffusion) {
				solveConcentrationEq(0, 5, 1);//cLevelIn=0; 5=1+4: Convection+DeltaVT; cLevelOut=1
				solveConcentrationEq(1, 6, 2);//cLevelIn=1; 6=2+4: Diffussion+DeltaVT; cLevelOut=2
			} else if(cd_c::phy.AlgorithmEvaporation==DiffusionConvection) {
				solveConcentrationEq(0, 6, 1);//cLevelIn=1; 6=2+4: Diffussion+DeltaVT; cLevelOut=1
				solveConcentrationEq(1, 5, 2);//cLevelIn=0; 5=1+4: Convection+DeltaVT; cLevelOut=2
			} else if(cd_c::phy.AlgorithmEvaporation==NonOperatorSplitting) {
				solveConcentrationEq(0, 7, 2);//cLevelIn=0; 7=1+2+4: Convection+Diffussion+DeltaVT; cLevelOut=2
			}
		}
	} else {
		//do nothing
	}

	if(update) {
		if(cd_c::phy.filmEffect) {
			for(size_t i=0; i<C.size(); ++i) {
				C[i].p[0]=C[i].p[2];
			}
			for(size_t i=0; i<P.size(); ++i) {
				P[i].f[0]=P[i].f[2];
				P[i].updatef2pr(1);//update the first copy
			}
		} else {
			for(size_t i=0; i<C.size(); ++i) {
				C[i].p[0]=C[i].p[2];
			}
			for(size_t i=0; i<P.size(); ++i) {
				P[i].p[0]=P[i].p[2];
			}
		}
	} else {
		//do not update
	}
//
//
//	for(size_t k=0; k<porousBlocks; ++k) {
//		if(Bk[k].B[n].T==PorousOpen) {
//			for(size_t i=0; i<=Bk[k].Nx; ++i) {
//				oflg<<C[P[Bk[k].PI[i][Bk[k].Ny]].CI[0]].p[2]<<"\t";
//			}
//			oflg<<endl;
//			for(size_t i=0; i<=Bk[k].Nx; ++i) {
//				oflg<<"----------\t";
//			}
//			oflg<<endl;
//		}
//		for(index_t j=Bk[k].Ny; j>=0; --j) {
//			if(Bk[k].B[w].T==PorousOpen) {
//				oflg<<C[P[Bk[k].PI[0][j]].CI[0]].p[2]<<"\t|\t";
//			}
//			for(size_t i=0; i<=Bk[k].Nx; ++i) {
//				oflg<<P[Bk[k].PI[i][j]].p[2]<<"\t";
//			}
//			if(Bk[k].B[e].T==PorousOpen) {
//				oflg<<"|\t"<<C[P[Bk[k].PI[Bk[k].Ny][j]].CI[0]].p[2];
//			}
//			oflg<<endl;
//		}
//		if(Bk[k].B[s].T==PorousOpen) {
//			for(size_t i=0; i<=Bk[k].Nx; ++i) {
//				oflg<<"----------\t";
//			}
//			oflg<<endl;
//			for(size_t i=0; i<=Bk[k].Nx; ++i) {
//				oflg<<C[P[Bk[k].PI[i][0]].CI[0]].p[2]<<"\t";
//			}
//			oflg<<endl;
//		}
//	}

	return true;
}

//----------------------------------------------------------------------------------------------------
bool solveConcentrationEq(const size_t & cli, const std::bitset<EqComponent> & ec, const size_t & clo) {
	if(cell_c::phy.ImplicitEvaporation) {
		VectorXd X(CPsize), B(CPsize);

//		MatrixXd A(CPsize, CPsize);
//		for(size_t i=0; i<CPsize; ++i) {
//			for(size_t j=0; j<CPsize; ++j) {
//				A(i, j)=0;
//			}
//			X(i)=0;
//			B(i)=0;
//		}
//		for(size_t i=0; i<C.size(); ++i) {
//			C[i].calculateCoefficient(0, cli, ec);//uvli=0
//			B(i)=C[i].b;
//			A(i, i)=C[i].aP;
//			for(size_t j=0; j<C[i].CI.size(); ++j) {
//				A(i, C[i].CI[j])=C[i].a[j];
//			}
//			for(size_t j=0; j<C[i].PI.size(); ++j) {
//				A(i, C[i].PI[j]+C.size())=C[i].ap[j];
//			}
//		}
//		for(size_t i=0, I=C.size(); i<P.size(); ++i, ++I) {
//			P[i].calculateCoefficient(cli, ec);//uvli=0; 2: there is only diffussion
//			B(I)=P[i].b;
//			A(I, I)=P[i].aP;
//			for(size_t j=0; j<P[i].CI.size(); ++j) {
//				A(I, P[i].CI[j])=P[i].a[j];
//			}
//			for(size_t j=0; j<P[i].PI.size(); ++j) {
//				A(I, P[i].PI[j]+C.size())=P[i].ap[j];
//			}
//		}

//		SparseMatrix<double> SA=A.sparseView();
//		SparseLU<SparseMatrix<double>> solver;
//		solver.analyzePattern(SA);
//		solver.factorize(SA);
//		X= solver.solve(B);

//		X=A.colPivHouseholderQr().solve(B);

		int EachSize(0);
		if(cd_c::phy.SchemeEvaporation==Hybrid     ) EachSize=5;
		if(cd_c::phy.SchemeEvaporation==HayaseQUICK) EachSize=9;
		typedef Eigen::Triplet<double> Tri;
		vector<Tri> TripletList;
		TripletList.reserve(CPsize*EachSize);

		for(size_t i=0; i<C.size(); ++i) {
			C[i].calculateCoefficient(0, cli, ec);//uvli=0
			B(i)=C[i].b;
			TripletList.push_back(Tri(i, i, C[i].aP));
			for(size_t j=0; j<C[i].CI.size(); ++j) {
				TripletList.push_back(Tri(i, C[i].CI[j], C[i].a[j]));
			}
			for(size_t j=0; j<C[i].PI.size(); ++j) {
				TripletList.push_back(Tri(i, C[i].PI[j]+C.size(), C[i].ap[j]));
			}
		}

		for(size_t i=0, I=C.size(); i<P.size(); ++i, ++I) {
			P[i].calculateCoefficient(cli, ec);//uvli=0; 2: there is only diffussion
			B(I)=P[i].b;
			TripletList.push_back(Tri(I, I, P[i].aP));
			for(size_t j=0; j<P[i].CI.size(); ++j) {
				TripletList.push_back(Tri(I, P[i].CI[j], P[i].a[j]));
			}
			for(size_t j=0; j<P[i].PI.size(); ++j) {
				TripletList.push_back(Tri(I, P[i].PI[j]+C.size(), P[i].ap[j]));
			}
		}

		SparseMatrix<double, Eigen::ColMajor> A(CPsize, CPsize);
		A.setFromTriplets(TripletList.begin(), TripletList.end());
		SparseLU<SparseMatrix<double>> solver;
		solver.analyzePattern(A);
		solver.factorize(A);
		X= solver.solve(B);

//		SparseMatrix<double> CA=A.sparseView();
//		SparseMatrix<double, Eigen::RowMajor> SA(CA);
////		SparseMatrix<double, Eigen::RowMajor> SA=A.sparseView();
//		BiCGSTAB<SparseMatrix<double> > solver;
//		solver.compute(SA);
//		X = solver.solve(B);


		for(size_t i=0; i<C.size(); ++i) {
			C[i].p[clo]=X[i];
		}
		for(size_t i=0, I=C.size(); i<P.size(); ++i, ++I) {
			P[i].p[clo]=X[I];
		}
	} else {
		for(size_t i=0; i<C.size(); ++i) {
			C[i].calculateCoefficient(0, cli, ec);
			for(size_t j=0; j<C[i].CI.size(); ++j) {
				C[i].b-=C[C[i].CI[j]].p[cli]*C[i].a [j];
			}
			for(size_t j=0; j<C[i].PI.size(); ++j) {
				C[i].b-=P[C[i].PI[j]].p[cli]*C[i].ap[j];
			}
			C[i].p[clo]=C[i].b/C[i].aP;
		}
		for(size_t i=0; i<P.size(); ++i) {
			P[i].calculateCoefficient(cli, ec);
			for(size_t j=0; j<P[i].CI.size(); ++j) {
				P[i].b-=C[P[i].CI[j]].p[cli]*P[i].a [j];
			}
			for(size_t j=0; j<P[i].PI.size(); ++j) {
				P[i].b-=P[P[i].PI[j]].p[cli]*P[i].ap[j];
			}
			P[i].p[clo]=P[i].b/P[i].aP;
		}
	}

	return true;
}

//----------------------------------------------------------------------------------------------------
bool solvePhiEq(const size_t & cli, const std::bitset<EqComponent> & ec, const size_t & clo) {
	if(cell_c::phy.ImplicitEvaporation) {
		VectorXd X(CPsize), B(CPsize);

//		MatrixXd A(CPsize, CPsize);
//		for(size_t i=0; i<CPsize; ++i) {
//			for(size_t j=0; j<CPsize; ++j) {
//				A(i, j)=0;
//			}
//			X(i)=0;
//			B(i)=0;
//		}
//		for(size_t i=0; i<C.size(); ++i) {
//			C[i].calculateCoefficient(0, cli, ec);//uvli=0
//			B(i)=C[i].b;
//			A(i, i)=C[i].aP;
//			for(size_t j=0; j<C[i].CI.size(); ++j) {
//				A(i, C[i].CI[j])=C[i].a[j];
//			}
//			for(size_t j=0; j<C[i].PI.size(); ++j) {
//				A(i, C[i].PI[j]+C.size())=C[i].ap[j];
//			}
//		}
//		for(size_t i=0, I=C.size(); i<P.size(); ++i, ++I) {
//			P[i].calculatePhiCoefficient(cli, ec);//uvli=0; 2: there is only diffussion
//			B(I)=P[i].b;
//			A(I, I)=P[i].aP;
//			for(size_t j=0; j<P[i].CI.size(); ++j) {
//				A(I, P[i].CI[j])=P[i].a[j];
//			}
//			for(size_t j=0; j<P[i].PI.size(); ++j) {
//				A(I, P[i].PI[j]+C.size())=P[i].ap[j];
//			}
//		}

//		SparseMatrix<double> SA=A.sparseView();
//		SparseLU<SparseMatrix<double>>   solver;
//		solver.analyzePattern(SA);
//		solver.factorize(SA);
//		X= solver.solve(B);

//		X=A.colPivHouseholderQr().solve(B);

		int EachSize(0);
		if(cd_c::phy.SchemeEvaporation==Hybrid     ) EachSize=5;
		if(cd_c::phy.SchemeEvaporation==HayaseQUICK) EachSize=9;
		typedef Eigen::Triplet<double> Tri;
		vector<Tri> TripletList;
		TripletList.reserve(CPsize*EachSize);
		for(size_t i=0; i<C.size(); ++i) {
			C[i].calculateCoefficient(0, cli, ec);//uvli=0
			B(i)=C[i].b;
			TripletList.push_back(Tri(i, i, C[i].aP));
			for(size_t j=0; j<C[i].CI.size(); ++j) {
				TripletList.push_back(Tri(i, C[i].CI[j], C[i].a[j]));
			}
			for(size_t j=0; j<C[i].PI.size(); ++j) {
				TripletList.push_back(Tri(i, C[i].PI[j]+C.size(), C[i].ap[j]));
			}
		}
		for(size_t i=0, I=C.size(); i<P.size(); ++i, ++I) {
			P[i].calculatePhiCoefficient(cli, ec);//uvli=0; 2: there is only diffussion
			B(I)=P[i].b;
			TripletList.push_back(Tri(I, I, P[i].aP));
			for(size_t j=0; j<P[i].CI.size(); ++j) {
				TripletList.push_back(Tri(I, P[i].CI[j], P[i].a[j]));
			}
			for(size_t j=0; j<P[i].PI.size(); ++j) {
				TripletList.push_back(Tri(I, P[i].PI[j]+C.size(), P[i].ap[j]));
			}
		}
		SparseMatrix<double, Eigen::ColMajor> A(CPsize, CPsize);
		A.setFromTriplets(TripletList.begin(), TripletList.end());
		SparseLU<SparseMatrix<double>> solver;
		solver.analyzePattern(A);
		solver.factorize(A);
		X= solver.solve(B);

//		SparseMatrix<double> CA=A.sparseView();
//		SparseMatrix<double, Eigen::RowMajor> SA(CA);
////		SparseMatrix<double, Eigen::RowMajor> SA=A.sparseView();
//		BiCGSTAB<SparseMatrix<double> > solver;
//		solver.compute(SA);
//		X = solver.solve(B);

		for(size_t i=0; i<C.size(); ++i) {
			C[i].p[clo]=X[i];
		}
		for(size_t i=0, I=C.size(); i<P.size(); ++i, ++I) {
			P[i].f[clo]=X[I];
		}
	} else {
		for(size_t i=0; i<C.size(); ++i) {
			C[i].calculateCoefficient(0, cli, ec);
			for(size_t j=0; j<C[i].CI.size(); ++j) {
				C[i].b-=C[C[i].CI[j]].p[cli]*C[i].a [j];
			}
			for(size_t j=0; j<C[i].PI.size(); ++j) {
				C[i].b-=P[C[i].PI[j]].f[cli]*C[i].ap[j];
			}
			C[i].p[clo]=C[i].b/C[i].aP;
		}
		for(size_t i=0; i<P.size(); ++i) {
			P[i].calculatePhiCoefficient(cli, ec);
			for(size_t j=0; j<P[i].CI.size(); ++j) {
				P[i].b-=C[P[i].CI[j]].p[cli]*P[i].a [j];
			}
			for(size_t j=0; j<P[i].PI.size(); ++j) {
				P[i].b-=P[P[i].PI[j]].f[cli]*P[i].ap[j];
			}
			P[i].f[clo]=P[i].b/P[i].aP;
		}
	}

	return true;
}

/*----------------------------------------------------------------------------------------------------
 *
 *                                        Function frequentOutput
 *
  ----------------------------------------------------------------------------------------------------*/
bool frequentOutput(const size_t & sc) {
	stringstream iss("");
	iss<<sc;
	string siss(iss.str());
	string   osa(caseName+".A"+siss+".vtk"); string   ose(caseName+".E"+siss+".vtk");// string   oss(caseName+".S"+siss);
	ofstream ofa(osa.c_str())              ; ofstream ofe(ose.c_str())              ;// ofstream ofs(oss.c_str())       ;
	plotAllField(ofa)                      ; plotExtField(ofe)                      ;// plotSaturation(ofs)             ;
	ofa.close()                            ; ofe.close()                            ;// ofs.close()                     ;

	return true;
}

//----------------------------------------------------------------------------------------------------
bool plotAllField(ofstream & ofs) {//plot real value
	ofs<<"# vtk DataFile Version 2.0"<<endl;
	ofs<<"All Field"<<endl;
	ofs<<"ASCII"<<endl;
	ofs<<"DATASET UNSTRUCTURED_GRID"<<endl;
	ofs<<"POINTS\t"<<  Pt.size()<<"\tdouble"<<endl;
	for(size_t i=0; i<Pt.size(); ++i) {
		ofs<<Pt[i].x*cd_c::phy.RefLength<<"\t"
		   <<Pt[i].y*cd_c::phy.RefLength<<"\t"<<" 0  "<<endl;
	}

	size_t pptc(0);
	for(size_t i=0; i<P.size(); ++i) {
		pptc++;
		pptc+=P[i].PtI.size();
	}

	ofs<<"CELLS\t"<<C.size()+P.size()+T.size()<<"\t"<<C.size()*5+pptc+T.size()*5<<endl;
	for(size_t i=0; i<C.size(); ++i) {
		ofs<<"4\t"<<          C[i].FI[ne]<<"\t"<<          C[i].FI[nw]<<"\t"<<          C[i].FI[sw]<<"\t"<<          C[i].FI[se]<<endl;
	}
	for(size_t i=0; i<P.size(); ++i) {
		ofs<<P[i].PtI.size();
		for(size_t j=0; j<P[i].PtI.size(); ++j) {
			ofs<<"\t"<<P[i].PtI[j];
		}
		ofs<<endl;
	}
	for(size_t i=0; i<T.size(); ++i) {
		ofs<<"4\t"<<          T[i].PtI[0]<<"\t"<<          T[i].PtI[1]<<"\t"<<          T[i].PtI[2]<<"\t"<<          T[i].PtI[3]<<endl;
	}

	ofs<<"CELL_TYPES\t"<<C.size()+P.size()+T.size()<<endl;
	for(size_t i=0; i<C.size(); ++i) {
		ofs<<" 9"<<endl;
	}
	for(size_t i=0; i<P.size(); ++i) {
		if(P[i].PtI.size()==3) {
			ofs<<" 5"<<endl;
		} else if(P[i].PtI.size()==4) {
			ofs<<" 9"<<endl;
		} else {
			ofs<<" 7"<<endl;
		}
	}
	for(size_t i=0; i<T.size(); ++i) {
		ofs<<" 9"<<endl;
	}

	ofs<<"CELL_DATA\t"<<C.size()+P.size()+T.size()<<endl;
	ofs<<"SCALARS\t"<<"Concentration\t"<<"double\t"<<"1"<<endl;
	ofs<<"LOOKUP_TABLE\t"<<"default"<<endl;
	for(size_t i=0; i<C.size(); ++i) {
		ofs<<C[i].p[0]*cd_c::phy.RefConcentration<<endl;
	}
	for(size_t i=0; i<P.size(); ++i) {
		ofs<<P[i].p[0]*cd_c::phy.RefConcentration<<endl;
	}
	for(size_t i=0; i<T.size(); ++i) {
		ofs<<0.5*(T[i].PP[0]->p[0]+T[i].PP[1]->p[0])*cd_c::phy.RefConcentration<<endl;
	}

	ofs<<"SCALARS\t"<<"Saturation\t"<<"int\t"<<"1"<<endl;
	ofs<<"LOOKUP_TABLE\t"<<"default"<<endl;
	for(size_t i=0; i<C.size(); ++i) {
		ofs<<empty<<endl;
	}
	for(size_t i=0; i<P.size(); ++i) {
		if(cd_c::phy.filmEffect) {
			if(P[i].sStatus) {
				ofs<<P[i].sStatus<<endl;
			} else {
				ofs<<(P[i].f[0]>=P[i].PhiFGI?1:0)<<endl;
			}
		} else {
			ofs<<P[i].sStatus<<endl;
		}
	}
	for(size_t i=0; i<T.size(); ++i) {
		if(cd_c::phy.filmEffect) {
			if(T[i].sStatus) {
				ofs<<T[i].sStatus<<endl;
			} else {
				ofs<<((P[T[i].PI[0]].f[0]>=P[T[i].PI[0]].PhiFGI && P[T[i].PI[1]].f[0]>=P[T[i].PI[1]].PhiFGI)?1:0)<<endl;
			}
		} else {
			ofs<<T[i].sStatus<<endl;
		}
	}

	if(cd_c::phy.filmEffect) {
		ofs<<"SCALARS\t"<<"FilmRadius\t"<<"double\t"<<"1"<<endl;
		ofs<<"LOOKUP_TABLE\t"<<"default"<<endl;
		for(size_t i=0; i<C.size(); ++i) {
			ofs<<0<<endl;
		}
		for(size_t i=0; i<P.size(); ++i) {
			ofs<<( (P[i].f[0]<=P[i].PhiFGI)?0:( pow( (P[i].f[0]-P[i].PhiFGI)/P[i].PhiF, 1./3. )*cd_c::phy.RefLength ) )<<endl;
		}
		for(size_t i=0; i<T.size(); ++i) {
			ofs<<(
					(
						((P[T[i].PI[0]].f[0]<=P[T[i].PI[0]].PhiFGI)?0:( pow( (P[T[i].PI[0]].f[0]-P[T[i].PI[0]].PhiFGI)/P[T[i].PI[0]].PhiF, 1./3. ) ))
				       +((P[T[i].PI[1]].f[0]<=P[T[i].PI[1]].PhiFGI)?0:( pow( (P[T[i].PI[1]].f[0]-P[T[i].PI[1]].PhiFGI)/P[T[i].PI[1]].PhiF, 1./3. ) ))
					)
				 /2*cd_c::phy.RefLength
				 )<<endl;
		}
	}

	ofs<<"SCALARS\t"<<"Pressure\t"<<"double\t"<<"1"<<endl;
	ofs<<"LOOKUP_TABLE\t"<<"default"<<endl;
	for(size_t i=0; i<C.size(); ++i) {
		ofs<<C[i].p[3]*cd_c::phy.RefPressure<<endl;
	}
	for(size_t i=0; i<P.size(); ++i) {
		ofs<<0<<endl;
	}
	for(size_t i=0; i<T.size(); ++i) {
		ofs<<0<<endl;
	}

	ofs<<"VECTORS\t"<<"Velocity\t"<<"double"<<endl;
	for(size_t i=0; i<C.size(); ++i) {
		ofs<<0.5*(U[C[i].FI[e]].p[0]+U[C[i].FI[w]].p[0])*cd_c::phy.RefVelocity<<"\t"
		   <<0.5*(V[C[i].FI[n]].p[0]+V[C[i].FI[s]].p[0])*cd_c::phy.RefVelocity<<"\t"<<"0"<<endl;
	}
	for(size_t i=0; i<P.size(); ++i) {
		ofs<<0<<"\t"
		   <<0<<"\t"<<"0"<<endl;
	}
	for(size_t i=0; i<T.size(); ++i) {
		ofs<<0<<"\t"
		   <<0<<"\t"<<"0"<<endl;
	}

	return true;
}

//----------------------------------------------------------------------------------------------------
bool plotExtField(ofstream & ofs) {
	ofs<<"# vtk DataFile Version 2.0"<<endl;
	ofs<<"External Field"<<endl;
	ofs<<"ASCII"<<endl;
	ofs<<"DATASET UNSTRUCTURED_GRID"<<endl;
	ofs<<"POINTS\t"<<ExtPtSize<<"\tdouble"<<endl;
	for(size_t i=0; i<ExtPtSize; ++i) {
		ofs<<Pt[i].x*cd_c::phy.RefLength<<"\t"
		   <<Pt[i].y*cd_c::phy.RefLength<<"\t"<<" 0  "<<endl;
	}

	ofs<<"CELLS\t"<<C.size()<<"\t"<<C.size()*5<<endl;
	for(size_t i=0; i<C.size(); ++i) {
		ofs<<"4\t"<<          C[i].FI[ne]<<"\t"<<          C[i].FI[nw]<<"\t"<<          C[i].FI[sw]<<"\t"<<          C[i].FI[se]<<endl;
	}

	ofs<<"CELL_TYPES\t"<<C.size()<<endl;
	for(size_t i=0; i<C.size(); ++i) {
		ofs<<" 9"<<endl;
	}

	ofs<<"CELL_DATA\t"<<C.size()<<endl;
	ofs<<"SCALARS\t"<<"Concentration\t"<<"double\t"<<"1"<<endl;
	ofs<<"LOOKUP_TABLE\t"<<"default"<<endl;
	for(size_t i=0; i<C.size(); ++i) {
		ofs<<C[i].p[0]*cd_c::phy.RefConcentration<<endl;
	}

	ofs<<"SCALARS\t"<<"Saturation\t"<<"int\t"<<"1"<<endl;
	ofs<<"LOOKUP_TABLE\t"<<"default"<<endl;
	for(size_t i=0; i<C.size(); ++i) {
		ofs<<empty<<endl;
	}

	ofs<<"SCALARS\t"<<"Pressure\t"<<"double\t"<<"1"<<endl;
	ofs<<"LOOKUP_TABLE\t"<<"default"<<endl;
	for(size_t i=0; i<C.size(); ++i) {
		ofs<<C[i].p[3]*cd_c::phy.RefPressure<<endl;
	}

	ofs<<"VECTORS\t"<<"Velocity\t"<<"double"<<endl;
	for(size_t i=0; i<C.size(); ++i) {
		ofs<<0.5*(U[C[i].FI[e]].p[0]+U[C[i].FI[w]].p[0])*cd_c::phy.RefVelocity<<"\t"
		   <<0.5*(V[C[i].FI[n]].p[0]+V[C[i].FI[s]].p[0])*cd_c::phy.RefVelocity<<"\t"<<"0"<<endl;
	}

	return true;
}

//----------------------------------------------------------------------------------------------------
bool plotSaturation(ofstream & ofs) {
	ofs<<"PoreCount=   "<<P.size()<<endl;
	ofs<<"ThroatCount= "<<T.size()<<endl;

	for(size_t i=0; i<P.size(); ++i) {
		ofs<<i<<"\t"<<P[i].x<<"\t"<<P[i].y<<"\t"<<(P[i].sStatus==empty?0:1)<<endl;
	}

	for(size_t i=0; i<T.size(); ++i) {
		ofs<<i<<"\t"<<T[i].PI[0]<<"\t"<<T[i].PI[1]<<"\t"<<T[i].A<<"\t"<<T[i].L<<"\t"<<(T[i].isEmpty()?0:1)<<endl;
	}

	return true;
}

/*----------------------------------------------------------------------------------------------------
 *
 *                                        Function emptyingOneThroat
 *
  ----------------------------------------------------------------------------------------------------*/
bool emptyingOneThroat() {



	return true;
}



/*----------------------------------------------------------------------------------------------------
 *
 *                                        Function plotEvaporateRateHistory
 *
 ----------------------------------------------------------------------------------------------------*/
bool plotEvaporateRateHistory(const ofstream & ofs) {

	return true;
}

/*----------------------------------------------------------------------------------------------------
 *
 *                                        Function plotDryingHistory
 *
 ----------------------------------------------------------------------------------------------------*/
bool plotDryingHistory(const ofstream & ofs) {

	return true;
}

bool liquidExistYet() {
	for(size_t i=0; i<T.size(); i++) {
		if(T[i].sStatus!=empty) return true;
	}

	return false;
}

bool wetYet() {
	for(size_t i=0; i<P.size(); i++) {
		if(P[i].p[0]>0) return true;
	}

	return false;
}
numeric_t TotMassFluxOut() {
	numeric_t tmfo(0);
	for(size_t i=0; i<pore_c::openPI.size(); ++i) {
		for(size_t j=0; j<P[pore_c::openPI[i]].CI.size(); ++j) {

		}
	}

	return tmfo;
}

/*----------------------------------------------------------------------------------------------------
 *
 *                                        Function calcVaporPressure
 *
 ----------------------------------------------------------------------------------------------------*/
bool calcVaporPressure() {
//	unknownP.clear();
//	for(size_t i=0; i!=P.size(); ++i) {
//		if(P[i].vpStatus==unknown) unknownP.push_back(i);
//	}
//	size_t unSize=unknownP.size()+C.size();
//
//	if(unSize>0) {
//		MatrixXd A(unSize, unSize);
//		VectorXd X(unSize);
//		VectorXd B(unSize);
//
//		for(size_t i=0; i<unSize; i++) {
//			X(i)=0.;
//			B(i)=0.;
//			for(size_t j=0; j<unSize; j++) {
//				A(i,j)=0.;
//			}
//		}
//
//		for(size_t i=0; i!=unSize; ++i) {
//
//
//			for(size_t j=0; j<P[unknownP[i]].PI.size(); j++) {
//				double g=T[P[unknownP[i]].TI[j]].g;
//
//				G(i,i)+=g;
//
//				if(P[P[unknownP[i]].PI[j]].vpStatus==unknown) {
//					G(i,unknownPi(P[unknownP[i]].PI[j]))=-g;
//				} else if(P[P[unknownP[i]].PI[j]].vpStatus==saturated || P[P[unknownP[i]].PI[j]].vpStatus==out) {
//					double lnp=P[P[unknownP[i]].PI[j]].lnp();
//					PA(i)+=g*lnp;
//				}
//			}
//		}
//
//		LNP=G.colPivHouseholderQr().solve(PA);
//
//		for(vector<index_t>::size_type k=0; k<unSize; k++) {
//			numeric_t x=LNP(k);
//
//			P[unknownP[k]].lnp2vp(x);
//
//			P[unknownP[k]].vpStatus=found;
//		}
//	}

	return true;
}

//----------------------------------------------------------------------------------------------------
size_t unknownPi(size_t & index) {
	for(size_t i=0; i!=unknownP.size(); ++i) {
		if(index==unknownP[i]) return i;
	}
	return unknownP.size();
}

/*----------------------------------------------------------------------------------------------------
 *
 *                                        cluster_c member functions
 *
 ----------------------------------------------------------------------------------------------------*/
bool cluster_c::clear() {
	PI.clear();
	TI.clear();

	maxThroat.clear();

	MassOutFlux=0;
	timeToEmptyMax=0;
	MassOut=0;

	return true;
}

numeric_t cluster_c::mass() {
    numeric_t m(0);
	for(size_t i=0; i<TI.size(); i++) {
		if(T[TI[i]].sStatus!=empty) m+=T[TI[i]].M;
	}
	return m;
}

bool cluster_c::calcMassOutFlux(const size_t & cli) {
	MassOutFlux=0;
	for(size_t i=0; i<TI.size(); i++) {
		if(T[TI[i]].isMeniscus()) {
//			oflg<<"T#"<<TI[i]<<" is meniscus, and the MassOutFlux="<<T[TI[i]].MassOutFlux(cli);
			MassOutFlux+=T[TI[i]].MassOutFlux(cli);
//			oflg<<"; so cluster MassOutFlux="<<MassOutFlux<<endl;
		}
	}
	if(MassOutFlux<=0) {
//		cerr<<"MassOutFlux="<<MassOutFlux<<endl;
//		MassOutFlux=1e-20;
		MassOutFlux=0;
	}
	return true;
}

bool cluster_c::findMaxThroat() {
	maxThroat.clear();
	numeric_t maxRI(0);
	for(size_t i=0; i<TI.size(); i++) {
		if(T[TI[i]].isMeniscus()) {
			if(T[TI[i]].RI>maxRI) {
				maxThroat.clear();
				maxThroat.push_back(TI[i]);
				maxRI=T[TI[i]].RI;
			} else if(T[TI[i]].RI==maxRI) {
				maxThroat.push_back(TI[i]);
			}
		}
	}

	return true;
}

inline numeric_t cluster_c::DividedMassOutFlux() {return MassOutFlux/maxThroat.size();}

numeric_t cluster_c::maxMass() {
	numeric_t minMass(T[maxThroat[0]].M);
	index_t temp(voidIndex);
	for(size_t i=1; i<maxThroat.size(); ++i) {
		if(T[maxThroat[i]].M<minMass) {
			temp=maxThroat[0];
			maxThroat[0]=maxThroat[i];
			maxThroat[i]=temp;

			minMass=T[maxThroat[0]].M;
		}
	}

	return minMass;
}



bool cluster_c::calctimeToEmptyMax() {
	timeToEmptyMax=(maxMass()-cd_c::phy.EnvironmentConcentration*T[maxThroat[0]].V)/DividedMassOutFlux();
	return true;
}

bool cluster_c::drying(const numeric_t & time, const bool GotEmptyItNow) {
	MassOut=mass();
	if(GotEmptyItNow) {
		for(size_t i=0; i<maxThroat.size(); ++i) {
			if(i==0) {
				T[maxThroat[i]].M =cd_c::phy.EnvironmentConcentration*T[maxThroat[i]].V;
				T[maxThroat[i]].sStatus=empty;
			} else {
				T[maxThroat[i]].M-=DividedMassOutFlux()*cd_c::EDT;
				if(T[maxThroat[i]].isEmpty()) {
					T[maxThroat[i]].sStatus=empty;
					T[maxThroat[i]].M=cd_c::phy.EnvironmentConcentration*T[maxThroat[i]].V;
				}
			}
		}
	} else {
		for(size_t i=0; i<maxThroat.size(); ++i) {
			T[maxThroat[i]].M-=DividedMassOutFlux()*time;
			if(T[maxThroat[i]].isEmpty()) {
				T[maxThroat[i]].sStatus=empty;
				T[maxThroat[i]].M=cd_c::phy.EnvironmentConcentration*T[maxThroat[i]].V;
			}
		}
	}

//	for(size_t i=0; i<maxThroat.size(); ++i) {
//		T[maxThroat[i]].M-=DividedMassOutFlux()*time;
//		if(T[maxThroat[i]].isEmpty()) {
//			T[maxThroat[i]].sStatus=empty;
//		}
//	}

	MassOut-=mass();
	return true;
}

