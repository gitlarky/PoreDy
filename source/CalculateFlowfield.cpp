/*
 * CalculateFlowfield.cpp
 *
 *  Created on: Apr 15, 2014
 *      Author: david
 */


#include "GlobalDeclaration.h"

using namespace std;
using namespace Eigen;

std::ofstream offh;

bool solveMomentum(const size_t & uvli, const size_t & pli, const size_t & uvlo);
bool solveContinuity(const size_t & uvli, const size_t & pli, const std::bitset<VariableCount> & vc
		           , const size_t & uvlo, const size_t & plo);
void updateMassRatio(const size_t & uvli);
void updateBCForNewMR();
void updateBUVForNewMR(const size_t & uvlio);
void calcFlowfieldRelaxFactor(const numeric_t & avgUVError, const numeric_t & avgUVPError, const numeric_t & avgError);

void outputFlowfield(const size_t & n);
void plotFlowfield(ofstream & ofs);


/*====================================================================================================

                                       function calculateFlowfield

====================================================================================================*/
bool calculateFlowfield() {
	cout<<"Start calculating flow field......"<<endl;

	string FFHistoryName(caseName+".fch");
	offh.open(FFHistoryName.c_str());
	offh<<scientific<<setprecision(3)<<setw(10);

	numeric_t avgUError(0), avgVError(0), avgPError(0), avgMError(0), avgUVError(0), avgUVPError(0), avgError(0);
	numeric_t avgUErr0r(0), avgVErr0r(0), avgPErr0r(0), avgMErr0r(0);

	vector<numeric_t> U0(U.size(), 0), V0(V.size(), 0), P0(C.size(), 0);
	size_t stepCount(0);
	size_t uvLevelIn(0), uvLevelOut(1), pLevelIn(3), pLevelOut(4);//u, v ,p values, which to use and where to put results

	cout<<"Start Iterating......"<<endl;
	do {
		if(stepCount%cd_c::phy.FlowOutputFrequency==0) outputFlowfield(stepCount);

        calcFlowfieldRelaxFactor(avgUVError, avgUVPError, avgError);

		if(cd_c::phy.AlgorithmFlowfield==SIMPLE) {
			for(size_t i=0; i<U.size(); ++i) {
				U[i].calculateCoefficient(uvLevelIn, pLevelIn, 15);//uvLevelIn=0, pLevelIn=0; 15=1+2+4+8: Convection+Diffussion+DeltaPA+r; uvLevelOut=1
			}
			for(size_t i=0; i<V.size(); ++i) {
				V[i].calculateCoefficient(uvLevelIn, pLevelIn, 15);//uvLevelIn=0, pLevelIn=0; 15=1+2+4+8: Convection+Diffussion+DeltaPA+r; uvLevelOut=1
			}
		} else if(cd_c::phy.AlgorithmFlowfield==SIMPLER) {
			for(size_t i=0; i<U.size(); ++i) {
				U[i].calculateCoefficient(uvLevelIn, pLevelIn, 11);//uvLevelIn=0, pLevelIn=0; 11=1+2+8: Convection+Diffussion+r; uvLevelOut=1

				U[i].p[uvLevelOut]=U[i].b;
				for(size_t j=0; j<U[i].CI.size(); ++j) {
					U[i].p[uvLevelOut]-=U[i].a[j]*U[U[i].CI[j]].p[uvLevelIn];
				}
				U[i].p[uvLevelOut]/=U[i].aP;//Notice: There is no DeltaPA item

//				oflg<<"U # "<<i<<"\tCoefficient: aP="<<U[i].aP<<"\tb="<<U[i].b<<"\tr="<<U[i].r;
//				for(size_t j=0; j<U[i].CI.size(); ++j) {
//					oflg<<"\ta["<<U[i].CI[j]<<"]="<<U[i].a[j];
//				}
//				oflg<<endl;
			}
			for(size_t i=0; i<V.size(); ++i) {
				V[i].calculateCoefficient(uvLevelIn, pLevelIn, 11);//uvLevelIn=0, pLevelIn=0; 11=1+2+8: Convection+Diffussion+r; uvLevelOut=1

				V[i].p[uvLevelOut]=V[i].b;
				for(size_t j=0; j<V[i].CI.size(); ++j) {
					V[i].p[uvLevelOut]-=V[i].a[j]*V[V[i].CI[j]].p[uvLevelIn];
				}
				V[i].p[uvLevelOut]/=V[i].aP;//Notice: There is no DeltaPA item

//				oflg<<"V # "<<i<<"\tCoefficient: aP="<<V[i].aP<<"\tb="<<V[i].b<<"\tr="<<V[i].r;
//				for(size_t j=0; j<V[i].CI.size(); ++j) {
//					oflg<<"\ta["<<V[i].CI[j]<<"]="<<V[i].a[j];
//				}
//				oflg<<endl;
			}

			uvLevelIn=uvLevelOut;
			for(size_t i=0; i<C.size(); ++i) {
				C[i].calculatePCoefficient(uvLevelIn, pLevelIn, 4);//4 only update p, that is calculate a new p not p'
			}

			solveContinuity(uvLevelIn, pLevelIn, 4, uvLevelOut, pLevelOut);//u v p c; 4 is for p, only update p, that is calculate a new p not p'

			uvLevelIn=0;
			pLevelIn=pLevelOut;
			for(size_t i=0; i<U.size(); ++i) {
				U[i].calculateCoefficient(uvLevelIn, pLevelIn, 4);//uvLevelIn=0, pLevelIn=4; 4: only calculate DeltaPA
			}
			for(size_t i=0; i<V.size(); ++i) {
				V[i].calculateCoefficient(uvLevelIn, pLevelIn, 4);//uvLevelIn=0, pLevelIn=4; 4: only calculate DeltaPA
			}
		} else {
			//Do nothing
		}

//		oflg<<"Solve momentum equations......"<<endl;
		solveMomentum(uvLevelIn, pLevelIn, uvLevelOut);//uvLevelIn=0, pLevelIn=3; uvLevelOut=1

//		for(size_t i=0; i<U.size(); ++i) {
//			oflg<<"U # "<<i<<" : U="<<U[i].p[uvLevelOut]<<endl;
//		}
//		for(size_t i=0; i<V.size(); ++i) {
//			oflg<<"V # "<<i<<" : V="<<V[i].p[uvLevelOut]<<endl;
//		}
//		for(size_t i=0; i<C.size(); ++i) {
//			oflg<<"C # "<<i<<" : P="<<C[i].p[ pLevelIn ]<<endl;
//		}

		uvLevelIn=uvLevelOut; uvLevelOut=2; pLevelOut=5;//uvLevelIn=1, pLevelIn=3; uvLevelOut=2

//		oflg<<"Calculate Pressure Correction Coefficient......"<<endl;

		for(size_t i=0; i<C.size(); ++i) {
			C[i].calculatePCoefficient(uvLevelIn, pLevelIn, 7);
		}

//		for(size_t i=0; i<C.size(); ++i) {
//			oflg<<"C # "<<i<<"\t Coefficient: aP="<<C[i].aP<<"\t b="<<C[i].b;
//			for(size_t j=0; j<C[i].CI.size(); ++j) {
//				oflg<<"\t a["<<C[i].CI[j]<<"]="<<C[i].a[j];
//			}
//			oflg<<endl;
//		}

		if(cd_c::phy.AlgorithmFlowfield==SIMPLE) {
			solveContinuity(uvLevelIn, pLevelIn, 7, uvLevelOut, pLevelOut);//uvLevelIn=1, pLevelIn=3/4; 7=1+2+4: update u v p; uvLevelOut=2, pLevelOut=5
		} else if(cd_c::phy.AlgorithmFlowfield==SIMPLER) {
			solveContinuity(uvLevelIn, pLevelIn, 3, uvLevelOut, pLevelOut);//uvLevelIn=1, pLevelIn=3/4; 3=1+2  : update u v  ; uvLevelOut=2, pLevelOut=5
		} else {
			//do nothing
		}


		updateMassRatio(uvLevelOut);
//		oflg<<"MassIn="<<uvCell_c::MassIn<<"\tMassOut="<<uvCell_c::MassOut<<"\tMassRatio="<<uvCell_c::MassRatio<<endl;
		updateBCForNewMR();
//		updateBUVForNewMR(uvLevelOut);

//		oflg<<"Calculating Errors & Recording Convergence Hisotry......"<<endl;
		avgUError=0; avgVError=0; avgPError=0; avgError=0;
		for(size_t i=0; i<U.size(); ++i) {
			avgUError+=abs(U[i].p[2]-U[i].p[0]);
//			oflg<<"U"<<i<<" Error: "<<U[i].p[2]-U[i].p[0]<<endl;
		}
		avgUError/=U.size();
		for(size_t i=0; i<V.size(); ++i) {
			avgVError+=abs(V[i].p[2]-V[i].p[0]);
//			oflg<<"V"<<i<<" Error: "<<V[i].p[2]-V[i].p[0]<<endl;
		}
		avgUError/=V.size();
		for(size_t i=0; i<C.size(); ++i) {
			avgPError+=abs(C[i].p[5]-C[i].p[3]);
//			oflg<<"C"<<i<<" Error: "<<C[i].p[5]-C[i].p[3]<<endl;
		}
		avgPError/=C.size();
		avgMError=fabs(uvCell_c::MassRatio-1);

		if(stepCount==0) {
			avgUErr0r=avgUError; avgVErr0r=avgVError; avgPErr0r=avgPError; avgMErr0r=avgMError;
			cout<<"AbsoluteValue: "<<"AvgUErr0r="<<avgUErr0r<<"\tAvgVErr0r="<<avgVErr0r<<"\tAvgPErr0r="<<avgPErr0r<<"\tAvgMErr0r="<<avgMErr0r<<endl;
			offh<<"AbsoluteValue: "<<"AvgUErr0r="<<avgUErr0r<<"\tAvgVErr0r="<<avgVErr0r<<"\tAvgPErr0r="<<avgPErr0r<<"\tAvgMErr0r="<<avgMErr0r<<endl;
		}
		avgUError/=avgUErr0r; avgVError/=avgVErr0r; avgPError/=avgPErr0r;
//		avgMError/=avgMErr0r;
		avgUVError =max(avgUError, avgVError  );
		avgUVPError=max(avgPError, avgUVError );

		avgError   =avgUVPError;
//		avgError   =max(avgMError, avgUVPError);

		if(stepCount%10==0) cout<<"Step # "<<"\tAvgUError"<<"\tAvgVError"<<"\tAvgPError"<<"\tAvgMError"<<endl;
		if(stepCount   ==0) offh<<"Step # "<<"\tAvgUError"<<"\tAvgVError"<<"\tAvgPError"<<"\tAvgMError"<<endl;
		stepCount++;
		cout<<stepCount<<"\t"<<avgUError<<"\t"<<avgVError<<"\t"<<avgPError<<"\t"<<avgMError<<endl;
		offh<<stepCount<<"\t"<<avgUError<<"\t"<<avgVError<<"\t"<<avgPError<<"\t"<<avgMError<<endl;

//		oflg<<"Updating the values in each cells with new one......"<<endl;
//		if(cd_c::phy.UVRelaxFactor==0 && cd_c::phy.PRelaxFactor==0) calcFlowfieldRelaxFactor(avgUVError, avgUVPError, avgError);
		for(size_t i=0; i<U.size(); ++i) {
//			U[i].p[0]=(1-cd_c::phy.UVRelaxFactor)*U[i].p[0]+cd_c::phy.UVRelaxFactor*U[i].p[2];
			U[i].p[0]=U[i].p[2];
		}
		for(size_t i=0; i<V.size(); ++i) {
//			V[i].p[0]=(1-cd_c::phy.UVRelaxFactor)*V[i].p[0]+cd_c::phy.UVRelaxFactor*V[i].p[2];
			V[i].p[0]=V[i].p[2];
		}

		for(size_t i=0; i<C.size(); ++i) {
//			C[i].p[3]=(1-cd_c::phy. PRelaxFactor)*C[i].p[3]+cd_c::phy. PRelaxFactor*C[i].p[5];
			C[i].p[3]=C[i].p[5];
		}

//		oflg<<"Reset LevelIn & LevelOut......"<<endl;
		uvLevelIn=0; pLevelIn=3; uvLevelOut=1; pLevelOut=4;

//		if(stepCount%cd_c::phy.DryOutputFrequency==0) outputFlowfield(stepCount);
	} while(avgError>cd_c::phy.FlowCriteria && stepCount<cd_c::phy.FlowfieldMaxStep);


	offh.close();

	if(avgError<=cd_c::phy.FlowCriteria) {
		cout<<"Flowfield Calculation Converged!"<<endl;
	} else {
		cout<<"Flowfield Calculation Reached Max Steps!"<<endl;
	}
	outputFlowfield(stepCount);

	return true;
}

/*----------------------------------------------------------------------------------------------------

                                       function solveMomentum

----------------------------------------------------------------------------------------------------*/
bool solveMomentum(const size_t & uvli, const size_t & pli, const size_t & uvlo) {
	if(cell_c::phy.ImplicitFlowfield) {
		VectorXd X(UVsize), B(UVsize);

//		MatrixXd A(UVsize, UVsize);
//
//		for(size_t i=0; i<UVsize; ++i) {
//			for(size_t j=0; j<UVsize; ++j) {
//				A(i, j)=0;
//			}
//			X(i)=0;
//			B(i)=0;
//		}
//
//		for(size_t i=0            ; i<U.size(); ++i     ) {
//			B(i)=U[i].b+U[i].DeltaPA;
//			A(i, i)=U[i].aP;
//			for(size_t j=0; j<U[i].CI.size(); ++j) {
//				A(i, U[i].CI[j])=U[i].a[j];
//			}
//		}
//		for(size_t i=0, I=U.size(); i<V.size(); ++i, ++I) {
//			B(I)=V[i].b+V[i].DeltaPA;
//			A(I, I)=V[i].aP;
//			for(size_t j=0; j<V[i].CI.size(); ++j) {
//				A(I, V[i].CI[j]+U.size())=V[i].a[j];
//			}
//		}

//		X=A.colPivHouseholderQr().solve(B);
//		X=A.fullPivLu().solve(B);
//		X=A.fullPivHouseholderQr().solve(B);

//		oflg<<"AAAAAAAAAAAAAAAAAAAAA"<<endl;
//		for(size_t i=0; i<UVsize; ++i) {
//			for(size_t j=0; j<UVsize; ++j) {
//				oflg<<A(i, j)<<"\t";
//			}
//			oflg<<"|\t"<<B(i)<<"\t|\t"<<X(i);
//			oflg<<endl;
//		}

//		SparseMatrix<double> SA=A.sparseView();



		int EachSize(0);
		if(cd_c::phy.SchemeFlowfield==Hybrid     ) EachSize=5;
		if(cd_c::phy.SchemeFlowfield==HayaseQUICK) EachSize=9;
		typedef Eigen::Triplet<double> Tri;
		vector<Tri> TripletList;
		TripletList.reserve(UVsize*EachSize);
		for(size_t i=0            ; i<U.size(); ++i     ) {
			B(i)=U[i].b+U[i].DeltaPA;
			TripletList.push_back(Tri(i, i, U[i].aP));
			for(size_t j=0; j<U[i].CI.size(); ++j) {
				TripletList.push_back(Tri(i, U[i].CI[j], U[i].a[j]));
			}
		}
		for(size_t i=0, I=U.size(); i<V.size(); ++i, ++I) {
			B(I)=V[i].b+V[i].DeltaPA;
			TripletList.push_back(Tri(I, I, V[i].aP));
			for(size_t j=0; j<V[i].CI.size(); ++j) {
				TripletList.push_back(Tri(I, V[i].CI[j]+U.size(), V[i].a[j]));
			}
		}

//		for(size_t i=0            ; i<U.size(); ++i     ) {
//			B(i)=U[i].b+U[i].DeltaPA;
//			TripletList.push_back(Tri(i, i, U[i].aP));
//			for(size_t j=0; j<U[i].CI.size(); ++j) {
//				TripletList.push_back(Tri(U[i].CI[j], i, U[i].a[j]));
//			}
//		}
//		for(size_t i=0, I=U.size(); i<V.size(); ++i, ++I) {
//			B(I)=V[i].b+V[i].DeltaPA;
//			TripletList.push_back(Tri(I, I, V[i].aP));
//			for(size_t j=0; j<V[i].CI.size(); ++j) {
//				TripletList.push_back(Tri(V[i].CI[j]+U.size(), I, V[i].a[j]));
//			}
//		}

		SparseMatrix<double, Eigen::ColMajor> A(UVsize, UVsize);
		A.setFromTriplets(TripletList.begin(), TripletList.end());
		SparseLU<SparseMatrix<double>> solver;
		solver.analyzePattern(A);
		solver.factorize(A);
		X= solver.solve(B);

//		SparseMatrix<double, Eigen::RowMajor> A(UVsize, UVsize);
//		A.setFromTriplets(TripletList.begin(), TripletList.end());
//		BiCGSTAB<SparseMatrix<double> > solver;
//		solver.compute(A);
//		X = solver.solve(B);



//		oflg<<"AAAAAAAAAAAAAAAAAAAAA"<<endl;
//		for(size_t i=0; i<UVsize; ++i) {
//			for(size_t j=0; j<UVsize; ++j) {
//				oflg<<A(i, j)<<"\t";
//			}
//			oflg<<"|\t"<<B(i)<<"\t|\t"<<X(i);
//			oflg<<endl;
//		}

		for(size_t i=0            ; i<U.size(); ++i      ) {
			U[i].p[uvlo]=X(i);
		}
		for(size_t i=0, I=U.size(); i<V.size(); ++i, ++I) {
			V[i].p[uvlo]=X(I);
		}
	} else {
		for(size_t i=0; i<U.size(); ++i) {
			U[i].p[uvlo]=U[i].b+U[i].DeltaPA;
			for(size_t i=0; i<U[i].CI.size(); ++i) {
				U[i].p[uvlo]-=U[i].a[i]*U[U[i].CI[i]].p[uvli];
			}
			U[i].p[uvlo]/=U[i].aP;
		}//Notice: There is every item
		for(size_t i=0; i<V.size(); ++i) {
			V[i].p[uvlo]=V[i].b+V[i].DeltaPA;
			for(size_t i=0; i<V[i].CI.size(); ++i) {
				V[i].p[uvlo]-=V[i].a[i]*V[V[i].CI[i]].p[uvli];
			}
			V[i].p[uvlo]/=V[i].aP;
		}//Notice: There is every item
	}

	return true;
}
/*----------------------------------------------------------------------------------------------------

                                       function solveContinuity

----------------------------------------------------------------------------------------------------*/
bool solveContinuity(const size_t & uvli, const size_t & pli, const std::bitset<VariableCount> & vc
		           , const size_t & uvlo, const size_t & plo) {
	VectorXd X(C.size()), B(C.size());

//	MatrixXd A(C.size(), C.size());
//	for(size_t i=0; i<C.size(); ++i) {
//		for(size_t j=0; j<C.size(); ++j) {
//			A(i, j)=0;
//		}
//		X(i)=0;
//		B(i)=0;
//	}
//
//	for(size_t i=0; i<C.size(); ++i) {
//		B(i)=C[i].b;
//		A(i, i)=C[i].aP;
//		for(size_t j=0; j<C[i].CI.size(); ++j) {
//			A(i, C[i].CI[j])=C[i].a[j];
//		}
//	}
//
//	X=A.colPivHouseholderQr().solve(B);
//	X=A.fullPivLu().solve(B);
//	X=A.fullPivHouseholderQr().solve(B);

	int EachSize(0);
	if(cd_c::phy.SchemeFlowfield==Hybrid     ) EachSize=5;
	if(cd_c::phy.SchemeFlowfield==HayaseQUICK) EachSize=9;
	typedef Eigen::Triplet<double> Tri;
	vector<Tri> TripletList;
	TripletList.reserve(C.size()*EachSize);
	for(size_t i=0; i<C.size(); ++i) {
		B(i)=C[i].b;
		TripletList.push_back(Tri(i, i, C[i].aP));
		for(size_t j=0; j<C[i].CI.size(); ++j) {
			TripletList.push_back(Tri(i, C[i].CI[j], C[i].a[j]));
		}
	}

//	for(size_t i=0; i<C.size(); ++i) {
//		B(i)=C[i].b;
//		TripletList.push_back(Tri(i, i, C[i].aP));
//		for(size_t j=0; j<C[i].CI.size(); ++j) {
//			TripletList.push_back(Tri(C[i].CI[j], i, C[i].a[j]));
//		}
//	}



//	oflg<<"ContinuityContinuityContinuityContinuityContinuityContinuityContinuityContinuityContinuity"<<endl;
//	for(size_t i=0; i<C.size(); ++i) {
//		for(size_t j=0; j<C.size(); ++j) {
//			oflg<<A(i, j)<<"\t";
//		}
//		oflg<<"|\t"<<B(i)<<"\t|\t"<<X(i);
//		oflg<<endl;
//	}

	SparseMatrix<double, Eigen::ColMajor> A(C.size(), C.size());
	A.setFromTriplets(TripletList.begin(), TripletList.end());
	SparseLU<SparseMatrix<double>> solver;
	solver.analyzePattern(A);
	solver.factorize(A);
	X= solver.solve(B);

//	SparseMatrix<double, Eigen::RowMajor> A(C.size(), C.size());
//	A.setFromTriplets(TripletList.begin(), TripletList.end());
//	BiCGSTAB<SparseMatrix<double> > solver;
//	solver.compute(A);
//	X = solver.solve(B);




//	for(size_t i=0; i<C.size(); ++i) {
//		oflg<<"C # "<<i<<" : Pressure Correction="<<X[i];
//		oflg<<endl;
//	}

	if(vc==7) {
		for(size_t i=0; i<U.size(); ++i) {
			if(U[i].FB.T==NoBC) {
				U[i].p[uvlo]=U[i].p[uvli]+U[i].r*(X(U[i].FI[w])-X(U[i].FI[e]));
			} else if(U[i].FB.T==Inlet || U[i].FB.T==VelocityInlet) {
				U[i].p[uvlo]=U[i].p[uvli];
			} else {
				U[i].p[uvlo]=U[i].p[uvli];
				//For Outlet/PressureOutlet/PeriodicDownstream, the correction should happen in calculatePcoefficient().
			}
		}
		for(size_t i=0; i<V.size(); ++i) {
			if(V[i].FB.T==NoBC) {
				V[i].p[uvlo]=V[i].p[uvli]+V[i].r*(X(V[i].FI[s])-X(V[i].FI[n]));
			} else if(V[i].FB.T==Inlet || V[i].FB.T==VelocityInlet) {
				V[i].p[uvlo]=V[i].p[uvli];
			} else {
				V[i].p[uvlo]=V[i].p[uvli];
				//For Outlet/PressureOutlet/PeriodicDownstream, the correction should happen in calculatePcoefficient().
			}
		}
		for(size_t i=0; i<C.size(); ++i) {
			if(C[i].FB.T==Inlet || C[i].FB.T==PressureOutlet || C[i].FB.T==PeriodicUpstream || C[i].FB.T==PeriodicDownstream) {
				C[i].p[plo]=C[i].p[pli];
			} else {
				C[i].p[plo]=C[i].p[pli]+X(i)*cd_c::phy.PRelaxFactor;
			}
		}
	} else if(vc==4) {
		for(size_t i=0; i<C.size(); ++i) {
			if(C[i].FB.T==Inlet || C[i].FB.T==PressureOutlet || C[i].FB.T==PeriodicUpstream || C[i].FB.T==PeriodicDownstream) {
				C[i].p[plo]=C[i].p[pli];
			} else {
				C[i].p[plo]=X(i);
			}
		}
	} else if(vc==3) {
		for(size_t i=0; i<U.size(); ++i) {
			if(U[i].FB.T==NoBC) {
				U[i].p[uvlo]=U[i].p[uvli]+U[i].r*(X(U[i].FI[w])-X(U[i].FI[e]));
			} else if(U[i].FB.T==Inlet || U[i].FB.T==VelocityInlet) {
				U[i].p[uvlo]=U[i].p[uvli];
			} else {
				U[i].p[uvlo]=U[i].p[uvli];
				//For Outlet/PressureOutlet/PeriodicDownstream, the correction should happen in calculatePcoefficient().
			}
		}
		for(size_t i=0; i<V.size(); ++i) {
			if(V[i].FB.T==NoBC) {
				V[i].p[uvlo]=V[i].p[uvli]+V[i].r*(X(V[i].FI[s])-X(V[i].FI[n]));
			} else if(V[i].FB.T==Inlet || V[i].FB.T==VelocityInlet) {
				V[i].p[uvlo]=V[i].p[uvli];
			} else {
				V[i].p[uvlo]=V[i].p[uvli];
				//For Outlet/PressureOutlet/PeriodicDownstream, the correction should happen in calculatePcoefficient().
			}
		}
		for(size_t i=0; i<C.size(); ++i) {
			C[i].p[plo]=C[i].p[pli];
		}
	} else {
		//no such case
	}

	return true;
}

/*----------------------------------------------------------------------------------------------------

                                       function updateMassRatio

----------------------------------------------------------------------------------------------------*/
void updateMassRatio(const size_t & uvli) {
	uvCell_c::MassIn=0; uvCell_c::MassOut=0;
	for(size_t i=0; i<uvCell_c::InletUI.size(); ++i) {
		uvCell_c::MassIn +=U[uvCell_c::InletUI[i]].p[uvli]*U[uvCell_c::InletUI[i]].NormalA;
	}
	for(size_t i=0; i<uvCell_c::InletVI.size(); ++i) {
		uvCell_c::MassIn +=V[uvCell_c::InletVI[i]].p[uvli]*V[uvCell_c::InletVI[i]].NormalA;
	}

	for(size_t i=0; i<uvCell_c::OutletUI.size(); ++i) {
		uvCell_c::MassOut+=U[uvCell_c::OutletUI[i]].p[uvli]*U[uvCell_c::OutletUI[i]].NormalA;
	}
	for(size_t i=0; i<uvCell_c::OutletVI.size(); ++i) {
		uvCell_c::MassOut+=V[uvCell_c::OutletVI[i]].p[uvli]*V[uvCell_c::OutletVI[i]].NormalA;
	}

//	uvCell_c::MassRatio=uvCell_c::MassOut!=0?fabs(uvCell_c::MassIn/uvCell_c::MassOut):1;
	uvCell_c::MassRatio=uvCell_c::MassOut!=0?    (uvCell_c::MassIn/uvCell_c::MassOut):1;
}

/*----------------------------------------------------------------------------------------------------

                                       function updateBCForNewMR

----------------------------------------------------------------------------------------------------*/
void updateBCForNewMR() {
	for(size_t i=0; i<uvCell_c::OutletUI.size(); ++i) {
		U[uvCell_c::OutletUI[i]].FB.a[0]=-uvCell_c::MassRatio;
	}

	for(size_t i=0; i<uvCell_c::OutletVI.size(); ++i) {
		V[uvCell_c::OutletVI[i]].FB.a[0]=-uvCell_c::MassRatio;
	}
}

/*----------------------------------------------------------------------------------------------------

                                       function updatePForNewMR

----------------------------------------------------------------------------------------------------*/
void updateBUVForNewMR(const size_t & uvlio) {
	for(size_t i=0; i<uvCell_c::OutletUI.size(); ++i) {
		U[uvCell_c::OutletUI[i]].p[uvlio]=(U[uvCell_c::OutletUI[i]].FB.a[0]*U[U[uvCell_c::OutletUI[i]].FB.I[0]].p[uvlio])
				                         /U[uvCell_c::OutletUI[i]].FB.aP;
	}

	for(size_t i=0; i<uvCell_c::OutletVI.size(); ++i) {
		V[uvCell_c::OutletVI[i]].p[uvlio]=(V[uvCell_c::OutletVI[i]].FB.a[0]*V[V[uvCell_c::OutletVI[i]].FB.I[0]].p[uvlio])
		                                 /V[uvCell_c::OutletVI[i]].FB.aP;
	}

	for(size_t i=0; i<uvCell_c::InletUI.size(); ++i) {
		if(U[uvCell_c::InletUI[i]].FB.T==PeriodicUpstream) {
			U[uvCell_c::InletUI[i]].p[uvlio]=(U[uvCell_c::InletUI[i]].FB.a[0]*U[U[uvCell_c::InletUI[i]].FB.I[0]].p[uvlio])
		                                    /U[uvCell_c::InletUI[i]].FB.aP;
		}
	}

	for(size_t i=0; i<uvCell_c::InletVI.size(); ++i) {
		if(V[uvCell_c::InletVI[i]].FB.T==PeriodicUpstream) {
			V[uvCell_c::InletVI[i]].p[uvlio]=(V[uvCell_c::InletVI[i]].FB.a[0]*V[V[uvCell_c::InletVI[i]].FB.I[0]].p[uvlio])
		                                    /U[uvCell_c::InletVI[i]].FB.aP;
		}
	}
}

/*----------------------------------------------------------------------------------------------------

                                       function calculateRelaxFactor

----------------------------------------------------------------------------------------------------*/
void calcFlowfieldRelaxFactor(const numeric_t & avgUVError, const numeric_t & avgUVPError, const numeric_t & avgError) {
	if(       avgError>=0.5                    ) {
		cd_c::phy.UVRelaxFactor=0.1; cd_c::phy.PRelaxFactor=0.2;
	} else if(avgError>=0.1  && avgError<1   ) {
		cd_c::phy.UVRelaxFactor=0.2; cd_c::phy.PRelaxFactor=0.4;
	} else if(avgError>=0.05 && avgError<0.1 ) {
		cd_c::phy.UVRelaxFactor=0.3; cd_c::phy.PRelaxFactor=0.6;
	} else if(avgError>=0.01 && avgError<0.05) {
		cd_c::phy.UVRelaxFactor=0.4; cd_c::phy.PRelaxFactor=0.8;
	} else if(avgError< 0.01                 ) {
		cd_c::phy.UVRelaxFactor=0.5; cd_c::phy.PRelaxFactor=1  ;
	}
}

/*----------------------------------------------------------------------------------------------------

                                       function initializeFlowfield

----------------------------------------------------------------------------------------------------*/
bool initializeFlowfield() {
	cout<<"Initialize the flow field......"<<endl;

	for(size_t i=0; i<U.size(); ++i) {
		for(size_t j=0; j<HowManyCopy; ++j) {
			U[i].p[j]=cd_c::phy.FeatureVelocity/cd_c::phy.RefVelocity;
		}
		if(U[i].FB.T==Inlet || U[i].FB.T==VelocityInlet || U[i].FB.T==Wall || U[i].FB.T==PorousOpen) {
			if((U[i].FB.T==Inlet || U[i].FB.T==VelocityInlet) && cd_c::phy.InletFlowDeveloped) {
				U[i].FB.b=1.5*(cd_c::phy.FeatureVelocity/cd_c::phy.RefVelocity)
						 *(1-4*pow(((U[i].y-0.5*(cd_c::phy.FeatureLength/cd_c::phy.RefLength))/(cd_c::phy.FeatureLength/cd_c::phy.RefLength)), 2));
			}
			U[i].p[0]=U[i].FB.b/U[i].FB.aP;
		}
	}
	for(size_t i=0; i<V.size(); ++i) {
		for(size_t j=0; j<HowManyCopy; ++j) {
			V[i].p[j]=0;
		}
		if(V[i].FB.T==Inlet || V[i].FB.T==VelocityInlet || V[i].FB.T==Wall || V[i].FB.T==PorousOpen) {
			V[i].p[0]=V[i].FB.b/V[i].FB.aP;
		}
	}
	for(size_t i=0; i<C.size(); ++i) {
		for(size_t j=3; j<2*HowManyCopy; ++j) {
			C[i].p[j]=0;
		}
		if(C[i].FB.T==PeriodicUpstream || C[i].FB.T==PeriodicDownstream || C[i].FB.T==Inlet || C[i].FB.T==PressureOutlet) {
			C[i].p[3]=C[i].FB.b/C[i].FB.aP;
		}
	}

	return true;
}

/*----------------------------------------------------------------------------------------------------

                                       function readFlowfield

----------------------------------------------------------------------------------------------------*/
bool readFlowfield() {

	return true;
}
/*----------------------------------------------------------------------------------------------------

                                       function readFlowfield

----------------------------------------------------------------------------------------------------*/
bool setFlowfield() {
	cCell_c::dpdd=12*(cd_c::phy.FeatureVelocity/cd_c::phy.RefVelocity)*(cd_c::phy.Dry.dynamicViscosity/cd_c::phy.RefViscosity)
			     /pow((cd_c::phy.FeatureLength/cd_c::phy.RefLength), 2)
			     /cd_c::phy.Re;

	for(size_t i=0; i<U.size(); ++i) {
		for(size_t j=0; j<HowManyCopy; ++j) {
			U[i].p[j]=1.5*(cd_c::phy.FeatureVelocity/cd_c::phy.RefVelocity)
				     *(1-4*pow(((U[i].y-0.5*(cd_c::phy.FeatureLength/cd_c::phy.RefLength))/(cd_c::phy.FeatureLength/cd_c::phy.RefLength)), 2));
//			U[i].p[j]=cd_c::phy.FeatureVelocity/cd_c::phy.RefVelocity;
		}
	}
	for(size_t i=0; i<V.size(); ++i) {
		for(size_t j=0; j<HowManyCopy; ++j) {
			V[i].p[j]=0;
		}
	}
	for(size_t i=0; i<C.size(); ++i) {
		for(size_t j=3; j<2*HowManyCopy; ++j) {
			C[i].p[j]=0-cCell_c::dpdd*C[i].x;
		}
	}
	return true;
}

/*----------------------------------------------------------------------------------------------------

                                       function outputFlowfield

----------------------------------------------------------------------------------------------------*/
void outputFlowfield(const size_t & n) {
	stringstream iss("");
	iss<<n;
	string siss(iss.str());
	string   osf(caseName+".F"+siss+".vtk");
	ofstream off(osf.c_str())              ;
	plotFlowfield(off)                     ;
	off.close()                            ;
}
/*----------------------------------------------------------------------------------------------------

                                       function plotFlowfield

----------------------------------------------------------------------------------------------------*/
void plotFlowfield(ofstream & ofs) {//Nondimensionalized by Feature Properties
	ofs<<"# vtk DataFile Version 2.0"<<endl;
	ofs<<"Flowfield"<<endl;
	ofs<<"ASCII"<<endl;
	ofs<<"DATASET UNSTRUCTURED_GRID"<<endl;

	ofs<<"POINTS\t"<<ExtPtSize<<"\tdouble"<<endl;
	for(size_t i=0; i<ExtPtSize; ++i) {
		ofs<<Pt[i].x*cd_c::phy.RefLength/cd_c::phy.FeatureLength<<"\t"
		   <<Pt[i].y*cd_c::phy.RefLength/cd_c::phy.FeatureLength<<"\t"<<" 0"<<endl;
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
	ofs<<"SCALARS\t"<<"Pressure\t"<<"double\t"<<"1"<<endl;
	ofs<<"LOOKUP_TABLE\t"<<"default"<<endl;
	for(size_t i=0; i<C.size(); ++i) {
		ofs<<C[i].p[3]*cd_c::phy.RefPressure/cd_c::phy.FeaturePressure<<endl;
	}
	ofs<<"VECTORS\t"<<"Velocity\t"<<"double"<<endl;
	for(size_t i=0; i<C.size(); ++i) {
		ofs<<0.5*(U[C[i].FI[e]].p[0]+U[C[i].FI[w]].p[0])*cd_c::phy.RefVelocity/cd_c::phy.FeatureVelocity<<"\t"
		   <<0.5*(V[C[i].FI[n]].p[0]+V[C[i].FI[s]].p[0])*cd_c::phy.RefVelocity/cd_c::phy.FeatureVelocity<<"\t"<<" 0"<<endl;
	}

	cout<<"Flowfield Ploted!"<<endl;
}
