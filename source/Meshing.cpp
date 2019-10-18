/*
 * Meshing.cpp
 *
 *  Created on: Apr 15, 2014
 *      Author: david
 */

#include "GlobalDeclaration.h"

using namespace std;
using namespace Eigen;

size_t porousBlocks(0), flowBlocks(0);

inline numeric_t Cxof(const point_c & p0=point_c(), const point_c & p1=point_c(), const point_c & p2=point_c(), const point_c & p3=point_c()) {
	return 0.25*(p0.x+p1.x+p2.x+p3.x);
}

inline numeric_t Cyof(const point_c & p0=point_c(), const point_c & p1=point_c(), const point_c & p2=point_c(), const point_c & p3=point_c()) {
	return 0.25*(p0.y+p1.y+p2.y+p3.y);
}

numeric_t RandomPercentage(const numeric_t & rseed);
bool BlocksBetweenIsNotPorous(const index_t & i, const index_t & j);
bool ThroatExist(const throat_c & t);

direction_e PoreRelation(const index_t & i, const index_t & j);
bool findPoreLocation(index_t & b, index_t & i, index_t & j, const index_t & pi);

point_c crossPoint(const point_c & p0, const point_c & p1, const numeric_t & os1, const point_c & p2, const numeric_t & os2);
point_c noCrossPoint(const point_c & p0, const point_c & p1, const numeric_t & os);

bool plotMesh(ofstream & ofs);

/*----------------------------------------------------------------------------------------------------

                                        function meshing

----------------------------------------------------------------------------------------------------*/
bool meshing() {
	cout<<"Start meshing......"<<endl;

	string meshName=caseName+".cd";
	ofstream meshOFS(meshName.c_str());

	cout<<"Sorting blocks by its type, PorousMedia Blocks goes first, and then follows FlowField Blocks......"<<endl;
	vector<block_c> blk(0);
	for(size_t i=0; i<Bk.size(); ++i) {
		if(Bk[i].T==porousmedia) blk.push_back(Bk[i]);
	}
	porousBlocks=blk.size();
	for(size_t i=0; i<Bk.size(); ++i) {
		if(Bk[i].T==flowfield) blk.push_back(Bk[i]);
	}
	flowBlocks=blk.size()-porousBlocks;
	Bk.clear();
	Bk=blk;
	blk.clear();
	cout<<"Sorting blocks finished, there are all together "<<porousBlocks<<"PorousMedia Blocks and "
			                                                <<flowBlocks  <<"FlowField Blocks!"<<endl;
	for(size_t k=0; k<Bk.size(); ++k) {
		oflg<<"Block # "<<k<<":"<<endl;
		oflg<<Bk[k];
	}

	cout<<"Looking for blocks' neighbors......"<<endl;
	for(size_t i=0; i<Bk.size(); ++i) {
		for(size_t j=i+1; j<Bk.size(); ++j) {
			if(Bk[i].CnI[NE]==Bk[j].CnI[NW] && Bk[i].CnI[SE]==Bk[j].CnI[SW]) {
				Bk[i].NI[e ]=j;
				Bk[j].NI[w ]=i;
			} else if(Bk[i].CnI[NE]==Bk[j].CnI[SW]) {
				Bk[i].NI[ne]=j;
				Bk[j].NI[sw]=i;
			} else if(Bk[i].CnI[NE]==Bk[j].CnI[SE] && Bk[i].CnI[NW]==Bk[j].CnI[SW]) {
				Bk[i].NI[n ]=j;
				Bk[j].NI[s ]=i;
			} else if(Bk[i].CnI[NW]==Bk[j].CnI[SE]) {
				Bk[i].NI[nw]=j;
				Bk[j].NI[se]=i;
			} else if(Bk[i].CnI[NW]==Bk[j].CnI[NE] && Bk[i].CnI[SW]==Bk[j].CnI[SE]) {
				Bk[i].NI[w ]=j;
				Bk[j].NI[e ]=i;
			} else if(Bk[i].CnI[SW]==Bk[j].CnI[NE]) {
				Bk[i].NI[sw]=j;
				Bk[j].NI[ne]=i;
			} else if(Bk[i].CnI[SE]==Bk[j].CnI[NE] && Bk[i].CnI[SW]==Bk[j].CnI[NW]) {
				Bk[i].NI[s ]=j;
				Bk[j].NI[n ]=i;
			} else if(Bk[i].CnI[SE]==Bk[j].CnI[NW]) {
				Bk[i].NI[se]=j;
				Bk[j].NI[nw]=i;
			} else {
				//Do nothing! Index remain as voidIndex.
			}
		}
	}

	cout<<"Calculating dz......"<<endl;
	if(porousBlocks!=0) {
		numeric_t WholeAvgThroatLength(0);
		for(size_t i=0; i<porousBlocks; ++i) {
			WholeAvgThroatLength+=Bk[i].avgTLength;
		}
		cd_c::dz=WholeAvgThroatLength/porousBlocks;
//		cd_c::dz=WholeAvgThroatLength/porousBlocks*2;

//		numeric_t WholeAvgThroatDiameter(0);
//		size_t WholeDifferentThroat(0);
//		for(size_t k=0; k<porousBlocks; ++k) {
//			for(size_t i=0; i<Bk[k].avgTDiameter.size(); ++i) {
//				WholeAvgThroatDiameter+=Bk[i].avgTDiameter[i];
//			}
//			WholeDifferentThroat+=Bk[k].avgTDiameter.size();
//		}
//		cd_c::dz=WholeAvgThroatDiameter/WholeDifferentThroat;
	} else {
		cd_c::dz=1;
	}
	oflg<<"Cells Thickness dz="<<cd_c::dz<<endl;

	oflg<<"Blocks' relationships are as follows: "<<endl;
	for(size_t i=0; i!=Bk.size(); ++i) {
		oflg<<"Block # "<<i<<"\tType:"<<Bk[i].T<<"\tCorner Points Index: ";
		oflg<<Bk[i].CnI[e ]<<"\t"<<Bk[i].CnI[n ]<<"\t"<<Bk[i].CnI[w ]<<"\t"<<Bk[i].CnI[s ]<<"\tNeighbor Blocks: ";
		oflg<<Bk[i]. NI[e ]<<"\t"<<Bk[i]. NI[n ]<<"\t"<<Bk[i]. NI[w ]<<"\t"<<Bk[i]. NI[s ]<<"\t"
		    <<Bk[i]. NI[ne]<<"\t"<<Bk[i]. NI[nw]<<"\t"<<Bk[i]. NI[sw]<<"\t"<<Bk[i]. NI[se]<<endl;
	}

	cout<<"Creating lines between blocks......"<<endl;
	for(size_t i=0; i<Bk.size(); ++i) {
		for(size_t j=0; j<sideNumber; ++j) {
			if(Bk[i].T==porousmedia && Bk[i].B[j].T==Wall && Bk[i].NI[j]==voidIndex) {
				//under this condition, the line no need to exist
			} else {
				line_c newline; //create a newline
				switch (j) {
				case e:
					newline=line_c(Pt[Bk[i].CnI[SE]], Pt[Bk[i].CnI[NE]], Bk[i].CnI[SE], Bk[i].CnI[NE], yAxis,
							       segment_c(Bk[i].Ny, Bk[i].G[s], Bk[i].G[n]), Bk[i].B[e]);
					break;
				case n:
					newline=line_c(Pt[Bk[i].CnI[NW]], Pt[Bk[i].CnI[NE]], Bk[i].CnI[NW], Bk[i].CnI[NE], xAxis,
							       segment_c(Bk[i].Nx, Bk[i].G[w], Bk[i].G[e]), Bk[i].B[n]);
					break;
				case w:
					newline=line_c(Pt[Bk[i].CnI[SW]], Pt[Bk[i].CnI[NW]], Bk[i].CnI[SW], Bk[i].CnI[NW], yAxis,
							       segment_c(Bk[i].Ny, Bk[i].G[s], Bk[i].G[n]), Bk[i].B[w]);
					break;
				case s:
					newline=line_c(Pt[Bk[i].CnI[SW]], Pt[Bk[i].CnI[SE]], Bk[i].CnI[SW], Bk[i].CnI[SE], xAxis,
							       segment_c(Bk[i].Nx, Bk[i].G[w], Bk[i].G[e]), Bk[i].B[s]);
					break;
				}

				bool foundSameLine(false); //check if the newline exist or not
				index_t SameLinePosition(voidIndex);
				for(size_t k=0; k<Ln.size(); ++k) {
					if(newline==Ln[k]) {
						foundSameLine=true;
						SameLinePosition=k;
						Ln[k]=newline;//update segment, in case the porousmedia blocks' Nx and Ny is not right
						break;
					}
				}

				if(foundSameLine) {//if exist already, then directly assign index
					Bk[i].LnI[j]=SameLinePosition;
				} else {//if not exist, then create it and assign index
					Bk[i].LnI[j]=Ln.size();
					Ln.push_back(newline);
				}
			}
		}
	}
	cout<<"Looking for parallel equal lines......"<<endl;
	for(size_t i=0; i<Ln.size(); ++i) {
		for(size_t j=i+1; j<Ln.size(); ++j) {
			if(Ln[i].ParaEqTo(Ln[j]) && BlocksBetweenIsNotPorous(i, j)) {
				Ln[i].PELI.push_back(j);
				Ln[j].PELI.push_back(i);
			}
		}
	}

	cout<<"Matching In & Out Boundary Conditions......"<<endl;
	for(size_t i=0; i<Ln.size(); ++i) {
		if(Ln[i].B[0].T==PeriodicUpstream) {
			for(size_t j=i+1; j<Ln.size(); ++j) {
				if(Ln[i].ParaEqTo(Ln[j]) && Ln[j].B[0].T==PeriodicDownstream) {
					Ln[i].B[0].DownstreamLI=j;
					Ln[j].B[0].  UpstreamLI=i;
				}
			}
		} else if(Ln[i].B[0].T==PeriodicDownstream) {
			for(size_t j=i+1; j<Ln.size(); ++j) {
				if(Ln[i].ParaEqTo(Ln[j]) && Ln[j].B[0].T==PeriodicUpstream) {
					Ln[i].B[0].  UpstreamLI=j;
					Ln[j].B[0].DownstreamLI=i;
				}
			}
//		} else if(Ln[i].B[0].T==Inlet) {
//			for(size_t j=i+1; j<Ln.size(); ++j) {
//				if(Ln[i].ParaEqTo(Ln[j]) && Ln[j].B[0].T==Outlet) {
//					Ln[i].B[0].DownstreamLI=j;
//					Ln[j].B[0].  UpstreamLI=i;
//				}
//			}
//		} else if(Ln[i].B[0].T==Outlet) {
//			for(size_t j=i+1; j<Ln.size(); ++j) {
//				if(Ln[i].ParaEqTo(Ln[j]) && Ln[j].B[0].T==Inlet) {
//					Ln[i].B[0].  UpstreamLI=j;
//					Ln[j].B[0].DownstreamLI=i;
//				}
//			}
//		} else if(Ln[i].B[0].T==VelocityInlet) {
//			for(size_t j=i+1; j<Ln.size(); ++j) {
//				if(Ln[i].ParaEqTo(Ln[j]) && Ln[j].B[0].T==PressureOutlet) {
//					Ln[i].B[0].DownstreamLI=j;
//					Ln[j].B[0].  UpstreamLI=i;
//				}
//			}
//		} else if(Ln[i].B[0].T==PressureOutlet) {
//			for(size_t j=i+1; j<Ln.size(); ++j) {
//				if(Ln[i].ParaEqTo(Ln[j]) && Ln[j].B[0].T==VelocityInlet) {
//					Ln[i].B[0].  UpstreamLI=j;
//					Ln[j].B[0].DownstreamLI=i;
//				}
//			}
		} else {
			//Do nothing
		}
	}
	for(size_t k=0; k<Bk.size(); ++k) {
		for(size_t d=e; d<=s; ++d) {
			if(Bk[k].B[d].T==PeriodicUpstream) {
				Bk[k].B[d].DownstreamLI=Ln[Bk[k].LnI[d]].B[0].DownstreamLI;
			} else if(Bk[k].B[d].T==PeriodicDownstream) {
				Bk[k].B[d].  UpstreamLI=Ln[Bk[k].LnI[d]].B[0].  UpstreamLI;
			} else {}//do nothing
		}
	}

	cout<<"Making sure lines between PeriodicUpstream and PeriodicDownstream are all ParaEq......"<<endl;
	for(size_t i=0; i<Ln.size(); ++i) {
		if(Ln[i].B[0].T==PeriodicUpstream) {
			bool AlreadyIn(false);
			for(size_t peli=0; peli<Ln[i].PELI.size(); ++peli) {
				if(Ln[i].B[0].DownstreamLI==Ln[i].PELI[peli]) AlreadyIn=true;
			}
			if(!AlreadyIn) Ln[i].PELI.push_back(Ln[i].B[0].DownstreamLI);

			for(size_t j=0; j<Ln[Ln[i].B[0].DownstreamLI].PELI.size(); ++j) {
				bool peliInpeli(false);
				for(size_t peli=0; peli<Ln[i].PELI.size(); ++peli) {
					if(Ln[Ln[i].B[0].DownstreamLI].PELI[j]==Ln[i].PELI[peli]) peliInpeli=true;
				}
				if(!peliInpeli && Ln[Ln[i].B[0].DownstreamLI].PELI[j]!=i) Ln[i].PELI.push_back(Ln[Ln[i].B[0].DownstreamLI].PELI[j]);
			}
		} else if(Ln[i].B[0].T==PeriodicDownstream) {
			bool AlreadyIn(false);
			for(size_t peli=0; peli<Ln[i].PELI.size(); ++peli) {
				if(Ln[i].B[0].UpstreamLI==Ln[i].PELI[peli]) AlreadyIn=true;
			}
			if(!AlreadyIn) Ln[i].PELI.push_back(Ln[i].B[0].UpstreamLI);

			for(size_t j=0; j<Ln[Ln[i].B[0].UpstreamLI].PELI.size(); ++j) {
				bool peliInpeli(false);
				for(size_t peli=0; peli<Ln[i].PELI.size(); ++peli) {
					if(Ln[Ln[i].B[0].UpstreamLI].PELI[j]==Ln[i].PELI[peli]) peliInpeli=true;
				}
				if(!peliInpeli && Ln[Ln[i].B[0].UpstreamLI].PELI[j]!=i) Ln[i].PELI.push_back(Ln[Ln[i].B[0].UpstreamLI].PELI[j]);
			}
		}
	}

	for(size_t i=0; i<Ln.size(); ++i) {
		if(Ln[i].B[0].T==PeriodicUpstream) {
			cCell_c::dpdd=(Ln[i].B[0].p-Ln[Ln[i].B[0].DownstreamLI].B[0].p)/Ln[i].distance(Ln[Ln[i].B[0].DownstreamLI]);
			oflg<<"Between Line # "<<i<<" and Line # "<<Ln[i].B[0].DownstreamLI<<" : dpdd= "<<cCell_c::dpdd<<endl;
		} else if(Ln[i].B[0].T==PeriodicDownstream) {
			cCell_c::dpdd=(Ln[Ln[i].B[0].UpstreamLI].B[0].p-Ln[i].B[0].p)/Ln[i].distance(Ln[Ln[i].B[0].UpstreamLI]);
			oflg<<"Between Line # "<<i<<" and Line # "<<Ln[i].B[0].  UpstreamLI<<" : dpdd= "<<cCell_c::dpdd<<endl;
		}
	}

	cout<<"Lines' initial parameters & relationships setting finished!"<<endl;
	for(size_t k=0; k<Bk.size(); ++k) {
		oflg<<"Surrounding Block# "<<k<<", there are following Lines: ";
		for(size_t d=e; d<=s; ++d) {
			oflg<<"\t"<<Bk[k].LnI[d];
		}
		oflg<<endl;
	}
	for(size_t i=0; i<Ln.size(); ++i) {
		oflg<<"Line # "<<i<<"\tBegin Point: "<<*Ln[i].PtI.begin()<<"\tEnd Point: "<<*(Ln[i].PtI.end()-1)<<"\tParallel Equal Lines: ";
		for(size_t j=0; j<Ln[i].PELI.size(); ++j) {
			oflg<<Ln[i].PELI[j]<<"\t";
		}
		oflg<<"\tBC: "<<Ln[i].B[0].T;
		if(Ln[i].B[0].T==PeriodicUpstream) oflg<<"DownstreamLI: "<<Ln[i].B[0].DownstreamLI<<endl;
		else if(Ln[i].B[0].T==PeriodicDownstream) oflg<<"UpstreamLI: "<<Ln[i].B[0].UpstreamLI<<endl;
		oflg<<endl;
	}

	cout<<"Starting to generate points according to initial parameters......"<<endl;
	cout<<"Meshing PorousMedia blocks and Set breakpoints for the side lines......"<<endl;
	for(size_t k=0; k<porousBlocks; ++k) {
		cout<<"Adding pores for Block #"<<k<<"......"<<endl;

		Bk[k].Lx=Pt[Bk[k].CnI[0]].x-Pt[Bk[k].CnI[2]].x;
		Bk[k].Ly=Pt[Bk[k].CnI[0]].y-Pt[Bk[k].CnI[2]].y;
		Bk[k].Ox=Pt[Bk[k].CnI[2]].x;
		Bk[k].Oy=Pt[Bk[k].CnI[2]].y;
		numeric_t orgx(Pt[Bk[k].CnI[2]].x+0.5*(Bk[k].Lx-Bk[k].Nx*Bk[k].avgTLength));
		numeric_t orgy(Pt[Bk[k].CnI[2]].y+0.5*(Bk[k].Ly-Bk[k].Ny*Bk[k].avgTLength)); // original point for basic regular network
		if(Bk[k].B[e].T==PorousOpen) orgx=Pt[Bk[k].CnI[0]].x-Bk[k].Nx*Bk[k].avgTLength;
		if(Bk[k].B[w].T==PorousOpen) orgx=Pt[Bk[k].CnI[2]].x;
		if(Bk[k].B[n].T==PorousOpen) orgy=Pt[Bk[k].CnI[0]].y-Bk[k].Ny*Bk[k].avgTLength;
		if(Bk[k].B[s].T==PorousOpen) orgy=Pt[Bk[k].CnI[2]].y;

		for(size_t i=0; i<=Bk[k].Nx; ++i) {
			vector<index_t> vertical;
			for(size_t j=0; j<=Bk[k].Ny; ++j) {
//				cout<<"    "<<i<<"     "<<j<<endl;
				if(   (i==Bk[k].Nx && j==Bk[k].Ny && Bk[k].B[e].T==PorousOpen && Bk[k].B[n].T==PorousOpen)
				   || (i==0        && j==Bk[k].Ny && Bk[k].B[n].T==PorousOpen && Bk[k].B[w].T==PorousOpen)
				   || (i==0        && j==0        && Bk[k].B[w].T==PorousOpen && Bk[k].B[s].T==PorousOpen)
				   || (i==Bk[k].Nx && j==0        && Bk[k].B[s].T==PorousOpen && Bk[k].B[e].T==PorousOpen)   ) {
					vertical.push_back(voidIndex);
				} else {
					numeric_t x(orgx+i*Bk[k].avgTLength), y(orgy+j*Bk[k].avgTLength);
					if(   (i==Bk[k].Nx && Bk[k].B[e].T==PorousOpen)
					   || (i==0        && Bk[k].B[w].T==PorousOpen)
					   || (j==Bk[k].Ny && Bk[k].B[n].T==PorousOpen)
					   || (j==0        && Bk[k].B[s].T==PorousOpen)   ) {
							vertical.push_back(P.size());
							P.push_back(pore_c(x, y));//**********Pores created here**********
							P[P.size()-1].setPosition(open);
					} else {
						unsigned long lengthSeed((P.size()+1)*Bk[k].RandomRealization) //**********Here make different realization different**********
						             , angleSeed((P.size()+1)*Bk[k].RandomRealization);//**********Here make different realization different**********
						numeric_t r(Bk[k].avgTLength*0.5*Bk[k].TLVarianceRate*RandomPercentage(lengthSeed))
						        , theta(2*pi*RandomPercentage(angleSeed));
						x+=r*cos(theta);
						y+=r*sin(theta);

						vertical.push_back(P.size());
						P.push_back(pore_c(x, y));//**********Pores created here**********
					}

					if((i>Bk[k].Nx/2 && Bk[k].porosityType==leftrightPorosity) ||
					   (j>Bk[k].Ny/2 && Bk[k].porosityType==downupPorosity)) {
						P[P.size()-1].RI=Bk[k].avgTDiameter[1]/2;
					} else {
						P[P.size()-1].RI=Bk[k].avgTDiameter[0]/2;
					}
				}
			}
			Bk[k].PI.push_back(vertical);
			vertical.clear();
		}

		oflg<<"Pores created for PorousMedia Block #"<<k<<", they are aligned as follows:"<<endl;
		for(index_t j=Bk[k].Ny; j>=0; --j) {
			for(size_t i=0; i<=Bk[k].Nx; ++i) {
				if(Bk[k].PI[i][j]!=voidIndex)
				oflg<<"P"<<Bk[k].PI[i][j]<<"("<<setprecision(4)<<P[Bk[k].PI[i][j]].x<<", "<<setprecision(4)<<P[Bk[k].PI[i][j]].y<<")\t";
			}
			oflg<<endl;
		}

		cout<<"Connecting selected Pores to form Throats......"<<endl;
		if(Bk[k].porosityType==assignPorosity) {
			cout<<"Creating Throats by directly assigning parameters for each Throat......"<<endl;

			cout<<"Reading Throats Assignment......"<<endl;
			string assignThroatsName=caseName+".at";
			ifstream atIFS(assignThroatsName.c_str());
			
			numeric_t totTLength(0), totWeighted(0); //total throat length, and total weighted throat length

			for(size_t indexj=Bk[k].Ny*2; indexj>=0; --indexj) {
				for(size_t indexi=0; indexi<=Bk[k].Nx*2; ++indexi) {
					numeric_t theDiameter(0); //Throat Diameter: Positive: Normal or J+I+; Negative: J+I-
					size_t thePolyN(0); //0: Circle; 3: Traiangle; 4: Square; 5: Pentagon; 6: Hexagon
					index_t i(0), I(0), j(0), J(0);

					atIFS>>theDiameter>>thePolyN;

					if((indexi%2==0 && indexj%2==0) || theDiameter==0) {
						//Do nothing, no throats
					} else {
						if(indexi%2==0) {
							i=indexi/2    ;
							I=i           ;
							j=(indexj-1)/2;
							J=(indexj+1)/2;
						} else if(indexj%2==0) {
							i=(indexi-1)/2;
							I=(indexi+1)/2;
							j=indexj/2    ;
							J=j           ;
						} else {
							if(theDiameter>0) {
								i=(indexi-1)/2;
								I=(indexi+1)/2;
								j=(indexj-1)/2;
								J=(indexj+1)/2;
							} else {
								i=(indexi+1)/2;
								I=(indexi-1)/2;
								j=(indexj-1)/2;
								J=(indexj+1)/2;
							}
						}

						Bk[k].TI.push_back(T.size());
						P[Bk[k].PI[i][j]].TI.push_back(T.size());
						P[Bk[k].PI[I][J]].TI.push_back(T.size());
						T.push_back(throat_c(Bk[k].PI[i][j], Bk[k].PI[I][J], thePolyN, abs(theDiameter)/2, P[Bk[k].PI[i][j]].distance(P[Bk[k].PI[I][J]])));

						totTLength +=T[T.size()-1].L;
						totWeighted+=T[T.size()-1].L*theDiameter*theDiameter;
					}
				}
			}

			Bk[k].avgTDiameter[0]=sqrt(totWeighted/totTLength);
		} else {
			cout<<"Creating Throats by connecting Pores using designated statistical parameters......"<<endl;

			for(size_t i=0; i<=Bk[k].Nx; ++i) {
				for(size_t j=0; j<=Bk[k].Ny; ++j) {
					if(Bk[k].PI[i][j]==voidIndex) {
						//do nothing
					} else {
						vector<bool> Connected(sideNumber, false);

						//if the pore is on the edge, then there is no need to connect out, so fake connected
						if(i==0       ) Connected[w]=true;
						if(i==Bk[k].Nx) Connected[e]=true;
						if(j==0       ) Connected[s]=true;
						if(j==Bk[k].Ny) Connected[n]=true;

						//if the pore is on the porousOpen border, then it is only allowed in one direction
						vector<bool> ConnectAllowed(nodeNumber, true);
						if(i==Bk[k].Nx && Bk[k].B[e].T==PorousOpen) {
							ConnectAllowed[n ]=false;
							ConnectAllowed[s ]=false;
							ConnectAllowed[nw]=false;
							ConnectAllowed[sw]=false;
							Connected[n]=true;
							Connected[s]=true;
						}
						if(i==0        && Bk[k].B[w].T==PorousOpen) {
							ConnectAllowed[n ]=false;
							ConnectAllowed[s ]=false;
							ConnectAllowed[ne]=false;
							ConnectAllowed[se]=false;
							Connected[n]=true;
							Connected[s]=true;
						}
						if(j==Bk[k].Ny && Bk[k].B[n].T==PorousOpen) {
							ConnectAllowed[e ]=false;
							ConnectAllowed[w ]=false;
							ConnectAllowed[se]=false;
							ConnectAllowed[sw]=false;
							Connected[e]=true;
							Connected[w]=true;
						}
						if(j==0        && Bk[k].B[s].T==PorousOpen) {
							ConnectAllowed[e ]=false;
							ConnectAllowed[w ]=false;
							ConnectAllowed[ne]=false;
							ConnectAllowed[nw]=false;
							Connected[e]=true;
							Connected[w]=true;
						}
						if((i==(Bk[k].Nx-1)) && Bk[k].B[e].T==PorousOpen) {
							ConnectAllowed[ne]=false;
							ConnectAllowed[se]=false;
						}
						if(i==1              && Bk[k].B[w].T==PorousOpen) {
							ConnectAllowed[nw]=false;
							ConnectAllowed[sw]=false;
						}
						if((j==(Bk[k].Ny-1)) && Bk[k].B[n].T==PorousOpen) {
							ConnectAllowed[ne]=false;
							ConnectAllowed[nw]=false;
						}
						if(j==1              && Bk[k].B[s].T==PorousOpen) {
							ConnectAllowed[se]=false;
							ConnectAllowed[sw]=false;
						}

						for(index_t d=se; d>=e; --d) {
		//					cout<<"Throat:"<<d<<"\t"<<T.size()<<endl;
							bool NoCross(true);
							index_t I(i), J(j);
							switch (d) {
							case ne: I=i+1; J=j+1; break;
							case nw: I=i-1; J=j+1; break;
							case sw: I=i-1; J=j-1; break;
							case se: I=i+1; J=j-1; break;
							case  e: I=i+1; J=j  ; break;
							case  n: I=i  ; J=j+1; break;
							case  w: I=i-1; J=j  ; break;
							case  s: I=i  ; J=j-1; break;
							}

		                    if(I>=0 && I<=Bk[k].Nx && J>=0 && J<=Bk[k].Ny && ConnectAllowed[d] && !ThroatExist(throat_c(Bk[k].PI[i][j], Bk[k].PI[I][J])) ) {//be sure no duplicates
								if(ThroatExist(throat_c(Bk[k].PI[i][J], Bk[k].PI[I][j]))) NoCross=false;

								size_t theSeed=(T.size()+1)*Bk[k].RandomRealization;//**********Here make different realization different**********
								numeric_t theDiameter(0);
								size_t thePolyN(0);
								if((i>=Bk[k].Nx/2 && Bk[k].porosityType==leftrightPorosity) ||
								   (j>=Bk[k].Ny/2 && Bk[k].porosityType==downupPorosity)) {
									theDiameter=Bk[k].avgTDiameter[1]*(1+Bk[k].TDVarianceRate[1]*(RandomPercentage(theSeed)-0.5));
									thePolyN=Bk[k].PolyN[1];
								} else {
									theDiameter=Bk[k].avgTDiameter[0]*(1+Bk[k].TDVarianceRate[0]*(RandomPercentage(theSeed)-0.5));
									thePolyN=Bk[k].PolyN[0];
								}

								switch (d) {
								case ne: case nw: case sw: case se:
									if(RandomPercentage(theSeed)<=Bk[k].PPCrossRate && NoCross) {
										Bk[k].TI.push_back(T.size());
										P[Bk[k].PI[i][j]].TI.push_back(T.size());
										P[Bk[k].PI[I][J]].TI.push_back(T.size());
										T.push_back(throat_c(Bk[k].PI[i][j], Bk[k].PI[I][J], thePolyN, theDiameter/2, P[Bk[k].PI[i][j]].distance(P[Bk[k].PI[I][J]])));
										Connected[d-4]=true;
										if(d==se) {
											Connected[e]=true;
										} else {
											Connected[d-3]=true;
										}
									}
									break;
								case e: case n: case w: case s:
									if(Connected[d]==false || RandomPercentage(theSeed)<=Bk[k].PPConnectionRate) {
										Bk[k].TI.push_back(T.size());
										P[Bk[k].PI[i][j]].TI.push_back(T.size());
										P[Bk[k].PI[I][J]].TI.push_back(T.size());
										T.push_back(throat_c(Bk[k].PI[i][j], Bk[k].PI[I][J], thePolyN, theDiameter/2, P[Bk[k].PI[i][j]].distance(P[Bk[k].PI[I][J]])));
										Connected[d]=true;
									}
									break;
								}
							} else {
								//if that pore doesn't exist or connection not allowed, then cannot create any throat
							}//finish if
						}//finish d
					}
				}// finish j
			}//finish i
		}//finish if, if not assignPorosity


		oflg<<"Throats added in Block #"<<k<<":"<<endl;
		for(size_t ti=0; ti!=Bk[k].TI.size(); ++ti) {
			oflg<<"T"<<Bk[k].TI[ti]<<"\t("<<T[Bk[k].TI[ti]].PI[0]<<"\t-\t"<<T[Bk[k].TI[ti]].PI[1]<<")\tLength: "<<T[Bk[k].TI[ti]].L<<"\tRadius: "<<T[Bk[k].TI[ti]].RI<<endl;
		}
//		oflg<<"Pores and their neighbor throats:"<<endl;
//		for(index_t j=Bk[k].Ny; j>=0; --j) {
//			for(size_t i=0; i<=Bk[k].Nx; ++i) {
//				oflg<<"P"<<Bk[k].PI[i][j]<<"(";
//				for(size_t PTI=0; PTI!=P[Bk[k].PI[i][j]].TI.size(); ++PTI) {
//					oflg<<"T"<<P[Bk[k].PI[i][j]].TI[PTI];
//					if(PTI<P[Bk[k].PI[i][j]].TI.size()-1) oflg<<", ";
//				}
//				oflg<<")\t";
//			}
//			oflg<<endl;
//		}

		cout<<"Creating interface with neighbor blocks, that is to override the interface original lines......"<<endl;
		for(size_t d=e; d<=s; ++d) {
			if(Bk[k].B[d].T==PorousOpen) {
				Ln[Bk[k].LnI[d]].clear();//If this is PorousOpen line, then the original line setting is erased

				numeric_t sinTheta(1), openDiameter(0);
				size_t i(0), j(0), I(0), J(0);

				switch (d) {
				case e: case w:
					if(d==e) {
						i=Bk[k].Nx;
						I=i-1;
					} else {
						i=0;
						I=i+1;
					}
					for(j=0; j<=Bk[k].Ny; ++j) {
						if(Bk[k].PI[i][j]==voidIndex) {
							//do nothing
						} else {
							if(Bk[k].NinDiameter!=0 && Bk[k].NinLength==0) {
								sinTheta=P[Bk[k].PI[i][j]].distanceX(P[Bk[k].PI[I][j]])/P[Bk[k].PI[i][j]].distance(P[Bk[k].PI[I][j]]);
								openDiameter=sqrt(4*T[P[Bk[k].PI[i][j]].TI[0]].A/pi)/sinTheta;

								segment_c s1, s3;
								if(j==0        || (Bk[k].PI[i][0       ]==voidIndex && j==1)) {
									s1=segment_c(0, 0, 1, 0, openDiameter/Bk[k].NinDiameter);
								} else {
									s1=segment_c(Bk[k].avgTLength/(openDiameter/Bk[k].NinDiameter));
								}
								if(j==Bk[k].Ny || (Bk[k].PI[i][Bk[k].Ny]==voidIndex && j==Bk[k].Ny-1)) s3=segment_c(0, 1, 0, openDiameter/Bk[k].NinDiameter, 0);

								Ln[Bk[k].LnI[d]].insertSegment(point_c(P[Bk[k].PI[i][j]].x, P[Bk[k].PI[i][j]].y-openDiameter/2),
										                       point_c(P[Bk[k].PI[i][j]].x, P[Bk[k].PI[i][j]].y+openDiameter/2),
										                       s1,
										                       segment_c(Bk[k].NinDiameter),
										                       s3,
										                       BC_c(Wall, 0.0, 0.0),
										                       BC_c(PorousOpen, Bk[k].PI[i][j], T[P[Bk[k].PI[i][j]].TI[0]].L),
										                       BC_c(Wall, 0.0, 0.0));
							} else if(Bk[k].NinDiameter==0 && Bk[k].NinLength!=0) {
								openDiameter=Bk[k].avgTLength/Bk[k].NinLength;

								if(j==0 || (Bk[k].PI[i][0       ]==voidIndex && j==1)) {
									if(P[Bk[k].PI[i][j]].y-openDiameter/2<=Ln[Bk[k].LnI[d]].Pt[0].y) {
										Ln[Bk[k].LnI[d]].insertPoint(point_c(P[Bk[k].PI[i][j]].x, P[Bk[k].PI[i][j]].y+openDiameter/2),
																	 segment_c(),
																	 segment_c(),
																	 BC_c(PorousOpen, Bk[k].PI[i][j], T[P[Bk[k].PI[i][j]].TI[0]].L),
																	 BC_c(Wall, 0.0, 0.0));
									} else {
										Ln[Bk[k].LnI[d]].insertSegment(point_c(P[Bk[k].PI[i][j]].x, P[Bk[k].PI[i][j]].y-openDiameter/2),
														               point_c(P[Bk[k].PI[i][j]].x, P[Bk[k].PI[i][j]].y+openDiameter/2),
														               segment_c(0, 0, 1, 0, openDiameter),
														               segment_c(),
														               segment_c(),
														               BC_c(Wall, 0.0, 0.0),
														               BC_c(PorousOpen, Bk[k].PI[i][j], T[P[Bk[k].PI[i][j]].TI[0]].L),
														               BC_c(Wall, 0.0, 0.0));
									}
								} else if(j==Bk[k].Ny || (Bk[k].PI[i][Bk[k].Ny]==voidIndex && j==Bk[k].Ny-1)) {
									if(Bk[k].NinLength==1) {
										if(P[Bk[k].PI[i][j]].y+openDiameter/2>=Ln[Bk[k].LnI[d]].Pt[Ln[Bk[k].LnI[d]].Pt.size()-1].y) {
											Ln[Bk[k].LnI[d]].B[Ln[Bk[k].LnI[d]].B.size()-1]=BC_c(PorousOpen, Bk[k].PI[i][j], T[P[Bk[k].PI[i][j]].TI[0]].L);
										} else {
											Ln[Bk[k].LnI[d]].insertPoint(point_c(P[Bk[k].PI[i][j]].x, P[Bk[k].PI[i][j]].y+openDiameter/2),
																		 segment_c(),
																		 segment_c(0, 1, 0, openDiameter, 0),
																		 BC_c(PorousOpen, Bk[k].PI[i][j], T[P[Bk[k].PI[i][j]].TI[0]].L),
																		 BC_c(Wall, 0.0, 0.0));
										}
									} else {
										if(P[Bk[k].PI[i][j]].y+openDiameter/2>=Ln[Bk[k].LnI[d]].Pt[Ln[Bk[k].LnI[d]].Pt.size()-1].y) {
											Ln[Bk[k].LnI[d]].insertPoint(point_c(P[Bk[k].PI[i][j]].x, P[Bk[k].PI[i][j]].y-openDiameter/2),
																		 segment_c(Bk[k].NinLength-1),
																		 segment_c(),
																		 BC_c(Wall, 0.0, 0.0),
																		 BC_c(PorousOpen, Bk[k].PI[i][j], T[P[Bk[k].PI[i][j]].TI[0]].L));
										} else {
											Ln[Bk[k].LnI[d]].insertSegment(point_c(P[Bk[k].PI[i][j]].x, P[Bk[k].PI[i][j]].y-openDiameter/2),
															               point_c(P[Bk[k].PI[i][j]].x, P[Bk[k].PI[i][j]].y+openDiameter/2),
															               segment_c(Bk[k].NinLength-1),
															               segment_c(),
																		   segment_c(0, 1, 0, openDiameter, 0),
															               BC_c(Wall, 0.0, 0.0),
															               BC_c(PorousOpen, Bk[k].PI[i][j], T[P[Bk[k].PI[i][j]].TI[0]].L),
															               BC_c(Wall, 0.0, 0.0));
										}
									}
								} else {
									if(Bk[k].NinLength==1) {
										Ln[Bk[k].LnI[d]].insertPoint(point_c(P[Bk[k].PI[i][j]].x, P[Bk[k].PI[i][j]].y+openDiameter/2),
																	 segment_c(),
																	 segment_c(),
																	 BC_c(PorousOpen, Bk[k].PI[i][j], T[P[Bk[k].PI[i][j]].TI[0]].L),
																	 BC_c(Wall, 0.0, 0.0));
									} else {
										Ln[Bk[k].LnI[d]].insertSegment(point_c(P[Bk[k].PI[i][j]].x, P[Bk[k].PI[i][j]].y-openDiameter/2),
														               point_c(P[Bk[k].PI[i][j]].x, P[Bk[k].PI[i][j]].y+openDiameter/2),
														               segment_c(Bk[k].NinLength-1),
														               segment_c(),
																	   segment_c(),
														               BC_c(Wall, 0.0, 0.0),
														               BC_c(PorousOpen, Bk[k].PI[i][j], T[P[Bk[k].PI[i][j]].TI[0]].L),
														               BC_c(Wall, 0.0, 0.0));
									}
								}
							}
						}
					}//jump out of loop j:0-Ny
					break;
				case n: case s:
					if(d==n) {
						j=Bk[k].Ny;
						J=j-1;
					} else {
						j=0;
						J=j+1;
					}

					for(i=0; i<=Bk[k].Nx; ++i) {
						if(Bk[k].PI[i][j]==voidIndex) {
							//do nothing
						} else {
							if(Bk[k].NinDiameter!=0 && Bk[k].NinLength==0) {
								sinTheta=P[Bk[k].PI[i][j]].distanceY(P[Bk[k].PI[i][J]])/P[Bk[k].PI[i][j]].distance(P[Bk[k].PI[i][J]]);
								openDiameter=sqrt(4*T[P[Bk[k].PI[i][j]].TI[0]].A/pi)/sinTheta;

								segment_c s1, s3;
								if(i==0        || (Bk[k].PI[0       ][j]==voidIndex && i==1)) {
									s1=segment_c(0, 0, 1, 0, openDiameter/Bk[k].NinDiameter);
							    } else {
							    	s1=segment_c(Bk[k].avgTLength/(openDiameter/Bk[k].NinDiameter));
							    }
								if(i==Bk[k].Nx || (Bk[k].PI[Bk[k].Nx][j]==voidIndex && i==Bk[k].Nx-1)) s3=segment_c(0, 1, 0, openDiameter/Bk[k].NinDiameter, 0);

								Ln[Bk[k].LnI[d]].insertSegment(point_c(P[Bk[k].PI[i][j]].x-openDiameter/2, P[Bk[k].PI[i][j]].y),
										                       point_c(P[Bk[k].PI[i][j]].x+openDiameter/2, P[Bk[k].PI[i][j]].y),
										                       s1,
										                       segment_c(Bk[k].NinDiameter),
										                       s3,
										                       BC_c(Wall, 0.0, 0.0),
										                       BC_c(PorousOpen, Bk[k].PI[i][j], T[P[Bk[k].PI[i][j]].TI[0]].L),
										                       BC_c(Wall, 0.0, 0.0));
							} else if(Bk[k].NinDiameter==0 && Bk[k].NinLength!=0) {
								openDiameter=Bk[k].avgTLength/Bk[k].NinLength;

								if(i==0 || (Bk[k].PI[0       ][j]==voidIndex && i==1)) {
									if(P[Bk[k].PI[i][j]].x-openDiameter/2<=Ln[Bk[k].LnI[d]].Pt[0].x) {
										Ln[Bk[k].LnI[d]].insertPoint(point_c(P[Bk[k].PI[i][j]].x+openDiameter/2, P[Bk[k].PI[i][j]].y),
																	 segment_c(),
																	 segment_c(),
																	 BC_c(PorousOpen, Bk[k].PI[i][j], T[P[Bk[k].PI[i][j]].TI[0]].L),
																	 BC_c(Wall, 0.0, 0.0));
									} else {
										Ln[Bk[k].LnI[d]].insertSegment(point_c(P[Bk[k].PI[i][j]].x-openDiameter/2, P[Bk[k].PI[i][j]].y),
														               point_c(P[Bk[k].PI[i][j]].x+openDiameter/2, P[Bk[k].PI[i][j]].y),
														               segment_c(0, 0, 1, 0, openDiameter),
														               segment_c(),
														               segment_c(),
														               BC_c(Wall, 0.0, 0.0),
														               BC_c(PorousOpen, Bk[k].PI[i][j], T[P[Bk[k].PI[i][j]].TI[0]].L),
														               BC_c(Wall, 0.0, 0.0));
									}
								} else if(i==Bk[k].Nx || (Bk[k].PI[Bk[k].Nx][j]==voidIndex && i==Bk[k].Nx-1)) {
									if(Bk[k].NinLength==1) {
										if(P[Bk[k].PI[i][j]].x+openDiameter/2>=Ln[Bk[k].LnI[d]].Pt[Ln[Bk[k].LnI[d]].Pt.size()-1].x) {
											Ln[Bk[k].LnI[d]].B[Ln[Bk[k].LnI[d]].B.size()-1]=BC_c(PorousOpen, Bk[k].PI[i][j], T[P[Bk[k].PI[i][j]].TI[0]].L);
										} else {
											Ln[Bk[k].LnI[d]].insertPoint(point_c(P[Bk[k].PI[i][j]].x+openDiameter/2, P[Bk[k].PI[i][j]].y),
																		 segment_c(),
																		 segment_c(0, 1, 0, openDiameter, 0),
																		 BC_c(PorousOpen, Bk[k].PI[i][j], T[P[Bk[k].PI[i][j]].TI[0]].L),
																		 BC_c(Wall, 0.0, 0.0));
										}
									} else {
										if(P[Bk[k].PI[i][j]].x+openDiameter/2>=Ln[Bk[k].LnI[d]].Pt[Ln[Bk[k].LnI[d]].Pt.size()-1].x) {
											Ln[Bk[k].LnI[d]].insertPoint(point_c(P[Bk[k].PI[i][j]].x-openDiameter/2, P[Bk[k].PI[i][j]].y),
																		 segment_c(Bk[k].NinLength-1),
																		 segment_c(),
																		 BC_c(Wall, 0.0, 0.0),
																		 BC_c(PorousOpen, Bk[k].PI[i][j], T[P[Bk[k].PI[i][j]].TI[0]].L));
										} else {
											Ln[Bk[k].LnI[d]].insertSegment(point_c(P[Bk[k].PI[i][j]].x-openDiameter/2, P[Bk[k].PI[i][j]].y),
															               point_c(P[Bk[k].PI[i][j]].x+openDiameter/2, P[Bk[k].PI[i][j]].y),
															               segment_c(Bk[k].NinLength-1),
															               segment_c(),
																		   segment_c(0, 1, 0, openDiameter, 0),
															               BC_c(Wall, 0.0, 0.0),
															               BC_c(PorousOpen, Bk[k].PI[i][j], T[P[Bk[k].PI[i][j]].TI[0]].L),
															               BC_c(Wall, 0.0, 0.0));
										}
									}

								} else {
									if(Bk[k].NinLength==1) {
										Ln[Bk[k].LnI[d]].insertPoint(point_c(P[Bk[k].PI[i][j]].x+openDiameter/2, P[Bk[k].PI[i][j]].y),
																	 segment_c(),
																	 segment_c(),
																	 BC_c(PorousOpen, Bk[k].PI[i][j], T[P[Bk[k].PI[i][j]].TI[0]].L),
																	 BC_c(Wall, 0.0, 0.0));
									} else {
										Ln[Bk[k].LnI[d]].insertSegment(point_c(P[Bk[k].PI[i][j]].x-openDiameter/2, P[Bk[k].PI[i][j]].y),
														               point_c(P[Bk[k].PI[i][j]].x+openDiameter/2, P[Bk[k].PI[i][j]].y),
														               segment_c(Bk[k].NinLength-1),
														               segment_c(),
														               segment_c(),
														               BC_c(Wall, 0.0, 0.0),
														               BC_c(PorousOpen, Bk[k].PI[i][j], T[P[Bk[k].PI[i][j]].TI[0]].L),
														               BC_c(Wall, 0.0, 0.0));
									}
								}
							}
						}
					}//jump out loop i:0-Nx
					break;
				}//jump out of switch

				//copy to parallel equal lines
				for(size_t peli=0; peli<Ln[Bk[k].LnI[d]].PELI.size(); ++peli) {
					Ln[Ln[Bk[k].LnI[d]].PELI[peli]].copyParaEq(Ln[Bk[k].LnI[d]]);
				}
			}//jump out of if(it is porous open)
		}//jump out of direction loop
	}//jump out of processing the porous blocks
	cd_c::phy.DryCriteria=cd_c::phy.RelativeConvergeCriteria/T.size(); //Modify criteria value according to throat size.

	oflg<<"Lines information before filling......"<<endl;
	for(size_t L=0; L!=Ln.size(); ++L) {
		oflg<<"Line # "<<L<<":\tPt Size="<<Ln[L].Pt.size()<<"\tPtI Size="<<Ln[L].PtI.size()<<"\tB Size="<<Ln[L].B.size()<<"\tS Size="<<Ln[L].S.size()
		                  <<"\tSCI Size="<<Ln[L].SCI.size()<<"\tPELI: ";
		for(size_t ls=0; ls< Ln[L].PELI.size(); ++ls) {
			oflg<<Ln[L].PELI[ls]<<"  ";
		}
		oflg<<endl;
		for(size_t ls=0; ls<=Ln[L].S.size(); ++ls) {
			oflg<<Ln[L].Pt[ls]<<"\t";
		}
		oflg<<endl;
		for(size_t ls=0; ls<=Ln[L].S.size(); ++ls) {
			oflg<<Ln[L].PtI[ls]<<"\t";
		}
		oflg<<endl;
		for(size_t ls=0; ls< Ln[L].SCI.size(); ++ls) {
			oflg<<Ln[L].SCI[ls]<<"\t";
		}
		oflg<<endl;
		for(size_t ls=0; ls< Ln[L].S.size(); ++ls) {
			oflg<<Ln[L].B[ls];
		}
		for(size_t ls=0; ls< Ln[L].S.size(); ++ls) {
			oflg<<"("<<Ln[L].S[ls].SN<<", "<<Ln[L].S[ls].GS<<", "<<Ln[L].S[ls].GL<<", "<<Ln[L].S[ls].LS<<", "<<Ln[L].S[ls].LL<<")\t";
		}
		oflg<<endl;
	}

	cout<<"Generating points for FlowField blocks......"<<endl;
	cout<<"Filling points on lines and Updating PtI on Lines......"<<endl;
	for(size_t L=0; L!=Ln.size(); ++L) {
		Ln[L].fill();
		for(size_t i=1; i<Ln[L].Pt.size()-1; ++i) {
			Ln[L].PtI[i]=Pt.size();
			Pt.push_back(Ln[L].Pt[i]);
		}
	}

	cout<<"Updating flow field Blocks' Nx Ny according to filled Lines......"<<endl;
	cout<<"Making PtI, CI, UI, VI containers in Blocks......"<<endl;
	for(size_t k=porousBlocks; k<Bk.size(); ++k) {
		Bk[k].Nx=Ln[Bk[k].LnI[s]].S.size();
		Bk[k].Ny=Ln[Bk[k].LnI[w]].S.size();

		for(size_t i=0; i<=Bk[k].Nx; ++i) {
			vector<index_t> vertical;
			for(size_t j=0; j<=Bk[k].Ny; ++j) {
				vertical.push_back(voidIndex);
			}
			Bk[k].PtI.push_back(vertical);
			vertical.clear();
		}
		for(size_t i=0; i<Bk[k].Nx; ++i) {
			vector<index_t> vertical;
			for(size_t j=0; j<Bk[k].Ny; ++j) {
				vertical.push_back(voidIndex);
			}
			Bk[k].CI.push_back(vertical);
			vertical.clear();
		}
		for(size_t i=0; i<=Bk[k].Nx; ++i) {
			vector<index_t> vertical;
			for(size_t j=0; j<Bk[k].Ny; ++j) {
				vertical.push_back(voidIndex);
			}
			Bk[k].UI.push_back(vertical);
			vertical.clear();
		}
		for(size_t i=0; i<Bk[k].Nx; ++i) {
			vector<index_t> vertical;
			for(size_t j=0; j<=Bk[k].Ny; ++j) {
				vertical.push_back(voidIndex);
			}
			Bk[k].VI.push_back(vertical);
			vertical.clear();
		}
	}

	cout<<"Filling Points and Updating PtI in Blocks......"<<endl;
	for(size_t k=porousBlocks; k<Bk.size(); ++k) {
		cout<<"Generating Points in FlowField Block # "<<k<<"......"<<endl;
		Bk[k].PtI[0       ][0       ]=Bk[k].CnI[SW];
		Bk[k].PtI[0       ][Bk[k].Ny]=Bk[k].CnI[NW];
		Bk[k].PtI[Bk[k].Nx][0       ]=Bk[k].CnI[SE];
		Bk[k].PtI[Bk[k].Nx][Bk[k].Ny]=Bk[k].CnI[NE];

		for(size_t i=1; i<Bk[k].Nx; ++i) {
			Bk[k].PtI[i][0       ]=Ln[Bk[k].LnI[s]].PtI[i];
			Bk[k].PtI[i][Bk[k].Ny]=Ln[Bk[k].LnI[n]].PtI[i];
		}
		for(size_t j=1; j<Bk[k].Ny; ++j) {
			Bk[k].PtI[0       ][j]=Ln[Bk[k].LnI[w]].PtI[j];
			Bk[k].PtI[Bk[k].Nx][j]=Ln[Bk[k].LnI[e]].PtI[j];
		}

		for(size_t i=1; i<Bk[k].Nx; ++i) {
			for(size_t j=1; j<Bk[k].Ny; ++j) {
			    Bk[k].PtI[i][j]=Pt.size();
				Pt.push_back( point_c(Ln[Bk[k].LnI[s]].Pt[i].x, Ln[Bk[k].LnI[w]].Pt[j].y) );
			}
		}

		oflg<<"Block # "<<k<<": PtI updated as follows:"<<endl;
		for(index_t j=Bk[k].Ny; j>=0; --j) {
			for(size_t i=0; i<=Bk[k].Nx; ++i) {
				oflg<<Bk[k].PtI[i][j]<<"\t";
			}
			oflg<<endl;
		}
	}

	cout<<"Generating cCell, uCell and vCell for FlowField Blocks......"<<endl;
	cout<<"Filling uCell and vCell on Lines......"<<endl;
	for(size_t L=0; L!=Ln.size(); ++L) {
		if(Ln[L].XorY==yAxis) {//Generating uCell_c for line
			for(size_t i=0; i<Ln[L].S.size(); ++i) {
				Ln[L].SCI[i]=U.size();
				U.push_back( uCell_c(0.5*(Pt[Ln[L].PtI[i]].x+Pt[Ln[L].PtI[i+1]].x),
			                         0.5*(Pt[Ln[L].PtI[i]].y+Pt[Ln[L].PtI[i+1]].y)) );
			}
		} else if(Ln[L].XorY==xAxis) {//Generating vCell_c for line
			for(size_t i=0; i<Ln[L].S.size(); ++i) {
				Ln[L].SCI[i]=V.size();
				V.push_back( vCell_c(0.5*(Pt[Ln[L].PtI[i]].x+Pt[Ln[L].PtI[i+1]].x),
			                         0.5*(Pt[Ln[L].PtI[i]].y+Pt[Ln[L].PtI[i+1]].y)) );
			}
		} else {}
	}

	oflg<<"Lines information after filling and updating......"<<endl;
	for(size_t L=0; L!=Ln.size(); ++L) {
		oflg<<"Line # "<<L<<":\tPt Size="<<Ln[L].Pt.size()<<"\tPtI Size="<<Ln[L].PtI.size()<<"\tB Size="<<Ln[L].B.size()<<"\tS Size="<<Ln[L].S.size()
		                  <<"\tSCI Size="<<Ln[L].SCI.size()<<"\tPELI: ";
		for(size_t ls=0; ls< Ln[L].PELI.size(); ++ls) {
			oflg<<Ln[L].PELI[ls]<<"  ";
		}
		oflg<<endl;
		for(size_t ls=0; ls<=Ln[L].S.size(); ++ls) {
			oflg<<Ln[L].Pt[ls]<<"\t";
		}
		oflg<<endl;
		for(size_t ls=0; ls<=Ln[L].S.size(); ++ls) {
			oflg<<Ln[L].PtI[ls]<<"\t";
		}
		oflg<<endl;
		for(size_t ls=0; ls< Ln[L].SCI.size(); ++ls) {
			oflg<<Ln[L].SCI[ls]<<"\t";
		}
		oflg<<endl;
		for(size_t ls=0; ls< Ln[L].S.size(); ++ls) {
			oflg<<Ln[L].B[ls];
		}
		for(size_t ls=0; ls< Ln[L].S.size(); ++ls) {
			oflg<<"("<<Ln[L].S[ls].SN<<", "<<Ln[L].S[ls].GS<<", "<<Ln[L].S[ls].GL<<", "<<Ln[L].S[ls].LS<<", "<<Ln[L].S[ls].LL<<")\t";
		}
		oflg<<endl;
	}

	cout<<"Filling cCell, uCell and vCell in Blocks......"<<endl;
	for(size_t k=porousBlocks; k<Bk.size(); ++k) {
		cout<<"Generating cCell in FlowField Block # "<<k<<"......"<<endl;
		for(size_t i=0; i!=Bk[k].Nx; ++i) {
			for(size_t j=0; j!=Bk[k].Ny; ++j) {
			    Bk[k].CI[i][j]=C.size();
				C.push_back( cCell_c(Cxof(Pt[Bk[k].PtI[i][j]], Pt[Bk[k].PtI[i+1][j]], Pt[Bk[k].PtI[i+1][j+1]], Pt[Bk[k].PtI[i][j+1]]),
						             Cyof(Pt[Bk[k].PtI[i][j]], Pt[Bk[k].PtI[i+1][j]], Pt[Bk[k].PtI[i+1][j+1]], Pt[Bk[k].PtI[i][j+1]])) );
			}
		}

		oflg<<"Block # "<<k<<": CI updated as follows:"<<endl;
		for(index_t j=Bk[k].Ny-1; j>=0; --j) {
			for(size_t i=0; i<Bk[k].Nx; ++i) {
				oflg<<Bk[k].CI[i][j]<<"\t";
			}
			oflg<<endl;
		}

		cout<<"Generating uCell in FlowField Block # "<<k<<"......"<<endl;
		for(size_t j=0; j<Bk[k].Ny; ++j) {
		    Bk[k].UI[0       ][j]=Ln[Bk[k].LnI[w]].SCI[j];
		    Bk[k].UI[Bk[k].Nx][j]=Ln[Bk[k].LnI[e]].SCI[j];
			for(size_t i=1; i<Bk[k].Nx; ++i) {
			    Bk[k].UI[i][j]=U.size();
				U.push_back( uCell_c(0.5*(Pt[Bk[k].PtI[i][j+1]].x+Pt[Bk[k].PtI[i][j]].x),
			                         0.5*(Pt[Bk[k].PtI[i][j+1]].y+Pt[Bk[k].PtI[i][j]].y)) );
			}
		}

		oflg<<"Block # "<<k<<": UI updated as follows:"<<endl;
		for(index_t j=Bk[k].Ny-1; j>=0; --j) {
			for(size_t i=0; i<=Bk[k].Nx; ++i) {
				oflg<<Bk[k].UI[i][j]<<"\t";
			}
			oflg<<endl;
		}

		cout<<"Generating vCell in FlowField Block # "<<k<<"......"<<endl;
		for(size_t i=0; i<Bk[k].Nx; ++i) {
		    Bk[k].VI[i][0       ]=Ln[Bk[k].LnI[s]].SCI[i];
		    Bk[k].VI[i][Bk[k].Ny]=Ln[Bk[k].LnI[n]].SCI[i];
			for(size_t j=1; j<Bk[k].Ny; ++j) {
			    Bk[k].VI[i][j]=V.size();
				V.push_back( vCell_c(0.5*(Pt[Bk[k].PtI[i+1][j]].x+Pt[Bk[k].PtI[i][j]].x),
			                         0.5*(Pt[Bk[k].PtI[i+1][j]].y+Pt[Bk[k].PtI[i][j]].y)) );
			}
		}

		oflg<<"Block # "<<k<<": VI updated as follows:"<<endl;
		for(index_t j=Bk[k].Ny; j>=0; --j) {
			for(size_t i=0; i<Bk[k].Nx; ++i) {
				oflg<<Bk[k].VI[i][j]<<"\t";
			}
			oflg<<endl;
		}
	}

	cout<<"Generate fCells according to Points......"<<endl;
	F.resize(Pt.size(), point_c());
	for(size_t i=0; i<Pt.size(); ++i) {
		F[i]=fCell_c(Pt[i]);
	}

	cout<<"Setting up relationships for all Cells......"<<endl;//**********Setting up relationships for all Cells**********
	cout<<"Looking for close neighbors and surrounding nodes for cCell......"<<endl;
	for(size_t k=porousBlocks; k<Bk.size(); ++k) {
		for(index_t i=0; i<Bk[k].Nx; ++i) {
			for(index_t j=0; j<Bk[k].Ny; ++j) {

				C[Bk[k].CI[i][j]].NI[E]=(i+1)< Bk[k].Nx?Bk[k].CI[i+1][j  ]:
						                                (Bk[k].B[E].T==NoBC?Bk[Bk[k].NI[E]].CI[0                   ][j                   ]:voidIndex);
				C[Bk[k].CI[i][j]].NI[N]=(j+1)< Bk[k].Ny?Bk[k].CI[i  ][j+1]:
						                                (Bk[k].B[N].T==NoBC?Bk[Bk[k].NI[N]].CI[i                   ][0                   ]:voidIndex);
				C[Bk[k].CI[i][j]].NI[W]=(i-1)>=0       ?Bk[k].CI[i-1][j  ]:
						                                (Bk[k].B[W].T==NoBC?Bk[Bk[k].NI[W]].CI[Bk[Bk[k].NI[W]].Nx-1][j                   ]:voidIndex);
				C[Bk[k].CI[i][j]].NI[S]=(j-1)>=0       ?Bk[k].CI[i  ][j-1]:
						                                (Bk[k].B[S].T==NoBC?Bk[Bk[k].NI[S]].CI[i                   ][Bk[Bk[k].NI[S]].Ny-1]:voidIndex);

				C[Bk[k].CI[i][j]].FI[e ]=Bk[k].UI[i+1][j  ];
				C[Bk[k].CI[i][j]].FI[n ]=Bk[k].VI[i  ][j+1];
				C[Bk[k].CI[i][j]].FI[w ]=Bk[k].UI[i  ][j  ];
				C[Bk[k].CI[i][j]].FI[s ]=Bk[k].VI[i  ][j  ];

				C[Bk[k].CI[i][j]].FI[ne]=Bk[k].PtI[i+1][j+1];
				C[Bk[k].CI[i][j]].FI[nw]=Bk[k].PtI[i  ][j+1];
				C[Bk[k].CI[i][j]].FI[sw]=Bk[k].PtI[i  ][j  ];
				C[Bk[k].CI[i][j]].FI[se]=Bk[k].PtI[i+1][j  ];
			}
		}
	}

	cout<<"Looking for close neighbors and side nodes for uCell......"<<endl;
	for(size_t k=porousBlocks; k<Bk.size(); ++k) {
		for(index_t j=0; j<Bk[k].Ny; ++j) {
			for(index_t i=0; i<=Bk[k].Nx; ++i) {
				U[Bk[k].UI[i][j]].NI[E]=(i+1)<=Bk[k].Nx?Bk[k].UI[i+1][j  ]:
						                                (Bk[k].B[E].T==NoBC?Bk[Bk[k].NI[E]].UI[1                   ][j                   ]:voidIndex);
				U[Bk[k].UI[i][j]].NI[W]=(i-1)>=0       ?Bk[k].UI[i-1][j  ]:
						                                (Bk[k].B[W].T==NoBC?Bk[Bk[k].NI[W]].UI[Bk[Bk[k].NI[W]].Nx-1][j                   ]:voidIndex);
				if(i==0) {
					U[Bk[k].UI[i][j]].NI[N]=(j+1)< Bk[k].Ny?Bk[k].UI[i  ][j+1]:
							                                (Bk[k].B[N].T==NoBC?Bk[Bk[k].NI[N]].UI[i                   ][0                   ]:
							                                ((Bk[k].B[W].T==NoBC && Bk[Bk[k].NI[W]].B[N].T==NoBC)?Bk[Bk[Bk[k].NI[W]].NI[N]].UI[Bk[Bk[Bk[k].NI[W]].NI[N]].Nx][0]:voidIndex));
					U[Bk[k].UI[i][j]].NI[S]=(j-1)>=0       ?Bk[k].UI[i  ][j-1]:
							                                (Bk[k].B[S].T==NoBC?Bk[Bk[k].NI[S]].UI[i                   ][Bk[Bk[k].NI[S]].Ny-1]:
							                                ((Bk[k].B[W].T==NoBC && Bk[Bk[k].NI[W]].B[S].T==NoBC)?
							                                		                 Bk[Bk[Bk[k].NI[W]].NI[S]].UI[Bk[Bk[Bk[k].NI[W]].NI[S]].Nx][Bk[Bk[Bk[k].NI[W]].NI[S]].Ny-1]:voidIndex));
				} else if(i==Bk[k].Nx) {
					U[Bk[k].UI[i][j]].NI[N]=(j+1)< Bk[k].Ny?Bk[k].UI[i  ][j+1]:
							                                (Bk[k].B[N].T==NoBC?Bk[Bk[k].NI[N]].UI[i                   ][0                   ]:
							                                ((Bk[k].B[E].T==NoBC && Bk[Bk[k].NI[E]].B[N].T==NoBC)?Bk[Bk[Bk[k].NI[E]].NI[N]].UI[0][0]:voidIndex));
					U[Bk[k].UI[i][j]].NI[S]=(j-1)>=0       ?Bk[k].UI[i  ][j-1]:
							                                (Bk[k].B[S].T==NoBC?Bk[Bk[k].NI[S]].UI[i                   ][Bk[Bk[k].NI[S]].Ny-1]:
							                                ((Bk[k].B[E].T==NoBC && Bk[Bk[k].NI[E]].B[S].T==NoBC)?
							                                		                 Bk[Bk[Bk[k].NI[E]].NI[S]].UI[0][Bk[Bk[Bk[k].NI[E]].NI[S]].Ny-1]:voidIndex));
				} else {
					U[Bk[k].UI[i][j]].NI[N]=(j+1)< Bk[k].Ny?Bk[k].UI[i  ][j+1]:
							                                (Bk[k].B[N].T==NoBC?Bk[Bk[k].NI[N]].UI[i                   ][0                   ]:voidIndex);
					U[Bk[k].UI[i][j]].NI[S]=(j-1)>=0       ?Bk[k].UI[i  ][j-1]:
							                                (Bk[k].B[S].T==NoBC?Bk[Bk[k].NI[S]].UI[i                   ][Bk[Bk[k].NI[S]].Ny-1]:voidIndex);
				}

				U[Bk[k].UI[i][j]].FI[e]=i<Bk[k].Nx?Bk[k].CI[i  ][j  ]:
                                                   (Bk[k].B[E].T==NoBC?Bk[Bk[k].NI[E]].CI[0                   ][j                   ]:voidIndex);
				U[Bk[k].UI[i][j]].FI[w]=i>0       ?Bk[k].CI[i-1][j  ]:
                                                   (Bk[k].B[W].T==NoBC?Bk[Bk[k].NI[W]].CI[Bk[Bk[k].NI[W]].Nx-1][j                   ]:voidIndex);
				U[Bk[k].UI[i][j]].FI[n]=Bk[k].PtI[i  ][j+1];
				U[Bk[k].UI[i][j]].FI[s]=Bk[k].PtI[i  ][j  ];
			}
		}
	}

	cout<<"Looking for close neighbors and side nodes for vCell......"<<endl;
	for(size_t k=porousBlocks; k<Bk.size(); ++k) {
		for(index_t i=0; i<Bk[k].Nx; ++i) {
			for(index_t j=0; j<=Bk[k].Ny; ++j) {
				V[Bk[k].VI[i][j]].NI[N]=(j+1)<=Bk[k].Ny?Bk[k].VI[i  ][j+1]:
						                                (Bk[k].B[N].T==NoBC?Bk[Bk[k].NI[N]].VI[i                   ][1                   ]:voidIndex);
				V[Bk[k].VI[i][j]].NI[S]=(j-1)>=0       ?Bk[k].VI[i  ][j-1]:
						                                (Bk[k].B[S].T==NoBC?Bk[Bk[k].NI[S]].VI[i                   ][Bk[Bk[k].NI[S]].Ny-1]:voidIndex);
				if(j==0) {
					V[Bk[k].VI[i][j]].NI[E]=(i+1)< Bk[k].Nx?Bk[k].VI[i+1][j  ]:
							                                (Bk[k].B[E].T==NoBC?Bk[Bk[k].NI[E]].VI[0                   ][j                   ]:
							                                ((Bk[k].B[S].T==NoBC && Bk[Bk[k].NI[S]].B[E].T==NoBC)?Bk[Bk[Bk[k].NI[S]].NI[E]].VI[0][Bk[Bk[Bk[k].NI[S]].NI[E]].Ny]:voidIndex));
					V[Bk[k].VI[i][j]].NI[W]=(i-1)>=0       ?Bk[k].VI[i-1][j  ]:
							                                (Bk[k].B[W].T==NoBC?Bk[Bk[k].NI[W]].VI[Bk[Bk[k].NI[W]].Nx-1][j                   ]:
									                        ((Bk[k].B[S].T==NoBC && Bk[Bk[k].NI[S]].B[W].T==NoBC)?
									                        		                 Bk[Bk[Bk[k].NI[S]].NI[W]].VI[Bk[Bk[Bk[k].NI[S]].NI[W]].Nx-1][Bk[Bk[Bk[k].NI[S]].NI[W]].Ny]:voidIndex));
				} else if(j==Bk[k].Ny) {
					V[Bk[k].VI[i][j]].NI[E]=(i+1)< Bk[k].Nx?Bk[k].VI[i+1][j  ]:
							                                (Bk[k].B[E].T==NoBC?Bk[Bk[k].NI[E]].VI[0                   ][j                   ]:
							                                ((Bk[k].B[N].T==NoBC && Bk[Bk[k].NI[N]].B[E].T==NoBC)?Bk[Bk[Bk[k].NI[N]].NI[E]].VI[0][0]:voidIndex));
					V[Bk[k].VI[i][j]].NI[W]=(i-1)>=0       ?Bk[k].VI[i-1][j  ]:
							                                (Bk[k].B[W].T==NoBC?Bk[Bk[k].NI[W]].VI[Bk[Bk[k].NI[W]].Nx-1][j                   ]:
									                        ((Bk[k].B[N].T==NoBC && Bk[Bk[k].NI[N]].B[W].T==NoBC)?
									                        		                 Bk[Bk[Bk[k].NI[N]].NI[W]].VI[Bk[Bk[Bk[k].NI[N]].NI[W]].Nx-1][0]:voidIndex));
				} else {
					V[Bk[k].VI[i][j]].NI[E]=(i+1)< Bk[k].Nx?Bk[k].VI[i+1][j  ]:
							                                (Bk[k].B[E].T==NoBC?Bk[Bk[k].NI[E]].VI[0                   ][j                   ]:voidIndex);
					V[Bk[k].VI[i][j]].NI[W]=(i-1)>=0       ?Bk[k].VI[i-1][j  ]:
							                                (Bk[k].B[W].T==NoBC?Bk[Bk[k].NI[W]].VI[Bk[Bk[k].NI[W]].Nx-1][j                   ]:voidIndex);
				}

				V[Bk[k].VI[i][j]].FI[e]=Bk[k].PtI[i+1][j  ];
				V[Bk[k].VI[i][j]].FI[w]=Bk[k].PtI[i  ][j  ];
				V[Bk[k].VI[i][j]].FI[n]=j<Bk[k].Ny?Bk[k].CI[i  ][j  ]:
                                                   (Bk[k].B[N].T==NoBC?Bk[Bk[k].NI[N]].CI[i                   ][0                   ]:voidIndex);
				V[Bk[k].VI[i][j]].FI[s]=j>0       ?Bk[k].CI[i  ][j-1]:
                                                   (Bk[k].B[S].T==NoBC?Bk[Bk[k].NI[S]].CI[i                   ][Bk[Bk[k].NI[S]].Ny-1]:voidIndex);
			}
		}
	}

	cout<<"Looking for close neighbors and side nodes for fCell......"<<endl;
	for(size_t k=porousBlocks; k<Bk.size(); ++k) {
		for(index_t i=0; i<=Bk[k].Nx; ++i) {
			for(index_t j=0; j<=Bk[k].Ny; ++j) {
				if(i==0) {
					if(j==0) {
						F[Bk[k].PtI[i][j]].NI[W ]=Bk[k].B[W].T==NoBC?Bk[Bk[k].NI[w]].PtI[Bk[Bk[k].NI[w]].Nx-1][j]:
						    ((Bk[k].B[S].T==NoBC && Bk[Bk[k].NI[s]].B[W].T==NoBC)?Bk[Bk[k].NI[sw]].PtI[Bk[Bk[k].NI[sw]].Nx-1][Bk[Bk[k].NI[sw]].Ny]:voidIndex);
						F[Bk[k].PtI[i][j]].FI[w ]=Bk[k].B[W].T==NoBC?Bk[Bk[k].NI[w]]. VI[Bk[Bk[k].NI[w]].Nx-1][j]:
				            ((Bk[k].B[S].T==NoBC && Bk[Bk[k].NI[s]].B[W].T==NoBC)?Bk[Bk[k].NI[sw]]. VI[Bk[Bk[k].NI[sw]].Nx-1][Bk[Bk[k].NI[sw]].Ny]:voidIndex);
					} else if(j==Bk[k].Ny) {
						F[Bk[k].PtI[i][j]].NI[W ]=Bk[k].B[W].T==NoBC?Bk[Bk[k].NI[w]].PtI[Bk[Bk[k].NI[w]].Nx-1][j]:
						    ((Bk[k].B[N].T==NoBC && Bk[Bk[k].NI[n]].B[W].T==NoBC)?Bk[Bk[k].NI[nw]].PtI[Bk[Bk[k].NI[nw]].Nx-1][0                  ]:voidIndex);
						F[Bk[k].PtI[i][j]].FI[w ]=Bk[k].B[W].T==NoBC?Bk[Bk[k].NI[w]]. VI[Bk[Bk[k].NI[w]].Nx-1][j]:
				            ((Bk[k].B[N].T==NoBC && Bk[Bk[k].NI[n]].B[W].T==NoBC)?Bk[Bk[k].NI[nw]]. VI[Bk[Bk[k].NI[nw]].Nx-1][0                  ]:voidIndex);
					} else {
						F[Bk[k].PtI[i][j]].NI[W ]=Bk[k].B[W].T==NoBC?Bk[Bk[k].NI[w]].PtI[Bk[Bk[k].NI[w]].Nx-1][j]:voidIndex;
						F[Bk[k].PtI[i][j]].FI[w ]=Bk[k].B[W].T==NoBC?Bk[Bk[k].NI[w]]. VI[Bk[Bk[k].NI[w]].Nx-1][j]:voidIndex;
					}
				} else {
					F[Bk[k].PtI[i][j]].NI[W ]=Bk[k].PtI[i-1][j  ];
					F[Bk[k].PtI[i][j]].FI[w ]=Bk[k]. VI[i-1][j  ];
				}

                if(i==Bk[k].Nx) {
					if(j==0) {
						F[Bk[k].PtI[i][j]].NI[E ]=Bk[k].B[E].T==NoBC?Bk[Bk[k].NI[e]].PtI[1                   ][j]:
						    ((Bk[k].B[S].T==NoBC && Bk[Bk[k].NI[s]].B[E].T==NoBC)?Bk[Bk[k].NI[se]].PtI[1                    ][Bk[Bk[k].NI[se]].Ny]:voidIndex);
						F[Bk[k].PtI[i][j]].FI[e ]=Bk[k].B[E].T==NoBC?Bk[Bk[k].NI[e]]. VI[0                   ][j]:
				            ((Bk[k].B[S].T==NoBC && Bk[Bk[k].NI[s]].B[E].T==NoBC)?Bk[Bk[k].NI[se]]. VI[0                    ][Bk[Bk[k].NI[se]].Ny]:voidIndex);
					} else if(j==Bk[k].Ny) {
						F[Bk[k].PtI[i][j]].NI[E ]=Bk[k].B[E].T==NoBC?Bk[Bk[k].NI[e]].PtI[1                   ][j]:
						    ((Bk[k].B[N].T==NoBC && Bk[Bk[k].NI[n]].B[E].T==NoBC)?Bk[Bk[k].NI[ne]].PtI[1                    ][0                  ]:voidIndex);
						F[Bk[k].PtI[i][j]].FI[e ]=Bk[k].B[E].T==NoBC?Bk[Bk[k].NI[e]]. VI[0                   ][j]:
				            ((Bk[k].B[N].T==NoBC && Bk[Bk[k].NI[n]].B[E].T==NoBC)?Bk[Bk[k].NI[ne]]. VI[0                    ][0                  ]:voidIndex);
					} else {
						F[Bk[k].PtI[i][j]].NI[E ]=Bk[k].B[E].T==NoBC?Bk[Bk[k].NI[e]].PtI[1                   ][j]:voidIndex;
						F[Bk[k].PtI[i][j]].FI[e ]=Bk[k].B[E].T==NoBC?Bk[Bk[k].NI[e]]. VI[0                   ][j]:voidIndex;
					}
				} else {
					F[Bk[k].PtI[i][j]].NI[E ]=Bk[k].PtI[i+1][j  ];
					F[Bk[k].PtI[i][j]].FI[e ]=Bk[k]. VI[i  ][j  ];
				}

				if(j==0) {
					if(i==0) {
						F[Bk[k].PtI[i][j]].NI[S ]=Bk[k].B[S].T==NoBC?Bk[Bk[k].NI[s]].PtI[i][Bk[Bk[k].NI[s]].Ny-1]:
						    ((Bk[k].B[W].T==NoBC && Bk[Bk[k].NI[w]].B[S].T==NoBC)?Bk[Bk[k].NI[sw]].PtI[Bk[Bk[k].NI[sw]].Nx][Bk[Bk[k].NI[sw]].Ny-1]:voidIndex);
						F[Bk[k].PtI[i][j]].FI[s ]=Bk[k].B[S].T==NoBC?Bk[Bk[k].NI[s]]. UI[i][Bk[Bk[k].NI[s]].Ny-1]:
				            ((Bk[k].B[W].T==NoBC && Bk[Bk[k].NI[w]].B[S].T==NoBC)?Bk[Bk[k].NI[sw]]. UI[Bk[Bk[k].NI[sw]].Nx][Bk[Bk[k].NI[sw]].Ny-1]:voidIndex);
					} else if(i==Bk[k].Nx) {
						F[Bk[k].PtI[i][j]].NI[S ]=Bk[k].B[S].T==NoBC?Bk[Bk[k].NI[s]].PtI[i][Bk[Bk[k].NI[s]].Ny-1]:
						    ((Bk[k].B[E].T==NoBC && Bk[Bk[k].NI[e]].B[S].T==NoBC)?Bk[Bk[k].NI[se]].PtI[0                  ][Bk[Bk[k].NI[se]].Ny-1]:voidIndex);
						F[Bk[k].PtI[i][j]].FI[s ]=Bk[k].B[S].T==NoBC?Bk[Bk[k].NI[s]]. UI[i][Bk[Bk[k].NI[s]].Ny-1]:
				            ((Bk[k].B[E].T==NoBC && Bk[Bk[k].NI[e]].B[S].T==NoBC)?Bk[Bk[k].NI[se]]. UI[0                  ][Bk[Bk[k].NI[se]].Ny-1]:voidIndex);
					} else {
						F[Bk[k].PtI[i][j]].NI[S ]=Bk[k].B[S].T==NoBC?Bk[Bk[k].NI[s]].PtI[i][Bk[Bk[k].NI[s]].Ny-1]:voidIndex;
						F[Bk[k].PtI[i][j]].FI[s ]=Bk[k].B[S].T==NoBC?Bk[Bk[k].NI[s]]. UI[i][Bk[Bk[k].NI[s]].Ny-1]:voidIndex;
					}
				} else {
					F[Bk[k].PtI[i][j]].NI[S ]=Bk[k].PtI[i  ][j-1];
					F[Bk[k].PtI[i][j]].FI[s ]=Bk[k]. UI[i  ][j-1];
				}

				if(j==Bk[k].Ny) {
					if(i==0) {
						F[Bk[k].PtI[i][j]].NI[N ]=Bk[k].B[N].T==NoBC?Bk[Bk[k].NI[n]].PtI[i][1                   ]:
						    ((Bk[k].B[W].T==NoBC && Bk[Bk[k].NI[w]].B[N].T==NoBC)?Bk[Bk[k].NI[nw]].PtI[Bk[Bk[k].NI[nw]].Nx][1                    ]:voidIndex);
						F[Bk[k].PtI[i][j]].FI[n ]=Bk[k].B[N].T==NoBC?Bk[Bk[k].NI[n]]. UI[i][0                   ]:
				            ((Bk[k].B[W].T==NoBC && Bk[Bk[k].NI[w]].B[N].T==NoBC)?Bk[Bk[k].NI[nw]]. UI[Bk[Bk[k].NI[nw]].Nx][0                    ]:voidIndex);
					} else if(i==Bk[k].Nx) {
						F[Bk[k].PtI[i][j]].NI[N ]=Bk[k].B[N].T==NoBC?Bk[Bk[k].NI[n]].PtI[i][1                   ]:
						    ((Bk[k].B[E].T==NoBC && Bk[Bk[k].NI[e]].B[N].T==NoBC)?Bk[Bk[k].NI[ne]].PtI[0                  ][1                    ]:voidIndex);
						F[Bk[k].PtI[i][j]].FI[n ]=Bk[k].B[N].T==NoBC?Bk[Bk[k].NI[n]]. UI[i][0                   ]:
				            ((Bk[k].B[E].T==NoBC && Bk[Bk[k].NI[e]].B[N].T==NoBC)?Bk[Bk[k].NI[ne]]. UI[0                  ][0                    ]:voidIndex);
					} else {
						F[Bk[k].PtI[i][j]].NI[N ]=Bk[k].B[N].T==NoBC?Bk[Bk[k].NI[n]].PtI[i][1                   ]:voidIndex;
						F[Bk[k].PtI[i][j]].FI[n ]=Bk[k].B[N].T==NoBC?Bk[Bk[k].NI[n]]. UI[i][0                   ]:voidIndex;
					}
				} else {
					F[Bk[k].PtI[i][j]].NI[N ]=Bk[k].PtI[i  ][j+1];
					F[Bk[k].PtI[i][j]].FI[n ]=Bk[k]. UI[i  ][j  ];
				}
			}
		}
	}

	cout<<"Looking for far neighbors for all cCells......"<<endl;
	for(size_t i=0; i<C.size(); ++i) {
		C[i].NI[EE]= C[i].NI[E]!=voidIndex? C[C[i].NI[E]].NI[E]:voidIndex;
		C[i].NI[NN]= C[i].NI[N]!=voidIndex? C[C[i].NI[N]].NI[N]:voidIndex;
		C[i].NI[WW]= C[i].NI[W]!=voidIndex? C[C[i].NI[W]].NI[W]:voidIndex;
		C[i].NI[SS]= C[i].NI[S]!=voidIndex? C[C[i].NI[S]].NI[S]:voidIndex;
	}

	cout<<"Looking for far neighbors and corner nodes for all uCells......"<<endl;
	for(size_t i=0; i<U.size(); ++i) {
		U[i].NI[EE]= U[i].NI[E]!=voidIndex? U[U[i].NI[E]].NI[E]:voidIndex;
		U[i].NI[NN]= U[i].NI[N]!=voidIndex? U[U[i].NI[N]].NI[N]:voidIndex;
		U[i].NI[WW]= U[i].NI[W]!=voidIndex? U[U[i].NI[W]].NI[W]:voidIndex;
		U[i].NI[SS]= U[i].NI[S]!=voidIndex? U[U[i].NI[S]].NI[S]:voidIndex;

		U[i].FI[ne]= U[i].FI[e]!=voidIndex? C[U[i].FI[e]].FI[n]:voidIndex;
		U[i].FI[nw]= U[i].FI[w]!=voidIndex? C[U[i].FI[w]].FI[n]:voidIndex;
		U[i].FI[sw]= U[i].FI[w]!=voidIndex? C[U[i].FI[w]].FI[s]:voidIndex;
		U[i].FI[se]= U[i].FI[e]!=voidIndex? C[U[i].FI[e]].FI[s]:voidIndex;
	}

	cout<<"Looking for far neighbors and corner nodes for all vCells......"<<endl;
	for(size_t i=0; i<V.size(); ++i) {
		V[i].NI[EE]= V[i].NI[E]!=voidIndex? V[V[i].NI[E]].NI[E]:voidIndex;
		V[i].NI[NN]= V[i].NI[N]!=voidIndex? V[V[i].NI[N]].NI[N]:voidIndex;
		V[i].NI[WW]= V[i].NI[W]!=voidIndex? V[V[i].NI[W]].NI[W]:voidIndex;
		V[i].NI[SS]= V[i].NI[S]!=voidIndex? V[V[i].NI[S]].NI[S]:voidIndex;

		V[i].FI[ne]= V[i].FI[n]!=voidIndex? C[V[i].FI[n]].FI[e]:voidIndex;
		V[i].FI[nw]= V[i].FI[n]!=voidIndex? C[V[i].FI[n]].FI[w]:voidIndex;
		V[i].FI[sw]= V[i].FI[s]!=voidIndex? C[V[i].FI[s]].FI[w]:voidIndex;
		V[i].FI[se]= V[i].FI[s]!=voidIndex? C[V[i].FI[s]].FI[e]:voidIndex;
	}

	cout<<"Looking for far neighbors and corner nodes for all fCells......"<<endl;
	for(size_t i=0; i<F.size(); ++i) {
		F[i].NI[EE]= F[i].NI[E]!=voidIndex? F[F[i].NI[E]].NI[E]:voidIndex;
		F[i].NI[NN]= F[i].NI[N]!=voidIndex? F[F[i].NI[N]].NI[N]:voidIndex;
		F[i].NI[WW]= F[i].NI[W]!=voidIndex? F[F[i].NI[W]].NI[W]:voidIndex;
		F[i].NI[SS]= F[i].NI[S]!=voidIndex? F[F[i].NI[S]].NI[S]:voidIndex;

		index_t fineu= F[i].FI[n]!=voidIndex? U[F[i].FI[n]].FI[e]:voidIndex;
		index_t finwu= F[i].FI[n]!=voidIndex? U[F[i].FI[n]].FI[w]:voidIndex;
		index_t fiswu= F[i].FI[s]!=voidIndex? U[F[i].FI[s]].FI[w]:voidIndex;
		index_t fiseu= F[i].FI[s]!=voidIndex? U[F[i].FI[s]].FI[e]:voidIndex;

		index_t finev= F[i].FI[e]!=voidIndex? V[F[i].FI[e]].FI[n]:voidIndex;
		index_t finwv= F[i].FI[w]!=voidIndex? V[F[i].FI[w]].FI[n]:voidIndex;
		index_t fiswv= F[i].FI[w]!=voidIndex? V[F[i].FI[w]].FI[s]:voidIndex;
		index_t fisev= F[i].FI[e]!=voidIndex? V[F[i].FI[e]].FI[s]:voidIndex;

		if(fineu==finev) {
			F[i].FI[ne]= fineu;
		} else {
			cout<<"FI[ne] error"<<endl;
		}

		if(finwu==finwv) {
			F[i].FI[ne]= finwu;
		} else {
			cout<<"FI[nw] error"<<endl;
		}

		if(fiswu==fiswv) {
			F[i].FI[ne]= fiswu;
		} else {
			cout<<"FI[sw] error"<<endl;
		}

		if(fiseu==fisev) {
			F[i].FI[ne]= fiseu;
		} else {
			cout<<"FI[se] error"<<endl;
		}
	}

	for(size_t k=porousBlocks; k<Bk.size(); ++k) {
		oflg<<"Block # "<<k<<", Checking Information:"<<endl;
		oflg<<"PtI(xCoordinate, yCoordinate):"<<endl;
		for(index_t j=Bk[k].Ny; j>=0; --j) {
			for(size_t i=0; i<=Bk[k].Nx; ++i) {
				oflg<<Bk[k].PtI[i][j]<<Pt[Bk[k].PtI[i][j]]<<"\t";
			}
			oflg<<endl;
		}

		oflg<<"CI(xCoordinate, yCoordinate):"<<endl;
		for(index_t j=Bk[k].Ny-1; j>=0; --j) {
			for(size_t i=0; i<Bk[k].Nx; ++i) {
				oflg<<Bk[k].CI[i][j]<<"("<<C[Bk[k].CI[i][j]].x<<", "<<C[Bk[k].CI[i][j]].y<<")\t";
			}
			oflg<<endl;
		}

		oflg<<"UI(xCoordinate, yCoordinate):"<<endl;
		for(index_t j=Bk[k].Ny-1; j>=0; --j) {
			for(size_t i=0; i<=Bk[k].Nx; ++i) {
				oflg<<Bk[k].UI[i][j]<<"("<<U[Bk[k].UI[i][j]].x<<", "<<U[Bk[k].UI[i][j]].y<<")\t";
			}
			oflg<<endl;
		}

		oflg<<"VI(xCoordinate, yCoordinate):"<<endl;
		for(index_t j=Bk[k].Ny; j>=0; --j) {
			for(size_t i=0; i<Bk[k].Nx; ++i) {
				oflg<<Bk[k].VI[i][j]<<"("<<V[Bk[k].VI[i][j]].x<<", "<<V[Bk[k].VI[i][j]].y<<")\t";
			}
			oflg<<endl;
		}
	}

	cout<<"Calculating dimensions of cells......"<<endl;
	for(size_t i=0; i<C.size(); ++i) {
		C[i].d[e]=C[i].distance( U[C[i].FI[e]]);
		C[i].d[n]=C[i].distance( V[C[i].FI[n]]);
		C[i].d[w]=C[i].distance( U[C[i].FI[w]]);
		C[i].d[s]=C[i].distance( V[C[i].FI[s]]);

		C[i].Ax  =C[i].dy()*cd_c::dz;
		C[i].Ay  =C[i].dx()*cd_c::dz;
		C[i].V   =C[i].dx()*C[i].dy()*cd_c::dz;

		C[i].D[E]=C[i].NI[E]!=voidIndex?C[i].distance( C[C[i].NI[E]]):0;
		C[i].D[N]=C[i].NI[N]!=voidIndex?C[i].distance( C[C[i].NI[N]]):0;
		C[i].D[W]=C[i].NI[W]!=voidIndex?C[i].distance( C[C[i].NI[W]]):0;
		C[i].D[S]=C[i].NI[S]!=voidIndex?C[i].distance( C[C[i].NI[S]]):0;
	}

	for(size_t i=0; i<U.size(); ++i) {
		U[i].d[e]=U[i].FI[e]!=voidIndex?U[i].distance( C[U[i].FI[e]]):0;
		U[i].d[n]=U[i].distance(Pt[U[i].FI[n]]);
		U[i].d[w]=U[i].FI[w]!=voidIndex?U[i].distance( C[U[i].FI[w]]):0;
		U[i].d[s]=U[i].distance(Pt[U[i].FI[s]]);

		U[i].Ax  =U[i].dy()*cd_c::dz;
		U[i].Ay  =U[i].dx()*cd_c::dz;
		U[i].V   =U[i].dx()*U[i].dy()*cd_c::dz;

		U[i].D[E]=U[i].NI[E]!=voidIndex?U[i].distance( U[U[i].NI[E]]):0;
		U[i].D[N]=U[i].NI[N]!=voidIndex?U[i].distance( U[U[i].NI[N]]):0;
		U[i].D[W]=U[i].NI[W]!=voidIndex?U[i].distance( U[U[i].NI[W]]):0;
		U[i].D[S]=U[i].NI[S]!=voidIndex?U[i].distance( U[U[i].NI[S]]):0;
	}

	for(size_t i=0; i<V.size(); ++i) {
		V[i].d[e]=V[i].distance(Pt[V[i].FI[e]]);
		V[i].d[n]=V[i].FI[n]!=voidIndex?V[i].distance( C[V[i].FI[n]]):0;
		V[i].d[w]=V[i].distance(Pt[V[i].FI[w]]);
		V[i].d[s]=V[i].FI[s]!=voidIndex?V[i].distance( C[V[i].FI[s]]):0;

		V[i].Ax  =V[i].dy()*cd_c::dz;
		V[i].Ay  =V[i].dx()*cd_c::dz;
		V[i].V   =V[i].dx()*V[i].dy()*cd_c::dz;

		V[i].D[E]=V[i].NI[E]!=voidIndex?V[i].distance( V[V[i].NI[E]]):0;
		V[i].D[N]=V[i].NI[N]!=voidIndex?V[i].distance( V[V[i].NI[N]]):0;
		V[i].D[W]=V[i].NI[W]!=voidIndex?V[i].distance( V[V[i].NI[W]]):0;
		V[i].D[S]=V[i].NI[S]!=voidIndex?V[i].distance( V[V[i].NI[S]]):0;
	}

	cout<<"Setting Boundary Conditions for cells......"<<endl;
	cout<<"Inlet, Outlet, PeriodicUpstream, PeriodicDownstream, VelocityInlet, PressureOutlet go first......"<<endl;
	for(size_t k=porousBlocks; k<Bk.size(); ++k) {
		size_t i(0), j(0), d(e);
		i=Bk[k].Nx; d=e;
		if(Bk[k].B[d].T==Inlet) {
			for(size_t j=0; j< Bk[k].Ny; ++j) {
				U[Bk[k]. UI[i  ][j]].FB=CBC_c(Inlet, 1, Bk[k].B[d].u);
				U[Bk[k]. UI[i  ][j]].NormalA=-U[Bk[k].UI[i  ][j]].Ax; uvCell_c::InletUI.push_back(Bk[k].UI[i  ][j]);
				C[Bk[k]. CI[i-1][j]].FB=CBC_c(Inlet, 1, Bk[k].B[d].p);
				C[Bk[k]. CI[i-1][j]].EB=CBC_c(Inlet, 1, Bk[k].B[d].c);
			}
			for(size_t j=0; j<=Bk[k].Ny; ++j) {
				V[Bk[k]. VI[i-1][j]].FB=CBC_c(Inlet, 1, Bk[k].B[d].v);
			}
		} else if(Bk[k].B[d].T==Outlet) {
			for(size_t j=0; j< Bk[k].Ny; ++j) {
				U[Bk[k]. UI[i  ][j]].FB=CBC_c(Outlet, 1, 0, Bk[k].UI[i-1][j], -uvCell_c::MassRatio);
				U[Bk[k]. UI[i  ][j]].NormalA=U[Bk[k].UI[i  ][j]].Ax; uvCell_c::OutletUI.push_back(Bk[k].UI[i  ][j]);
//						C[Bk[k]. CI[i-1][j]].FB=CBC_c(Outlet, 1, 0, Bk[k].CI[i-2][j], -1);
				C[Bk[k]. CI[i-1][j]].EB=CBC_c(Outlet, 1, 0, Bk[k].CI[i-2][j], -1);
			}
			for(size_t j=0; j<=Bk[k].Ny; ++j) {
				V[Bk[k]. VI[i-1][j]].FB=CBC_c(Outlet, 1, 0, Bk[k].VI[i-2][j], -1);
			}
		} else if(Bk[k].B[d].T==PeriodicDownstream) {
			for(size_t j=0; j< Bk[k].Ny; ++j) {
				U[Bk[k]. UI[i  ][j]].FB=CBC_c(PeriodicDownstream, 1, 0, Bk[k].UI[i-1][j], -uvCell_c::MassRatio);
				U[Bk[k]. UI[i  ][j]].NormalA=U[Bk[k].UI[i  ][j]].Ax; uvCell_c::OutletUI.push_back(Bk[k].UI[i  ][j]);
				C[Bk[k]. CI[i-1][j]].FB=CBC_c(PeriodicDownstream, 1, Bk[k].B[d].p+cCell_c::dpdd*C[Bk[k].CI[i-1][j]].d[e]);
				C[Bk[k]. CI[i-1][j]].EB=CBC_c(PeriodicDownstream, 1, 0, Bk[k].CI[i-2][j], -1);
			}
			for(size_t j=0; j<=Bk[k].Ny; ++j) {
				V[Bk[k]. VI[i-1][j]].FB=CBC_c(PeriodicDownstream, 1, 0, Bk[k].VI[i-2][j], -1);
			}
		} else if(Bk[k].B[d].T==PeriodicUpstream) {
			for(size_t j=0; j< Bk[k].Ny; ++j) {
				U[Bk[k]. UI[i  ][j]].FB=CBC_c(PeriodicUpstream, 1, 0, Ln[Bk[k].B[d].DownstreamLI].SCI[j], -1);
				U[Bk[k]. UI[i  ][j]].NormalA=-U[Bk[k].UI[i  ][j]].Ax; uvCell_c::InletUI.push_back(Bk[k].UI[i  ][j]);
				C[Bk[k]. CI[i-1][j]].FB=CBC_c(PeriodicUpstream, 1, Bk[k].B[d].p-cCell_c::dpdd*C[Bk[k].CI[i-1][j]].d[e]);
				C[Bk[k]. CI[i-1][j]].EB=CBC_c(PeriodicUpstream, 1, Bk[k].B[d].c);
				V[Bk[k]. VI[i-1][j]].FB=CBC_c(PeriodicUpstream, 1, 0, U[Ln[Bk[k].B[d].DownstreamLI].SCI[j  ]].FI[se], -1);
			}
			    size_t j=Bk[k].Ny;
			    V[Bk[k]. VI[i-1][j]].FB=CBC_c(PeriodicUpstream, 1, 0, U[Ln[Bk[k].B[d].DownstreamLI].SCI[j-1]].FI[ne], -1);
		} else if(Bk[k].B[d].T==VelocityInlet) {
			for(size_t j=0; j< Bk[k].Ny; ++j) {
				U[Bk[k]. UI[i  ][j]].FB=CBC_c(VelocityInlet, 1, Bk[k].B[d].u);
				U[Bk[k]. UI[i  ][j]].NormalA=-U[Bk[k].UI[i  ][j]].Ax; uvCell_c::InletUI.push_back(Bk[k].UI[i  ][j]);
//						C[Bk[k]. CI[i-1][j]].FB=CBC_c(VelocityInlet, 1, 0, Bk[k].CI[i-2][j], -1);
				C[Bk[k]. CI[i-1][j]].EB=CBC_c(VelocityInlet, 1, Bk[k].B[d].c);
			}
			for(size_t j=0; j<=Bk[k].Ny; ++j) {
				V[Bk[k]. VI[i-1][j]].FB=CBC_c(VelocityInlet, 1, Bk[k].B[d].v);
			}
		} else if(Bk[k].B[d].T==PressureOutlet) {
			for(size_t j=0; j< Bk[k].Ny; ++j) {
				U[Bk[k]. UI[i  ][j]].FB=CBC_c(PressureOutlet, 1, 0, Bk[k].UI[i-1][j], -uvCell_c::MassRatio);
				U[Bk[k]. UI[i  ][j]].NormalA=U[Bk[k].UI[i  ][j]].Ax; uvCell_c::OutletUI.push_back(Bk[k].UI[i  ][j]);
				C[Bk[k]. CI[i-1][j]].FB=CBC_c(PressureOutlet, 1, Bk[k].B[d].p); //
				C[Bk[k]. CI[i-1][j]].EB=CBC_c(PressureOutlet, 1, 0, Bk[k].CI[i-2][j], -1);
			}
			for(size_t j=0; j<=Bk[k].Ny; ++j) {
				V[Bk[k]. VI[i-1][j]].FB=CBC_c(PressureOutlet, 1, 0, Bk[k].VI[i-2][j], -1);
			}
		} else {
			//Do nothing when the boundary condition is NoBC; Wall or PorousOpen will be treated next
		}

		i=0; d=w;
		if(Bk[k].B[d].T==Inlet) {
			for(size_t j=0; j< Bk[k].Ny; ++j) {
				U[Bk[k]. UI[i][j]].FB=CBC_c(Inlet, 1, Bk[k].B[d].u);
				U[Bk[k]. UI[i][j]].NormalA=U[Bk[k].UI[i][j]].Ax; uvCell_c::InletUI.push_back(Bk[k].UI[i][j]);
				C[Bk[k]. CI[i][j]].FB=CBC_c(Inlet, 1, Bk[k].B[d].p);
				C[Bk[k]. CI[i][j]].EB=CBC_c(Inlet, 1, Bk[k].B[d].c);
			}
			for(size_t j=0; j<=Bk[k].Ny; ++j) {
				V[Bk[k]. VI[i][j]].FB=CBC_c(Inlet, 1, Bk[k].B[d].v);
			}
		} else if(Bk[k].B[d].T==Outlet) {
			for(size_t j=0; j< Bk[k].Ny; ++j) {
				U[Bk[k]. UI[i][j]].FB=CBC_c(Outlet, 1, 0, Bk[k].UI[i+1][j], -uvCell_c::MassRatio);
				U[Bk[k]. UI[i][j]].NormalA=-U[Bk[k].UI[i][j]].Ax; uvCell_c::OutletUI.push_back(Bk[k].UI[i][j]);
//						C[Bk[k]. CI[i][j]].FB=CBC_c(Outlet, 1, 0, Bk[k].CI[i+1][j], -1);
				C[Bk[k]. CI[i][j]].EB=CBC_c(Outlet, 1, 0, Bk[k].CI[i+1][j], -1);
			}
			for(size_t j=0; j<=Bk[k].Ny; ++j) {
				V[Bk[k]. VI[i][j]].FB=CBC_c(Outlet, 1, 0, Bk[k].VI[i+1][j], -1);
			}
		} else if(Bk[k].B[d].T==PeriodicDownstream) {
			for(size_t j=0; j< Bk[k].Ny; ++j) {
				U[Bk[k]. UI[i][j]].FB=CBC_c(PeriodicDownstream, 1, 0, Bk[k].UI[i+1][j], -uvCell_c::MassRatio);
				U[Bk[k]. UI[i][j]].NormalA=-U[Bk[k].UI[i][j]].Ax; uvCell_c::OutletUI.push_back(Bk[k].UI[i][j]);
				C[Bk[k]. CI[i][j]].FB=CBC_c(PeriodicDownstream, 1, Bk[k].B[d].p+cCell_c::dpdd*C[Bk[k].CI[i][j]].d[w]);
				C[Bk[k]. CI[i][j]].EB=CBC_c(PeriodicDownstream, 1, 0, Bk[k].CI[i+1][j], -1);
			}
			for(size_t j=0; j<=Bk[k].Ny; ++j) {
				V[Bk[k]. VI[i][j]].FB=CBC_c(PeriodicDownstream, 1, 0, Bk[k].VI[i+1][j], -1);
			}
		} else if(Bk[k].B[d].T==PeriodicUpstream) {
			for(size_t j=0; j<Bk[k].Ny; ++j) {
				U[Bk[k]. UI[i][j]].FB=CBC_c(PeriodicUpstream, 1, 0, Ln[Bk[k].B[d].DownstreamLI].SCI[j], -1);
				U[Bk[k]. UI[i][j]].NormalA=U[Bk[k].UI[i][j]].Ax; uvCell_c::InletUI.push_back(Bk[k].UI[i][j]);
				C[Bk[k]. CI[i][j]].FB=CBC_c(PeriodicUpstream, 1, Bk[k].B[d].p-cCell_c::dpdd*C[Bk[k].CI[i][j]].d[w]);
				C[Bk[k]. CI[i][j]].EB=CBC_c(PeriodicUpstream, 1, Bk[k].B[d].c);
				V[Bk[k]. VI[i][j]].FB=CBC_c(PeriodicUpstream, 1, 0, U[Ln[Bk[k].B[d].DownstreamLI].SCI[j  ]].FI[sw], -1);
			}
			    size_t j=Bk[k].Ny;
			    V[Bk[k]. VI[i][j]].FB=CBC_c(PeriodicUpstream, 1, 0, U[Ln[Bk[k].B[d].DownstreamLI].SCI[j-1]].FI[nw], -1);
		} else if(Bk[k].B[d].T==VelocityInlet) {
			for(size_t j=0; j< Bk[k].Ny; ++j) {
				U[Bk[k]. UI[i][j]].FB=CBC_c(VelocityInlet, 1, Bk[k].B[d].u);
				U[Bk[k]. UI[i][j]].NormalA=U[Bk[k].UI[i][j]].Ax; uvCell_c::InletUI.push_back(Bk[k].UI[i][j]);
//						C[Bk[k]. CI[i][j]].FB=CBC_c(VelocityInlet, 1, 0, Bk[k].CI[i+1][j], -1);
				C[Bk[k]. CI[i][j]].EB=CBC_c(VelocityInlet, 1, Bk[k].B[d].c);
			}
			for(size_t j=0; j<=Bk[k].Ny; ++j) {
				V[Bk[k]. VI[i][j]].FB=CBC_c(VelocityInlet, 1, Bk[k].B[d].v);
			}
		} else if(Bk[k].B[d].T==PressureOutlet) {
			for(size_t j=0; j< Bk[k].Ny; ++j) {
				U[Bk[k]. UI[i][j]].FB=CBC_c(PressureOutlet, 1, 0, Bk[k].UI[i+1][j], -uvCell_c::MassRatio);
				U[Bk[k]. UI[i][j]].NormalA=-U[Bk[k].UI[i][j]].Ax; uvCell_c::OutletUI.push_back(Bk[k].UI[i][j]);
				C[Bk[k]. CI[i][j]].FB=CBC_c(PressureOutlet, 1, Bk[k].B[d].p);
				C[Bk[k]. CI[i][j]].EB=CBC_c(PressureOutlet, 1, 0, Bk[k].CI[i+1][j], -1);
			}
			for(size_t j=0; j<=Bk[k].Ny; ++j) {
				V[Bk[k]. VI[i][j]].FB=CBC_c(PressureOutlet, 1, 0, Bk[k].VI[i+1][j], -1);
			}
		} else {
			//Do nothing when the boundary condition is NoBC; Wall or PorousOpen will be treated next
		}

		j=Bk[k].Ny; d=n;
		if(Bk[k].B[d].T==Inlet) {
			for(size_t i=0; i< Bk[k].Nx; ++i) {
				V[Bk[k]. VI[i][j  ]].FB=CBC_c(Inlet, 1, Bk[k].B[d].v);
				V[Bk[k]. VI[i][j  ]].NormalA=-V[Bk[k].VI[i][j  ]].Ay; uvCell_c::InletVI.push_back(Bk[k].VI[i][j  ]);
				C[Bk[k]. CI[i][j-1]].FB=CBC_c(Inlet, 1, Bk[k].B[d].p);
				C[Bk[k]. CI[i][j-1]].EB=CBC_c(Inlet, 1, Bk[k].B[d].c);
			}
			for(size_t i=0; i<=Bk[k].Nx; ++i) {
				U[Bk[k]. UI[i][j]].FB=CBC_c(Inlet, 1, Bk[k].B[d].u);
			}
		} else if(Bk[k].B[d].T==Outlet) {
			for(size_t i=0; i< Bk[k].Nx; ++i) {
				V[Bk[k]. VI[i][j  ]].FB=CBC_c(Outlet, 1, 0, Bk[k].VI[i][j-1], -uvCell_c::MassRatio);
				V[Bk[k]. VI[i][j  ]].NormalA=V[Bk[k].VI[i][j  ]].Ay; uvCell_c::OutletVI.push_back(Bk[k].VI[i][j  ]);
//						C[Bk[k]. CI[i][j-1]].FB=CBC_c(Outlet, 1, 0, Bk[k].CI[i][j-2], -1);
				C[Bk[k]. CI[i][j-1]].EB=CBC_c(Outlet, 1, 0, Bk[k].CI[i][j-2], -1);
			}
			for(size_t i=0; i<=Bk[k].Nx; ++i) {
				U[Bk[k]. UI[i][j-1]].FB=CBC_c(Outlet, 1, 0, Bk[k].UI[i][j-2], -1);
			}
		} else if(Bk[k].B[d].T==PeriodicDownstream) {
			for(size_t i=0; i< Bk[k].Nx; ++i) {
				V[Bk[k]. VI[i][j  ]].FB=CBC_c(PeriodicDownstream, 1, 0, Bk[k].VI[i][j-1], -uvCell_c::MassRatio);
				V[Bk[k]. VI[i][j  ]].NormalA=V[Bk[k].VI[i][j  ]].Ay; uvCell_c::OutletVI.push_back(Bk[k].VI[i][j  ]);
				C[Bk[k]. CI[i][j-1]].FB=CBC_c(PeriodicDownstream, 1, Bk[k].B[d].p+cCell_c::dpdd*C[Bk[k].CI[i][j-1]].d[n]);
				C[Bk[k]. CI[i][j-1]].EB=CBC_c(PeriodicDownstream, 1, 0, Bk[k].CI[i][j-2], -1);
			}
			for(size_t i=0; i<=Bk[k].Nx; ++i) {
				U[Bk[k]. UI[i][j-1]].FB=CBC_c(PeriodicDownstream, 1, 0, Bk[k].UI[i][j-2], -1);
			}
		} else if(Bk[k].B[d].T==PeriodicUpstream) {
			for(size_t i=0; i< Bk[k].Nx; ++i) {
				V[Bk[k]. VI[i][j  ]].FB=CBC_c(PeriodicUpstream, 1, 0, Ln[Bk[k].B[d].DownstreamLI].SCI[i], -1);
				V[Bk[k]. VI[i][j  ]].NormalA=-V[Bk[k].VI[i][j  ]].Ay; uvCell_c::InletVI.push_back(Bk[k].VI[i][j  ]);
				C[Bk[k]. CI[i][j-1]].FB=CBC_c(PeriodicUpstream, 1, Bk[k].B[d].p-cCell_c::dpdd*C[Bk[k].CI[i][j-1]].d[n]);
				C[Bk[k]. CI[i][j-1]].EB=CBC_c(PeriodicUpstream, 1, Bk[k].B[d].c);
				U[Bk[k]. UI[i][j-1]].FB=CBC_c(PeriodicUpstream, 1, 0, V[Ln[Bk[k].B[d].DownstreamLI].SCI[i  ]].FI[nw], -1);
			}
			    size_t i=Bk[k].Nx;
				U[Bk[k]. UI[i][j-1]].FB=CBC_c(PeriodicUpstream, 1, 0, V[Ln[Bk[k].B[d].DownstreamLI].SCI[i-1]].FI[ne], -1);
		} else if(Bk[k].B[d].T==VelocityInlet) {
			for(size_t i=0; i< Bk[k].Nx; ++i) {
				V[Bk[k]. VI[i][j  ]].FB=CBC_c(VelocityInlet, 1, Bk[k].B[d].v);
				V[Bk[k]. VI[i][j  ]].NormalA=-V[Bk[k].VI[i][j  ]].Ay; uvCell_c::InletVI.push_back(Bk[k].VI[i][j  ]);
//						C[Bk[k]. CI[i][j-1]].FB=CBC_c(VelocityInlet, 1, 0, Bk[k].CI[i][j-1], -1);
				C[Bk[k]. CI[i][j-1]].EB=CBC_c(VelocityInlet, 1, Bk[k].B[d].c);
			}
			for(size_t i=0; i<=Bk[k].Nx; ++i) {
				U[Bk[k]. UI[i][j-1]].FB=CBC_c(VelocityInlet, 1, Bk[k].B[d].u);
			}
		} else if(Bk[k].B[d].T==PressureOutlet) {
			for(size_t i=0; i< Bk[k].Nx; ++i) {
				V[Bk[k]. VI[i][j  ]].FB=CBC_c(PressureOutlet, 1, 0, Bk[k].VI[i][j-1], -uvCell_c::MassRatio);
				V[Bk[k]. VI[i][j  ]].NormalA=V[Bk[k].VI[i][j  ]].Ay; uvCell_c::OutletVI.push_back(Bk[k].VI[i][j  ]);
				C[Bk[k]. CI[i][j-1]].FB=CBC_c(PressureOutlet, 1, Bk[k].B[d].p);
				C[Bk[k]. CI[i][j-1]].EB=CBC_c(PressureOutlet, 1, 0, Bk[k].CI[i][j-2], -1);
			}
			for(size_t i=0; i<=Bk[k].Nx; ++i) {
				U[Bk[k]. UI[i][j-1]].FB=CBC_c(PressureOutlet, 1, 0, Bk[k].UI[i][j-2], -1);
			}
		} else {
			//Do nothing when the boundary condition is NoBC; Wall or PorousOpen will be treated next
		}

		j=0; d=s;
		if(Bk[k].B[d].T==Inlet) {
			for(size_t i=0; i< Bk[k].Nx; ++i) {
				V[Bk[k]. VI[i][j]].FB=CBC_c(Inlet, 1, Bk[k].B[d].v);
				V[Bk[k]. VI[i][j]].NormalA=V[Bk[k].VI[i][j]].Ay; uvCell_c::InletVI.push_back(Bk[k].VI[i][j]);
				C[Bk[k]. CI[i][j]].FB=CBC_c(Inlet, 1, Bk[k].B[d].p);
				C[Bk[k]. CI[i][j]].EB=CBC_c(Inlet, 1, Bk[k].B[d].c);
			}
			for(size_t i=0; i<=Bk[k].Nx; ++i) {
				U[Bk[k]. UI[i][j]].FB=CBC_c(Inlet, 1, Bk[k].B[d].u);
			}
		} else if(Bk[k].B[d].T==Outlet) {
			for(size_t i=0; i< Bk[k].Nx; ++i) {
				V[Bk[k]. VI[i][j]].FB=CBC_c(Outlet, 1, 0, Bk[k].VI[i][j+1], -uvCell_c::MassRatio);
				V[Bk[k]. VI[i][j]].NormalA=-V[Bk[k].VI[i][j]].Ay; uvCell_c::OutletVI.push_back(Bk[k].VI[i][j]);
//						C[Bk[k]. CI[i][j]].FB=CBC_c(Outlet, 1, 0, Bk[k].CI[i][j+1], -1);
				C[Bk[k]. CI[i][j]].EB=CBC_c(Outlet, 1, 0, Bk[k].CI[i][j+1], -1);
			}
			for(size_t i=0; i<=Bk[k].Nx; ++i) {
				U[Bk[k]. UI[i][j]].FB=CBC_c(Outlet, 1, 0, Bk[k].UI[i][j+1], -1);
			}
		} else if(Bk[k].B[d].T==PeriodicDownstream) {
			for(size_t i=0; i< Bk[k].Nx; ++i) {
				V[Bk[k]. VI[i][j]].FB=CBC_c(PeriodicDownstream, 1, 0, Bk[k].UI[i][j+1], -uvCell_c::MassRatio);
				V[Bk[k]. VI[i][j]].NormalA=-V[Bk[k].VI[i][j]].Ay; uvCell_c::OutletVI.push_back(Bk[k].VI[i][j]);
				C[Bk[k]. CI[i][j]].FB=CBC_c(PeriodicDownstream, 1, Bk[k].B[d].p+cCell_c::dpdd*C[Bk[k].CI[i][j]].d[s]);
				C[Bk[k]. CI[i][j]].EB=CBC_c(PeriodicDownstream, 1, 0, Bk[k].CI[i][j+1], -1);
			}
			for(size_t i=0; i<=Bk[k].Nx; ++i) {
				U[Bk[k]. UI[i][j]].FB=CBC_c(PeriodicDownstream, 1, 0, Bk[k].UI[i][j+1], -1);
			}
		} else if(Bk[k].B[d].T==PeriodicUpstream) {
			for(size_t i=0; i< Bk[k].Nx; ++i) {
				V[Bk[k]. VI[i][j]].FB=CBC_c(PeriodicUpstream, 1, 0, Ln[Bk[k].B[d].DownstreamLI].SCI[j], -1);
				V[Bk[k]. VI[i][j]].NormalA=V[Bk[k].VI[i][j]].Ay; uvCell_c::InletVI.push_back(Bk[k].VI[i][j]);
				C[Bk[k]. CI[i][j]].FB=CBC_c(PeriodicUpstream, 1, Bk[k].B[d].p-cCell_c::dpdd*C[Bk[k].CI[i][j]].d[s]);
				C[Bk[k]. CI[i][j]].EB=CBC_c(PeriodicUpstream, 1, Bk[k].B[d].c);
				U[Bk[k]. UI[i][j]].FB=CBC_c(PeriodicUpstream, 1, 0, V[Ln[Bk[k].B[d].DownstreamLI].SCI[j]].FI[sw], -1);
			}
			    size_t i=Bk[k].Nx;
				U[Bk[k]. UI[i][j]].FB=CBC_c(PeriodicUpstream, 1, 0, V[Ln[Bk[k].B[d].DownstreamLI].SCI[j]].FI[se], -1);
		} else if(Bk[k].B[d].T==VelocityInlet) {
			for(size_t i=0; i< Bk[k].Nx; ++i) {
				V[Bk[k]. VI[i][j]].FB=CBC_c(VelocityInlet, 1, Bk[k].B[d].v);
				V[Bk[k]. VI[i][j]].NormalA=V[Bk[k].VI[i][j]].Ay; uvCell_c::InletVI.push_back(Bk[k].VI[i][j]);
//						C[Bk[k]. CI[i][j]].FB=CBC_c(VelocityInlet, 1, 0, Bk[k].CI[i][j+1], -1);
				C[Bk[k]. CI[i][j]].EB=CBC_c(VelocityInlet, 1, Bk[k].B[d].c);
			}
			for(size_t i=0; i<=Bk[k].Nx; ++i) {
				U[Bk[k]. UI[i][j]].FB=CBC_c(VelocityInlet, 1, Bk[k].B[d].u);
			}
		} else if(Bk[k].B[d].T==PressureOutlet) {
			for(size_t i=0; i< Bk[k].Nx; ++i) {
				V[Bk[k]. VI[i][j]].FB=CBC_c(PressureOutlet, 1, 0, Bk[k].VI[i][j+1], -uvCell_c::MassRatio);
				V[Bk[k]. VI[i][j]].NormalA=-V[Bk[k].VI[i][j]].Ay; uvCell_c::OutletVI.push_back(Bk[k].VI[i][j]);
				C[Bk[k]. CI[i][j]].FB=CBC_c(PressureOutlet, 1, Bk[k].B[d].p);
				C[Bk[k]. CI[i][j]].EB=CBC_c(PressureOutlet, 1, 0, Bk[k].CI[i][j+1], -1);
			}
			for(size_t i=0; i<=Bk[k].Nx; ++i) {
				U[Bk[k]. UI[i][j]].FB=CBC_c(PressureOutlet, 1, 0, Bk[k].UI[i][j+1], -1);
			}
		} else {
			//Do nothing when the boundary condition is NoBC; Wall or PorousOpen will be treated next
		}
	}

	cout<<"Processing Wall and PorousOpen BCs......"<<endl;
	for(size_t k=porousBlocks; k<Bk.size(); ++k) {
		size_t i(0), j(0), d(e);

		i=Bk[k].Nx; d=e;
		if(Bk[k].B[d].T==Wall) {
			for(size_t j=0; j< Bk[k].Ny; ++j) {
				U[Bk[k].UI[i  ][j]        ].FB=CBC_c(Wall, 1, Bk[k].B[d].u);
			}
		} else if(Bk[k].B[d].T==PorousOpen) {
			for(size_t j=0; j< Bk[k].Ny; ++j) {
				U[Bk[k].UI[i  ][j]        ].FB=CBC_c(Ln[Bk[k].LnI[d]].B[j].T, 1, 0);
				if(Ln[Bk[k].LnI[d]].B[j].PI!=voidIndex) {
					C[Bk[k].CI[i-1][j]        ].PI.push_back(  Ln[Bk[k].LnI[d]].B  [j].PI    );
					P[Ln[Bk[k].LnI[d]].B[j].PI].CI.push_back(U[Ln[Bk[k].LnI[d]].SCI[j]].FI[w]);
//					P[Ln[Bk[k].LnI[d]].B[j].PI].CD.push_back(C[U[Ln[Bk[k].LnI[d]].SCI[j]].FI[w]].d[e]);
//					P[Ln[Bk[k].LnI[d]].B[j].PI].CD.push_back(C[U[Ln[Bk[k].LnI[d]].SCI[j]].FI[w]].d[e]*2);
					P[Ln[Bk[k].LnI[d]].B[j].PI].CD.push_back(cd_c::dz/8);
//					P[Ln[Bk[k].LnI[d]].B[j].PI].CD.push_back(cd_c::dz/16);
					if(Bk[Bk[k].NI[d]].NinDiameter!=0 && Bk[Bk[k].NI[d]].NinLength==0) {
						P[Ln[Bk[k].LnI[d]].B[j].PI].CA.push_back(T[P[Ln[Bk[k].LnI[d]].B[j].PI].TI[0]].A/Bk[Bk[k].NI[d]].NinDiameter);
					} else if(Bk[Bk[k].NI[d]].NinDiameter==0 && Bk[Bk[k].NI[d]].NinLength!=0) {
						P[Ln[Bk[k].LnI[d]].B[j].PI].CA.push_back(T[P[Ln[Bk[k].LnI[d]].B[j].PI].TI[0]].A);
					} else {
						//do nothing
					}
				}
			}
		} else {
			//Do nothing when the boundary condition is NoBC;
		}

		i=0; d=w;
		if(Bk[k].B[d].T==Wall) {
			for(size_t j=0; j< Bk[k].Ny; ++j) {
				U[Bk[k].UI[i][j]          ].FB=CBC_c(Wall, 1, Bk[k].B[d].u);
			}
		} else if(Bk[k].B[d].T==PorousOpen) {
			for(size_t j=0; j<Bk[k].Ny; ++j) {
				U[Bk[k].UI[i][j]          ].FB=CBC_c(Ln[Bk[k].LnI[d]].B[j].T, 1, 0);
				if(Ln[Bk[k].LnI[d]].B[j].PI!=voidIndex) {
					C[Bk[k].CI[i][j]          ].PI.push_back(  Ln[Bk[k].LnI[d]].B  [j].PI    );
					P[Ln[Bk[k].LnI[d]].B[j].PI].CI.push_back(U[Ln[Bk[k].LnI[d]].SCI[j]].FI[e]);
//					P[Ln[Bk[k].LnI[d]].B[j].PI].CD.push_back(C[U[Ln[Bk[k].LnI[d]].SCI[j]].FI[e]].d[w]);
//					P[Ln[Bk[k].LnI[d]].B[j].PI].CD.push_back(C[U[Ln[Bk[k].LnI[d]].SCI[j]].FI[e]].d[w]*2);
					P[Ln[Bk[k].LnI[d]].B[j].PI].CD.push_back(cd_c::dz/8);
//					P[Ln[Bk[k].LnI[d]].B[j].PI].CD.push_back(cd_c::dz/16);
					if(Bk[Bk[k].NI[d]].NinDiameter!=0 && Bk[Bk[k].NI[d]].NinLength==0) {
						P[Ln[Bk[k].LnI[d]].B[j].PI].CA.push_back(T[P[Ln[Bk[k].LnI[d]].B[j].PI].TI[0]].A/Bk[Bk[k].NI[d]].NinDiameter);
					} else if(Bk[Bk[k].NI[d]].NinDiameter==0 && Bk[Bk[k].NI[d]].NinLength!=0) {
						P[Ln[Bk[k].LnI[d]].B[j].PI].CA.push_back(T[P[Ln[Bk[k].LnI[d]].B[j].PI].TI[0]].A);
					} else {
						//do nothing
					}
				}
			}
		} else {
			//Do nothing when the boundary condition is NoBC; Wall or PorousOpen will be treated next
		}

		j=Bk[k].Ny; d=n;
		if(Bk[k].B[d].T==Wall) {
			for(size_t i=0; i< Bk[k].Nx; ++i) {
				V[Bk[k].VI[i][j  ]        ].FB=CBC_c(Wall, 1, Bk[k].B[d].v);
			}
		} else if(Bk[k].B[d].T==PorousOpen) {
			for(size_t i=0; i< Bk[k].Nx; ++i) {
				V[Bk[k].VI[i][j  ]        ].FB=CBC_c(Ln[Bk[k].LnI[d]].B[i].T, 1, 0);
				if(Ln[Bk[k].LnI[d]].B[i].PI!=voidIndex) {
					C[Bk[k].CI[i][j-1]        ].PI.push_back(  Ln[Bk[k].LnI[d]].B  [i].PI    );
					P[Ln[Bk[k].LnI[d]].B[i].PI].CI.push_back(V[Ln[Bk[k].LnI[d]].SCI[i]].FI[s]);
//					P[Ln[Bk[k].LnI[d]].B[i].PI].CD.push_back(C[V[Ln[Bk[k].LnI[d]].SCI[i]].FI[s]].d[n]);
//					P[Ln[Bk[k].LnI[d]].B[i].PI].CD.push_back(C[V[Ln[Bk[k].LnI[d]].SCI[i]].FI[s]].d[n]*2);
					P[Ln[Bk[k].LnI[d]].B[i].PI].CD.push_back(cd_c::dz/8);
//					P[Ln[Bk[k].LnI[d]].B[i].PI].CD.push_back(cd_c::dz/16);
					if(Bk[Bk[k].NI[d]].NinDiameter!=0 && Bk[Bk[k].NI[d]].NinLength==0) {
						P[Ln[Bk[k].LnI[d]].B[i].PI].CA.push_back(T[P[Ln[Bk[k].LnI[d]].B[i].PI].TI[0]].A/Bk[Bk[k].NI[d]].NinDiameter);
					} else if(Bk[Bk[k].NI[d]].NinDiameter==0 && Bk[Bk[k].NI[d]].NinLength!=0) {
						P[Ln[Bk[k].LnI[d]].B[i].PI].CA.push_back(T[P[Ln[Bk[k].LnI[d]].B[i].PI].TI[0]].A);
					} else {
						//do nothing
					}
				}
			}
		} else {
			//Do nothing when the boundary condition is NoBC; Wall or PorousOpen will be treated next
		}

		j=0; d=s;
		if(Bk[k].B[d].T==Wall) {
			for(size_t i=0; i< Bk[k].Nx; ++i) {
				V[Bk[k].VI[i][j]          ].FB=CBC_c(Wall, 1, Bk[k].B[d].v);
			}
		} else if(Bk[k].B[d].T==PorousOpen) {
			for(size_t i=0; i< Bk[k].Nx; ++i) {
				V[Bk[k].VI[i][j]          ].FB=CBC_c(Ln[Bk[k].LnI[d]].B[i].T, 1, 0);
				if(Ln[Bk[k].LnI[d]].B[i].PI!=voidIndex) {
					C[Bk[k].CI[i][j]          ].PI.push_back(  Ln[Bk[k].LnI[d]].B  [i].PI    );
					P[Ln[Bk[k].LnI[d]].B[i].PI].CI.push_back(V[Ln[Bk[k].LnI[d]].SCI[i]].FI[n]);
//					P[Ln[Bk[k].LnI[d]].B[i].PI].CD.push_back(C[V[Ln[Bk[k].LnI[d]].SCI[i]].FI[n]].d[s]);
//					P[Ln[Bk[k].LnI[d]].B[i].PI].CD.push_back(C[V[Ln[Bk[k].LnI[d]].SCI[i]].FI[n]].d[s]*2);
					P[Ln[Bk[k].LnI[d]].B[i].PI].CD.push_back(cd_c::dz/8);
//					P[Ln[Bk[k].LnI[d]].B[i].PI].CD.push_back(cd_c::dz/16);
					if(Bk[Bk[k].NI[d]].NinDiameter!=0 && Bk[Bk[k].NI[d]].NinLength==0) {
						P[Ln[Bk[k].LnI[d]].B[i].PI].CA.push_back(T[P[Ln[Bk[k].LnI[d]].B[i].PI].TI[0]].A/Bk[Bk[k].NI[d]].NinDiameter);
					} else if(Bk[Bk[k].NI[d]].NinDiameter==0 && Bk[Bk[k].NI[d]].NinLength!=0) {
						P[Ln[Bk[k].LnI[d]].B[i].PI].CA.push_back(T[P[Ln[Bk[k].LnI[d]].B[i].PI].TI[0]].A);
					} else {
						//do nothing
					}
				}
			}
		} else {
			//Do nothing when the boundary condition is NoBC; Wall or PorousOpen will be treated next
		}
	}//finish block loop

	cout<<"Processing Constant Concentration BCs......"<<endl;
	for(size_t k=porousBlocks; k<Bk.size(); ++k) {
		size_t i(0), j(0), d(e);

		i=Bk[k].Nx-1; d=e;
		if(Bk[k].B[d].T==ConstConcentration) {
			for(size_t j=0; j< Bk[k].Ny; ++j) {
				C[Bk[k].CI[i][j]].EB=CBC_c(ConstConcentration, 1, Bk[k].B[d].c);
			}
		}

		i=0; d=w;
		if(Bk[k].B[d].T==ConstConcentration) {
			for(size_t j=0; j< Bk[k].Ny; ++j) {
				C[Bk[k].CI[i][j]].EB=CBC_c(ConstConcentration, 1, Bk[k].B[d].c);
			}
		}

		j=Bk[k].Ny-1; d=n;
		if(Bk[k].B[d].T==ConstConcentration) {
			for(size_t i=0; i< Bk[k].Nx; ++i) {
				C[Bk[k].CI[i][j]].EB=CBC_c(ConstConcentration, 1, Bk[k].B[d].c);
			}
		}

		j=0; d=s;
		if(Bk[k].B[d].T==ConstConcentration) {
			for(size_t i=0; i< Bk[k].Nx; ++i) {
				C[Bk[k].CI[i][j]].EB=CBC_c(ConstConcentration, 1, Bk[k].B[d].c);
			}
		}
	}//finish block loop

	for(size_t i=0; i<C.size(); ++i) {
		C[i].ap.resize(C[i].PI.size(), 0);
		for(size_t j=0; j<C[i].PI.size(); ++j) {
			C[i].TI.push_back(P[C[i].PI[j]].TI[0]);
		}
	}

	for(size_t i=0; i<P.size(); ++i) {
		for(size_t j=0; j<P[i].TI.size(); ++j) {
			P[i].PI.push_back(T[P[i].TI[j]].otherP(i));
		}
		P[i].a .resize(P[i].CI.size(), 0);
		P[i].ap.resize(P[i].PI.size(), 0);
	}
	cout<<"Find out cells near Boundary Conditions......"<<endl;
	for(size_t k=porousBlocks; k<Bk.size(); ++k) {
		if(Bk[k].B[e].T!=NoBC && Bk[k].Nx>=2) {
			size_t i=Bk[k].Nx;
			for(size_t j=0; j< Bk[k].Ny; ++j) {
				U[Bk[k].UI[i-1][j]].nearBC=true;
				C[Bk[k].CI[i-2][j]].nearBC=true;
			}
			for(size_t j=0; j<=Bk[k].Ny; ++j) {
				V[Bk[k].VI[i-2][j]].nearBC=true;
			}
		}
		if(Bk[k].B[w].T!=NoBC && Bk[k].Nx>=2) {
			size_t i=0;
			for(size_t j=0; j< Bk[k].Ny; ++j) {
				U[Bk[k].UI[i+1][j]].nearBC=true;
				C[Bk[k].CI[i+2][j]].nearBC=true;
			}
			for(size_t j=0; j<=Bk[k].Ny; ++j) {
				V[Bk[k].VI[i+2][j]].nearBC=true;
			}
		}
		if(Bk[k].B[n].T!=NoBC && Bk[k].Ny>=2) {
			size_t j=Bk[k].Ny;
			for(size_t i=0; i< Bk[k].Nx; ++i) {
				V[Bk[k].VI[i][j-1]].nearBC=true;
				C[Bk[k].CI[i][j-2]].nearBC=true;
			}
			for(size_t i=0; i<=Bk[k].Nx; ++i) {
				U[Bk[k].UI[i][j-2]].nearBC=true;
			}
		}
		if(Bk[k].B[n].T!=NoBC && Bk[k].Ny>=2) {
			size_t j=0;
			for(size_t i=0; i< Bk[k].Nx; ++i) {
				V[Bk[k].VI[i][j+1]].nearBC=true;
				C[Bk[k].CI[i][j+2]].nearBC=true;
			}
			for(size_t i=0; i<=Bk[k].Nx; ++i) {
				U[Bk[k].UI[i][j+2]].nearBC=true;
			}
		}
	}

	cout<<"Set isWall for all cells......"<<endl;
	for(size_t k=porousBlocks; k<Bk.size(); ++k) {
		if(Bk[k].B[e].T==Wall || Bk[k].B[e].T==PorousOpen) {
			size_t i=Bk[k].Nx;
			for(size_t j=0; j< Bk[k].Ny; ++j) {
				C[Bk[k].CI[i-1][j]].isWall.set(e);
			}
			for(size_t j=0; j<=Bk[k].Ny; ++j) {
				V[Bk[k].VI[i-1][j]].isWall.set(e);
			}
		}
		if(Bk[k].B[w].T==Wall || Bk[k].B[w].T==PorousOpen) {
			size_t i=0;
			for(size_t j=0; j< Bk[k].Ny; ++j) {
				C[Bk[k].CI[i  ][j]].isWall.set(w);
			}
			for(size_t j=0; j<=Bk[k].Ny; ++j) {
				V[Bk[k].VI[i  ][j]].isWall.set(w);
			}
		}
		if(Bk[k].B[n].T==Wall || Bk[k].B[n].T==PorousOpen) {
			size_t j=Bk[k].Ny;
			for(size_t i=0; i< Bk[k].Nx; ++i) {
				C[Bk[k].CI[i][j-1]].isWall.set(n);
			}
			for(size_t i=0; i<=Bk[k].Nx; ++i) {
				U[Bk[k].UI[i][j-1]].isWall.set(n);
			}
		}
		if(Bk[k].B[s].T==Wall || Bk[k].B[s].T==PorousOpen) {
			size_t j=0;
			for(size_t i=0; i< Bk[k].Nx; ++i) {
				C[Bk[k].CI[i][j  ]].isWall.set(s);
			}
			for(size_t i=0; i<=Bk[k].Nx; ++i) {
				U[Bk[k].UI[i][j  ]].isWall.set(s);
			}
		}
	}
	oflg<<"Cells side is wall: "<<endl;
	for(size_t i=0; i<U.size(); ++i) {
		oflg<<"U # "<<i;
		for(size_t d=e; d<=s; ++d) {
			oflg<<"\t"<<U[i].isWall[d];
		}
		oflg<<endl;
	}
	for(size_t i=0; i<V.size(); ++i) {
		oflg<<"V # "<<i;
		for(size_t d=e; d<=s; ++d) {
			oflg<<"\t"<<V[i].isWall[d];
		}
		oflg<<endl;
	}
	for(size_t i=0; i<C.size(); ++i) {
		oflg<<"C # "<<i;
		for(size_t d=e; d<=s; ++d) {
			oflg<<"\t"<<C[i].isWall[d];
		}
		oflg<<endl;
	}
	oflg<<"cCell, Checking relationships:";
	for(size_t i=0; i<C.size(); ++i) {
		oflg<<endl<<"C#"<<i<<":\t";
		for(size_t n=0; n<neighborNumber; ++n) {
			oflg<<C[i].NI[n]<<"\t";
		}
		oflg<<"|";
		for(size_t n=0; n<nodeNumber; ++n) {
			oflg<<"\t"<<C[i].FI[n];
		}
		oflg<<"\t|";
		for(size_t n=0; n<C[i].PI.size(); ++n) {
			oflg<<"\tP"<<C[i].PI[n];
		}
		oflg<<"\t|";
		for(size_t n=0; n<C[i].TI.size(); ++n) {
			oflg<<"\tT"<<C[i].TI[n];
		}
	}

	oflg<<endl<<"uCell, Checking relationships:";
	for(size_t i=0; i<U.size(); ++i) {
		oflg<<endl<<"U#"<<i<<":\t";
		for(size_t n=0; n<neighborNumber; ++n) {
			oflg<<U[i].NI[n]<<"\t";
		}
		oflg<<"|";
		for(size_t n=0; n<nodeNumber; ++n) {
			oflg<<"\t"<<U[i].FI[n];
		}
	}

	oflg<<endl<<"vCell, Checking relationships:";
	for(size_t i=0; i<V.size(); ++i) {
		oflg<<endl<<"V#"<<i<<":\t";
		for(size_t n=0; n<neighborNumber; ++n) {
			oflg<<V[i].NI[n]<<"\t";
		}
		oflg<<"|";
		for(size_t n=0; n<nodeNumber; ++n) {
			oflg<<"\t"<<V[i].FI[n];
		}
	}

	oflg<<endl<<"fCell, Checking relationships:";
	for(size_t i=0; i<F.size(); ++i) {
		oflg<<endl<<"F#"<<i<<":\t";
		for(size_t n=0; n<neighborNumber; ++n) {
			oflg<<F[i].NI[n]<<"\t";
		}
		oflg<<"|";
		for(size_t n=0; n<nodeNumber; ++n) {
			oflg<<"\t"<<F[i].FI[n];
		}
	}
	oflg<<endl;

	cout<<"Setting Points for Pore-Network......"<<endl;
	ExtPtSize=Pt.size();
	NetworkPtSize=0;

	for(size_t i=0; i<P.size(); ++i) {
		if(P[i].CI.size()!=0) {
			P[i].PtI.resize(P[i].CI.size()+3, voidIndex);
			P[i].PtP.resize(P[i].CI.size()+3, nullptr);
		}
	}
	for(size_t i=0; i<T.size(); ++i) {
		T[i].PtI.resize(cornerNumber, voidIndex);
		T[i].PtP.resize(cornerNumber, nullptr);
	}

	//set throats radius for display, this need to be reversed later
	numeric_t MinTRadius(0), MaxTRadius(0);
	for(size_t k=0; k<porousBlocks; ++k) {
		for(size_t i=0; i<Bk[k].avgTDiameter.size(); ++i) {
			if(k==0 && i==0) {
				MinTRadius=cd_c::phy.RefLength*Bk[k].avgTDiameter[i]/2;
				MaxTRadius=cd_c::phy.RefLength*Bk[k].avgTDiameter[i]/2;
			} else {
				if(cd_c::phy.RefLength*Bk[k].avgTDiameter[i]/2<MinTRadius) MinTRadius=cd_c::phy.RefLength*Bk[k].avgTDiameter[i]/2;
				if(cd_c::phy.RefLength*Bk[k].avgTDiameter[i]/2>MaxTRadius) MaxTRadius=cd_c::phy.RefLength*Bk[k].avgTDiameter[i]/2;
			}
		}
	}
	for(size_t k=0; k<porousBlocks; ++k) {
		for(size_t i=0; i<Bk[k].avgTDiameter.size(); ++i) {
			if(MinTRadius==MaxTRadius && MinTRadius<MinDisplayThroatRadius) {
				Bk[k].TDisplayScale[i]=MinDisplayThroatRadius/MinTRadius;
			} else if(MinTRadius<MaxTRadius) {
//				numeric_t ourradius(cd_c::phy.RefLength*Bk[k].avgTDiameter[i]/2);
//				Bk[k].TDisplayScale[i]=(ourradius-MinTRadius)*(MaxDisplayThroatRadius-MinDisplayThroatRadius)/(MaxTRadius-MinTRadius)/(ourradius)
//						              +MinDisplayThroatRadius/(ourradius);
				Bk[k].TDisplayScale[i]=(MaxDisplayThroatRadius+MinDisplayThroatRadius)/(MaxTRadius+MinTRadius);
			}
		}
	}
	for(size_t k=0; k<porousBlocks; ++k) {
		for(size_t i=0; i<Bk[k].TI.size(); ++i) {
			T[Bk[k].TI[i]].RI*=Bk[k].TDisplayScale[0];
		}
	}

	for(size_t k=0; k<porousBlocks; ++k) {
		if(Bk[k].B[e].T==PorousOpen) {
			size_t i=Bk[k].Nx;
			for(size_t j=0; j<=Bk[k].Ny; ++j) {
				if(Bk[k].PI[i][j]==voidIndex) {
					//do nothing
				} else {
					numeric_t X(0), Y(0);
					numeric_t y1(P[P[Bk[k].PI[i][j]].PI[0]].y), x1(P[P[Bk[k].PI[i][j]].PI[0]].x)
							, y0(P[Bk[k].PI[i][j]].y)         , x0(P[Bk[k].PI[i][j]].x);

					X=P[Bk[k].PI[i][j]].x-T[P[Bk[k].PI[i][j]].TI[0]].RI;

					Y=(y1-y0)*(X-x0)/(x1-x0)+y0+T[P[Bk[k].PI[i][j]].TI[0]].RI*sqrt(pow((y1-y0)/(x1-x0), 2)+1);
					P[Bk[k].PI[i][j]].PtI[0]=Pt.size(); T[P[Bk[k].PI[i][j]].TI[0]].PtI[0]=Pt.size();
					Pt.push_back(point_c(X, Y));
					NetworkPtSize++;

					Y=(y1-y0)*(X-x0)/(x1-x0)+y0-T[P[Bk[k].PI[i][j]].TI[0]].RI*sqrt(pow((y1-y0)/(x1-x0), 2)+1);
					P[Bk[k].PI[i][j]].PtI[1]=Pt.size(); T[P[Bk[k].PI[i][j]].TI[0]].PtI[3]=Pt.size();
					Pt.push_back(point_c(X, Y));
					NetworkPtSize++;

					if(       Bk[k].NinLength!=0 || Bk[k].NinDiameter==1) {
						P[Bk[k].PI[i][j]].PtI[2]=C[P[Bk[k].PI[i][j]].CI[0]].FI[sw];
						P[Bk[k].PI[i][j]].PtI[3]=C[P[Bk[k].PI[i][j]].CI[0]].FI[nw];
					} else if(Bk[k].NinDiameter>1) {
						vector<index_t> CopyCI(P[Bk[k].PI[i][j]].CI), SortedCI(0);
						while(CopyCI.size()>0) {
							numeric_t min(C[P[Bk[k].PI[i][j]].CI[0]].y);
							size_t minindex(0);

							for(size_t l=1; l<CopyCI.size(); ++l) {
								if(C[P[Bk[k].PI[i][j]].CI[l]].y<min) {
									min=C[P[Bk[k].PI[i][j]].CI[l]].y;
									minindex=l;
								}
							}

							SortedCI.push_back(CopyCI[minindex]);
							CopyCI.erase(CopyCI.begin()+minindex);
						}

						for(size_t l=0; l<SortedCI.size(); ++l) {
							if(l==0) {
								P[Bk[k].PI[i][j]].PtI[2  ]=C[SortedCI[l]].FI[sw];
								P[Bk[k].PI[i][j]].PtI[3  ]=C[SortedCI[l]].FI[nw];
							} else {
								P[Bk[k].PI[i][j]].PtI[3+l]=C[SortedCI[l]].FI[nw];
							}
						}
					}
				}//jump out of if void
			}
		}//jump out of if B[e] is open
		if(Bk[k].B[w].T==PorousOpen) {
			size_t i=0;
			for(size_t j=0; j<=Bk[k].Ny; ++j) {
				if(Bk[k].PI[i][j]==voidIndex) {
					//do nothing
				} else {
					numeric_t X(0), Y(0);
					numeric_t y1(P[P[Bk[k].PI[i][j]].PI[0]].y), x1(P[P[Bk[k].PI[i][j]].PI[0]].x)
							, y0(P[Bk[k].PI[i][j]].y)         , x0(P[Bk[k].PI[i][j]].x);

					X=P[Bk[k].PI[i][j]].x+T[P[Bk[k].PI[i][j]].TI[0]].RI;

					Y=(y1-y0)*(X-x0)/(x1-x0)+y0-T[P[Bk[k].PI[i][j]].TI[0]].RI*sqrt(pow((y1-y0)/(x1-x0), 2)+1);
					P[Bk[k].PI[i][j]].PtI[0]=Pt.size(); T[P[Bk[k].PI[i][j]].TI[0]].PtI[2]=Pt.size();
					Pt.push_back(point_c(X, Y));
					NetworkPtSize++;

					Y=(y1-y0)*(X-x0)/(x1-x0)+y0+T[P[Bk[k].PI[i][j]].TI[0]].RI*sqrt(pow((y1-y0)/(x1-x0), 2)+1);
					P[Bk[k].PI[i][j]].PtI[1]=Pt.size(); T[P[Bk[k].PI[i][j]].TI[0]].PtI[1]=Pt.size();
					Pt.push_back(point_c(X, Y));
					NetworkPtSize++;

					if(       Bk[k].NinLength!=0 || Bk[k].NinDiameter==1) {
						P[Bk[k].PI[i][j]].PtI[2]=C[P[Bk[k].PI[i][j]].CI[0]].FI[ne];
						P[Bk[k].PI[i][j]].PtI[3]=C[P[Bk[k].PI[i][j]].CI[0]].FI[se];
					} else if(Bk[k].NinDiameter>1) {
						vector<index_t> CopyCI(P[Bk[k].PI[i][j]].CI), SortedCI(0);
						while(CopyCI.size()>0) {
							numeric_t max(C[P[Bk[k].PI[i][j]].CI[0]].y);
							size_t maxindex(0);

							for(size_t l=1; l<CopyCI.size(); ++l) {
								if(C[P[Bk[k].PI[i][j]].CI[l]].y>max) {
									max=C[P[Bk[k].PI[i][j]].CI[l]].y;
									maxindex=l;
								}
							}

							SortedCI.push_back(CopyCI[maxindex]);
							CopyCI.erase(CopyCI.begin()+maxindex);
						}

						for(size_t l=0; l<SortedCI.size(); ++l) {
							if(l==0) {
								P[Bk[k].PI[i][j]].PtI[2  ]=C[SortedCI[l]].FI[ne];
								P[Bk[k].PI[i][j]].PtI[3  ]=C[SortedCI[l]].FI[se];
							} else {
								P[Bk[k].PI[i][j]].PtI[3+l]=C[SortedCI[l]].FI[se];
							}
						}
					}
				}//jump out of if void
			}
		}//jump out of B[w] is open
		if(Bk[k].B[n].T==PorousOpen) {
			size_t j=Bk[k].Ny;
			for(size_t i=0; i<=Bk[k].Nx; ++i) {
				if(Bk[k].PI[i][j]==voidIndex) {
					//do nothing
				} else {
					numeric_t X(0), Y(0);
					numeric_t y1(P[P[Bk[k].PI[i][j]].PI[0]].y), x1(P[P[Bk[k].PI[i][j]].PI[0]].x)
							, y0(P[Bk[k].PI[i][j]].y)         , x0(P[Bk[k].PI[i][j]].x);

					Y=P[Bk[k].PI[i][j]].y-T[P[Bk[k].PI[i][j]].TI[0]].RI;

					X=(x1-x0)*(Y-y0)/(y1-y0)+x0-T[P[Bk[k].PI[i][j]].TI[0]].RI*sqrt(pow((x1-x0)/(y1-y0), 2)+1);
					P[Bk[k].PI[i][j]].PtI[0]=Pt.size(); T[P[Bk[k].PI[i][j]].TI[0]].PtI[1]=Pt.size();
					Pt.push_back(point_c(X, Y));
					NetworkPtSize++;

					X=(x1-x0)*(Y-y0)/(y1-y0)+x0+T[P[Bk[k].PI[i][j]].TI[0]].RI*sqrt(pow((x1-x0)/(y1-y0), 2)+1);
					P[Bk[k].PI[i][j]].PtI[1]=Pt.size(); T[P[Bk[k].PI[i][j]].TI[0]].PtI[0]=Pt.size();
					Pt.push_back(point_c(X, Y));
					NetworkPtSize++;

					if(       Bk[k].NinLength!=0 || Bk[k].NinDiameter==1) {
						P[Bk[k].PI[i][j]].PtI[2]=C[P[Bk[k].PI[i][j]].CI[0]].FI[se];
						P[Bk[k].PI[i][j]].PtI[3]=C[P[Bk[k].PI[i][j]].CI[0]].FI[sw];
					} else if(Bk[k].NinDiameter>1) {
						vector<index_t> CopyCI(P[Bk[k].PI[i][j]].CI), SortedCI(0);
						while(CopyCI.size()>0) {
							numeric_t max(C[P[Bk[k].PI[i][j]].CI[0]].x);
							size_t maxindex(0);

							for(size_t l=1; l<CopyCI.size(); ++l) {
								if(C[P[Bk[k].PI[i][j]].CI[l]].x>max) {
									max=C[P[Bk[k].PI[i][j]].CI[l]].x;
									maxindex=l;
								}
							}

							SortedCI.push_back(CopyCI[maxindex]);
							CopyCI.erase(CopyCI.begin()+maxindex);
						}

						for(size_t l=0; l<SortedCI.size(); ++l) {
							if(l==0) {
								P[Bk[k].PI[i][j]].PtI[2  ]=C[SortedCI[l]].FI[se];
								P[Bk[k].PI[i][j]].PtI[3  ]=C[SortedCI[l]].FI[sw];
							} else {
								P[Bk[k].PI[i][j]].PtI[3+l]=C[SortedCI[l]].FI[sw];
							}
						}
					}
				}//jump out of if void
			}
		}
		if(Bk[k].B[s].T==PorousOpen) {
			size_t j=0;
			for(size_t i=0; i<=Bk[k].Nx; ++i) {
				if(Bk[k].PI[i][j]==voidIndex) {
					//do nothing
				} else {
					numeric_t X(0), Y(0);
					numeric_t y1(P[P[Bk[k].PI[i][j]].PI[0]].y), x1(P[P[Bk[k].PI[i][j]].PI[0]].x)
							, y0(P[Bk[k].PI[i][j]].y)         , x0(P[Bk[k].PI[i][j]].x);

					Y=P[Bk[k].PI[i][j]].y+T[P[Bk[k].PI[i][j]].TI[0]].RI;

					X=(x1-x0)*(Y-y0)/(y1-y0)+x0+T[P[Bk[k].PI[i][j]].TI[0]].RI*sqrt(pow((x1-x0)/(y1-y0), 2)+1);
					P[Bk[k].PI[i][j]].PtI[0]=Pt.size(); T[P[Bk[k].PI[i][j]].TI[0]].PtI[3]=Pt.size();
					Pt.push_back(point_c(X, Y));
					NetworkPtSize++;

					X=(x1-x0)*(Y-y0)/(y1-y0)+x0-T[P[Bk[k].PI[i][j]].TI[0]].RI*sqrt(pow((x1-x0)/(y1-y0), 2)+1);
					P[Bk[k].PI[i][j]].PtI[1]=Pt.size(); T[P[Bk[k].PI[i][j]].TI[0]].PtI[2]=Pt.size();
					Pt.push_back(point_c(X, Y));
					NetworkPtSize++;

					if(       Bk[k].NinLength!=0 || Bk[k].NinDiameter==1) {
						P[Bk[k].PI[i][j]].PtI[2]=C[P[Bk[k].PI[i][j]].CI[0]].FI[nw];
						P[Bk[k].PI[i][j]].PtI[3]=C[P[Bk[k].PI[i][j]].CI[0]].FI[ne];
					} else if(Bk[k].NinDiameter>1) {
						vector<index_t> CopyCI(P[Bk[k].PI[i][j]].CI), SortedCI(0);
						while(CopyCI.size()>0) {
							numeric_t min(C[P[Bk[k].PI[i][j]].CI[0]].x);
							size_t minindex(0);

							for(size_t l=1; l<CopyCI.size(); ++l) {
								if(C[P[Bk[k].PI[i][j]].CI[l]].x<min) {
									min=C[P[Bk[k].PI[i][j]].CI[l]].x;
									minindex=l;
								}
							}

							SortedCI.push_back(CopyCI[minindex]);
							CopyCI.erase(CopyCI.begin()+minindex);
						}

						for(size_t l=0; l<SortedCI.size(); ++l) {
							if(l==0) {
								P[Bk[k].PI[i][j]].PtI[2  ]=C[SortedCI[l]].FI[nw];
								P[Bk[k].PI[i][j]].PtI[3  ]=C[SortedCI[l]].FI[ne];
							} else {
								P[Bk[k].PI[i][j]].PtI[3+l]=C[SortedCI[l]].FI[ne];
							}
						}
					}
				}//jump out if void
			}
		}
	}

	for(size_t i=0; i<P.size(); ++i) {
		if(P[i].CI.size()==0) {
			vector<direction_e> PD(0);//Neighbor Pore Direction
			for(size_t j=0; j<P[i].PI.size(); ++j) {
				PD.push_back(PoreRelation(i, P[i].PI[j]));
			}

			vector<index_t> fullPI (nodeNumber, voidIndex);
			vector<index_t> fullTI (nodeNumber, voidIndex);
			vector<index_t> whichPT(nodeNumber, voidIndex);

			for(size_t j=0; j<P[i].PI.size(); ++j) {
				if(PD[j]== e) {fullPI[0]=P[i].PI[j]; fullTI[0]=P[i].TI[j]; whichPT[0]=j;}
				if(PD[j]==ne) {fullPI[1]=P[i].PI[j]; fullTI[1]=P[i].TI[j]; whichPT[1]=j;}
				if(PD[j]==n ) {fullPI[2]=P[i].PI[j]; fullTI[2]=P[i].TI[j]; whichPT[2]=j;}
				if(PD[j]==nw) {fullPI[3]=P[i].PI[j]; fullTI[3]=P[i].TI[j]; whichPT[3]=j;}
				if(PD[j]== w) {fullPI[4]=P[i].PI[j]; fullTI[4]=P[i].TI[j]; whichPT[4]=j;}
				if(PD[j]==sw) {fullPI[5]=P[i].PI[j]; fullTI[5]=P[i].TI[j]; whichPT[5]=j;}
				if(PD[j]==s ) {fullPI[6]=P[i].PI[j]; fullTI[6]=P[i].TI[j]; whichPT[6]=j;}
				if(PD[j]==se) {fullPI[7]=P[i].PI[j]; fullTI[7]=P[i].TI[j]; whichPT[7]=j;}
			}

			for(size_t d=0; d<=7; ++d) {
				if(whichPT[d]!=voidIndex) {
					size_t nextd=whichPT[(d+1)%8]!=voidIndex?(d+1)%8:(
							     whichPT[(d+2)%8]!=voidIndex?(d+2)%8:(
							     whichPT[(d+3)%8]!=voidIndex?(d+3)%8:(
							     whichPT[(d+4)%8]!=voidIndex?(d+4)%8:(
								 whichPT[(d+5)%8]!=voidIndex?(d+5)%8:(
								 whichPT[(d+6)%8]!=voidIndex?(d+6)%8:(
								 whichPT[(d+7)%8]!=voidIndex?(d+7)%8:(d+8)%8))))));
					size_t GapBT=whichPT[(d+1)%8]!=voidIndex?1:(
							     whichPT[(d+2)%8]!=voidIndex?2:(
							     whichPT[(d+3)%8]!=voidIndex?3:(
							     whichPT[(d+4)%8]!=voidIndex?4:(
								 whichPT[(d+5)%8]!=voidIndex?5:(
								 whichPT[(d+6)%8]!=voidIndex?6:(
								 whichPT[(d+7)%8]!=voidIndex?7:0))))));

					switch (GapBT) {
					case 1: case 2: case 3:
							P[i].PtI.push_back(Pt.size());
							switch (    d) {
							case 0: case 1:
								T[fullTI[    d]].PtI[1]=Pt.size();
								break;
							case 2: case 3:
								T[fullTI[    d]].PtI[2]=Pt.size();
								break;
							case 4: case 5:
								T[fullTI[    d]].PtI[3]=Pt.size();
								break;
							case 6: case 7:
								T[fullTI[    d]].PtI[0]=Pt.size();
								break;
							}
							switch (nextd) {
							case 0: case 1:
								T[fullTI[nextd]].PtI[2]=Pt.size();
								break;
							case 2: case 3:
								T[fullTI[nextd]].PtI[3]=Pt.size();
								break;
							case 4: case 5:
								T[fullTI[nextd]].PtI[0]=Pt.size();
								break;
							case 6: case 7:
								T[fullTI[nextd]].PtI[1]=Pt.size();
								break;
							}
							Pt.push_back(crossPoint(P[i], P[fullPI[d]], T[fullTI[d]].RI, P[fullPI[nextd]], -T[fullTI[nextd]].RI));
							NetworkPtSize++;
							break;
					case 4:
							P[i].PtI.push_back(Pt.size());
							switch (    d) {
							case 0: case 1:
								T[fullTI[    d]].PtI[1]=Pt.size();
								break;
							case 2: case 3:
								T[fullTI[    d]].PtI[2]=Pt.size();
								break;
							case 4: case 5:
								T[fullTI[    d]].PtI[3]=Pt.size();
								break;
							case 6: case 7:
								T[fullTI[    d]].PtI[0]=Pt.size();
								break;
							}
							Pt.push_back(noCrossPoint(P[i], P[fullPI[d]], T[fullTI[d]].RI));
							NetworkPtSize++;

							P[i].PtI.push_back(Pt.size());
							switch (nextd) {
							case 0: case 1:
								T[fullTI[nextd]].PtI[2]=Pt.size();
								break;
							case 2: case 3:
								T[fullTI[nextd]].PtI[3]=Pt.size();
								break;
							case 4: case 5:
								T[fullTI[nextd]].PtI[0]=Pt.size();
								break;
							case 6: case 7:
								T[fullTI[nextd]].PtI[1]=Pt.size();
								break;
							}
							Pt.push_back(noCrossPoint(P[i], P[fullPI[nextd]], -T[fullTI[nextd]].RI));
							NetworkPtSize++;
							break;
					case 5: case 6: case 7:
							P[i].PtI.push_back(Pt.size());
							switch (    d) {
							case 0: case 1:
								T[fullTI[    d]].PtI[1]=Pt.size();
								break;
							case 2: case 3:
								T[fullTI[    d]].PtI[2]=Pt.size();
								break;
							case 4: case 5:
								T[fullTI[    d]].PtI[3]=Pt.size();
								break;
							case 6: case 7:
								T[fullTI[    d]].PtI[0]=Pt.size();
								break;
							}
							Pt.push_back(noCrossPoint(P[i], P[fullPI[d]], T[fullTI[d]].RI));
							NetworkPtSize++;

							P[i].PtI.push_back(Pt.size());
							Pt.push_back(crossPoint(P[i], P[fullPI[d]], T[fullTI[d]].RI, P[fullPI[nextd]], -T[fullTI[nextd]].RI));
							NetworkPtSize++;

							P[i].PtI.push_back(Pt.size());
							switch (nextd) {
							case 0: case 1:
								T[fullTI[nextd]].PtI[2]=Pt.size();
								break;
							case 2: case 3:
								T[fullTI[nextd]].PtI[3]=Pt.size();
								break;
							case 4: case 5:
								T[fullTI[nextd]].PtI[0]=Pt.size();
								break;
							case 6: case 7:
								T[fullTI[nextd]].PtI[1]=Pt.size();
								break;
							}
							Pt.push_back(noCrossPoint(P[i], P[fullPI[nextd]], -T[fullTI[nextd]].RI));
							NetworkPtSize++;
							break;
					case 0:
							P[i].PtI.push_back(Pt.size());
							switch (    d) {
							case 0: case 1:
								T[fullTI[    d]].PtI[1]=Pt.size();
								break;
							case 2: case 3:
								T[fullTI[    d]].PtI[2]=Pt.size();
								break;
							case 4: case 5:
								T[fullTI[    d]].PtI[3]=Pt.size();
								break;
							case 6: case 7:
								T[fullTI[    d]].PtI[0]=Pt.size();
								break;
							}
							Pt.push_back(noCrossPoint(P[i], P[fullPI[d]], T[fullTI[d]].RI));
							NetworkPtSize++;

							P[i].PtI.push_back(Pt.size());
							Pt.push_back(noCrossPoint(P[i], point_c(-P[fullPI[d]].x, -P[fullPI[d]].y), -T[fullTI[d]].RI));
							NetworkPtSize++;

							P[i].PtI.push_back(Pt.size());
							Pt.push_back(noCrossPoint(P[i], point_c(-P[fullPI[d]].x, -P[fullPI[d]].y),  T[fullTI[d]].RI));
							NetworkPtSize++;

							P[i].PtI.push_back(Pt.size());
							switch (nextd) {
							case 0: case 1:
								T[fullTI[nextd]].PtI[2]=Pt.size();
								break;
							case 2: case 3:
								T[fullTI[nextd]].PtI[3]=Pt.size();
								break;
							case 4: case 5:
								T[fullTI[nextd]].PtI[0]=Pt.size();
								break;
							case 6: case 7:
								T[fullTI[nextd]].PtI[1]=Pt.size();
								break;
							}
							Pt.push_back(noCrossPoint(P[i], P[fullPI[nextd]], -T[fullTI[nextd]].RI));
							NetworkPtSize++;

							break;
					}
				}
			}//get out direction loop
		}
	}

	for(size_t k=0; k<porousBlocks; ++k) {
		for(size_t i=0; i<Bk[k].TI.size(); ++i) {
			T[Bk[k].TI[i]].RI/=Bk[k].TDisplayScale[0];
		}
	}
	if(NetworkPtSize==Pt.size()-ExtPtSize) cout<<"Points added for pore net work!"<<endl;
	else cerr<<"Something is wrong with the points in pore network!"<<endl;

	cout<<"Updating global pointers......"<<endl;
	cout<<"Resizing NP, FP, PP......"<<endl;
	for(size_t i=0; i<C.size(); ++i) {
		C[i].NP.resize(C[i].NI.size(), nullptr);
		C[i].FP.resize(C[i].FI.size(), nullptr);
		C[i].PP.resize(C[i].PI.size(), nullptr);
	}

	for(size_t i=0; i<U.size(); ++i) {
		U[i].NP.resize(U[i].NI.size(), nullptr);
		U[i].FP.resize(U[i].FI.size(), nullptr);
	}

	for(size_t i=0; i<V.size(); ++i) {
		V[i].NP.resize(V[i].NI.size(), nullptr);
		V[i].FP.resize(V[i].FI.size(), nullptr);
	}

	for(size_t i=0; i<P.size(); ++i) {
		P[i].NP.resize(P[i].CI.size(), nullptr);
		P[i].PP.resize(P[i].PI.size(), nullptr);
		P[i].TP.resize(P[i].TI.size(), nullptr);
		P[i].PtP.resize(P[i].PtI.size(), nullptr);
	}

	for(size_t i=0; i<T.size(); ++i) {
		T[i].PP.resize(T[i].PI.size(), nullptr);
		T[i].PtP.resize(T[i].PtI.size(), nullptr);
	}

	cout<<"Set neighbors and corner nodes pointer for cCells......"<<endl;
	for(size_t i=0; i<C.size(); ++i) {
		C[i].NP[E ]=C[i].NI[E ]==voidIndex?nullptr:& C[C[i].NI[E ]];
		C[i].NP[N ]=C[i].NI[N ]==voidIndex?nullptr:& C[C[i].NI[N ]];
		C[i].NP[W ]=C[i].NI[W ]==voidIndex?nullptr:& C[C[i].NI[W ]];
		C[i].NP[S ]=C[i].NI[S ]==voidIndex?nullptr:& C[C[i].NI[S ]];
		C[i].NP[EE]=C[i].NI[EE]==voidIndex?nullptr:& C[C[i].NI[EE]];
		C[i].NP[NN]=C[i].NI[NN]==voidIndex?nullptr:& C[C[i].NI[NN]];
		C[i].NP[WW]=C[i].NI[WW]==voidIndex?nullptr:& C[C[i].NI[WW]];
		C[i].NP[SS]=C[i].NI[SS]==voidIndex?nullptr:& C[C[i].NI[SS]];

		C[i].FP[e ]=C[i].FI[e ]==voidIndex?nullptr:& U[C[i].FI[e ]];
		C[i].FP[n ]=C[i].FI[n ]==voidIndex?nullptr:& V[C[i].FI[n ]];
		C[i].FP[w ]=C[i].FI[w ]==voidIndex?nullptr:& U[C[i].FI[w ]];
		C[i].FP[s ]=C[i].FI[s ]==voidIndex?nullptr:& V[C[i].FI[s ]];
		C[i].FP[ne]=C[i].FI[ne]==voidIndex?nullptr:& F[C[i].FI[ne]];
		C[i].FP[nw]=C[i].FI[nw]==voidIndex?nullptr:& F[C[i].FI[nw]];
		C[i].FP[sw]=C[i].FI[sw]==voidIndex?nullptr:& F[C[i].FI[sw]];
		C[i].FP[se]=C[i].FI[se]==voidIndex?nullptr:& F[C[i].FI[se]];

		for(size_t j=0; j<C[i].PP.size(); ++j) {
			C[i].PP[j]=& P[C[i].PI[j]];
		}
	}

	cout<<"Set neighbors and corner nodes pointer for uCells......"<<endl;
	for(size_t i=0; i<U.size(); ++i) {
		U[i].NP[E ]=U[i].NI[E ]==voidIndex?nullptr:& U[U[i].NI[E ]];
		U[i].NP[N ]=U[i].NI[N ]==voidIndex?nullptr:& U[U[i].NI[N ]];
		U[i].NP[W ]=U[i].NI[W ]==voidIndex?nullptr:& U[U[i].NI[W ]];
		U[i].NP[S ]=U[i].NI[S ]==voidIndex?nullptr:& U[U[i].NI[S ]];
		U[i].NP[EE]=U[i].NI[EE]==voidIndex?nullptr:& U[U[i].NI[EE]];
		U[i].NP[NN]=U[i].NI[NN]==voidIndex?nullptr:& U[U[i].NI[NN]];
		U[i].NP[WW]=U[i].NI[WW]==voidIndex?nullptr:& U[U[i].NI[WW]];
		U[i].NP[SS]=U[i].NI[SS]==voidIndex?nullptr:& U[U[i].NI[SS]];

		U[i].FP[e ]=U[i].FI[e ]==voidIndex?nullptr:& C[U[i].FI[e ]];
		U[i].FP[n ]=U[i].FI[n ]==voidIndex?nullptr:& F[U[i].FI[n ]];
		U[i].FP[w ]=U[i].FI[w ]==voidIndex?nullptr:& C[U[i].FI[w ]];
		U[i].FP[s ]=U[i].FI[s ]==voidIndex?nullptr:& F[U[i].FI[s ]];
		U[i].FP[ne]=U[i].FI[ne]==voidIndex?nullptr:& V[U[i].FI[ne]];
		U[i].FP[nw]=U[i].FI[nw]==voidIndex?nullptr:& V[U[i].FI[nw]];
		U[i].FP[sw]=U[i].FI[sw]==voidIndex?nullptr:& V[U[i].FI[sw]];
		U[i].FP[se]=U[i].FI[se]==voidIndex?nullptr:& V[U[i].FI[se]];
	}

	cout<<"Set neighbors and corner nodes pointer for vCells......"<<endl;
	for(size_t i=0; i<V.size(); ++i) {
		V[i].NP[E ]=V[i].NI[E ]==voidIndex?nullptr:& V[V[i].NI[E ]];
		V[i].NP[N ]=V[i].NI[N ]==voidIndex?nullptr:& V[V[i].NI[N ]];
		V[i].NP[W ]=V[i].NI[W ]==voidIndex?nullptr:& V[V[i].NI[W ]];
		V[i].NP[S ]=V[i].NI[S ]==voidIndex?nullptr:& V[V[i].NI[S ]];
		V[i].NP[EE]=V[i].NI[EE]==voidIndex?nullptr:& V[V[i].NI[EE]];
		V[i].NP[NN]=V[i].NI[NN]==voidIndex?nullptr:& V[V[i].NI[NN]];
		V[i].NP[WW]=V[i].NI[WW]==voidIndex?nullptr:& V[V[i].NI[WW]];
		V[i].NP[SS]=V[i].NI[SS]==voidIndex?nullptr:& V[V[i].NI[SS]];

		V[i].FP[e ]=V[i].FI[e ]==voidIndex?nullptr:& F[V[i].FI[e ]];
		V[i].FP[n ]=V[i].FI[n ]==voidIndex?nullptr:& C[V[i].FI[n ]];
		V[i].FP[w ]=V[i].FI[w ]==voidIndex?nullptr:& F[V[i].FI[w ]];
		V[i].FP[s ]=V[i].FI[s ]==voidIndex?nullptr:& C[V[i].FI[s ]];
		V[i].FP[ne]=V[i].FI[ne]==voidIndex?nullptr:& U[V[i].FI[ne]];
		V[i].FP[nw]=V[i].FI[nw]==voidIndex?nullptr:& U[V[i].FI[nw]];
		V[i].FP[sw]=V[i].FI[sw]==voidIndex?nullptr:& U[V[i].FI[sw]];
		V[i].FP[se]=V[i].FI[se]==voidIndex?nullptr:& U[V[i].FI[se]];
	}

	cout<<"Set member pointers for pores & throats......"<<endl;
	for(size_t i=0; i<P.size(); ++i) {
		for(size_t j=0; j<P[i].CI.size(); ++j) {
			P[i].NP[j]=& C[P[i].CI[j]];
		}
		for(size_t j=0; j<P[i].PI.size(); ++j) {
			P[i].PP[j]=& P[P[i].PI[j]];
		}
		for(size_t j=0; j<P[i].TI.size(); ++j) {
			P[i].TP[j]=& T[P[i].TI[j]];
		}
		for(size_t j=0; j<P[i].PtI.size(); ++j) {
			P[i].PtP[j]=& Pt[P[i].PtI[j]];
		}
	}
	for(size_t i=0; i<T.size(); ++i) {
		for(size_t j=0; j<T[i].PI.size(); ++j) {
			T[i].PP[j]=& P[T[i].PI[j]];
		}
		for(size_t j=0; j<T[i].PtI.size(); ++j) {
			T[i].PtP[j]=& Pt[T[i].PtI[j]];
		}
	}

	oflg<<endl<<"PorousOpen connection information:"<<endl;
	for(size_t i=0; i<C.size(); ++i) {
		if(C[i].PI.size()!=0) {
			oflg<<"C # "<<i<<":";
			for(size_t j=0; j<C[i].PI.size(); ++j) {
				oflg<<"\tP"<<C[i].PI[j];
			}
			oflg<<endl;
		}
	}
	for(size_t i=0; i<P.size(); ++i) {
		if(P[i].CI.size()!=0) {
			oflg<<"P # "<<i<<":";
			for(size_t j=0; j<P[i].CI.size(); ++j) {
				oflg<<"\tC"<<P[i].CI[j]<<"\tL="<<P[i].CD[j]<<"\tA="<<P[i].CA[j];
			}
			oflg<<endl;
		}
	}

	oflg<<"Cells dimensions information:";
	oflg<<endl<<"cCell, Checking dimensions:";
	for(size_t i=0; i<C.size(); ++i) {
		oflg<<endl<<"C#"<<i<<":\t";
		for(size_t n=0; n<sideNumber; ++n) {
			oflg<<C[i].d[n]<<"\t";
		}
		oflg<<"|";
		for(size_t n=0; n<sideNumber; ++n) {
			oflg<<"\t"<<C[i].D[n];
		}
	}
	oflg<<endl<<"uCell, Checking dimensions:";
	for(size_t i=0; i<C.size(); ++i) {
		oflg<<endl<<"U#"<<i<<":\t";
		for(size_t n=0; n<sideNumber; ++n) {
			oflg<<U[i].d[n]<<"\t";
		}
		oflg<<"|";
		for(size_t n=0; n<sideNumber; ++n) {
			oflg<<"\t"<<U[i].D[n];
		}
	}
	oflg<<endl<<"vCell, Checking dimensions:";
	for(size_t i=0; i<V.size(); ++i) {
		oflg<<endl<<"V#"<<i<<":\t";
		for(size_t n=0; n<sideNumber; ++n) {
			oflg<<V[i].d[n]<<"\t";
		}
		oflg<<"|";
		for(size_t n=0; n<sideNumber; ++n) {
			oflg<<"\t"<<V[i].D[n];
		}
	}
	oflg<<endl;

	oflg<<"Inlet cells:";
	for(size_t i=0; i<uvCell_c::InletUI.size(); ++i) {
		oflg<<"\tU#"<<uvCell_c::InletUI[i];
	}
	for(size_t i=0; i<uvCell_c::InletVI.size(); ++i) {
		oflg<<"\tV#"<<uvCell_c::InletVI[i];
	}
	oflg<<endl;

	oflg<<"Outlet cells:";
	for(size_t i=0; i<uvCell_c::OutletUI.size(); ++i) {
		oflg<<"\tU#"<<uvCell_c::OutletUI[i];
	}
	for(size_t i=0; i<uvCell_c::OutletVI.size(); ++i) {
		oflg<<"\tV#"<<uvCell_c::OutletVI[i];
	}
	oflg<<endl;

	UVsize=U.size()+V.size(); CPsize=C.size()+P.size();

	cout<<"Meshing Finished!"<<endl;
	return true;
}

//----------------------------------------------------------------------------------------------------
direction_e PoreRelation(const index_t & i, const index_t & j) {
	index_t ib(voidIndex), ii(voidIndex), ij(voidIndex);
	findPoreLocation(ib, ii, ij, i);
	index_t jb(voidIndex), ji(voidIndex), jj(voidIndex);
	findPoreLocation(jb, ji, jj, j);

	if(     ib==jb && ji-ii== 1 && jj-ij== 0) return e;
	else if(ib==jb && ji-ii== 0 && jj-ij== 1) return n;
	else if(ib==jb && ji-ii==-1 && jj-ij== 0) return w;
	else if(ib==jb && ji-ii== 0 && jj-ij==-1) return s;
	else if(ib==jb && ji-ii== 1 && jj-ij== 1) return ne;
	else if(ib==jb && ji-ii==-1 && jj-ij== 1) return nw;
	else if(ib==jb && ji-ii==-1 && jj-ij==-1) return sw;
	else if(ib==jb && ji-ii== 1 && jj-ij==-1) return se;
	else cerr<<"Shit, these two pores are not connected at all!"<<endl;

	return e;
}
//----------------------------------------------------------------------------------------------------
bool findPoreLocation(index_t & b, index_t & i, index_t & j, const index_t & pi) {
	for(size_t bi=0; bi<porousBlocks; ++bi) {
		for(size_t ii=0; ii<=Bk[bi].Nx; ++ii) {
			for(size_t jj=0; jj<=Bk[bi].Ny; ++jj) {
				if(pi==Bk[bi].PI[ii][jj]) {
					b=bi; i=ii; j=jj;
					return true;
				}
			}
		}
	}

	return false;
}
//----------------------------------------------------------------------------------------------------
point_c crossPoint(const point_c & p0, const point_c & p1, const numeric_t & os1
		                             , const point_c & p2, const numeric_t & os2) {
	numeric_t d1(sqrt(pow(p1.x-p0.x, 2)+pow(p1.y-p0.y, 2))), d2(sqrt(pow(p2.x-p0.x, 2)+pow(p2.y-p0.y, 2)));
	numeric_t st1((p1.y-p0.y)/d1), ct1((p1.x-p0.x)/d1), st2((p2.y-p0.y)/d2), ct2((p2.x-p0.x)/d2);
	numeric_t rs1(p0.x*st1-p0.y*ct1-os1), rs2(p0.x*st2-p0.y*ct2-os2);

	point_c cp;

	cp.x=(rs2*ct1-rs1*ct2)/(ct1*st2-ct2*st1);
	cp.y=(rs2*st1-rs1*st2)/(ct1*st2-ct2*st1);

	return cp;
}

point_c noCrossPoint(const point_c & p0, const point_c & p1, const numeric_t & os) {
	numeric_t distance(sqrt(pow(p1.x-p0.x, 2)+pow(p1.y-p0.y, 2)));

	numeric_t sintheta((p1.y-p0.y)/distance), costheta((p1.x-p0.x)/distance);

	point_c ncp(point_c(abs(os)*costheta-os*sintheta+p0.x, abs(os)*sintheta+os*costheta+p0.y));

	return ncp;
}
//----------------------------------------------------------------------------------------------------
numeric_t RandomPercentage(const numeric_t & rseed) {
	numeric_t R(0);
//	srand(rseed);
//	R=rand();
	srandom(rseed);
	R=random();
	R/=RAND_MAX;
	return R;
}

//----------------------------------------------------------------------------------------------------
inline numeric_t CxofBk(const index_t & i) {return Cxof(Pt[Bk[i].CnI[0]], Pt[Bk[i].CnI[1]], Pt[Bk[i].CnI[2]], Pt[Bk[i].CnI[3]]);}
inline numeric_t CyofBk(const index_t & i) {return Cyof(Pt[Bk[i].CnI[0]], Pt[Bk[i].CnI[1]], Pt[Bk[i].CnI[2]], Pt[Bk[i].CnI[3]]);}

bool BlocksBetweenIsNotPorous(const index_t & i, const index_t & j) {
	numeric_t xsmall(0), xlarge(0), ysmall(0), ylarge(0);
	if(Ln[i].XorY==xAxis && Ln[j].XorY==xAxis) {
		xsmall=Ln[i].Pt[0                ].x;
		xlarge=Ln[i].Pt[Ln[i].Pt.size()-1].x;
		if(Ln[i].Pt[0].y<Ln[j].Pt[0].y) {
			ysmall=Ln[i].Pt.begin()->y;
			ylarge=Ln[j].Pt.begin()->y;
		} else if((Ln[i].Pt[0].y>Ln[j].Pt[0].y)) {
			ysmall=Ln[j].Pt.begin()->y;
			ylarge=Ln[i].Pt.begin()->y;
		} else {
			cerr<<"The two lines are overlapping......"<<endl;
		}
	} else if(Ln[i].XorY==yAxis && Ln[j].XorY==yAxis) {
		ysmall=Ln[i].Pt[0                ].y;
		ylarge=Ln[i].Pt[Ln[i].Pt.size()-1].y;
		if(Ln[i].Pt.begin()->x<Ln[j].Pt.begin()->x) {
			xsmall=Ln[i].Pt.begin()->x;
			xlarge=Ln[j].Pt.begin()->x;
		} else if(Ln[i].Pt.begin()->x>Ln[j].Pt.begin()->x) {
			xsmall=Ln[j].Pt.begin()->x;
			xlarge=Ln[i].Pt.begin()->x;
		} else {
			cerr<<"The two lines are overlapping......"<<endl;
		}
	} else {
		cerr<<"The two lines are not parallel, blocks between are meaningless!"<<endl;
		return false;
	}
	for(size_t k=0; k<Bk.size(); ++k) {
		numeric_t cx(CxofBk(k)), cy(CyofBk(k));
		if(cx>=xsmall && cx<=xlarge && cy>=ysmall && cy<=ylarge
		&& Bk[k].T==porousmedia) {
			return false;
		}
	}

	return true;
}

//----------------------------------------------------------------------------------------------------
bool ThroatExist(const throat_c & t=throat_c()) {
	for(size_t i=0; i<T.size(); ++i) {
		if(t==T[i]) return true;
	}
	return false;
}
//----------------------------------------------------------------------------------------------------
//
//                                        function readMesh
//
//----------------------------------------------------------------------------------------------------
bool readMesh() {
	return true;
}

/*----------------------------------------------------------------------------------------------------

                                       function plotMesh

----------------------------------------------------------------------------------------------------*/
bool plotMesh(ofstream & ofs) {//Real Size
	ofs<<"# vtk DataFile Version 2.0"<<endl;
	ofs<<"Mesh"<<endl;
	ofs<<"ASCII"<<endl;
	ofs<<"DATASET UNSTRUCTURED_GRID"<<endl;

	ofs<<"POINTS\t"<<  Pt.size()<<"\tdouble"<<endl;
	for(size_t i=0; i<Pt.size(); ++i) {
		ofs<<Pt[i].x*cd_c::phy.RefLength<<"\t"<<Pt[i].y*cd_c::phy.RefLength<<"\t"<<" 0  "<<endl;
	}

	ofs<<"CELLS\t"<<C.size()<<"\t"<<C.size()*5<<endl;
	for(size_t i=0; i<C.size(); ++i) {
		ofs<<"4\t"<<          C[i].FI[ne]<<"\t"<<          C[i].FI[nw]<<"\t"<<          C[i].FI[sw]<<"\t"<<          C[i].FI[se]<<endl;
	}

	ofs<<"CELL_TYPES\t"<<C.size()<<endl;
	for(size_t i=0; i<C.size(); ++i) {
		ofs<<" 9"<<endl;
	}

	cout<<"Flowfield Ploted!"<<endl;
	return true;
}
