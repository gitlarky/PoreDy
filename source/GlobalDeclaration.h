/*
 * GlobalDeclaration.h
 *
 *  Created on: Apr 15, 2014
 *      Author: david
 */

#ifndef GLOBALDECLARATION_H_
#define GLOBALDECLARATION_H_

#include "Basic.h"
#include "Physics.h"
#include "BoundaryCondition.h"
#include "Geometry.h"
#include "CalculationDomain.h"
#include "FiniteVolume.h"
#include "PoreThroat.h"

extern std::string caseName;
extern std::ofstream oflg;
extern std::ofstream offh;

extern std::vector<  point_c> Pt;
extern std::vector<   line_c> Ln;
extern std::vector<  block_c> Bk;

extern std::vector<  uCell_c> U;
extern std::vector<  vCell_c> V;
extern std::vector<  cCell_c> C;
extern std::vector<   pore_c> P;
extern std::vector< throat_c> T;
extern std::vector<cluster_c> Ct;

extern std::vector<  fCell_c> F;



extern size_t porousBlocks, flowBlocks;
extern size_t UVsize, CPsize, ExtPtSize, NetworkPtSize;


void Testing();

bool readCase();

bool meshing();
bool readMesh();

bool initializeFlowfield();
bool readFlowfield();
bool setFlowfield();

bool calculateFlowfield();

bool calculateEvaporation();

#endif /* GLOBALDECLARATION_H_ */
