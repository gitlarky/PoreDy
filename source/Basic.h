/*
 * Basic.h
 *
 *  Created on: Apr 15, 2014
 *      Author: david
 */

#ifndef BASIC_H_
#define BASIC_H_
/*----------------------------------------------------------------------------------------------------
 *
 *                                        include what?
 *
  ----------------------------------------------------------------------------------------------------*/
#include <iostream>
#include <fstream>
#include <sstream>

#include <vector>
#include <map>
#include <string>
#include <bitset>

#include <iomanip>
#include <cstdlib>
#include <cmath>

#include <limits>

#include <algorithm>

#include <Eigen/Dense>
#include <Eigen/Sparse>

/*----------------------------------------------------------------------------------------------------
 *
 *                                        type definition
 *
  ----------------------------------------------------------------------------------------------------*/
typedef double numeric_t;//typedef float numeric_t;
typedef int index_t; //Range: [-2147483648	, 2147483647]//typedef long index_t; //Range: [-9223372036854775808	, 9223372036854775807]

/*----------------------------------------------------------------------------------------------------
 *
 *                                        enumeration definition
 *
  ----------------------------------------------------------------------------------------------------*/
enum    direction_e {e, n, w, s, ne, nw, sw, se};
enum       corner_e {NE, NW, SW, SE};
enum     neighbor_e {E, N, W, S, EE, NN, WW, SS, EEE, NNN, WWW, SSS, EEEE, NNNN, WWWW, SSSS};
enum   coordinate_e {origin, xAxis, yAxis};
enum parallelLine_e {notParallel, allHorizontal, allVertical};

enum    blockType_e {ghost, flowfield, porousmedia};
enum physicsModel_e {incompressibleFLow, isothermalEvaporation, poreNetwork};
enum       eqItem_e {convectionItem, diffussionItem, deltaPAVT, RorD};

enum     cellType_e {fCell, uvCell, uCell, vCell, cCell, pore};

enum porosityType_e {singlePorosity, downupPorosity, leftrightPorosity, assignPorosity};
enum     position_e {liquidSide, inside, meniscus, open};
enum     vpStatus_e {none, unknown, found, out, saturated};
//enum      sStatus_e {empty, film, partial, liquid};
enum      sStatus_e {empty, film, liquid};

enum    PDEBCType_e {NoPDEBC, Dirichlet, Neumann, Cauchy};
enum       BCType_e {NoBC, Inlet, Outlet, Wall, PeriodicUpstream, PeriodicDownstream, PorousOpen, VelocityInlet, PressureOutlet, ConstConcentration, Symmetry
	                     , FixedValueBC, Relation0BC, Relation1BC, Relation2BC};

enum    algorithm_e {SIMPLE, SIMPLER, SIMPLEC, PISO, OperatorSplitting, NonOperatorSplitting, ConvectionDiffusion, DiffusionConvection};
enum       scheme_e {Upwind, Hybrid, PowerLaw, QUICK, HayaseQUICK, TVD};
enum  schemeOrder_e {firstOrder=1, secondOrder=2};

enum        shape_e {circle=0, triangle=3, square, pentagon, hexagon};

/*----------------------------------------------------------------------------------------------------
 *
 *                                        constant definition
 *
  ----------------------------------------------------------------------------------------------------*/
const    size_t outputPrecision(4);

const   index_t voidIndex(-1);

const    size_t HowManyCopy(3);
const    size_t HowManyVariable(4);//u v p c

const numeric_t defaultGrowth(1);
const    size_t defaultSegment(1);
const    size_t neighborNumber(8), nodeNumber(8), sideNumber(4), cornerNumber(4);
const    size_t EqComponent(4);
const    size_t VariableCount(4);

const numeric_t pi(3.1415926);
const numeric_t UniversalR(8.3144621);
const numeric_t atm(1.01325e5);

const numeric_t ThroatAverageLength(5.0e-4), ThroatAverageDiameter(9.0e-5);
const numeric_t Psat(2.333e3), Csat(0.9592), Cliquid(55394.01);
const numeric_t Pa(1.01325e5), Pout(0);
const numeric_t dryCriteria(1e-6);

const numeric_t MinDisplayThroatRadius(25e-6), MaxDisplayThroatRadius(50e-6);

//----------------------------------------------------------------------------------------------------
#endif /* BASIC_H_ */
