/*
 * CalculationDomain.h
 *
 *  Created on: Apr 29, 2015
 *      Author: david
 */

#ifndef CALCULATIONDOMAIN_H_
#define CALCULATIONDOMAIN_H_

#include "Basic.h"
#include "Physics.h"
#include "BoundaryCondition.h"
#include "Geometry.h"

class    cell_c;
class  uvCell_c;
class   uCell_c;
class   vCell_c;
class   cCell_c;
class   fCell_c;
class    pore_c;
class  throat_c;

/*----------------------------------------------------------------------------------------------------
 *
 *                                        CD class definition
 *
  ----------------------------------------------------------------------------------------------------*/
class cd_c {
public:
	static physics_c phy;
	static numeric_t TotLiquid, TotLiquid0;//TotLiquid0: Original Total Liquid
	static numeric_t TimeToEmptyMax, EDT;//EDT: Evaporation Delta Time
	static   index_t ESS;//ESS: Evaporation Sub Step
	static numeric_t dz;

	static std::vector<index_t> ClusterToEmpty;
//	static std::vector<index_t> BulkToEmpty, FilmToEmpty;


	cd_c() {}
	~cd_c()=default;

};

/*----------------------------------------------------------------------------------------------------
 *
 *                                        relationPoint_c class definition
 *
  ----------------------------------------------------------------------------------------------------*/
class cPoint_c: public cd_c, public point_c {
public:
	cellType_e T;

	cPoint_c(const numeric_t & xV=0, const numeric_t & yV=0
		   , const cellType_e & tv=fCell)
	: point_c(xV, yV)
	, T(tv) {}

	virtual ~cPoint_c()=default;
};

#endif /* CALCULATIONDOMAIN_H_ */

