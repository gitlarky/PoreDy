/*
 * CalculationDomain.cpp
 *
 *  Created on: May 11, 2015
 *      Author: david
 */

#include "CalculationDomain.h"

/*----------------------------------------------------------------------------------------------------
 *
 *                                        Initialize Static Members
 *
  ----------------------------------------------------------------------------------------------------*/
physics_c cd_c::phy           =physics_c();
numeric_t cd_c::TotLiquid     =0          ;
numeric_t cd_c::TotLiquid0    =0          ;
numeric_t cd_c::TimeToEmptyMax=0          ;
numeric_t cd_c::EDT           =0          ;
index_t   cd_c::ESS           =0          ;
numeric_t cd_c::dz            =1          ;//******************************attention********************

std::vector<index_t> cd_c::ClusterToEmpty(0, voidIndex);
//std::vector<index_t> cd_c::BulkToEmpty()
//index_t   cd_c::BulkToEmpty   =0          ;


