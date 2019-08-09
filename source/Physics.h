/*
 * Physics.h
 *
 *  Created on: Apr 15, 2014
 *      Author: david
 */

#ifndef PHYSICS_H_
#define PHYSICS_H_

#include "Basic.h"

class    material_c;
class dryMaterial_c;
class wetMaterial_c;
class     physics_c;

/*----------------------------------------------------------------------------------------------------
 *
 *                                        Physics-Material
 *
  ----------------------------------------------------------------------------------------------------*/
class material_c {
public:
	std::string material;
	numeric_t MolarDensity;//[Kg/mol]
	numeric_t density; //[kg/m^3]


	material_c(const std::string & s="Air", const numeric_t & md=0.02897, const numeric_t & d=1.204)
	:material(s), MolarDensity(md), density(d) {}

	friend std::istream & operator >> (std::istream &is , material_c & m);
	friend std::ostream & operator << (std::ostream & os, const material_c & m);
};

//----------------------------------------------------------------------------------------------------
class dryMaterial_c:public material_c {
public:
	numeric_t dynamicViscosity;

	dryMaterial_c(const std::string & s="Air", const numeric_t & md=0.02897, const numeric_t & d=1.204
			    , const numeric_t & dv=1.813e-5)
	: material_c(s, md, d)
	, dynamicViscosity(dv) {}

	friend std::istream & operator >> (std::istream & is, dryMaterial_c & m);
	friend std::ostream & operator << (std::ostream & os, const dryMaterial_c & m);
};

//----------------------------------------------------------------------------------------------------
class wetMaterial_c:public dryMaterial_c {
public:
	numeric_t massDiffusivity; //[m^2/s]
	numeric_t contactAngle; //[Radians]
	numeric_t surfaceTension; //[N/m]
	numeric_t saturatedConcentration; //[mol/m^3]

	wetMaterial_c(const std::string & s="Water", const numeric_t & md=0.01802, const numeric_t & d=998.2
			    , const numeric_t & dv=1.002e-3
			    , const numeric_t & mdi=2.119e-5, const numeric_t & cav=0, const numeric_t & stv=7.266e-2, const numeric_t & scv=Csat)
	: dryMaterial_c(s, md, d, dv)
	, massDiffusivity(mdi), contactAngle(cav), surfaceTension(stv), saturatedConcentration(scv){}

	friend std::istream & operator >> (std::istream & is, wetMaterial_c & m);
	friend std::ostream & operator << (std::ostream & os, const wetMaterial_c & m);
};

//----------------------------------------------------------------------------------------------------
class physics_c {
public:
	physicsModel_e phyModel;

	bool filmEffect, viscosityEffect;

	dryMaterial_c Dry;
	wetMaterial_c Wet;

	numeric_t RC; //Radius of Corner, decided by manufacturing technology

	numeric_t RefLength, RefM, RefDensity, RefViscosity, RefMassDiffusivity, RefSurfaceTension, RefConcentration;
	numeric_t RefTime, RefVelocity, RefPressure, RefMass, RefMassFlux;
	numeric_t Re, Pe, Ca; //used for calculation

	numeric_t FeatureLength, FeatureDensity, FeatureVelocity, FeaturePressure, FeatureViscosity;
	numeric_t ReFeature, PeFeature, CaFeature; //used for show the feathers

	numeric_t EnvironmentConcentration;
	numeric_t SaturatedConcentration, LiquidConcentration;

	algorithm_e AlgorithmFlowfield, AlgorithmEvaporation;
	scheme_e SchemeFlowfield, SchemeEvaporation;
	bool ImplicitFlowfield, ImplicitEvaporation;

	size_t FlowfieldMaxStep, EvaporationMaxStep;
	numeric_t RelativeConvergeCriteria, FlowCriteria, DryCriteria;

	size_t EvaporationSubStep;
	size_t FlowOutputFrequency, DryOutputFrequency;

	numeric_t UVRelaxFactor, PRelaxFactor, CRelaxFactor;

	size_t FilmApproximationType;
	numeric_t FilmApproximationCoefficient;

	bool InletFlowDeveloped;

	physics_c(const physicsModel_e & pm=isothermalEvaporation
			, const bool & feff=false, const bool & veff=false
			, const dryMaterial_c & dry=dryMaterial_c(), const wetMaterial_c & wet=wetMaterial_c()
			, const numeric_t & rcv=0
			, const numeric_t & rl=1
			, const numeric_t & flv=1, const numeric_t & fvv=0, const numeric_t & fpv=0
			, const numeric_t & ec=0
			, const algorithm_e & afl=SIMPLER          , const scheme_e & schfl=Hybrid, const bool & imf=true
			, const algorithm_e &  ae=OperatorSplitting, const scheme_e &  sche=Hybrid, const bool & ime=true
			, const size_t & fms=1000, const size_t & ems=1000
			, const numeric_t & rcc=1e-4, const numeric_t & fc=1e-4, const numeric_t & dc=1e-4, const numeric_t & emt=1
			, const size_t & ess=1
			, const size_t & fof=1, const size_t & dof=1
			, const numeric_t & uvrf=1, const numeric_t & prf=1, const numeric_t & crf=1
			, const size_t & fatv=1, const numeric_t & facv=3
			, const bool & ifdv=true)
	: phyModel(pm), filmEffect(feff), viscosityEffect(veff)
	, Dry(dry), Wet(wet)
	, RC(rcv)
    , RefLength(rl)
	, RefM(wet.MolarDensity), RefDensity(wet.MolarDensity*wet.saturatedConcentration), RefViscosity(wet.dynamicViscosity), RefMassDiffusivity(wet.massDiffusivity), RefSurfaceTension(wet.surfaceTension)
	, RefConcentration(wet.saturatedConcentration)
	, RefTime(0), RefVelocity(0), RefPressure(0)
	, RefMass(0), RefMassFlux(0)
	, Re(0), Pe(0), Ca(0)
	, FeatureLength(flv), FeatureDensity(dry.density), FeatureVelocity(fvv), FeaturePressure(fpv), FeatureViscosity(dry.dynamicViscosity)
	, ReFeature(0), PeFeature(0), CaFeature(0)
	, EnvironmentConcentration(ec)
	, SaturatedConcentration(Csat), LiquidConcentration(Cliquid)
	, AlgorithmFlowfield(afl), AlgorithmEvaporation(ae), SchemeFlowfield(schfl), SchemeEvaporation(sche), 	ImplicitFlowfield(imf), ImplicitEvaporation(ime)
	, FlowfieldMaxStep(fms), EvaporationMaxStep(ems), RelativeConvergeCriteria(rcc), FlowCriteria(fc), DryCriteria(dc)
	, EvaporationSubStep(ess)
	, FlowOutputFrequency(fof), DryOutputFrequency(dof)
	, UVRelaxFactor(uvrf), PRelaxFactor(prf), CRelaxFactor(crf)
	, FilmApproximationType(fatv), FilmApproximationCoefficient(facv)
	, InletFlowDeveloped(ifdv) {
		RefTime=RefDensity*RefLength*RefLength/RefMassDiffusivity/RefM/RefConcentration;
		RefVelocity=RefLength/RefTime;
		RefPressure=RefDensity*RefVelocity*RefVelocity;

		RefMass=RefDensity*pow(RefLength, 3);
		RefMassFlux=RefMass/RefTime;

		Re=RefDensity*RefVelocity*RefLength/RefViscosity;
		Pe=RefVelocity*RefLength/RefMassDiffusivity;
		Ca=RefViscosity*RefVelocity/RefSurfaceTension;

		if(FeatureVelocity!=0 && FeaturePressure==0) {
			FeaturePressure=FeatureDensity*FeatureVelocity*FeatureVelocity;
		} else if(FeatureVelocity==0 && FeaturePressure!=0) {
			FeatureVelocity=sqrt(FeaturePressure/FeatureDensity);
		} else {
			//do nothing
		}
		ReFeature=FeatureDensity*FeatureVelocity*FeatureLength/FeatureViscosity;//All these Feature Properties are used to calculate the External Flow Field
		PeFeature=FeatureVelocity*FeatureLength/RefMassDiffusivity;
		CaFeature=RefViscosity*FeatureVelocity/RefSurfaceTension;//Almost no use, does not reflect any phenomena
	}

	friend std::istream & operator >> (std::istream & is, physics_c & p);
	friend std::ostream & operator << (std::ostream & os, const physics_c & p);

	bool nondimensionalize();
};

#endif /* PHYSICS_H_ */
