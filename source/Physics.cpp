/*
 * Physics.cpp
 *
 *  Created on: Apr 15, 2014
 *      Author: david
 */

#include "Physics.h"

using namespace std;

//----------------------------------------------------------------------------------------------------
istream & operator >> (istream & is , material_c & m) {
	is>>m.material>>m.MolarDensity>>m.density;

	return is;
}
ostream & operator << (ostream & os, const material_c & m){
	os<<m.material<<":\tMolar Density:"<<m.MolarDensity<<"[kg/mol]\tDensity:"<<m.density<<"[kg/m^3]"<<endl;

	return os;
}

//----------------------------------------------------------------------------------------------------
istream & operator >> (istream & is, dryMaterial_c & m) {
	is>>m.material>>m.MolarDensity>>m.density
	  >>m.dynamicViscosity;

	return is;
}
ostream & operator << (ostream & os, const dryMaterial_c & m) {
	os<<m.material<<":\tMolar Density:"<<m.MolarDensity<<"[kg/mol]\tDensity:"<<m.density<<"[kg/m^3]"
	  <<"\tDynamic Viscosity:"<<m.dynamicViscosity<<"[Pa*s]"<<endl;

	return os;
}

//----------------------------------------------------------------------------------------------------
istream & operator >> (istream & is, wetMaterial_c & m) {
	is>>m.material>>m.MolarDensity>>m.density
	  >>m.dynamicViscosity
	  >>m.massDiffusivity>>m.contactAngle>>m.surfaceTension>>m.saturatedConcentration;

	return is;
}

ostream & operator << (ostream & os, const wetMaterial_c & m) {
	os<<m.material<<":\tMolar Density:"<<m.MolarDensity<<"[kg/mol]\tDensity:"<<m.density<<"[kg/m^3]"
	  <<"\tDynamic Viscosity:"<<m.dynamicViscosity<<"[Pa*s]"
	  <<"\tMass Diffusivity:"<<m.massDiffusivity<<"[m^2/s]\tContact Angle:"<<m.contactAngle
	  <<"[Radians]\tSurface Tension:"<<m.surfaceTension<<"[N/m]\tSaturated Concnetration:"<<m.saturatedConcentration<<"[mol/m^3]"<<endl;

	return os;
}

//----------------------------------------------------------------------------------------------------
istream & operator >> (istream & is, physics_c & p) {
	string str;

	is>>str; //Read In
	if(str=="incompressibleFLow") p.phyModel=incompressibleFLow;
	else if(str=="isothermalEvaporation") p.phyModel=isothermalEvaporation;
	else if(str=="poreNetwork") p.phyModel=poreNetwork;

	is>>p.filmEffect>>p.viscosityEffect;

	is>>p.Dry>>p.Wet>>p.RC>>p.RefLength; //Read In
	p.RefM=p.Wet.MolarDensity;
	p.RefDensity=p.Wet.MolarDensity*p.Wet.saturatedConcentration;
	p.RefViscosity=p.Wet.dynamicViscosity;
	p.RefMassDiffusivity=p.Wet.massDiffusivity;
	p.RefSurfaceTension=p.Wet.surfaceTension;
	p.RefConcentration=p.Wet.saturatedConcentration;
	p.RefTime=p.RefDensity*p.RefLength*p.RefLength/p.RefMassDiffusivity/p.RefM/p.RefConcentration;
	p.RefVelocity=p.RefLength/p.RefTime;
	p.RefPressure=p.RefDensity*p.RefVelocity*p.RefVelocity;
	p.RefMass=p.RefDensity*pow(p.RefLength, 3);
	p.RefMassFlux=p.RefMass/p.RefTime;
	p.Re=p.RefDensity*p.RefVelocity*p.RefLength/p.RefViscosity;
	p.Pe=p.RefVelocity*p.RefLength/p.RefMassDiffusivity;
	p.Ca=p.RefViscosity*p.RefVelocity/p.RefSurfaceTension;

	is>>p.FeatureLength>>p.FeatureVelocity>>p.FeaturePressure; //Read In
	p.FeatureDensity=p.Dry.density;
	p.FeatureViscosity=p.Dry.dynamicViscosity;
	if(p.FeatureVelocity!=0 && p.FeaturePressure==0) {
		p.FeaturePressure=p.FeatureDensity*p.FeatureVelocity*p.FeatureVelocity;
	} else if(p.FeatureVelocity==0 && p.FeaturePressure!=0) {
		p.FeatureVelocity=sqrt(p.FeaturePressure/p.FeatureDensity);
	} else {
		std::cout<<"Inputs about FeatureVelocity or FeaturePressure are wrong!"<<std::endl;
	}
	p.ReFeature=p.FeatureDensity*p.FeatureVelocity*p.FeatureLength/p.FeatureViscosity;
	p.PeFeature=p.FeatureVelocity*p.FeatureLength/p.RefMassDiffusivity;
	p.CaFeature=p.RefViscosity*p.FeatureVelocity/p.RefSurfaceTension;

	is>>p.EnvironmentConcentration; //Read In

	p.SaturatedConcentration=p.Wet.saturatedConcentration;
	p.LiquidConcentration=p.Wet.density/p.Wet.MolarDensity;

	is>>str; //Read In
	if(str=="SIMPLE") p.AlgorithmFlowfield=SIMPLE;
	else if(str=="SIMPLER") p.AlgorithmFlowfield=SIMPLER;
	else if(str=="SIMPLEC") p.AlgorithmFlowfield=SIMPLEC;
	else if(str=="OperatorSplitting") p.AlgorithmFlowfield=OperatorSplitting;
	else if(str=="PISO") p.AlgorithmFlowfield=PISO;
	else if(str=="NonOperatorSplitting") p.AlgorithmFlowfield=NonOperatorSplitting;
	else;
	is>>str; //Read In
	if(str=="Upwind") p.SchemeFlowfield=Upwind;
	else if(str=="Hybrid") p.SchemeFlowfield=Hybrid;
	else if(str=="PowerLaw") p.SchemeFlowfield=PowerLaw;
	else if(str=="QUICK") p.SchemeFlowfield=QUICK;
	else if(str=="HayaseQUICK") p.SchemeFlowfield=HayaseQUICK;
	else if(str=="TVD") p.SchemeFlowfield=TVD;
	else;
	is>>str; //Read In
	if(str=="Implicit") p.ImplicitFlowfield=true;
	else p.ImplicitFlowfield=false;

	is>>str; //Read In
	if(str=="SIMPLE") p.AlgorithmEvaporation=SIMPLE;
	else if(str=="SIMPLER") p.AlgorithmEvaporation=SIMPLER;
	else if(str=="SIMPLEC") p.AlgorithmEvaporation=SIMPLEC;
	else if(str=="OperatorSplitting") p.AlgorithmEvaporation=OperatorSplitting;
	else if(str=="PISO") p.AlgorithmEvaporation=PISO;
	else if(str=="NonOperatorSplitting") p.AlgorithmEvaporation=NonOperatorSplitting;
	else if(str=="ConvectionDiffusion") p.AlgorithmEvaporation=ConvectionDiffusion;
	else if(str=="DiffusionConvection") p.AlgorithmEvaporation=DiffusionConvection;
	else;
	is>>str; //Read In
	if(str=="Upwind") p.SchemeEvaporation=Upwind;
	else if(str=="Hybrid") p.SchemeEvaporation=Hybrid;
	else if(str=="PowerLaw") p.SchemeEvaporation=PowerLaw;
	else if(str=="QUICK") p.SchemeEvaporation=QUICK;
	else if(str=="HayaseQUICK") p.SchemeEvaporation=HayaseQUICK;
	else if(str=="TVD") p.SchemeEvaporation=TVD;
	else;
	is>>str; //Read In
	if(str=="Implicit") p.ImplicitEvaporation=true;
	else p.ImplicitEvaporation=false;

	is>>p.FlowfieldMaxStep>>p.EvaporationMaxStep>>p.RelativeConvergeCriteria>>p.FlowCriteria>>p.DryCriteria; //Read In

	is>>p.EvaporationSubStep;
	is>>p.FlowOutputFrequency>>p.DryOutputFrequency; //Read In

	is>>p.UVRelaxFactor>>p.PRelaxFactor>>p.CRelaxFactor; //Read In

	is>>p.FilmApproximationType>>p.FilmApproximationCoefficient; //Read In

	is>>p.InletFlowDeveloped; //Read In

	return is;
}

ostream & operator << (ostream & os, const physics_c & p) {
	os<<"Physics Model:\t";
	if(p.phyModel==incompressibleFLow) os<<"incompressibleFLow";
	else if(p.phyModel==isothermalEvaporation) os<<"isothermalEvaporation";
	else if(p.phyModel==poreNetwork) os<<"poreNetwork";

	os<<"Film Effect: "<<p.filmEffect<<"\tViscosity Effect: "<<p.viscosityEffect;

	os<<endl<<p.Dry<<p.Wet;

	os<<"Round Corner Radius: "<<p.RC<<endl;

	os<<"RefLength: "         <<p.RefLength         <<"[m]\t"
      <<"RefM: "              <<p.RefM              <<"[kg/mol]\t"
	  <<"RefDensity: "        <<p.RefDensity        <<"[kg/m^3]\t"
	  <<"RefViscosity: "      <<p.RefViscosity      <<"[kg/m-s]\t"
	  <<"RefMassDiffusivity: "<<p.RefMassDiffusivity<<"[m^2/s]\t"
	  <<"RefSurfaceTension: " <<p.RefSurfaceTension <<"[N/m]\t"
	  <<"RefConcentration: "  <<p.RefConcentration  <<"[mol/m^3]\t"<<endl;
	os<<"RefTime: "           <<p.RefTime           <<"[s]\t"
	  <<"RefVelocity: "       <<p.RefVelocity       <<"[m/s]\t"
	  <<"RefPressure: "       <<p.RefPressure       <<"[Pa]"
	  <<"RefMass: "           <<p.RefMass           <<"[kg]\t"
	  <<"RefMassFlux: "       <<p.RefMassFlux       <<"[kg/s]"<<endl;
	os<<"Reynolds Number: "<<p.Re<<"\tPeclet Number: "<<p.Pe<<"\tCapillary Number: "<<p.Ca<<endl;

	os<<"Feature Length: "            <<p.FeatureLength
	  <<"[m]\tFeature Density: "      <<p.FeatureDensity
	  <<"[kg/m^3]\tFeature Velocity: "<<p.FeatureVelocity
	  <<"[m/s]\tFeature Pressure: "   <<p.FeaturePressure<<"[Pa]"<<endl;
	os<<"Re Feature: "<<p.ReFeature<<"\tPe Feature: "<<p.PeFeature<<"\tCa Feature: "<<p.CaFeature<<endl;

	os<<"Environment Concentration: "<<p.EnvironmentConcentration<<endl;

	os<<"Saturated Concentration: "<<p.SaturatedConcentration<<endl;
	os<<"Liquid Concentration: "   <<p.LiquidConcentration   <<endl;

	os<<"Flowfield:\t";
	if(p.AlgorithmFlowfield==SIMPLE) os<<"SIMPLE";
	else if(p.AlgorithmFlowfield==SIMPLER) os<<"SIMPLER";
	else if(p.AlgorithmFlowfield==SIMPLEC) os<<"SIMPLEC";
	else if(p.AlgorithmFlowfield==OperatorSplitting) os<<"OperatorSplitting";
	else if(p.AlgorithmFlowfield==PISO) os<<"PISO";
	else if(p.AlgorithmFlowfield==NonOperatorSplitting) os<<"NonOperatorSplitting";
	else;
	os<<" + ";
	if(p.SchemeFlowfield==Upwind) os<<"Upwind";
	else if(p.SchemeFlowfield==Hybrid) os<<"Hybrid";
	else if(p.SchemeFlowfield==PowerLaw) os<<"PowerLaw";
	else if(p.SchemeFlowfield==QUICK) os<<"QUICK";
	else if(p.SchemeFlowfield==HayaseQUICK) os<<"HayaseQUICK";
	else if(p.SchemeFlowfield==TVD) os<<"TVD";
	else;
	os<<" + ";
	if(p.ImplicitFlowfield) os<<"Implicit";
	else os<<"Explicit";

	os<<endl<<"Evaporation:\t";
	if(p.AlgorithmEvaporation==SIMPLE) os<<"SIMPLE";
	else if(p.AlgorithmEvaporation==SIMPLER) os<<"SIMPLER";
	else if(p.AlgorithmEvaporation==SIMPLEC) os<<"SIMPLEC";
	else if(p.AlgorithmEvaporation==OperatorSplitting) os<<"OperatorSplitting";
	else if(p.AlgorithmEvaporation==PISO) os<<"PISO";
	else if(p.AlgorithmEvaporation==NonOperatorSplitting) os<<"NonOperatorSplitting";
	else if(p.AlgorithmEvaporation==ConvectionDiffusion) os<<"ConvectionDiffusion";
	else if(p.AlgorithmEvaporation==DiffusionConvection) os<<"DiffusionConvection";
	else;
	os<<" + ";
	if(p.SchemeEvaporation==Upwind) os<<"Upwind";
	else if(p.SchemeEvaporation==Hybrid) os<<"Hybrid";
	else if(p.SchemeEvaporation==PowerLaw) os<<"PowerLaw";
	else if(p.SchemeEvaporation==QUICK) os<<"QUICK";
	else if(p.SchemeEvaporation==HayaseQUICK) os<<"HayaseQUICK";
	else if(p.SchemeEvaporation==TVD) os<<"TVD";
	else;
	os<<" + ";
	if(p.ImplicitEvaporation) os<<"Implicit";
	else os<<"Explicit";
	os<<endl;

	os<<"FlowfieldMaxStep: "<<p.FlowfieldMaxStep<<"\tEvaporationMaxStep: "<<p.EvaporationMaxStep<<endl;
    os<<"FlowfieldConvergeCriteria: "<<p.RelativeConvergeCriteria<<"\tDryCriteria: "<<p.DryCriteria<<endl;

	os<<"EvaporationSubStep: "<<p.EvaporationSubStep<<endl;
	os<<"DryOutputFrequency: "<<p.DryOutputFrequency<<endl;
	os<<"UVRelaxFactor: "<<p.UVRelaxFactor<<"\tPRelaxFactor: "<<p.PRelaxFactor<<"\tCRelaxFactor: "<<p.CRelaxFactor<<endl;

	os<<"FilmApproximationType: "<<p.FilmApproximationType<<"\tFilmApproximationCoefficient: "<<p.FilmApproximationCoefficient<<endl;

	os<<"InletFlowDeveloped: "<<p.InletFlowDeveloped<<endl;

	return os;
}

bool physics_c::nondimensionalize() {
	RC/=RefLength;

	EnvironmentConcentration/=RefConcentration;

	SaturatedConcentration/=RefConcentration;
	LiquidConcentration/=RefConcentration;

	return true;
}
