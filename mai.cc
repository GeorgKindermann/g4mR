#include <vector>
#include <limits>
#include <cmath>
#include <valarray>

#include <Rcpp.h>

using namespace Rcpp;

namespace g4m {

  class mai {
  public:
    mai(const std::valarray<double>& c = {1.95866e-11,16.6774,5.30899,300,7.00127,0.0828,0.180158,50,0.002,0.85,20, 2, 0.5, 0.841, -22.131}  //Coefficients
	, const std::valarray<double>& c0 = {1.054684,1.099922,1.075196,0.980570,1.002155,1.044522,1.134524,1.073864,1.000548,1.070339,1.068615,1.086483,1.054495,1.036821,1.095323,1.008207,1.094867,1.031270,0.987843,1.035130,0.950606,1.074587,1.008381} //c0 Coefficients
	, const std::valarray<double>& t = std::valarray<double>(10., 24)  //Temperature of each month [deg C]
	, const std::valarray<double>& p = std::valarray<double>(70., 24)  //Precipitation [mm/30 Days]
	, const std::valarray<double>& r = std::valarray<double>(180., 24) //Radiation [Watt/m2]
	, double whc=100.    //Water holding capacity [mm]
	, double swr=0.  //Soil Water Regime, additional water from soil or irrigation [mm/month]
	, double co2=0.038   //CO2 concentration [volume %]
	, int soilType=0   //Soil type code
	, double altitude=0.    //Altitude [m]
	, unsigned char fNpp2mai=0 //Function type to convert npp to mai
	 , const std::valarray<double>& cNpp2mai = {0.35} //Coefficients to convert npp to mai
	//the followign temperatures set boundaries where this species can exist
	 , const std::valarray<double>& tMinJ = {0.} //Average annual temperature - minimum
	 , const std::valarray<double>& tMaxJ = {std::numeric_limits<double>::infinity()} //Average annual temperature - maximum
	 , const std::valarray<double>& pMinJ = {0.} //Annual precipitation - minimum
	 , const std::valarray<double>& pMaxJ = {std::numeric_limits<double>::infinity()} //Annual precipitation - maximum
	 , const std::valarray<double>& tMinM = {0.} //Month temperature - minimum
	 , const std::valarray<double>& tMaxM = {std::numeric_limits<double>::infinity()} //Month temperature - maximum
	 , const std::valarray<double>& pMinM = {0.} //Month precipitation - minimum
	 , const std::valarray<double>& pMaxM = {std::numeric_limits<double>::infinity()} //Month precipitation - maximum
	 , const std::valarray<double>& minNpp = {0.} //Minimum NPP
	, bool weatherIsDynamic = false
	);

    mai(const std::valarray<double>& c14t17 // = {2, 0.5, 0.841, -22.131}  //Coefficients c[14] to c[17]
	, const std::valarray<double>& c = {1.29460152,-0.09012495,0.17332495, 35, -1, 0.66170523, 2.8, 0.4138984, -54.4741443, -1.4, 1.155907e+00, 2.154028e-04, -3.733458e-01, 2.335792e-05}  //Coefficients
	, const std::valarray<double>& c0 = {0.06179233,0.06448844,0.07000044,0.07867775,0.06506758,0.08137664,0.06192571,0.07169721,0.07110523,0.06381677,0.05441309,0.06347873,0.07584091,0.07330926,0.05766713,0.07205265,0.05055277,0.06077571,0.07759581,0.05685617,0.06527024,0.05558023,0.06699292} //c0 Coefficients
	 , const std::valarray<double>& t = std::valarray<double>(10., 24)  //Temperature of each month [deg C]
	 , const std::valarray<double>& p = std::valarray<double>(70., 24)  //Precipitation [mm/month]
	, double whc=100.    //Water holding capacity [mm]
	, double swr=0.  //Soil Water Regime, additional water from soil or irrigation [mm/month]
	, double co2=0.038   //CO2 concentration [volume %]
	, int soilType=0   //Soil type code
	, double altitude=0.    //Altitude [m]
	, double latitude=0. //latitude of location [rad]
	, double soilWaterDecayRate=0.8  //lake of soil water
	, unsigned char fNpp2mai=0 //Function type to convert npp to mai
	 , const std::valarray<double>& cNpp2mai = {0.35} //Coefficients to convert npp to mai
	//the followign temperatures set boundaries where this species can exist
	 , const std::valarray<double>& tMinJ = {0.} //Average annual temperature - minimum
	 , const std::valarray<double>& tMaxJ = {std::numeric_limits<double>::infinity()} //Average annual temperature - maximum
	 , const std::valarray<double>& pMinJ = {0.} //Annual precipitation - minimum
	 , const std::valarray<double>& pMaxJ = {std::numeric_limits<double>::infinity()} //Annual precipitation - maximum
	 , const std::valarray<double>& tMinM = {0.} //Month temperature - minimum
	 , const std::valarray<double>& tMaxM = {std::numeric_limits<double>::infinity()} //Month temperature - maximum
	 , const std::valarray<double>& pMinM = {0.} //Month precipitation - minimum
	 , const std::valarray<double>& pMaxM = {std::numeric_limits<double>::infinity()} //Month precipitation - maximum
	 , const std::valarray<double>& minNpp = {0.} //Minimum NPP
	, bool weatherIsDynamic = false
	);

    double getNpp(unsigned int type=0, bool useMinNpp=false);
    std::valarray<double> getNpp(bool minNpp=false);
    std::valarray<double> getNpp(std::valarray<bool> dontNeed, bool useMinNpp=false);
    double getMai(unsigned int type=0, bool minNpp=false);
    std::valarray<double> getMai(bool minNpp=false);
    std::valarray<double> getMai(std::valarray<bool> dontNeed, bool useMinNpp=false);
    std::size_t setCoef(unsigned int type, const std::valarray<double>& c);
    std::size_t setCoefC0(unsigned int type, const std::valarray<double>& c);
    std::size_t setCoefC0(unsigned int type, std::valarray<double> c, double lo, double hi);
    bool setBoundaries(unsigned int type, const double& tMinJ, const double& tMaxJ, const double& pMinJ, const double& pMaxJ, const double& tMinM, const double& tMaxM, const double& pMinM, const double& pMaxM, const double& minNpp); 
    void setTemperature(const std::valarray<double>& t);
    void setPrecipitation(const std::valarray<double>& p);
    void setRadiation(const std::valarray<double>& r);
    void setCo2(const double& co2);
    void setSwr(const double& swr);
    void setWhc(const double& whc);
    void setAltitude(const double& altitude);
    void setSoilType(const int& soilType);
    void setLatitude(double latitude);
    double setSoilWaterDecayRate(double soilWaterDecayRate);
    unsigned int setcNpp2mai(const std::valarray<double>& acNpp2mai);
    bool setWeatherAsDynamic(bool weatherIsDynamic);
    bool testBoundaries(unsigned int type);
    std::valarray<bool> testBoundaries();
  private:
    int version;
    bool inputWasChanged;
    int updateSoilWater();
    void calcWalterLieth();
 
    std::valarray<double> c14t17;
    unsigned int numberOfTypes;
    const unsigned int nc; //Number of coefficients (new:15, old:14)
    const unsigned int nc0; //Number of c0 coefficients (23)
    std::valarray<double> c;  //18 Coefficients
    std::valarray<double> c0;  //23 c0 Coefficients
    std::valarray<double> tMinJ;
    std::valarray<double> tMaxJ;
    std::valarray<double> pMinJ;
    std::valarray<double> pMaxJ;
    std::valarray<double> tMinM;
    std::valarray<double> tMaxM;
    std::valarray<double> pMinM;
    std::valarray<double> pMaxM;
    std::valarray<double> minNpp;
    std::valarray<bool> outOfBoundaries;
    void resizeAll(unsigned int types);

    double walterLieth; //Walter lieth coefficient
    double whc;   //Water holding capacity 
    double swr;   //Soild water regime (Grundwassereinfluss)
    double co2;   //Co2 concentration
    double altitude;         //altitude
    int soilType; //Soil type
    std::valarray<double> t; //Temperature for each month of this year [12]
    std::valarray<double> p; //Precipitation for each month of this year [12]
    std::valarray<double> r; //Radiation for each month of this year [12]
    std::valarray<double> sw; //Soil water content [12]
    std::valarray<double> tp; //Temperature for each month of previous year [12]
    std::valarray<double> pp; //Precipitation for each month of previous year [12]
    std::valarray<double> rp; //Radiation for each month of previous year [12]
    bool weatherIsDynamic;
    double soilWaterDecayRate;
    double latitude;
    unsigned char fNpp2mai;
    std::valarray<double> cNpp2mai;
    
    double getNpp1(unsigned int type=0, bool useMinNpp=false);
    double getNpp2(unsigned int type=0, bool useMinNpp=false);
  };

    mai::mai(const std::valarray<double>& ac     //Coefficients
	   , const std::valarray<double>& ac0  //c0 Coefficients (soilType)
	   , const std::valarray<double>& at   //Temperature of each month [deg C]
	   , const std::valarray<double>& ap   //Precipitation [mm/30 Day]
	   , const std::valarray<double>& ar   //Radiation [Watt/m2]
	   , double awhc    //Water holding capacity [mm]
	   , double aswr    //Soil Water Regime, additional water from soil or irrigation [mm/month]
	   , double aco2     //CO2 concentration [volume %] (e.g. 0.038)
	   , int asoilType   //Soil type code
	   , double aaltitude    //Altitude [m]
	   , unsigned char afNpp2mai  //Function type to convert npp to mai
	   , const std::valarray<double>& acNpp2mai //Coefficients to convert npp to mai
	   //the followign temperatures set boundaries where this species can exist
	   , const std::valarray<double>& atMinJ //Average annual temperature - minimum
	   , const std::valarray<double>& atMaxJ //Average annual temperature - maximum
	   , const std::valarray<double>& apMinJ //Annual precipitation - minimum
	   , const std::valarray<double>& apMaxJ //Annual precipitation - maximum
	   , const std::valarray<double>& atMinM //Month temperature - minimum
	   , const std::valarray<double>& atMaxM //Month temperature - maximum
	   , const std::valarray<double>& apMinM //Month precipitation - minimum
	   , const std::valarray<double>& apMaxM //Month precipitation - maximum
	   , const std::valarray<double>& aminNpp //Minimum NPP
	   , bool aweatherIsDynamic
	   ) :
    version(2)
    , inputWasChanged(true)
    , nc(15)
    , nc0(23)
    , whc(awhc)
    , swr(aswr)
    , co2(aco2)
    , altitude(aaltitude)
    , soilType(asoilType)
    , weatherIsDynamic(aweatherIsDynamic)
    , fNpp2mai(afNpp2mai)
  {
    numberOfTypes = ac.size() / nc;
    c = ac;
    c0 = ac0;
    t = at[std::slice(at.size()-12,12,1)];
    p = ap[std::slice(ap.size()-12,12,1)];
    r = ar[std::slice(ar.size()-12,12,1)];
    tp = at[std::slice(0,12,1)];
    pp = ap[std::slice(0,12,1)];
    rp = ar[std::slice(0,12,1)];
    cNpp2mai = acNpp2mai;
    tMinJ = atMinJ;
    tMaxJ = atMaxJ;
    pMinJ = apMinJ;
    pMaxJ = apMaxJ;
    tMinM = atMinM;
    tMaxM = atMaxM;
    pMinM = apMinM;
    pMaxM = apMaxM;
    minNpp = aminNpp;
    sw.resize(12);
    outOfBoundaries.resize(numberOfTypes);
  }

    mai::mai(const std::valarray<double>& ac14t17  //Coefficients c[14] to c[17]
	   , const std::valarray<double>& ac     //Coefficients
	   , const std::valarray<double>& ac0  //c0 Coefficients
	   , const std::valarray<double>& at   //Temperature of each month [deg C]
	   , const std::valarray<double>& ap   //Precipitation [mm/month]
	   , double awhc    //Water holding capacity [mm]
	   , double aswr    //Soil Water Regime, additional water from soil or irrigation [mm/month]
	   , double aco2     //CO2 concentration [volume %] (e.g. 0.038)
	   , int asoilType   //Soil type code
	   , double aaltitude    //Altitude [m]
	   , double alatitude //latitude of location
	   , double asoilWaterDecayRate  //lake of soil water
	   , unsigned char afNpp2mai  //Function type to convert npp to mai
	   , const std::valarray<double>& acNpp2mai //Coefficients to convert npp to mai
	   //the followign temperatures set boundaries where this species can exist
	   , const std::valarray<double>& atMinJ //Average annual temperature - minimum
	   , const std::valarray<double>& atMaxJ //Average annual temperature - maximum
	   , const std::valarray<double>& apMinJ //Annual precipitation - minimum
	   , const std::valarray<double>& apMaxJ //Annual precipitation - maximum
	   , const std::valarray<double>& atMinM //Month temperature - minimum
	   , const std::valarray<double>& atMaxM //Month temperature - maximum
	   , const std::valarray<double>& apMinM //Month precipitation - minimum
	   , const std::valarray<double>& apMaxM //Month precipitation - maximum
	   , const std::valarray<double>& aminNpp //Minimum NPP
	   , bool aweatherIsDynamic
	   ) :
    version(1)
    , inputWasChanged(true)
    , nc(14)
    , nc0(23)
    , whc(awhc)
    , swr(aswr)
    , co2(aco2)
    , altitude(aaltitude)
    , soilType(asoilType)
    , weatherIsDynamic(aweatherIsDynamic)
    , soilWaterDecayRate(asoilWaterDecayRate)
    , latitude(alatitude*2.)
    , fNpp2mai(afNpp2mai)

  {
    c14t17 = ac14t17;
    numberOfTypes =ac.size() / nc;
    c = ac;
    c0 = ac0;
    t = at[std::slice(at.size()-12,12,1)];
    p = ap[std::slice(ap.size()-12,12,1)];
    tp = at[std::slice(0,12,1)];
    pp = ap[std::slice(0,12,1)];
    cNpp2mai = acNpp2mai;
    tMinJ = atMinJ;
    tMaxJ = atMaxJ;
    pMinJ = apMinJ;
    pMaxJ = apMaxJ;
    tMinM = atMinM;
    tMaxM = atMaxM;
    pMinM = apMinM;
    pMaxM = apMaxM;
    minNpp = aminNpp;
    sw.resize(12);
    outOfBoundaries.resize(numberOfTypes);
  }

  
  void mai::resizeAll(unsigned int types) {
    if(numberOfTypes < types) {
      numberOfTypes = types;
      std::valarray<double> tmp = c;
      c.resize((numberOfTypes)*nc);
      c[std::slice(0,tmp.size(),1)] = tmp;
      tmp = c0;
      c0.resize((numberOfTypes)*nc0);
      c0[std::slice(0,tmp.size(),1)] = tmp;
  tmp=tMinJ; tMinJ.resize(numberOfTypes); tMinJ[std::slice(0,tmp.size(),1)]=tmp;
  tmp=tMaxJ; tMaxJ.resize(numberOfTypes); tMaxJ[std::slice(0,tmp.size(),1)]=tmp;
  tmp=pMinJ; pMinJ.resize(numberOfTypes); pMinJ[std::slice(0,tmp.size(),1)]=tmp;
  tmp=pMaxJ; pMaxJ.resize(numberOfTypes); pMaxJ[std::slice(0,tmp.size(),1)]=tmp;
  tmp=tMinM; tMinM.resize(numberOfTypes); tMinM[std::slice(0,tmp.size(),1)]=tmp;
  tmp=tMaxM; tMaxM.resize(numberOfTypes); tMaxM[std::slice(0,tmp.size(),1)]=tmp;
  tmp=pMinM; pMinM.resize(numberOfTypes); pMinM[std::slice(0,tmp.size(),1)]=tmp;
  tmp=pMaxM; pMaxM.resize(numberOfTypes); pMaxM[std::slice(0,tmp.size(),1)]=tmp;
      tmp=minNpp; minNpp.resize(numberOfTypes);
      minNpp[std::slice(0,tmp.size(),1)]=tmp;
      std::valarray<bool> tmp2=outOfBoundaries;
      outOfBoundaries.resize(numberOfTypes);
      outOfBoundaries[std::slice(0,tmp2.size(),1)]=tmp2;
    }
  }

  std::size_t mai::setCoef(unsigned int type, const std::valarray<double>& ac) {
    resizeAll(type+1);
    c[std::slice(type*nc,nc,1)] = ac;
    return(c.size());
  }

  std::size_t mai::setCoefC0(unsigned int type, const std::valarray<double>& ac) {
    resizeAll(type+1);
    c0[std::slice(type*nc0,nc0,1)] = ac;
    return(c0.size());
  }

  std::size_t mai::setCoefC0(unsigned int type, std::valarray<double> ac, double lo, double hi) {
    for(unsigned int i=0; i<ac.size(); ++i) {
      if(ac[i] < lo) {ac[i] = lo;}
      if(ac[i] > hi) {ac[i] = hi;}
    }
    return(mai::setCoefC0(type, ac));
  }

  bool mai::setBoundaries(unsigned int type, const double& atMinJ, const double& atMaxJ, const double& apMinJ, const double& apMaxJ, const double& atMinM, const double& atMaxM, const double& apMinM, const double& apMaxM, const double& aminNpp) {
    resizeAll(type+1);
    tMinJ[type] = atMinJ;
    tMaxJ[type] = atMaxJ;
    pMinJ[type] = apMinJ;
    pMaxJ[type] = apMaxJ;
    tMinM[type] = atMinM;
    tMaxM[type] = atMaxM;
    pMinM[type] = apMinM;
    pMaxM[type] = apMaxM;
    minNpp[type] = aminNpp;
    return(testBoundaries(type));
  }

  void mai::calcWalterLieth() {
    walterLieth = c14t17[0] / (c14t17[1] + 1./(1. + exp(c14t17[2] + c14t17[3]*co2)));
  }

  int mai::updateSoilWater() {
    double tsw[15];
    for(int i=0; i<4; ++i) { //The previous Years
      double tmp;
      if(weatherIsDynamic) {
	tmp = tp[i+8];
	if(tmp < 0) {tmp=0;}
	tmp *= walterLieth;
	tmp = pp[i+8] - tmp;
      } else {
	tmp = t[i+8];
	if(tmp < 0) {tmp=0;}
	tmp *= walterLieth;
	tmp = p[i+8] - tmp;
      }
      tsw[i] = tmp;
    }
    for(int i=0; i<11; ++i) {
      double tmp = t[i];
      if(tmp < 0) {tmp=0;}
      tmp *= walterLieth;
      tmp = p[i] - tmp;
      tsw[i+4] = tmp;
    }
    for(int i=0; i<12; ++i) {
      double tmp = 0.;
      for(int j=0; j<4; ++j) {
	tmp += tsw[i+j];
	if(tmp < 0) {tmp=0;}
	tmp *= soilWaterDecayRate;
	if(tmp > whc) {tmp = whc;}
      }
      sw[i] = tmp + swr;
    }
    return(0);
  }

  double mai::getNpp1(unsigned int type, bool useMinNpp) {
    double ret=0.;
    if(inputWasChanged) {
      calcWalterLieth();
      updateSoilWater();
      inputWasChanged = false;
    }
    double t3 = c[6+type*nc] / (1. + exp(c[7+type*nc] + c[8+type*nc]*co2)) + c[9+type*nc];
    double days[12] = {31.,28.25,31.,30.,31.,30.,31.,31.,30.,31.,30.,31.};
    for(int month=0; month<12; ++month) {
      double t1 = c0[soilType+type*nc0]/(1. + std::exp(c[0+type*nc] + c[1+type*nc]*t[month])) - c0[soilType+type*nc0]*c[2+type*nc] - c0[soilType+type*nc0]/(1. + std::exp(c[3+type*nc]+c[4+type*nc]*t[month]));
      if(t1<0.) {t1=0.;}
      double t2 = 1. - 2./(1. + exp( (sw[month]+p[month]-std::max(0.,t[month])*walterLieth)/(std::max(1.,t[month]*c[5+type*nc]))));
      if(t2<0.) {t2=0.;}
      ret += t1 * t2 * t3 * days[month];
    }
    ret *= c[10+type*nc] + c[11+type*nc]*altitude + c[12+type*nc]*cos(latitude) + c[13+type*nc]*altitude*cos(latitude);
    if(ret<0.) {ret = 0.;}
    if(useMinNpp && ret < minNpp[type]) {ret = 0.;}
    return(ret);
   }
      
  double mai::getNpp2(unsigned int type, bool useMinNpp) {
    double ret=0.;
    static const double days[12] = {31.,28.25,31.,30.,31.,30.,31.,31.,30.,31.,30.,31.};
    double bodenwasser = 0.;
    for(int month=0; month<12; ++month) {
      double niederschlag = p[month];
      double temperatur = t[month];
      double verdunstungPot = 30. * exp(17.62*temperatur/(243.12 + temperatur));
      double evapotrans = tanh(c[6+type*nc]) * verdunstungPot;
      double interceptionPot = abs(c[10+type*nc]) + fmin(abs(c[7+type*nc]) * tanh(abs(c[8+type*nc]) * verdunstungPot), evapotrans);
      double interception = fmin(niederschlag, interceptionPot);
      double niederschlagBoden = niederschlag - interception;
      double verfuegbaresBodenwasser = bodenwasser * pow((bodenwasser / whc), c[9+type*nc]);
      double nutzbaresWasser = niederschlagBoden + verfuegbaresBodenwasser;
      double genutztesWasser = fmin(nutzbaresWasser, evapotrans - interception);
      bodenwasser += (niederschlagBoden - genutztesWasser) * days[month]/30.;
      if(bodenwasser > whc) {bodenwasser = whc;}
      if(bodenwasser < 0.) {bodenwasser = 0.;}
    }
    for(int month=0; month<12; ++month) {
      double niederschlag = p[month];
      double temperatur = t[month];
      double verdunstungPot = 30. * exp(17.62*temperatur/(243.12 + temperatur));
      double evapotrans = tanh(c[6+type*nc]) * verdunstungPot;
      double interceptionPot = abs(c[10+type*nc]) + fmin(abs(c[7+type*nc]) * tanh(abs(c[8+type*nc]) * verdunstungPot), evapotrans);
      double interception = fmin(niederschlag, interceptionPot);
      double niederschlagBoden = niederschlag - interception;
      double verfuegbaresBodenwasser = bodenwasser * pow((bodenwasser / whc), c[9+type*nc]);
      double nutzbaresWasser = niederschlagBoden + verfuegbaresBodenwasser;
      double genutztesWasser = fmin(nutzbaresWasser, evapotrans - interception);
      bodenwasser += (niederschlagBoden - genutztesWasser) * days[month]/30.;
      if(bodenwasser > whc) {bodenwasser = whc;}
      if(bodenwasser < 0.) {bodenwasser = 0.;}
      ret += days[month]
	* c0[soilType] * c[0+type*nc]*exp(-altitude/7990.)
	* r[month]
	* pow(fmax(0., (c[1+type*nc]+temperatur)), c[2+type*nc])
	* pow(fmax(0, 1. - verdunstungPot/tanh(c[5+type*nc]*nutzbaresWasser)/c[3+type*nc]), c[4+type*nc]);
    }
    if(ret<0.) {ret = 0.;}
    if(useMinNpp && ret < minNpp[type]) {ret = 0.;}
    return(ret);
  }

  double mai::getNpp(unsigned int type, bool useMinNpp) {
    double ret=0.;
    switch ( version ) {
    case 2:
      ret = getNpp2(type, useMinNpp);
      break;
    case 1:
      ret = getNpp1(type, useMinNpp);
      break;
    default:
      ret = getNpp2(type, useMinNpp);
      break;
    }
    return(ret);
  }

  
  std::valarray<double> mai::getNpp(bool useMinNpp) {
    std::valarray<double> ret(outOfBoundaries.size());
    for(unsigned int i=0; i<ret.size(); ++i) {ret[i] = getNpp(i, useMinNpp);}
    return(ret);
  }

  std::valarray<double> mai::getNpp(std::valarray<bool> dontNeed, bool useMinNpp) {
    std::valarray<double> ret(dontNeed.size());
    for(unsigned int i=0; i<ret.size(); ++i) {
      if(dontNeed[i]) {ret[i] = 0.;
      } else {ret[i] = getNpp(i, useMinNpp);}
    }
    return(ret);
  }

  bool mai::testBoundaries(unsigned int type) {
    outOfBoundaries[type] = false;
    if(t.min() < tMinM[type]) {outOfBoundaries[type] = true; return(outOfBoundaries[type]);}
    if(t.max() > tMaxM[type]) {outOfBoundaries[type] = true; return(outOfBoundaries[type]);}
    if(p.min() < pMinM[type]) {outOfBoundaries[type] = true; return(outOfBoundaries[type]);}
    if(p.max() > pMaxM[type]) {outOfBoundaries[type] = true; return(outOfBoundaries[type]);}
    double avg=t.sum()/12.;
    if(avg < tMinJ[type]) {outOfBoundaries[type] = true; return(outOfBoundaries[type]);}
    if(avg > tMaxJ[type]) {outOfBoundaries[type] = true; return(outOfBoundaries[type]);}
    avg=p.sum();
    if(avg < pMinJ[type]) {outOfBoundaries[type] = true; return(outOfBoundaries[type]);}
    if(avg > pMaxJ[type]) {outOfBoundaries[type] = true; return(outOfBoundaries[type]);}
    if(weatherIsDynamic) {
      if(tp.min() < tMinM[type]) {outOfBoundaries[type] = true; return(outOfBoundaries[type]);}
      if(tp.max() > tMaxM[type]) {outOfBoundaries[type] = true; return(outOfBoundaries[type]);}
      if(pp.min() < pMinM[type]) {outOfBoundaries[type] = true; return(outOfBoundaries[type]);}
      if(pp.max() > pMaxM[type]) {outOfBoundaries[type] = true; return(outOfBoundaries[type]);}
      double avg=tp.sum()/12.;
      if(avg < tMinJ[type]) {outOfBoundaries[type] = true; return(outOfBoundaries[type]);}
      if(avg > tMaxJ[type]) {outOfBoundaries[type] = true; return(outOfBoundaries[type]);}
      avg=pp.sum();
      if(avg < pMinJ[type]) {outOfBoundaries[type] = true; return(outOfBoundaries[type]);}
      if(avg > pMaxJ[type]) {outOfBoundaries[type] = true; return(outOfBoundaries[type]);}
    }
    return(outOfBoundaries[type]);
  }
  
  std::valarray<bool> mai::testBoundaries() {
    for(unsigned int i=0; i<outOfBoundaries.size(); ++i) {testBoundaries(i);}
    return(outOfBoundaries);
  }

  void mai::setTemperature(const std::valarray<double>& at) {
    inputWasChanged = true;
    if(weatherIsDynamic) {tp = t;}
    t = at;
  }

  void mai::setPrecipitation(const std::valarray<double>& ap) {
    inputWasChanged = true;
    if(weatherIsDynamic) {pp = p;}
    p = ap;
  }

  void mai::setRadiation(const std::valarray<double>& ar) {
    inputWasChanged = true;
    if(weatherIsDynamic) {rp = r;}
    r = ar;
  }

  void mai::setCo2(const double& aco2) {
    co2 = aco2;
    inputWasChanged = true;
  }

  void mai::setSwr(const double& aswr) {
    swr = aswr;
    inputWasChanged = true;
  }

  void mai::setWhc(const double& awhc) {
    whc = awhc;
    inputWasChanged = true;
  }

  void mai::setAltitude(const double& aaltitude) {
    altitude = aaltitude;
    inputWasChanged = true;
  }

  void mai::setSoilType(const int& asoilType) {
    soilType = asoilType;
  }
  
  void mai::setLatitude(double alatitude) {
    latitude = alatitude*2.;
    inputWasChanged = true;
  }

  double mai::getMai(unsigned int type, bool minNpp) {
    double npp = 0.;
    if(!testBoundaries(type)) {npp = getNpp(type, minNpp);}
    double mai = npp;
    if(fNpp2mai == 0) {
      mai *= cNpp2mai[0];
    } else {
      mai *= 0.35;
    }
    return(mai);
  }

  std::valarray<double> mai::getMai(bool minNpp) {
    std::valarray<double> ret(outOfBoundaries.size());
    for(unsigned int i=0; i<ret.size(); ++i) {ret[i] = getMai(i, minNpp);}
    return(ret);
  }

  std::valarray<double> mai::getMai(std::valarray<bool> dontNeed, bool minNpp) {
    std::valarray<double> ret(dontNeed.size());
    for(unsigned int i=0; i<ret.size(); ++i) {
      if(dontNeed[i]) {ret[i] = 0.;
      } else {ret[i] = getMai(i, minNpp);}
    }
    return(ret);
  }

  double mai::setSoilWaterDecayRate(double asoilWaterDecayRate) {
    soilWaterDecayRate = asoilWaterDecayRate;
    return(soilWaterDecayRate);
  }

  unsigned int mai::setcNpp2mai(const std::valarray<double>& acNpp2mai) {
    cNpp2mai = acNpp2mai;
    return(cNpp2mai.size());
  }

  bool mai::setWeatherAsDynamic(bool aweatherIsDynamic) {
    weatherIsDynamic = aweatherIsDynamic;
    return(weatherIsDynamic);
  }

}

g4m::mai mai;

// [[Rcpp::export]]
void maiInit() {
  mai.setSwr(0.8);
  mai.setcNpp2mai({1./3.});
  mai.setCo2(0.038);

  double BITT = 99.; //Boundary Increase Tropical Temperature

  //lc1er1 Nadel-Evergreen-Tropical
  mai.setCoef(0, {6.94357e-06,8,1.15427,300,2.38,0.025,0.8,50,0.002,0.85,20,0,0,0,0});
  mai.setCoefC0(0, {1.167273,1.181957,1.148051,1.144072,1.057347,1.183415,0.839456,0.998758,1.066749,0.998820,0.888920,1.132417,1.070195,1.146386,1.444292,1.338450,0.926202,1.164307,1.482320,1.161192,1.119969,1.010927,1.126816}, 0.7, 1.3);
  mai.setBoundaries(0,4.3 + 13.8,29.1 + 2.2 + BITT,200. - 0,6900. + 4000.,-3.3 + 14.6,33.5 + 2.3 + BITT,0,1300 + 2000,0.64);

  //lc2er1 Laub-Evergreen-Tropical
  mai.setCoef(1, {6.97374e-07,8,2.08263,300,2.36,0.025,0.6,50,0.002,0.85,20,0,0,0,0});
  mai.setCoefC0(1, {1.027236,0.930810,1.026614,1.058107,1.044301,1.088828,1.284784,0.997171,1.030403,0.930521,0.917611,1.071603,1.167654,1.063758,1.260789,1.139207,0.999577,0.966556,1.161625,0.990053,1.105658,1.056113,1.107417}, 0.7, 1.3);
  mai.setBoundaries(1,10.1 + 8.,28.3 + 3. + BITT,600. - 0,6900. + 4000.,7.2 + 4.1,33.6 + 2.2 + BITT,0,1300 + 2000,4.6);

  //lc3er1 Nadel-Deciduous-Tropical
mai.setCoef(2, {4.98896e-09,4,3.3326,300,3,0.04,0.8,50,0.002,0.85,20,0,0,0,0});
mai.setCoefC0(2, {1.25067,1.30666,1.42140,1.37749,1.19130,1.15079,0.42130,1.41935,1.35999,1.15866,1.22622,1.17439,1.46707,1.43959,1.89361,1.66058,1.01811,1.55586,1.65774,1.17218,1.36042,1.39463,1.66056}, 0.7, 1.3);
 mai.setBoundaries(2,8.3 + 9.8 - 3,28.5 + 1. + BITT,225. - 0,6900. + 4000.,1.9 + 9.4 - 3,35.2 + 1. + BITT,0,1300 + 2000,1.13);
 
  //lc4er1 Laub-Deciduous-Tropical
mai.setCoef(3, {3.17454e-07,8,1.98803,300,1.82,0.025,0.42,50,0.002,0.85,20,0,0,0,0});
mai.setCoefC0(3, {1.014203,0.655291,0.932896,1.681138,0.909431,0.724058,1.170861,0.730898,0.584047,0.946718,0.868518,0.866705,0.973624,1.075178,2.262815,0.867673,0.700617,0.794007,1.981404,1.001099,1.352024,1.284358,1.225052}, 0.7, 1.3);
 mai.setBoundaries(3,8.3 + 9.8 - 3,28.5 + 1. + BITT,225. - 0,6900. + 4000.,1.9 + 9.4 - 3,35.2 + 1. + BITT,0,1300 + 2000,1.13);
  
  //lc1er2 Nadel-Evergreen-Subtropical
mai.setCoef(4, {2.98529e-07,8,2.52942,300,5,0.025,0.8,50,0.002,0.85,20,0,0,0,0});
mai.setCoefC0(4, {1.084565,1.258374,1.544269,1.208804,1.048047,1.152335,0.798060,1.058165,1.172564,1.078057,1.483471,1.160240,1.351295,1.353240,1.784506,1.224317,1.034672,1.308062,1.714973,1.010345,1.210056,1.215814,0.994749}, 0.7, 1.3);
mai.setBoundaries(4,-3.1 + 15.5,22.9 - 4.2,280. - 0,6900. + 4000.,-11. + 12.2,30.0 - 0.1,0,1300 + 2000,1.2);

  //lc2er2 Laub-Evergreen-Subtropical
mai.setCoef(5, {.5*3.68661e-07,8,2.4796,300,3.4,0.025,0.6,50,0.002,0.85,20,0,0,0,0});
mai.setCoefC0(5, {1.0392105,0.9697758,1.1182245,1.0074064,0.9719184,0.9856556,0.8268068,0.9808688,1.0001164,0.9635926,1.0681479,1.0059678,1.2294432,1.0322276,1.0641367,0.9296392,1.0646439,1.1399712,1.0061634,1.1424225,0.9144326,1.1622985,1.0073817}, 0.7, 1.3);
mai.setBoundaries(5,3.4 + 9,23.7 - 5. + 2.5,290. - 0,6900. + 4000.,-5.8 + 7.,29.9 + 2.,0,1300 + 2000,0.85);

  //lc3er2 Nadel-Deciduous-Subtropical
mai.setCoef(6, {2.78094e-06,4,1.46148,300,3,0.04,0.8,50,0.002,0.85,20,0,0,0,0});
mai.setCoefC0(6, {1.33731,1.21174,2.09552,1.45462,1.11848,1.78982,0.44737,1.14455,1.33598,1.16703,1.59218,1.89613,1.54779,1.33450,1.61385,1.23386,1.15095,1.48598,1.50806,1.13044,1.04776,1.27060,1.15392}, 0.7, 1.3);
mai.setBoundaries(6,-3.1 + 15.5,22.9 - 4.2,280. - 0,6900. + 4000.,-11. + 12.2,30.0 - 0.1,0,1300 + 2000,1.2);
 
  //lc4er2 Laub-Deciduous-Subtropical
mai.setCoef(7, {5.22354e-08,3.81171,2.94545,300,3.18175,0.0232828,0.809172,50,0.002,0.851749,20,0,0,0,0});
mai.setCoefC0(7, {0.983435,1.216336,1.208377,1.164035,1.127213,1.427934,0.427887,1.041255,0.899436,0.951187,1.505357,1.606139,1.512237,1.196615,1.284658,1.063652,1.121445,1.179302,1.375756,0.836196,1.106857,1.427041,1.055856}, 0.7, 1.3);
mai.setBoundaries(7,1.9 + 10.5,23.8 - 5.1,190. - 0,6900. + 4000.,-10.5 + 11.7,29.9 - 0.,1,1300 + 2000,0.15);

  //lc1er3 Nadel-Evergreen-Temperate
mai.setCoef(8, {1.95866e-11,16.6774,5.30899,300,7.00127,0.0828,0.180158,50,0.002,0.85,20,0,0,0,0});
mai.setCoefC0(8, {1.054684,1.099922,1.075196,0.980570,1.002155,1.044522,1.134524,1.073864,1.000548,1.070339,1.068615,1.086483,1.054495,1.036821,1.095323,1.008207,1.094867,1.031270,0.987843,1.035130,0.950606,1.074587,1.008381}, 0.7, 1.3);
 mai.setBoundaries(8,-4.9 + 5,14.6 - 2. + 3,220. - 0,6900. + 4000.,-27.2 + 0.,25.1 + 4.5,0,1300 + 2000,1.18);

  //lc2er3 Laub-Evergreen-Temperate
mai.setCoef(9, {3.15e-09,9.35918,3.66485,300,4.64592,0.06,0.26,50,0.002,0.85,20,0,0,0,0});
mai.setCoefC0(9, {1.0118956,1.0615462,1.0779899,1.0120382,0.9627487,0.9898199,1.0791497,1.0310722,0.9754276,0.8438813,1.1233655,1.1068356,1.0901310,1.0101304,1.0565290,1.0779955,1.0621784,1.0007144,0.9295641,1.0472222,1.0161400,1.0173837,1.0171552}, 0.7, 1.3);
mai.setBoundaries(9,-5.8 + 5.9,15.9 - 3.3 + 2.,330. - 0,6900. + 4000.,-25.7 + 0.,25.1 - 0.,1,1300 + 2000,3.18);

  //lc3er3 Nadel-Deciduous-Temperate
mai.setCoef(10, {1.50349e-07,7.07274,2.54163,300,4.00133,0.06,0.25,50,0.002,0.85,20,0,0,0,0});
mai.setCoefC0(10, {0.862190,1.145883,1.135764,0.950486,1.102457,1.030473,0.956814,0.975896,1.170444,1.003162,1.169161,1.053181,1.126635,1.232963,0.864202,1.179097,0.774995,0.995233,1.321486,1.032954,1.041349,1.029158,1.056834}, 0.7, 1.3);
mai.setBoundaries(10,-9.0 + 9.1,15.4 - 2.8,100. + 150,6900. + 4000.,-32.6 + 0,27.9 - 0.,1,1300 + 2000,0.19);

  //lc4er3 Laub-Deciduous-Temperate
mai.setCoef(11, {6.30082e-09,9.35918,3.66485,300,4.64592,0.06,0.26,50,0.002,0.85,20,0,0,0,0});
mai.setCoefC0(11, {0.8426849,1.1336508,0.9219486,1.0428545,1.1223109,0.8339236,1.0641686,0.5452793,1.0681361,1.0280533,1.0013517,0.6224940,1.2335694,1.1491972,1.0465842,0.9879109,1.1978397,0.9430111,1.0909037,1.1027924,1.0124663,1.2082930,1.0756028}, 0.7, 1.3);
mai.setBoundaries(11,-4.9 + 5,14.6 - 2. + 3,305. - 0,6900. + 4000.,-29.7 + 0.,27.1 + 2.5,2,1300 + 2000,1.74);

  //lc1er4 Nadel-Evergreen-Boreal
mai.setCoef(12, {1.36116e-07,6.7432,2.82955,300,4.57307,0.06,0.50,50,0.002,0.85,20,0,0,0,0});
mai.setCoefC0(12, {1.035657,1.036384,1.032846,1.030142,1.055238,1.052846,1.045243,1.024852,1.002326,0.998242,1.081536,1.108310,1.013428,1.036258,1.054235,1.041855,1.031564,1.027304,1.115472,1.011239,1.016438,1.006963,1.013799}, 0.7, 1.3);
mai.setBoundaries(12,-10.3 - 8,6.6 - 6.,175. - 0,6900. + 4000.,-39.9 - 8.,19.0 + 2.,2-2,1300 + 2000,0.53);
 
  //lc2er4 Laub-Evergreen-Boreal
mai.setCoef(13, {6.94485e-11,10.3926,5.18179,300,8.21249,0.06,0.2,50,0.002,0.85,20,0,0,0,0});
mai.setCoefC0(13, {0.93566,1.04938,0.99955,1.17027,0.92923,0.90336,1.00938,1.02759,0.95376,0.90328,1.00736,0.99617,1.14578,1.09325,1.25368,1.11511,0.92586,0.99376,1.10367,0.90914,0.96348,0.94270,1.41604}, 0.7, 1.3);
mai.setBoundaries(13,-11.4 - 2.9,7.9 - 7.3,190. - 0,6900. + 4000.,-40.9 + 1.,21.4 + 1.,2,1300 + 2000,0.47);
 
  //lc3er4 Nadel-Deciduous-Boreal
mai.setCoef(14, {7.54653e-07,3.39758,2.39893,300,5.10519,0.06,0.59,50,0.002,0.85,20,0,0,0,0});
mai.setCoefC0(14, {0.997382,1.164807,1.137406,1.157535,1.025785,1.016309,0.922417,1.078434,1.033414,0.960146,1.018151,1.064526,1.083873,1.061930,1.108677,1.053499,1.025711,0.942647,1.277620,1.014716,1.114794,1.009307,1.115872}, 0.7, 1.3);
mai.setBoundaries(14,-10.3 - 8,6.6 - 6.,215. - 0,6900. + 4000.,-39.9 - 8.,19.6 + 2.,2-2,1300 + 2000,1.17);

  //lc3er4 Laub-Deciduous-Boreal
mai.setCoef(15, {1.38897e-10,10.3926,5.18179,300,8.21249,0.06,0.2,50,0.002,0.85,20,0,0,0,0});
mai.setCoefC0(15, {0.816807,1.171890,1.106970,1.007987,0.960848,0.923770,1.029693,1.114420,1.059554,1.006887,1.168652,0.913977,1.082629,1.532929,1.053006,1.072970,0.989127,0.994642,1.297751,0.959316,1.233419,0.953433,1.041302}, 0.7, 1.3);
mai.setBoundaries(15,-11.4 - 2.9,7.9 - 7.3,190. - 0,6900. + 4000.,-40.9 + 1.,21.4 + 1.,2,1300 + 2000,0.47);
}

// [[Rcpp::export]]
void maiSetSwr(double x) {mai.setSwr(x);}

// [[Rcpp::export]]
void maiSetCo2(double x) {mai.setCo2(x);}

// [[Rcpp::export]]
void maiSetTemperature(std::vector<double> x) {
  mai.setTemperature(std::valarray<double> (x.data(), x.size()));
}

// [[Rcpp::export]]
void maiSetPrecipitation(std::vector<double> x) {
  mai.setPrecipitation(std::valarray<double> (x.data(), x.size()));
}

// [[Rcpp::export]]
void maiSetRadiation(std::vector<double> x) {
  mai.setRadiation(std::valarray<double> (x.data(), x.size()));
}

// [[Rcpp::export]]
void maiSetWhc(double x) {mai.setWhc(x);}

// [[Rcpp::export]]
void maiSetAltitude(double x) {mai.setAltitude(x);}

// [[Rcpp::export]]
void maiSetSoilType(int x) {mai.setSoilType(x);}

// [[Rcpp::export]]
std::vector<double> maiGet() {
  std::valarray<double> rmai = mai.getMai(mai.testBoundaries(), true);
  return(std::vector<double> (std::begin(rmai), std::end(rmai)));
}
