/* RENYI ENTROPY -an addition to Nick Collins's SCMIRUGen, SpectralEntropy 
 * hayakzan 2021
 */


//#define SC_DARWIN
#include "SC_PlugIn.h"
#include "FFT_UGens.h"


InterfaceTable *ft;


//convert FFT
struct RenyiEntropy : public Unit
{
  //float entropy_;

  int numbands;
  int fftsize;
  int alpha;
  int* bandindices;
  float* intensities;
  float* entropies;

};

extern "C" {  

static void RenyiEntropy_next_k(RenyiEntropy *unit, int inNumSamples);
static void RenyiEntropy_Ctor(RenyiEntropy* unit);
static void RenyiEntropy_Dtor(RenyiEntropy* unit);

}



void RenyiEntropy_Ctor( RenyiEntropy* unit ) {

  int i, j;

  unit->fftsize = IN0(1);
  unit->numbands = IN0(2);
  unit->alpha = IN0(3);

  int numbins = unit->fftsize/2; //won't use actual Nyquist bin in this UGen-but this is Nyquist bin anyway...

  int split = numbins/(unit->numbands);

  if(split<1) {

    split = 1;
    unit->numbands = numbins;
  }


  //will include guard element at top
  unit->bandindices = (int *)RTAlloc(unit->mWorld, sizeof(int)*(unit->numbands+1));
  unit->entropies = (float *)RTAlloc(unit->mWorld, sizeof(float)*unit->numbands);
  unit->intensities = (float *)RTAlloc(unit->mWorld, sizeof(float)*numbins);

  for (i=0; i<unit->numbands; ++i) {

    unit->entropies[i] = 0.0f;

    unit->bandindices[i] = split*i;

  }

  //guard can be one above actual final array slot index since always use less than in code below
  unit->bandindices[unit->numbands] = numbins; //Nyquist position

  SETCALC(RenyiEntropy_next_k);


}


void RenyiEntropy_Dtor(RenyiEntropy *unit)
{

  RTFree(unit->mWorld, unit->bandindices);
  RTFree(unit->mWorld, unit->entropies);
  RTFree(unit->mWorld, unit->intensities);

}



void RenyiEntropy_next_k( RenyiEntropy *unit, int inNumSamples ) {

  int i,j;

  int numbands = unit->numbands;
  int alpha = unit->alpha;
  int * bandindices = unit->bandindices;
  float * entropies = unit->entropies;
  float * intensities = unit->intensities;

  //if input is legitimate buffer number:
  float fbufnum = IN0(0);

  //next FFT bufffer ready, update
  //assuming at this point that buffer precalculated for any resampling
  if (fbufnum > -0.01f) {

    int ibufnum = (uint32)fbufnum;

    World *world = unit->mWorld;
    SndBuf *buf;

    if (ibufnum >= world->mNumSndBufs) {
      int localBufNum = ibufnum - world->mNumSndBufs;
      Graph *parent = unit->mParent;
      if(localBufNum <= parent->localBufNum) {
	buf = parent->mLocalSndBufs + localBufNum;
      } else {
	buf = world->mSndBufs;
      }
    } else {
      buf = world->mSndBufs + ibufnum;
    }

    if(unit->fftsize == buf->frames) {

      //make sure in real and imag form
      //SCComplexBuf * complexbuf = ToComplexApx(buf);

      float * data = (float *)ToComplexApx(buf);

      //float * data= buf->data;

      float real, imag;
      float intensity;

      data[1] = 0.0f; //avoid issues with dc, nyquist packed together, just want dc here

      for  (j=0; j<numbands; ++j) {

	int start = bandindices[j];
	int end = bandindices[j+1];

	float max = 0.0f;
	float entropysum = 0.0f;
	float entropyval = 0.0f;

		//less than because of guard elements ???
		for  (i=start; i<end; ++i) {

		  int index = 2*i;
		  real= data[index];
		  imag= data[index+1];
		  intensity = (real*real) + (imag*imag);

		  intensities[i] = intensity;
		}

	for  (i=start; i<end; ++i) {

	    float prob = intensities[i]; //removed *max

	    //again, "intensities" is the combination of real and imaginary numbers.
	    //so prob must be normalized intensities...let's check that out.
	    
	    //negative worked in via -=
	    //Shannon:
          //	    if(prob>0.0f) entropysum -= prob * log2(prob); //-p(i)*log2(p(i))

	    
	    //Rényi: 
		if(prob>0.0f) { 
	 		entropysum += pow(prob, alpha); 
		}	   
	}
	
	entropyval = log2(entropysum);
	entropyval /= (1.0f-alpha);
	
	entropies[j] = entropyval; 

      }


    }

  }

    for (i=0; i<numbands; ++i)
		OUT0(i) = entropies[i];
    
}




PluginLoad(RenyiEntropy) {

  init_SCComplex(inTable);

  ft = inTable;

  DefineDtorCantAliasUnit(RenyiEntropy);
  //DefineDtorUnit(RenyiEntropy)
}




