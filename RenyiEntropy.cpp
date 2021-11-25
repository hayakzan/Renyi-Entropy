/* RENYI ENTROPY -adds a feature to Nick Collins's Sensory Dissonance SCMIR library for SuperCollider
 * https://github.com/sicklincoln/SCMIR
 * based on: https://digitalassets.lib.berkeley.edu/math/ucb/text/math_s4_v1_article-27.pdf
 */


//#define SC_DARWIN
#include "SC_PlugIn.h"
#include "FFT_UGens.h"


InterfaceTable *ft;


struct SensoryDissonance : public Unit
{
	int fftsize;
	int topbin;
	int frequencyperbin;

	float dissonance;

	int maxnumpeaks;
	float peakthreshold;

	float * peakfreqs;
	float * peakamps;

	float norm;
	int clamp;

	int initfftsize;
};


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

void SensoryDissonance_next_k(SensoryDissonance *unit, int inNumSamples);
void SensoryDissonance_Ctor(SensoryDissonance* unit);
void SensoryDissonance_Dtor(SensoryDissonance* unit);

static void RenyiEntropy_next_k(RenyiEntropy *unit, int inNumSamples);
static void RenyiEntropy_Ctor(RenyiEntropy* unit);
static void RenyiEntropy_Dtor(RenyiEntropy* unit);

}

void SensoryDissonance_Ctor( SensoryDissonance* unit ) {

	//int i, j;

	unit->initfftsize = 0; //must defer this till have a buffer to check

	unit->maxnumpeaks = IN0(1); //100;
	unit->peakthreshold = IN0(2);
	unit->peakfreqs = (float *)RTAlloc(unit->mWorld, sizeof(float)*unit->maxnumpeaks);
	unit->peakamps = (float *)RTAlloc(unit->mWorld, sizeof(float)*unit->maxnumpeaks);

	unit->norm = IN0(3); //0.01/unit->maxnumpeaks; //unit->fftsize;

	unit->clamp = IN0(4);

	SETCALC(SensoryDissonance_next_k);


}

void SensoryDissonance_Dtor(SensoryDissonance *unit)
{

	RTFree(unit->mWorld, unit->peakfreqs);
	RTFree(unit->mWorld, unit->peakamps);

}

//NEXT: destructor then next function using octaves and divisions slots, with appropriate power calc from fft data as go


void SensoryDissonance_next_k( SensoryDissonance *unit, int inNumSamples ) {

	//int i, j;

	//float *input = IN(0);

	//int numSamples = unit->mWorld->mFullRate.mBufLength;

	//if input is legitimate buffer number:
	float fbufnum = IN0(0);

	//next FFT buffer ready, update
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


		if(unit->initfftsize ==0) {

			double sr = unit->mWorld->mFullRate.mSampleRate; //never trust SAMPLERATE, gives UGens output rate, not audio rate
			//float nyquist = sr*0.5;

			unit->fftsize = buf->frames; //IN0(1);
			//printf("check fftsize %d \n",unit->fftsize);
			unit->topbin= unit->fftsize*0.25;

			unit->frequencyperbin = sr / unit->fftsize;

			unit->initfftsize = 1;
		}

		//make sure in real and imag form
		//SCComplexBuf * complexbuf = ToComplexApx(buf);

		float * data= (float *)ToComplexApx(buf);

		//float * data= buf->data;

		//int numindices= unit->numindices;

		float * peakfreqs= unit->peakfreqs;
		float * peakamps= unit->peakamps;

		float real, imag;
		int index;

		int numpeaks = 0;
		int maxnumpeaks = unit->maxnumpeaks;

		float intensity;
		float position;

		float threshold = unit->peakthreshold;

		//create powerspectrum

		float prev=0.0, now=0.0, next=0.0;

		float frequencyperbin = unit->frequencyperbin;

		//float totalpeakpower = 0.0f;
		float temp1, refinement;

		for (int j=1; j<=unit->topbin; ++j) {

				index = 2*j;
				real= data[index];
				imag= data[index+1];
				intensity = (real*real) + (imag*imag);
//

			next= intensity;

			if(j>=3) {

				//hunt for peaks

				//look for peak by scoring within +-3
				//assume peak must be centrally greater than 60dB say

				//powertest
				//minpeakdB was 60

				if (now>threshold)  {

					//y1= powerspectrum[i-1];
					//				//y2= valuenow;
					//				y3= powerspectrum[i+1];
					//
					if ((now>prev) && (now>next)) {

						//second peak condition; sum of second differences must be positive
						//NCfloat testsum= (valuenow - powerspectrum[i-2]) + (valuenow - powerspectrum[i+2]);

						//if (testsum>0.0) {

						//refine estimate of peak using quadratic function
						//see workbook 28th Jan 2010

						temp1= prev+next-(2*now);

						if (fabs(temp1)>0.00001) {
							position=(prev-next)/(2*temp1);

							//running quadratic formula
							refinement = (0.5*temp1*(position*position)) + (0.5*(next-prev)*position) + now;
							//refinement= y2 -  (((y3-y1)^2)/(8*temp1));

						} else {
							//degenerate straight line case; shouldn't occur
							//since require greater than for peak, not equality

							position=0.0; //may as well take centre

							//bettervalue= max([y1,y2,y3]); %straight line through them, find max

							refinement= now; //must be max for else would have picked another one in previous calculation! %max([y1,y2,y3]);

						}

						//correct??????????????????????????????
						peakfreqs[numpeaks] = (j-1+position)*frequencyperbin; //frequencyconversion;
						//printf("peakfrequencies %d is %f from i %d position %f freqperbin %f \n", numpeaks,peakfrequencies[numpeaks],i, position, frequencyperbin);

						peakamps[numpeaks] = sqrt(refinement); //Sethares formula requires amplitudes
						//totalpeakpower += refinement;

						//cout << " peak " << numpeaks << " " << peakfrequencies[numpeaks] << " " << refinement << " " ;

						++numpeaks;

						//}

					}

				}

				//test against maxnumberpeaks
				if ( numpeaks == maxnumpeaks )
					break;

			}

			prev = now; now=next;

		}


		//now have list of peaks: calculate total dissonance:

		//iterate through peaks, matching each to min of next 10, and no more than octave, using Sethares p. 58 CMJ article

		float dissonancesum = 0.0;

		float f1, v1, f2, v2;
		float d;
		float diff; //, minf;
		float s, a, b;
		float octave;

		for (int i=0; i<(numpeaks-1); ++i) {

			f1 = peakfreqs[i];
			v1 = peakamps[i];
			s = 0.24f/(0.21f*f1+19.f); //constant needed as denominator in formula
			a = -3.5f*s;
			b= -5.75f*s;

			octave = 2.0f*f1;

			for (int k=i+1; k<sc_min(i+20,numpeaks); ++k) {

				f2 = peakfreqs[k];
				v2 = peakamps[k];

				if(f2>octave) break; //shortcut escape if separated by more than an octave

				diff = f2-f1; //no need for fabs, f2>f1
				//minf =  //always f1 lower


				d = v1*v2*(exp(a*diff) - exp(b*diff));

				dissonancesum += d;
			}

		}

		unit->dissonance = sc_min(unit->clamp,dissonancesum*unit->norm); //numpeaks; //dissonancesum;  //divide by fftsize as compensation for amplitudes via FFT

	}

	//OUT0(i) = unit->dissonance;
	OUT0(0) = unit->dissonance;

}



// Renyi Entropy
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

  //next FFT buffer ready, update
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

      float * data= (float *)ToComplexApx(buf);

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

	  if(intensity>max) {

	    max = intensity;

	  }

	}

	if(max>0.0f) {

	  max = 1.0f/max; //will be used as straight multiplier in calculation now

	  for  (i=start; i<end; ++i) {

	    float prob = intensities[i] * max;

	    //again, "intensities" is the combination of real and imaginary numbers.
	    //so prob must be normalized intensities...let's check that out.

	    //Shannon:
          //	    if(prob>0.0f) entropysum -= prob * log2(prob); //-p(i)*log2(p(i))


	    //Rényi:
		entropysum += pow(prob, alpha); //base prob, exp alpha

	  }
	  entropyval =  log2(entropysum) / (1.0f-alpha);
	  entropies[i] = entropyval; 

	} else
	  entropies[i] = 0.0f;

      }

    }

  }

    for (i=0; i<numbands; ++i)
		OUT0(i) = entropies[i];

}

PluginLoad(RenyiEntropy) {

  init_SCComplex(inTable);

  ft = inTable;

DefineDtorCantAliasUnit(SensoryDissonance);
  DefineDtorCantAliasUnit(RenyiEntropy);
  //DefineDtorUnit(RenyiEntropy)
}
