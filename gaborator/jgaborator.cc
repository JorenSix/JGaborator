#include "jgaborator.h"
#include <math.h>


gaborator::parameters* paramsRef= NULL;
gaborator::analyzer<float>* analyzerRef = NULL;
gaborator::coefs<float>* coefsRef = NULL;

int64_t t_in = 0;
int min_band = 0;
int sample_rate = 0;
size_t analysis_support = 0;
//analysis output array, to be reused;
//allocation should depend on blocksize, now a fixed
//'big enough' number is picked FIXME:
jfloat output_array[131072*32];
jfloatArray outJNIArray;

JNIEXPORT jint JNICALL Java_be_ugent_jgaborator_JGaborator_initialize(JNIEnv * env, jobject object, jint blocksize, jdouble fs, jint bands_per_octave, jdouble ff_min , jdouble ff_ref, jdouble ff_max, jdouble overlap, jdouble max_error){

	paramsRef = new  gaborator::parameters(bands_per_octave, ff_min / fs, ff_ref / fs, overlap, max_error);
	analyzerRef = new gaborator::analyzer<float>(*paramsRef);
	coefsRef = new gaborator::coefs<float>(*analyzerRef);

	analysis_support = ceil(analyzerRef->analysis_support());

   //converts frequency (ff_max) in hertz to the number of bands above the min frequency
   //the ceil is used to end up at a full band 
   int interesting_bands = ceil(bands_per_octave * log(ff_max/ff_min)/log(2.0f));

   //since bands are ordered from high to low we are only interested in lower bands:
   //fs/2.0 is the nyquist frequency
   int total_bands = ceil(bands_per_octave * log(fs/2.0/ff_min)/log(2.0f));

   min_band = total_bands - interesting_bands;

   sample_rate = fs;


	int64_t t_in = 0;

	return analysis_support;
}



/*
 * Class:     be_ugent_jgaborator_JGaborator
 * Method:    analyse
 * Signature: ([F)[F
 */
JNIEXPORT jfloatArray JNICALL Java_be_ugent_jgaborator_JGaborator_analyse(JNIEnv * env, jobject, jfloatArray audio_block ){

   // Step 1: Convert the incoming JNI jintarray to C's jfloat[]
   jfloat *audio_block_c_array = env->GetFloatArrayElements(audio_block, NULL);
   if (NULL == audio_block_c_array) return NULL;
   jsize blocksize = env->GetArrayLength(audio_block);

   std::vector<float> buf(audio_block_c_array,audio_block_c_array + blocksize);

   //printf("%d\n", (int) blocksize);
   //printf("%d\n", (int) paramsRef->bands_per_octave);

   int output_index = 0;
   
   analyzerRef->analyze(buf.data(), t_in, t_in + blocksize, *coefsRef);
   apply(
            *analyzerRef, *coefsRef,
            [&](std::complex<float> coef, int band, int64_t audioSampleIndex ) {
                 //coef.norm = -coef;
               //ignores everything above the max_band
               if(band >= min_band){
                 //printf("%f %d %ld\n",std::abs(coef),band,audioSampleIndex);
                  output_array[output_index++] = band;
                  output_array[output_index++] = audioSampleIndex;
                  output_array[output_index++] = std::abs(coef);
               }
            },
            t_in - analysis_support,
            t_in - analysis_support + blocksize);

   t_in += blocksize;

   int64_t t_out = t_in - analysis_support;
   
   forget_before(*analyzerRef, *coefsRef, t_out - blocksize );
  
   //release audio block memory resources
   env->ReleaseFloatArrayElements(audio_block, audio_block_c_array, 0);

	// Convert the C's Native jfloat[] to JNI jfloatarray, and return
   outJNIArray = env->NewFloatArray(output_index);  // allocate
   //if (NULL == outJNIArray) return NULL;
   env->SetFloatArrayRegion(outJNIArray, 0 , output_index, output_array);  // copy
   return outJNIArray;
}

JNIEXPORT jfloatArray JNICALL Java_be_ugent_jgaborator_JGaborator_bandcenters(JNIEnv * env, jobject){
   int max_band = analyzerRef->bandpass_bands_end();
   float band_centers[max_band+1];

   for(int i = 0 ; i < max_band ; i++){
      if(i<min_band){
         band_centers[i]=-1;
      }else{
         band_centers[i]=analyzerRef->band_ff(i) * sample_rate;
      }
   }

   // Convert the C's Native jfloat[] to JNI jfloatarray, and return
   jfloatArray outJNIArray = env->NewFloatArray(max_band+1);  // allocate
   //if (NULL == outJNIArray) return NULL;
   env->SetFloatArrayRegion(outJNIArray, 0 , max_band+1, band_centers);  // copy
   return outJNIArray;
}

/*
 * Class:     be_ugent_jgaborator_JGaborator
 * Method:    release
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_be_ugent_jgaborator_JGaborator_release(JNIEnv *, jobject){
	//cleanup memory
	delete analyzerRef;
	delete coefsRef;
	delete paramsRef;

	//reset variables for reuse
	t_in = 0;
	analysis_support = 0;
	min_band = 0;
 	sample_rate = 0;
}

int main(){
   return 0;
}
