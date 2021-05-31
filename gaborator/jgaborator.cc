#include "jgaborator.h"
#include "gaborator-1.2/gaborator/gaborator.h"

#include <unordered_map>
#include <math.h>

const int C_ARRAY_SIZE = 300000 * 2;
struct GaboratorState {   
	gaborator::parameters* paramsRef;
	gaborator::analyzer<float>* analyzerRef;
	gaborator::coefs<float>* coefsRef;
	
	int64_t t_in;
	int min_band;
	int sample_rate;
	int64_t anal_support;

   jfloat *cArray;
};


std::unordered_map<uintptr_t, uintptr_t> stateMap;
std::mutex stateMutex;


JNIEXPORT jint JNICALL Java_be_ugent_jgaborator_JGaborator_initialize(JNIEnv * env, jobject object, jint blocksize, jdouble fs, jint bands_per_octave, jdouble ff_min , jdouble ff_ref, jdouble ff_max, jdouble overlap, jdouble max_error){
   
   std::unique_lock<std::mutex> lck (stateMutex);
   GaboratorState * state =  new GaboratorState();
   uintptr_t env_addresss = reinterpret_cast<uintptr_t>(env);

   state->paramsRef = new gaborator::parameters(bands_per_octave, ff_min / fs, ff_ref / fs, overlap, max_error);
   state->analyzerRef = new gaborator::analyzer<float>(*(state->paramsRef));
   state->coefsRef = new gaborator::coefs<float>(*(state->analyzerRef));

   //converts frequency (ff_max) in hertz to the number of bands above the min frequency
   //the ceil is used to end up at a full band 
   int interesting_bands = ceil(bands_per_octave * log(ff_max/ff_min)/log(2.0f));

   //since bands are ordered from high to low we are only interested in lower bands:
   //fs/2.0 is the nyquist frequency
   int total_bands = ceil(bands_per_octave * log(fs/2.0/ff_min)/log(2.0f));

   state->anal_support = (int64_t) ceil(state->analyzerRef->analysis_support());
   state->min_band = total_bands - interesting_bands;
   state->sample_rate = (int) fs;
   state->t_in = 0;
   state->cArray = new jfloat[C_ARRAY_SIZE];

   
   uintptr_t state_addresss = reinterpret_cast<uintptr_t>(state);
   assert(stateMap.count(env_addresss)==0);
   stateMap[env_addresss] = state_addresss;
   assert(stateMap.count(env_addresss)==1);
 

   assert(state->t_in == 0);

   return (int) state->anal_support;
}

/*
 * Class:     be_ugent_jgaborator_JGaborator
 * Method:    analyse
 * Signature: ([F)[F
 */
JNIEXPORT jfloatArray JNICALL Java_be_ugent_jgaborator_JGaborator_analyse(JNIEnv * env, jobject obj, jfloatArray audio_block ){

	//get a ref to the state pointer
   uintptr_t env_addresss = reinterpret_cast<uintptr_t>(env);
   assert(stateMap.count(env_addresss)==1);
   if(stateMap.count(env_addresss)==0){
      return NULL;
   }
   GaboratorState * state = reinterpret_cast<GaboratorState *>(stateMap[env_addresss]);

   //
   // Step 1: Convert the incoming JNI jintarray to C's jfloat[]
   jfloat *audio_block_c_array = env->GetFloatArrayElements(audio_block, NULL);
   if (NULL == audio_block_c_array) return NULL;
   jsize blocksize = env->GetArrayLength(audio_block);

   std::vector<float> buf(audio_block_c_array,audio_block_c_array + blocksize);

   //printf("Data analyis support %lld\n", state->anal_support);
   //printf("Audio block size: %d\n", (int) blocksize);
   //printf("Bands per octave:  %d\n", (int) state->paramsRef->bands_per_octave);
   //printf("t_in:  %lld\n", state->t_in);

   int output_index = 0;
   
   state->analyzerRef->analyze(buf.data(), state->t_in, state->t_in + blocksize, *(state->coefsRef));
   
   int64_t st0 = state->t_in - state->anal_support;
   int64_t st1 = state->t_in - state->anal_support + blocksize;

   
   
   apply(
            *state->analyzerRef, 
            *state->coefsRef,
            [&](std::complex<float> coef, int band, int64_t audioSampleIndex ) {
               //ignores everything above the max_band
               if(band >= state->min_band){
                  //printf("%f %d %ld\n",std::abs(coef),band,audioSampleIndex);
                  state->cArray[output_index++] = band;
                  state->cArray[output_index++] = audioSampleIndex;
                  state->cArray[output_index++] = std::abs(coef);
                  //printf("output_index:  %d\n", output_index++);
                  //output_index++;
               }
            },st0,
            st1);
   
   state->t_in += (int64_t) blocksize;
   
   //printf("After apply output_index %d\n", output_index);
   //printf("after apply Data analyis support %lld\n", state->anal_support);
   //printf("after apply Audio block size: %d\n", (int) blocksize);
   //printf("after apply Bands per octave:  %d\n", (int) state->paramsRef->bands_per_octave);
   //printf("after apply t_in:  %lld\n", state->t_in);

   int64_t t_out = state->t_in - state->anal_support;
   
   forget_before(*state->analyzerRef, *state->coefsRef, t_out - blocksize);
  
   //release audio block memory resources
   env->ReleaseFloatArrayElements(audio_block, audio_block_c_array, 0);
   //env->DeleteLocalRef(audio_block);
   //printf("out size: %d\n",output_index);

   jfloatArray outputArray = env->NewFloatArray(output_index);
   env->SetFloatArrayRegion(outputArray, 0 , output_index, state->cArray);  // copy


   //

   return outputArray;
}

JNIEXPORT jfloatArray JNICALL Java_be_ugent_jgaborator_JGaborator_bandcenters(JNIEnv * env, jobject obj){

   uintptr_t env_addresss = reinterpret_cast<uintptr_t>(env);
   assert(stateMap.count(env_addresss)==1);

   if(stateMap.count(env_addresss)==0){
      return NULL;
   }
   GaboratorState * state = reinterpret_cast<GaboratorState *>(stateMap[env_addresss]);
   
   int max_band = state->analyzerRef->bandpass_bands_end();
   float band_centers[max_band+1];

   for(int i = 0 ; i < max_band ; i++){
      if(i<state->min_band){
         band_centers[i]=-1;
      }else{
         band_centers[i]=state->analyzerRef->band_ff(i) * state->sample_rate;
      }
   }

   jfloatArray outputArray = env->NewFloatArray(max_band+1);
   env->SetFloatArrayRegion(outputArray, 0 , max_band+1, band_centers);  // copy
   return outputArray;
}

/*
 * Class:     be_ugent_jgaborator_JGaborator
 * Method:    release
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_be_ugent_jgaborator_JGaborator_release(JNIEnv *env , jobject obj){
   
   std::unique_lock<std::mutex> lck (stateMutex);

   uintptr_t env_addresss = reinterpret_cast<uintptr_t>(env);
   assert(stateMap.count(env_addresss)==1);
   if(stateMap.count(env_addresss)==0){
      return;
   }
   GaboratorState * state = reinterpret_cast<GaboratorState *>(stateMap[env_addresss]);
   
   assert(stateMap.count(env_addresss)==1);
   stateMap.erase(env_addresss); 
   assert(stateMap.count(env_addresss)==0);

	//cleanup memory
   //delete state->jniArray;
   //delete cArray;
	delete state->analyzerRef;
	delete state->coefsRef;
	delete state->paramsRef;
   delete [] state->cArray;

	delete state;
}

int main(){
   return 0;
}
