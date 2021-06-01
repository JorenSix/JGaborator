#include "jgaborator.h"
#include <unordered_map>

#include <chrono>
#include <thread>

struct GaboratorState { 	
	int64_t t_in;
	int min_band;
	int sample_rate;
	int64_t anal_support;

   int block;
};


std::unordered_map<uintptr_t, uintptr_t> stateMap;

std::mutex stateMutex;

const int C_ARRAY_SIZE = 300000;

JNIEXPORT jint JNICALL Java_be_ugent_jgaborator_JGaborator_initialize(JNIEnv * env, jobject object, jint blocksize, jdouble fs, jint bands_per_octave, jdouble ff_min , jdouble ff_ref, jdouble ff_max, jdouble overlap, jdouble max_error){
	
   //Makes sure only one thread writes to the stateMap
   std::unique_lock<std::mutex> lck (stateMutex);

   GaboratorState * state =  new GaboratorState();
   uintptr_t env_addresss = reinterpret_cast<uintptr_t>(env);
   
   state->anal_support = (int64_t) 16000;
	state->min_band = 64;
	state->sample_rate = (int) fs;
	state->t_in = 0;
   state->block= 0;

   std::chrono::milliseconds timespan((long)rand() % 1000); // or whatever
   std::this_thread::sleep_for(timespan);

   uintptr_t state_addresss = reinterpret_cast<uintptr_t>(state);
   assert(stateMap.count(env_addresss)==0);
   stateMap[env_addresss] = state_addresss;
   assert(stateMap.count(env_addresss)==1);
   
   assert(state->block == 0);

   printf("INITIAL : env: %p\n",env);
	return (int) 64;
}

/*
 * Class:     be_ugent_jgaborator_JGaborator
 * Method:    analyse
 * Signature: ([F)[F
 */
JNIEXPORT jfloatArray JNICALL Java_be_ugent_jgaborator_JGaborator_analyse(JNIEnv * env, jobject obj, jfloatArray input_jni_arrray){


	//get a ref to the state pointer
   uintptr_t env_addresss = reinterpret_cast<uintptr_t>(env);
   assert(stateMap.count(env_addresss)==1);
   if(stateMap.count(env_addresss)==0){
      return NULL;
   }
   GaboratorState * state = reinterpret_cast<GaboratorState *>(stateMap[env_addresss]);
 
   // Convert the incoming JNI jintarray to C's jfloat[]
   jfloat *input_array = env->GetFloatArrayElements(input_jni_arrray, NULL);
   if (NULL == input_array) return NULL;
   jsize array_length = env->GetArrayLength(input_jni_arrray);

   //create a new array of the same length
   
   //copy the contents

   jfloat *output_array = new jfloat[C_ARRAY_SIZE];

   for(size_t i = 0 ; i < C_ARRAY_SIZE ; i +=3 ){
      output_array[i+0]=i+0;
      output_array[i+1]=i+1;
      output_array[i+2]=i+2;
   }

   std::chrono::milliseconds timespan((long)rand() % 300); // or whatever
   std::this_thread::sleep_for(timespan);
  
   //release audio block memory resources
   env->ReleaseFloatArrayElements(input_jni_arrray, input_array, 0);

   //printf("ANALYSIS: env: %p data: %p data->block:%d \n",env,data,data->block);
   printf("ANALYSIS: env: %p\n",env);
   //change the state
   state->block++;

   printf("ANALYSIS: data: %p data->block:%d \n",state,state->block);

   //copy output to JNI array
   jfloatArray outputArray = env->NewFloatArray(C_ARRAY_SIZE);
   env->SetFloatArrayRegion(outputArray, 0 , C_ARRAY_SIZE, output_array);

   //delete [] output_array;

   return outputArray;
}

JNIEXPORT jfloatArray JNICALL Java_be_ugent_jgaborator_JGaborator_bandcenters(JNIEnv * env, jobject obj){

   

   uintptr_t env_addresss = reinterpret_cast<uintptr_t>(env);
   assert(stateMap.count(env_addresss)==1);
   if(stateMap.count(env_addresss)==0){
      return NULL;
   }
   GaboratorState * state = reinterpret_cast<GaboratorState *>(stateMap[env_addresss]);

   int max_band = 65;
   float band_centers[max_band+1];

   for(int i = 0 ; i< max_band+1  ; i++){
      band_centers[i]=i;
   }

   std::chrono::milliseconds timespan((long)rand() % 300); // or whatever
   std::this_thread::sleep_for(timespan);

   //printf("CENTERS : env: %p data: %p data->block:%d \n",env,data,data->block);

   printf("CENTERS : data: %p data->block:%d \n",state,state->block);
   
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
   //Makes sure only one thread writes to the stateMap
   std::unique_lock<std::mutex> lck (stateMutex);
   
   uintptr_t env_addresss = reinterpret_cast<uintptr_t>(env);
   assert(stateMap.count(env_addresss)==1);
   if(stateMap.count(env_addresss)==0){
      return;
   }
   GaboratorState * state = reinterpret_cast<GaboratorState *>(stateMap[env_addresss]);

   std::chrono::milliseconds timespan((long)rand() % 300); // or whatever
   std::this_thread::sleep_for(timespan);
   
   assert(stateMap.count(env_addresss)==1);
   stateMap.erase(env_addresss); 
   assert(stateMap.count(env_addresss)==0);

   printf("RELEASE : env: %p\n",env);
   printf("RELEASE : data: %p data->block:%d \n",state,state->block);
   delete state;
}

int main(){
   return 0;
}
