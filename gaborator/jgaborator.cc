#include "jgaborator.h"
#include <math.h>

typedef struct GaboratorData {   
	gaborator::parameters* paramsRef;
	gaborator::analyzer<float>* analyzerRef;
	gaborator::coefs<float>* coefsRef;
	
	int64_t t_in;
	int min_band;
	int sample_rate;
	int64_t anal_support;

   float* cArray;

   jfloatArray jniArray;
   
} gaborator_data;

jfieldID getPtrFieldId(JNIEnv * env, jobject obj){
    static jfieldID ptrFieldId = 0;
    if (!ptrFieldId){
        jclass c = env->GetObjectClass(obj);
        ptrFieldId = env->GetFieldID(c, "ptr", "J");
        env->DeleteLocalRef(c);
    }
    return ptrFieldId;
}


JNIEXPORT jint JNICALL Java_be_ugent_jgaborator_JGaborator_initialize(JNIEnv * env, jobject object, jint blocksize, jdouble fs, jint bands_per_octave, jdouble ff_min , jdouble ff_ref, jdouble ff_max, jdouble overlap, jdouble max_error){
	gaborator_data * data =  (gaborator_data *) malloc (sizeof(gaborator_data *));
	
	//https://stackoverflow.com/questions/24491965/store-a-c-object-instance-inside-jni-jobject-and-retrieve-later
	
	data->paramsRef = new gaborator::parameters(bands_per_octave, ff_min / fs, ff_ref / fs, overlap, max_error);
	data->analyzerRef = new gaborator::analyzer<float>(*(data->paramsRef));
	data->coefsRef = new gaborator::coefs<float>(*(data->analyzerRef));

   data->cArray =  (float *) malloc (sizeof(float) * 130000 * 32);
	
	//converts frequency (ff_max) in hertz to the number of bands above the min frequency
	//the ceil is used to end up at a full band 
	int interesting_bands = ceil(bands_per_octave * log(ff_max/ff_min)/log(2.0f));
	
	//since bands are ordered from high to low we are only interested in lower bands:
	//fs/2.0 is the nyquist frequency
	int total_bands = ceil(bands_per_octave * log(fs/2.0/ff_min)/log(2.0f));
	
   data->anal_support = (int64_t) ceil(data->analyzerRef->analysis_support());
	data->min_band = total_bands - interesting_bands;
	data->sample_rate = (int) fs;
	data->t_in = 0;

   data->jniArray = NULL;
	//save the pointer and set it as an object instance variable
	env->SetLongField(object, getPtrFieldId(env, object), (jlong) data);
	
	return (int) data->anal_support;
}

/*
 * Class:     be_ugent_jgaborator_JGaborator
 * Method:    analyse
 * Signature: ([F)[F
 */
JNIEXPORT jfloatArray JNICALL Java_be_ugent_jgaborator_JGaborator_analyse(JNIEnv * env, jobject obj, jfloatArray audio_block ){

	//get a ref to the data pointer
   gaborator_data * data = (gaborator_data *) env->GetLongField(obj, getPtrFieldId(env, obj));

   //

   // Step 1: Convert the incoming JNI jintarray to C's jfloat[]
   jfloat *audio_block_c_array = env->GetFloatArrayElements(audio_block, NULL);
   if (NULL == audio_block_c_array) return NULL;
   jsize blocksize = env->GetArrayLength(audio_block);

   std::vector<float> buf(audio_block_c_array,audio_block_c_array + blocksize);

   //printf("Data analyis support %lld\n", data->anal_support);
   //printf("Audio block size: %d\n", (int) blocksize);
   //printf("Bands per octave:  %d\n", (int) data->paramsRef->bands_per_octave);
   //printf("t_in:  %lld\n", data->t_in);

   int output_index = 0;
   
   data->analyzerRef->analyze(buf.data(), data->t_in, data->t_in + blocksize, *(data->coefsRef));
   
   int64_t st0 = data->t_in - data->anal_support;
   int64_t st1 = data->t_in - data->anal_support + blocksize;
 
   apply(
            *data->analyzerRef, 
            *data->coefsRef,
            [&](std::complex<float> coef, int band, int64_t audioSampleIndex ) {
               //ignores everything above the max_band
               if(band >= data->min_band){
                  //printf("%f %d %ld\n",std::abs(coef),band,audioSampleIndex);
                  data->cArray[output_index++] = band;
                  data->cArray[output_index++] = audioSampleIndex;
                  data->cArray[output_index++] = std::abs(coef);
                  //printf("output_index:  %d\n", output_index++);
                  //output_index++;
               }
            },st0,
            st1);
  
   data->t_in += (int64_t) blocksize;
   

   //printf("after apply output_index %d\n", output_index);
   //printf("after apply Data analyis support %lld\n", data->anal_support);
   //printf("after apply Audio block size: %d\n", (int) blocksize);
   //printf("after apply Bands per octave:  %d\n", (int) data->paramsRef->bands_per_octave);
   //printf("after apply t_in:  %lld\n", data->t_in);

   int64_t t_out = data->t_in - data->anal_support;
   
   forget_before(*data->analyzerRef, *data->coefsRef, t_out - blocksize );
  
   //release audio block memory resources
   env->ReleaseFloatArrayElements(audio_block, audio_block_c_array, 0);

   //see here: https://www3.ntu.edu.sg/home/ehchua/programming/java/JavaNativeInterface.html
   if(data->jniArray != NULL){
      //cleanup
       env->DeleteGlobalRef(data->jniArray);
   }
	jfloatArray tempOutputArray = env->NewFloatArray(output_index);
   data->jniArray = (jfloatArray)env->NewGlobalRef(tempOutputArray);
   env->SetFloatArrayRegion(data->jniArray, 0 , output_index, data->cArray);  // copy

   return data->jniArray;
}

JNIEXPORT jfloatArray JNICALL Java_be_ugent_jgaborator_JGaborator_bandcenters(JNIEnv * env, jobject obj){

   //get a ref to the data pointer
   gaborator_data * data = (gaborator_data *) env->GetLongField(obj, getPtrFieldId(env, obj));
   
   int max_band = data->analyzerRef->bandpass_bands_end();
   float band_centers[max_band+1];

   for(int i = 0 ; i < max_band ; i++){
      if(i<data->min_band){
         band_centers[i]=-1;
      }else{
         band_centers[i]=data->analyzerRef->band_ff(i) * data->sample_rate;
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
JNIEXPORT void JNICALL Java_be_ugent_jgaborator_JGaborator_release(JNIEnv *env , jobject obj){
   //get a ref to the data pointer
   gaborator_data * data = (gaborator_data *) env->GetLongField(obj, getPtrFieldId(env, obj));


   if(data->jniArray != NULL){
      //cleanup
      env->DeleteGlobalRef(data->jniArray);
   }
   

	//cleanup memory
   delete data->cArray;
	delete data->analyzerRef;
	delete data->coefsRef;
	delete data->paramsRef;
	delete data;
}

int main(){
   return 0;
}
