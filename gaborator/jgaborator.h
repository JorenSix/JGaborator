
#include <jni.h>
#include <stdio.h>

#include <memory.h>
#include <iostream>



#ifndef _Included_be_ugent_jgaborator_JGaborator
#define _Included_be_ugent_jgaborator_JGaborator
#ifdef __cplusplus
extern "C" {
#endif


/*
 * Class:     be_ugent_jgaborator_JGaboratorTest
 * Method:    initialize
 * Signature: (IDIDDDDD)V
 */
JNIEXPORT jint JNICALL Java_be_ugent_jgaborator_JGaborator_initialize
  (JNIEnv *, jobject, jint, jdouble, jint, jdouble, jdouble, jdouble, jdouble, jdouble);

/*
 * Class:     be_ugent_jgaborator_JGaborator
 * Method:    analyse
 * Signature: ([F)[F
 */
JNIEXPORT jfloatArray JNICALL Java_be_ugent_jgaborator_JGaborator_analyse
  (JNIEnv *, jobject, jfloatArray);

JNIEXPORT jfloatArray JNICALL Java_be_ugent_jgaborator_JGaborator_bandcenters
  (JNIEnv *, jobject);  

/*
 * Class:     be_ugent_jgaborator_JGaborator
 * Method:    release
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_be_ugent_jgaborator_JGaborator_release
  (JNIEnv *, jobject);

#ifdef __cplusplus
}
#endif
#endif
