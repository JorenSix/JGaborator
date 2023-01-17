package be.ugent.jgaborator;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.FloatBuffer;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import be.tarsos.dsp.AudioEvent;
import be.tarsos.dsp.AudioProcessor;
import be.ugent.jgaborator.util.ZigNativeUtils;

/**
 * The TarsosDSP audio processor to calculate the Gabor transform.
 */
public class JGaborator implements AudioProcessor{
	
	private final float[] bandcenterCache;
	private final int firstBandCache;

	private final double overlap = 1.0;//Gaborator library default, exposed here if change is needed
	private final double minError = 1e-5;//Gaborator library default, exposed here if change is needed

	private final double sampleRate;
	private final int audioBlockSize;
	private final int frequencyBinTimeStepSize;
	private final int bandsPerOctave;
	private final int latency;//processing latency in audio samples

	private final float[] audioDataToTransform;
	
	private final float[][] coefficients;//circular buffer with current coefficents
	private int coefficientIndexOffset;//The offset (in steps)
	private int mostRecentCoefficentIndex;//The index of the most recent coefficent (in steps)
	
	//a history with calculated coefficents
	private final List<float[]> fixedCoefficents;
	
	// Loads the jgaborator library
	static {
		// Load native library at runtime
		try {
			System.loadLibrary("jgaborator");
			System.err.println("Loaded jgaborator library from " +  System.getProperty("java.library.path"));
		}catch (UnsatisfiedLinkError e ){
			System.err.println("Could not load 'jgaborator' JNI library. \n Will attempt to use a precompiled version packed in the JAR archive\n" + e.getMessage());
			boolean libraryLoaded = ZigNativeUtils.loadLibraryFromJarWithOSDetection("/jni/" + System.mapLibraryName("jgaborator"));
		}

	}



	/**
	 * @param blocksize The size of a block of audio
	 * @param samplerate The sample rate of the audio in Hz. For example 8000Hz or 44100Hz
	 * @param bandsPerOctave The number of bands per octave. 12 is reasonable for eurogenetic music.
	 * @param minimumFrequency The minimum frequency to return in Hz. E.g. 110Hz.
	 * @param maximumFrequency The maximum frequency in Hz. If a value above the Nyquist frequency is given the maximum frequency is equal to the nyquist. 
	 * @param referenceFrequency To center the bins on a certain frequency given in Hz. For example 440Hz is common for eurogenetic music.
	 * @param stepSize The resolution of the resulting transform in audio samples. It needs to be a power of two and not too low otherwise memory use explodes. For example 32, 64 or 128.
	 * Memory usage is equal to samplerate/stepsize  * audio duration in seconds * 4 bytes per coefficient (float) * number of bands. 
	 */
	public JGaborator(int blocksize, double samplerate, int bandsPerOctave, double minimumFrequency,
			double maximumFrequency, double referenceFrequency,int stepSize) {


		//System.out.println("init " + Thread.currentThread());
		latency = doInitialization(blocksize, samplerate, bandsPerOctave, minimumFrequency, referenceFrequency, maximumFrequency,
				overlap, minError);
		//System.out.println("after init ");
		this.sampleRate = samplerate;
		this.audioBlockSize = blocksize;
		
		//System.out.println("Band centers ");
		bandcenterCache = getBandcenters();
		//System.out.println("After Band centers " + bandcenterCache.length);
		firstBandCache = firstBand();
		
		audioDataToTransform = new float[blocksize];
		
		fixedCoefficents=new ArrayList<float[]>();
				
		this.frequencyBinTimeStepSize = stepSize;
		this.bandsPerOctave = bandsPerOctave;
		
		int coefficentsSize =  (latency + 2*blocksize) / stepSize;
		coefficients = new float[coefficentsSize][numberOfBands()];
		coefficientIndexOffset = 0;
	}

	// Initializes the constant-q transform according to these parameters
	private native int initialize(int blockSize, double samplerate, int bandsPerOctave, double fMin, double fRef,
			double fmax, double overlap, double minError);
	
	private synchronized int doInitialization(int blockSize, double samplerate, int bandsPerOctave, double fMin, double fRef,
			double fmax, double overlap, double minError) {
		return initialize(blockSize,samplerate,bandsPerOctave,fMin,fRef,fmax,overlap,minError);
	}

	// Analyze a block of audio and get constant-q transform of previous blocks
	private native float[] analyse(float[] block);

	// Return the center frequencies in Hz for each band
	private native float[] bandcenters();
	
	private synchronized float[] getBandcenters(){
		return bandcenters();
	}

	// Release occupied resources
	private native void release();
	
	private synchronized void doRelease() {
		release();
	}

	private int numberOfBands() {
		float[] bandcenters = bandcenterCache;
		int numberOfBands = 0;
		for (int i = 0; i < bandcenters.length; i++) {
			if (bandcenters[i] > 0)
				numberOfBands++;
		}
		return numberOfBands;
	}

	private int firstBand() {
		float[] bandcenters = bandcenterCache;
		for (int i = 0; i < bandcenters.length; i++) {
			if (bandcenters[i] > 0)
				return i;
		}
		return -1;
	}

	/**
	 * The center for the frequency band
	 * @param bandIndex The index of the frequency band.
	 * @return The frequency of the center place of the band.
	 */
	public float bandCenters(int bandIndex) {
		return bandcenterCache[bandIndex + firstBandCache];
	}
	
	private synchronized void gaborTransform(float[] audioData) {
		// Analyze an audio block
		//System.out.println("Before analyse " + Thread.currentThread());
		float[] analysisResult = analyse(audioData);
		//System.out.println("After analyse");

		//The analysis result consists of a float array with three values:
		// a frequency band index [i] (always an integer)
		// an audio sample index [i+1] (expressed in audio samples)
		// a magnitude value [i+2] (the magnitude value)
		for (int i = 0; i < analysisResult.length; i += 3) {

			int band = (int) analysisResult[i];
			int audioSample = (int) analysisResult[i + 1];
			
			float coeficcient = analysisResult[i + 2];

			int coefficientIndex = audioSample / frequencyBinTimeStepSize - coefficientIndexOffset;
			int bandIndex = band - firstBandCache;
			
			int circularIndex = coefficientIndex % coefficients.length;
			
			// The first results have a negative audio sample index
			// ignore these
			if(coefficientIndex > 0 && bandIndex < coefficients[circularIndex].length ) {
				
				// If a new index is reached, save the old (fixed) coefficents in the history
				// Fill the array with zeros to get the max
				if(coefficientIndex > mostRecentCoefficentIndex && coefficientIndex > coefficients.length) {
					// keep the new maximum
					mostRecentCoefficentIndex = coefficientIndex;
					// copy the oldest data to the history
					fixedCoefficents.add(coefficients[circularIndex].clone());
					// fill the oldest with zeros
					Arrays.fill(coefficients[circularIndex], 0.0f);
				}
				// due to reduction in precision (from audio sample accuracy to steps) multiple 
				// magnitudes could be placed in the same stepIndex, bandIndex pair.
				// We take the maximum magnitudes value. 
				coefficients[circularIndex][bandIndex] = Math.max(coefficients[circularIndex][bandIndex], coeficcient);
			}

		}
	}

	/**
	 * Return the spectral coefficents.
	 * @return Return the spectral coefficents.
	 */
	public List<float[]> getCoefficents() {
		return fixedCoefficents;
		//ArrayList<float[]> fixedCoefficentsCopy = new ArrayList<>(fixedCoefficents);
		
		/*Goes from a sparse array to a filled array 
		for(int i = 1 ; i < fixedCoefficentsCopy.size() ; i++) {
			for(int j = 0 ; j < fixedCoefficentsCopy.get(i).length; j++) {
				if(fixedCoefficentsCopy.get(i)[j] == 0) {
					fixedCoefficentsCopy.get(i)[j] = fixedCoefficentsCopy.get(i-1)[j];
				}
			}
		}
		*/
		
		//fixedCoefficents.clear();
		//return fixedCoefficentsCopy;
	}

	/**
	 * Transform a raw audio file with float 32 bit little-endian encoded audio
	 * samples without a header.
	 * 
	 * @param fileName
	 *            The name of the file to transform
	 * @return The resulting Gabor coefficients and timing
	 * @throws IOException
	 *             if the file is not readable
	 */
	public float[][] gaborTransform(String fileName) throws IOException {
		Path path = Paths.get(fileName);
		byte[] fileContents = Files.readAllBytes(path);

		ByteBuffer audioFileDataBytes = ByteBuffer.wrap(fileContents);
		audioFileDataBytes.order(ByteOrder.LITTLE_ENDIAN);
		FloatBuffer audioFileData = audioFileDataBytes.asFloatBuffer();

		for (int floatIndex = 0; floatIndex < fileContents.length / 4 - audioBlockSize; floatIndex += audioBlockSize) {
			audioFileData.get(audioDataToTransform);
			gaborTransform(audioDataToTransform);
		}
		processingFinished();
		return coefficients;
	}
	
	@Override
	public boolean process(AudioEvent audioEvent) {
		float[] audioData = audioEvent.getFloatBuffer();
		gaborTransform(audioData);
		return true;
	}

	@Override
	public void processingFinished() {
		//System.out.println("Release "  + Thread.currentThread());
		//TODO: copy the final coefficents to the history arrayList
		doRelease();
		//System.out.println("after Release ");
	}

	/**
	 * The step size between in samples.
	 * @return  The step size between in samples.
	 */
	public int getStepSize() {
		return frequencyBinTimeStepSize;
	}

	/**
	 * The audio sample rate
	 * @return The audio sample rate in Hz
	 */
	public float getSampleRate() {
		return (float)sampleRate;
	}

	/**
	 * The with of each frequency band in cents.
	 * @return The width of each frequency band in cents.
	 */
	public float getBandWidth() {
		return  1200.0f/(float) bandsPerOctave;
	}

	/**
	 * The algorithmic latency in samples
	 * @return The algorithmic latency in samples
	 */
	public int getLatency() {
		return latency;
		
	}
	
}