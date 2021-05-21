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
import be.ugent.jgaborator.util.NativeUtils;

public class JGaborator implements AudioProcessor{

	final float[] bandcenterCache;
	final int firstBandCache;

	final double overlap = 1.0;//Gaborator library default, exposed here if change is needed
	final double minError = 1e-5;//Gaborator library default, exposed here if change is needed

	final double sampleRate;
	final int audioBlockSize;
	final int frequencyBinTimeStepSize;
	final int bandsPerOctave;
	final int latency;//processing latency in audio samples
	
	final float[] audioDataToTransform;
	
	private final float[][] coefficients;//circular buffer with current coefficents
	int coefficientIndexOffset;//The offset (in steps)
	int mostRecentCoefficentIndex;//The index of the most recent coefficent (in steps)
	
	//a history with calculated coefficents
	private final List<float[]> fixedCoefficents;
	
	// Loads the jgaborator library
	static {
		// Load native library at runtime
		try {
			System.loadLibrary("jgaborator"); 
		}catch (UnsatisfiedLinkError e){
			System.err.println("Could not load jgaborator JNI library. Will attempt to use a version packed in the JAR archive");
			System.err.println("  info : " + e.getMessage());
			try {
				NativeUtils.loadLibraryFromJar("/jni/" + System.mapLibraryName("jgaborator"));
				System.err.println("Loaded JNI jgaborator library from JAR archive.");
			} catch (IOException e1) {
				
				e1.printStackTrace();
			}
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
		
		latency = initialize(blocksize, samplerate, bandsPerOctave, minimumFrequency, referenceFrequency, maximumFrequency,
				overlap, minError);
		this.sampleRate = samplerate;
		this.audioBlockSize = blocksize;
		bandcenterCache = bandcenters();
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

	// Analyze a block of audio and get constant-q transform of previous blocks
	private native float[] analyse(float[] block);

	// Return the center frequencies in Hz for each band
	private native float[] bandcenters();

	// Release occupied resources
	private native void release();

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

	public float bandCenters(int bandIndex) {
		return bandcenterCache[bandIndex + firstBandCache];
	}
	
	private void gaborTransform(float[] audioData) {
		// Analyze an audio block
		float[] analysisResult = analyse(audioData);

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
		//TODO: copy the final coefficents to the history arrayList
	
		release();
	}
	

	public int getStepSize() {
		return frequencyBinTimeStepSize;
	}

	public float getSampleRate() {
		return (float)sampleRate;
	}

	/**
	 * @return The width of each band in cents.
	 */
	public float getBandWidth() {
		return  1200.0f/(float) bandsPerOctave;
	}

	public int getLatency() {
		return latency;
		
	}
	
}