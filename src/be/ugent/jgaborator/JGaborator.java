package be.ugent.jgaborator;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.FloatBuffer;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import be.tarsos.dsp.AudioEvent;
import be.tarsos.dsp.AudioProcessor;

public class JGaborator implements AudioProcessor{

	final float[] bandcenterCache;
	final int firstBandCache;

	final double overlap = 1.0;//Gaborator library default, exposed here if change is needed
	final double minError = 1e-5;//Gaborator library default, exposed here if change is needed

	final double sampleRate;
	final int blocksize;
	final int stepSize;
	final int bandsPerOctave;

	final float[] audioDataToTransform;
	
	float[][] coefficients;
	int coefficientOffset;

	// Loads the jgaborator library
	static {
		// Load native library at runtime
		System.loadLibrary("jgaborator"); 
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
		
		initialize(blocksize, samplerate, bandsPerOctave, minimumFrequency, referenceFrequency, maximumFrequency,
				overlap, minError);
		this.sampleRate = samplerate;
		this.blocksize = blocksize;
		bandcenterCache = bandcenters();
		firstBandCache = firstBand();
		
		audioDataToTransform = new float[blocksize];
		this.stepSize = stepSize;
		this.bandsPerOctave = bandsPerOctave;
		coefficients = new float[(int)(60*samplerate/stepSize)][numberOfBands()];//start with 60seconds
		coefficientOffset = 0;
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

		// System.out.println(analysisResult.length);
		for (int i = 0; i < analysisResult.length; i += 3) {

			int audioSample = (int) analysisResult[i + 1];
			int band = (int) analysisResult[i];
			float coeficcient = analysisResult[i + 2];

			int coefficientIndex = audioSample / stepSize - coefficientOffset;
			int bandIndex = band - firstBandCache;
			
			// Expand the coefficent array if no more place is available (double the size)
			if(coefficientIndex >= coefficients.length) {				
				float[][] moreCoefficents = new float[coefficients.length + coefficients.length][coefficients[0].length];
				System.arraycopy(coefficients, 0, moreCoefficents, 0, coefficients.length);
				coefficients = moreCoefficents;
			}
			
			// The first results have a negative audio sample index
						// ignore these
			if (coefficientIndex > 0 &&  bandIndex < coefficients[coefficientIndex].length )
				coefficients[coefficientIndex][bandIndex] = Math.max(coefficients[coefficientIndex][bandIndex], coeficcient);
		}
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

		for (int floatIndex = 0; floatIndex < fileContents.length / 4 - blocksize; floatIndex += blocksize) {
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
		//Goes from a sparse array to a filled array 
		for(int i = 1 ; i < coefficients.length ; i++) {
			for(int j = 0 ; j < coefficients[i].length; j++) {
				if(coefficients[i][j] == 0) {
					coefficients[i][j] = coefficients[i-1][j];
				}
			}
		}
		release();
	}
	

	public int getStepSize() {
		return stepSize;
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

	public float[][] getCoefficents() {
		return coefficients;
	}

	
}