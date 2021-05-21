package be.ugent.jgaborator.ui;

import java.awt.Color;
import java.awt.Graphics2D;
import java.util.List;

import be.tarsos.dsp.ui.Axis;
import be.tarsos.dsp.ui.CoordinateSystem;
import be.tarsos.dsp.ui.layers.Layer;
import be.tarsos.dsp.util.PitchConverter;
import be.ugent.jgaborator.JGaborator;

public class GaborLayer implements Layer {
	
	private final CoordinateSystem cs;
	
	List<float[]> coefficents;
	JGaborator jgaborator;
	
	
	public GaborLayer(CoordinateSystem cs, JGaborator jgaborator) {
		
		this.coefficents = jgaborator.getCoefficents();
		this.jgaborator = jgaborator;
		this.cs = cs;
	}
	
	private float maxInRegion(int timeIndex) {
		float currentMax = -10;
		for(int i = Math.max(timeIndex - 5,0) ; i < Math.min(timeIndex +5, coefficents.size()) ; i++) {
			float[] magnitudes = coefficents.get(i);
    		// draw the pixels
			for (int j = 0; j < magnitudes.length; j++) {
				currentMax = Math.max(currentMax, magnitudes[j]);
			}
		}
		return currentMax;
	}

	@Override
	public void draw(Graphics2D graphics) {
	
		int stepSize = jgaborator.getStepSize();
		float sampleRate = jgaborator.getSampleRate();
		float bandWidthInCents = jgaborator.getBandWidth();
		float frameDurationInMS = stepSize/sampleRate * 1000; 
		
		float startTimeInSeconds = cs.getMin(Axis.X) / 1000.0f; 
	    float stopTimeInSeconds = cs.getMax(Axis.X) / 1000.0f;
	    
	    int startTimeInSamples = ((int) Math.floor(startTimeInSeconds * sampleRate / stepSize)) * stepSize;
	    int stopTimeInSamples = ((int) Math.floor(stopTimeInSeconds * sampleRate / stepSize)) * stepSize;
	    
	    for(int timeInSamples = startTimeInSamples ; timeInSamples <= stopTimeInSamples ; timeInSamples += stepSize) {
	    	int timeIndex = timeInSamples/stepSize;
	    	if(timeIndex > 0 && timeIndex < coefficents.size()) {
	    		float[] magnitudes = coefficents.get(timeIndex);
	    		float maxInRegion = maxInRegion(timeIndex);
	    		
	    		// draw the pixels
				for (int i = 0; i < magnitudes.length; i++) {
					Color color = Color.black;
					
					//actual energy at frame.frequencyEstimates[i];
					
					//remove 5 cents to make sure bands overlap (and no white stripes show)
					float centsStartingPoint = (float) PitchConverter.hertzToAbsoluteCent(jgaborator.bandCenters(i)) - bandWidthInCents/2 - 5;
					
					// only draw the visible frequency range
					if (centsStartingPoint >= cs.getMin(Axis.Y)
							&& centsStartingPoint <= cs.getMax(Axis.Y)) {
						
						if(magnitudes[i]==0) {
							for(int prevIndex = 0 ; prevIndex < 10 * stepSize ; prevIndex += 1) {
								int prevTimeIndex = timeIndex-prevIndex;
								if( prevTimeIndex >= 0 && prevTimeIndex < coefficents.size() && coefficents.get(prevTimeIndex)[i] != 0) {
									magnitudes[i]=coefficents.get(prevTimeIndex)[i] * ((100-prevIndex)/100.0f);
									break;
								}
							}
						}
					
						int greyValue = 255 - (int) (magnitudes[i] / maxInRegion  * 255);
						greyValue = Math.max(0, greyValue);
						greyValue = Math.min(255, greyValue);
						color = new Color(greyValue, greyValue, greyValue);
						graphics.setColor(color);
						if(greyValue < 253)
						graphics.fillRect((int) Math.round(timeInSamples / sampleRate * 1000.0 ),
								Math.round(centsStartingPoint),
								(int) Math.round(frameDurationInMS),
								(int) bandWidthInCents + 5); //remove 5 cents to make sure bands overlap (and no white stripes show)
					}
				}
	    	}	
	    }
	}
	

	@Override
	public String getName() {		
		return "Gabor Layer";
	}
}
